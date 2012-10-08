/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/


// C++ includes
#include <numeric>
#include <vector>

#include "asinh.hpp" // for asinh

// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "resistance_region.h"
#include "insulator_region.h"
#include "boundary_condition_resistance_ohmic.h"
#include "parallel.h"
#include "mathfunc.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


void IF_Metal_OhmicBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag )
{
  std::complex <PetscScalar> j(0.0, 1.0);
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *> ( _r2 );

  // we set Jacobian entries for ohmic boundary condition
  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;


    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data = resistance_node->node_data();

    // find the node variable offset
    unsigned int semiconductor_n_node_var   = semiconductor_region->ebm_n_variables();
    unsigned int resistance_n_node_var      = resistance_region->ebm_n_variables();


    unsigned int node_psi_offset = semiconductor_region->ebm_variable_offset(POTENTIAL);
    unsigned int node_n_offset   = semiconductor_region->ebm_variable_offset(ELECTRON);
    unsigned int node_p_offset   = semiconductor_region->ebm_variable_offset(HOLE);
    unsigned int node_Tl_offset  = semiconductor_region->ebm_variable_offset(TEMPERATURE);
    unsigned int node_Tn_offset  = semiconductor_region->ebm_variable_offset(E_TEMP);
    unsigned int node_Tp_offset  = semiconductor_region->ebm_variable_offset(H_TEMP);

    //the indepedent variable number
    adtl::AutoDScalar::numdir = semiconductor_n_node_var + resistance_n_node_var;
    //synchronize with material database
    semiconductor_region->material()->set_ad_num(adtl::AutoDScalar::numdir);
    semiconductor_region->material()->mapping(semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock);

    PetscInt ad_index = 0;
    AutoDScalar V = semiconductor_node_data->psi();  V.setADValue(ad_index++, 1.0);  // psi of this node
    AutoDScalar n = semiconductor_node_data->n();    n.setADValue(ad_index++, 1.0);  // electron density
    AutoDScalar p = semiconductor_node_data->p();    p.setADValue(ad_index++, 1.0);  // hole density

    AutoDScalar T  =  T_external();
    AutoDScalar Tn =  T_external();
    AutoDScalar Tp =  T_external();

    // lattice temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tl())
    {
      T =  semiconductor_node_data->T();
      T.setADValue(ad_index++, 1.0);
    }

    // electron temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tn())
    {
      AutoDScalar nTn = semiconductor_node_data->n()*semiconductor_node_data->Tn();
      nTn.setADValue(ad_index++, 1.0);
      Tn = nTn/n;
    }

    // hole temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tp())
    {
      AutoDScalar pTp = semiconductor_node_data->p()*semiconductor_node_data->Tp();
      pTp.setADValue(ad_index++, 1.0);
      Tp = pTp/p;
    }

    // the electrode potential in current iteration
    AutoDScalar V_resistance = resistance_node_data->psi();        V_resistance.setADValue(ad_index++, 1.0);
    AutoDScalar T_resistance = T_external();
    // lattice temperature if required
    if(resistance_region->get_advanced_model()->enable_Tl())
    {
      T_resistance =  resistance_node_data->T();
      T_resistance.setADValue(ad_index++, 1.0);
    }

    // the insert position
    std::vector<PetscInt> row_re, col_re;
    std::vector<PetscInt> row_im, col_im;
    for(unsigned int nv=0; nv<semiconductor_n_node_var; ++nv)
    {
      row_re.push_back(semiconductor_node->global_offset() + nv);
      row_im.push_back(semiconductor_node->global_offset() + semiconductor_n_node_var + nv);
    }

    col_re = row_re;
    col_im = row_im;

    for(unsigned int nv=0; nv<resistance_n_node_var; ++nv)
    {
      col_re.push_back(resistance_node->global_offset() + nv);
      col_im.push_back(resistance_node->global_offset() + resistance_n_node_var + nv);
    }

    AutoDScalar nie = semiconductor_region->material()->band->nie(p, n, T);
    AutoDScalar Nc  = semiconductor_region->material()->band->Nc(T);
    AutoDScalar Nv  = semiconductor_region->material()->band->Nv(T);
    AutoDScalar Eg  = semiconductor_region->material()->band->Eg(T);


    //governing equation of pis/n/p for Ohmic contact boundary
    AutoDScalar ff1,ff2,ff3;
    if(semiconductor_region->get_advanced_model()->Fermi) //Fermi
    {
      AutoDScalar Ec =  -(e*V + semiconductor_node_data->affinity() );
      AutoDScalar Ev =  -(e*V + semiconductor_node_data->affinity() + Eg);

      // the quasi-fermi potential equals to electrode Vapp
      AutoDScalar phin = V_resistance + resistance_node_data->affinity()/e;
      AutoDScalar phip = V_resistance + resistance_node_data->affinity()/e;

      AutoDScalar etan = (-e*phin-Ec)/kb/T;
      AutoDScalar etap = (Ev+e*phip)/kb/T;

      ff1 =  Nc*fermi_half(etan) - Nv*fermi_half(etap) -semiconductor_node_data->Net_doping();
      ff2 =  n - Nc*fermi_half(etan);
      ff3 =  p - Nv*fermi_half(etap);
    }
    else //Boltzmann
    {

      ff1 = V - kb*T/e*adtl::asinh(semiconductor_node_data->Net_doping()/(2*nie))
            + Eg/(2*e)
            + kb*T*log(Nc/Nv)/(2*e)
            + semiconductor_node_data->affinity()/e
            - (V_resistance + resistance_node_data->affinity()/e);

      AutoDScalar  electron_density;
      AutoDScalar  hole_density;
      PetscScalar  net_dpoing = semiconductor_node_data->Net_doping();
      if( net_dpoing <0 )   //p-type
      {
        hole_density = (-net_dpoing + sqrt(net_dpoing*net_dpoing + 4*nie*nie))/2.0;
        electron_density = nie*nie/hole_density;
      }
      else                               //n-type
      {
        electron_density = (net_dpoing + sqrt(net_dpoing*net_dpoing + 4*nie*nie))/2.0;
        hole_density = nie*nie/electron_density;
      }

      ff2 =  n - electron_density;  //governing equation for electron density
      ff3 =  p - hole_density;      //governing equation for hole density
    }

    // set Jacobian of governing equations of psi/n/p
    MatSetValues(A, 1, &row_re[node_psi_offset], col_re.size(), &col_re[0], ff1.getADValue(), ADD_VALUES);
    MatSetValues(A, 1, &row_im[node_psi_offset], col_im.size(), &col_im[0], ff1.getADValue(), ADD_VALUES);

    MatSetValues(A, 1, &row_re[node_n_offset],   col_re.size(), &col_re[0], ff2.getADValue(), ADD_VALUES);
    MatSetValues(A, 1, &row_im[node_n_offset],   col_im.size(), &col_im[0], ff2.getADValue(), ADD_VALUES);

    MatSetValues(A, 1, &row_re[node_p_offset],   col_re.size(), &col_re[0], ff3.getADValue(), ADD_VALUES);
    MatSetValues(A, 1, &row_im[node_p_offset],   col_im.size(), &col_im[0], ff3.getADValue(), ADD_VALUES);


    //governing equation of T
    if(semiconductor_region->get_advanced_model()->enable_Tl())
    {
      // fill A with J
      semiconductor_region->DDMAC_Fill_Nodal_Matrix_Vector(semiconductor_node, TEMPERATURE, A, b, J, omega, add_value_flag, resistance_region, resistance_node);
      semiconductor_region->DDMAC_Force_equal(semiconductor_node, TEMPERATURE, A, add_value_flag, resistance_region, resistance_node);
    }

    // electron temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tn())
    {
      AutoDScalar fTn = (n*(Tn - T));
      MatSetValues(A, 1, &row_re[node_Tn_offset], col_re.size(), &col_re[0], fTn.getADValue(),  ADD_VALUES);
      MatSetValues(A, 1, &row_im[node_Tn_offset], col_im.size(), &col_im[0], fTn.getADValue(),  ADD_VALUES);
    }

    // hole temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tp())
    {
      AutoDScalar fTp = (p*(Tp - T));
      MatSetValues(A, 1, &row_re[node_Tp_offset], col_re.size(), &col_re[0], fTp.getADValue(),  ADD_VALUES);
      MatSetValues(A, 1, &row_im[node_Tp_offset], col_im.size(), &col_im[0], fTp.getADValue(),  ADD_VALUES);
    }




    // for resistance region
    PetscInt resistance_node_psi_offset_re = resistance_node->global_offset() + resistance_region->ebm_variable_offset(POTENTIAL);
    PetscInt resistance_node_psi_offset_im = resistance_node_psi_offset_re + resistance_region->ebm_n_variables();

    resistance_region->DDMAC_Fill_Nodal_Matrix_Vector(resistance_node, POTENTIAL, A, b, J, omega, add_value_flag);
    if(resistance_region->get_advanced_model()->enable_Tl())
      resistance_region->DDMAC_Fill_Nodal_Matrix_Vector(resistance_node, TEMPERATURE, A, b, J, omega, add_value_flag);

    // add ohmic current to resistance region
    {
      std::vector< std::vector<PetscInt> >    cols_re, cols_im;
      std::vector< std::vector<PetscScalar> > current_jacobian_re, current_jacobian_im;

      std::vector<PetscScalar> A1(semiconductor_n_node_var), A2(semiconductor_n_node_var);
      std::vector<PetscScalar> BJ1(semiconductor_n_node_var), BJ2(semiconductor_n_node_var), BJ3(semiconductor_n_node_var), BJ4(semiconductor_n_node_var);
      // get derivative item of obmic electrode current from Jacobian matrx
      MatGetValues(J, 1, &row_re[node_n_offset], semiconductor_n_node_var, &row_re[0], &A1[0]);
      MatGetValues(J, 1, &row_re[node_p_offset], semiconductor_n_node_var, &row_re[0], &A2[0]);

      // build AC matrix items
      for(unsigned int nv=0; nv<semiconductor_n_node_var; ++nv)
      {
        BJ1[nv] =  -(A1[nv]-A2[nv]); // (real , real )
        BJ2[nv] =  0; // (real , image)
        BJ3[nv] =  0; // (image, real )
        BJ4[nv] =  -(A1[nv]-A2[nv]); // (image, image)
      }

      // the real equation location of ohmic bc in matrix entry
      cols_re.push_back(row_re);
      current_jacobian_re.push_back(BJ1);
      cols_re.push_back(row_im);
      current_jacobian_re.push_back(BJ2);

      // the image equation location of ohmic bc in matrix entry
      cols_im.push_back(row_re);
      current_jacobian_im.push_back(BJ3);
      cols_im.push_back(row_im);
      current_jacobian_im.push_back(BJ4);

      FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
      for(; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it)
      {
        const FVM_Node *nb_node = (*nb_it).first;
        std::vector<PetscInt>    nb_col_re(semiconductor_n_node_var), nb_col_im(semiconductor_n_node_var);
        for(unsigned int i=0; i<semiconductor_n_node_var; ++i)
        {
          nb_col_re[i] = nb_node->global_offset() + i;
          nb_col_im[i] = nb_node->global_offset() + semiconductor_n_node_var + i;
        }

        MatGetValues(J, 1, &row_re[node_n_offset], semiconductor_n_node_var, &nb_col_re[0], &A1[0]);
        MatGetValues(J, 1, &row_re[node_p_offset], semiconductor_n_node_var, &nb_col_re[0], &A2[0]);

        for(unsigned int nv=0; nv<semiconductor_n_node_var; ++nv)
        {
          BJ1[nv] =  -(A1[nv]-A2[nv]); // (real , real )
          BJ2[nv] =  0; // (real , image)
          BJ3[nv] =  0; // (image, real )
          BJ4[nv] =  -(A1[nv]-A2[nv]); // (image, image)
        }

        // the real equation location of ohmic bc in matrix entry
        cols_re.push_back(nb_col_re);
        current_jacobian_re.push_back(BJ1);
        cols_re.push_back(nb_col_im);
        current_jacobian_re.push_back(BJ2);

        // the image equation location of ohmic bc in matrix entry
        cols_im.push_back(nb_col_re);
        current_jacobian_im.push_back(BJ3);
        cols_im.push_back(nb_col_im);
        current_jacobian_im.push_back(BJ4);
      }

      for(unsigned int n=0; n<cols_re.size(); ++n)
        MatSetValues(A, 1, &resistance_node_psi_offset_re, cols_re[n].size(), &(cols_re[n])[0], &(current_jacobian_re[n])[0], ADD_VALUES);

      for(unsigned int n=0; n<cols_im.size(); ++n)
        MatSetValues(A, 1, &resistance_node_psi_offset_im, cols_im[n].size(), &(cols_im[n])[0], &(current_jacobian_im[n])[0], ADD_VALUES);
    }

    // compute displacement current
    {
      FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
      for(; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it)
      {
        const FVM_Node *nb_node = (*nb_it).first;

        // the psi of this node
        AutoDScalar  V = semiconductor_node_data->psi(); V.setADValue(0, 1.0);
        // the psi of neighbor node
        AutoDScalar V_nb = nb_node->node_data()->psi(); V_nb.setADValue(1, 1.0);

        // distance from nb node to this node
        PetscScalar distance = (*(semiconductor_node->root_node()) - *(nb_node->root_node())).size();

        // area of out surface of control volume related with neighbor node
        PetscScalar cv_boundary = semiconductor_node->cv_surface_area(nb_node);
        AutoDScalar D = semiconductor_node_data->eps()*(V-V_nb)/distance;

        // the 1/dt is replaced by j*omega.
        std::complex <PetscScalar> dJdisp_dV  = -cv_boundary*D.getADValue(0)*j*omega;
        std::complex <PetscScalar> dJdisp_dVn = -cv_boundary*D.getADValue(1)*j*omega;

        // V
        MatSetValue(A, resistance_node_psi_offset_re, semiconductor_node->global_offset(), dJdisp_dV.real(), ADD_VALUES);
        MatSetValue(A, resistance_node_psi_offset_re, semiconductor_node->global_offset()+semiconductor_region->ebm_n_variables(), -dJdisp_dV.imag(), ADD_VALUES);

        MatSetValue(A, resistance_node_psi_offset_im, semiconductor_node->global_offset(), dJdisp_dV.imag(), ADD_VALUES);
        MatSetValue(A, resistance_node_psi_offset_im, semiconductor_node->global_offset()+semiconductor_region->ebm_n_variables(),  dJdisp_dV.real(), ADD_VALUES);

        // V_nb
        MatSetValue(A, resistance_node_psi_offset_re, nb_node->global_offset(), dJdisp_dVn.real(), ADD_VALUES);
        MatSetValue(A, resistance_node_psi_offset_re, nb_node->global_offset()+semiconductor_region->ebm_n_variables(), -dJdisp_dVn.imag(), ADD_VALUES);

        MatSetValue(A, resistance_node_psi_offset_im, nb_node->global_offset(), dJdisp_dVn.imag(), ADD_VALUES);
        MatSetValue(A, resistance_node_psi_offset_im, nb_node->global_offset()+semiconductor_region->ebm_n_variables(),  dJdisp_dVn.real(), ADD_VALUES);
      }
    }

    // if we have extra regions associated with this node?

    std::vector<SimulationRegion *> node_extra_regions = this->extra_regions(*node_it);
    for(unsigned int n=0; n<node_extra_regions.size(); n++)
    {
      const SimulationRegion * region = node_extra_regions[n];
      const FVM_Node * fvm_node  = get_region_fvm_node(*node_it, region);

      // the psi of insulator region equals to metal region
      region->DDMAC_Force_equal(fvm_node, POTENTIAL, A, add_value_flag, resistance_region, resistance_node);

      if(region->get_advanced_model()->enable_Tl())
      {
        // load Jacobian entry of this node from J and fill into AC matrix A
        region->DDMAC_Fill_Nodal_Matrix_Vector(fvm_node, TEMPERATURE, A, b, J, omega, add_value_flag, resistance_region, resistance_node);
        region->DDMAC_Force_equal(fvm_node, TEMPERATURE, A, add_value_flag, resistance_region, resistance_node);
      }
    }

  }

}

void IF_Metal_OhmicBC::DDMAC_Update_Solution(const PetscScalar * lxx, const Mat, const double omega)
{}


