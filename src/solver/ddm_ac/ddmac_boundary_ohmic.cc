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

#include "asinh.hpp" // for asinh

// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_ohmic.h"
#include "parallel.h"
#include "mathfunc.h"
#include "petsc_utils.h"
#include "TNT/tnt_cmat.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;





/*---------------------------------------------------------------------
 * build AC matrix for DDMAC solver
 */
void OhmicContactBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const PetscScalar omega, InsertMode & add_value_flag )
{
  // the DDM AC matrix of Ohmic boundary condition is processed here

  // for AC scan, the electrode should always be voltage driven
  assert (this->ext_circuit()->is_voltage_driven());

  PetscInt bc_global_offset_re = this->global_offset();
  PetscInt bc_global_offset_im = this->global_offset()+1;

  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");

  // impedance at frequency omega
  std::complex <PetscScalar> j(0.0, 1.0);

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();


  // we set Jacobian entries for ohmic boundary condition
  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

      switch ( regions[i]->type() )
      {
          // Semiconductor Region of course owns OhmicContactBC
          case SemiconductorRegion:
          {
            // semiconductor region should be the first region

            const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(regions[i]);

            // find the node variable offset
            unsigned int n_node_var      = semi_region->ebm_n_variables();
            unsigned int node_psi_offset = semi_region->ebm_variable_offset(POTENTIAL);
            unsigned int node_n_offset   = semi_region->ebm_variable_offset(ELECTRON);
            unsigned int node_p_offset   = semi_region->ebm_variable_offset(HOLE);
            unsigned int node_Tl_offset  = semi_region->ebm_variable_offset(TEMPERATURE);
            unsigned int node_Tn_offset  = semi_region->ebm_variable_offset(E_TEMP);
            unsigned int node_Tp_offset  = semi_region->ebm_variable_offset(H_TEMP);

            //the indepedent variable number, we only need n_node_var plus electrode potential here.
            adtl::AutoDScalar::numdir = n_node_var + 1;
            //synchronize with material database
            semi_region->material()->set_ad_num(adtl::AutoDScalar::numdir);
            semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

            AutoDScalar V = node_data->psi();  V.setADValue(node_psi_offset, 1.0);  // psi of this node
            AutoDScalar n = node_data->n();    n.setADValue(node_n_offset, 1.0);  // electron density
            AutoDScalar p = node_data->p();    p.setADValue(node_p_offset, 1.0);  // hole density

            AutoDScalar T  =  T_external();
            AutoDScalar Tn =  T_external();
            AutoDScalar Tp =  T_external();

            // lattice temperature if required
            if(semi_region->get_advanced_model()->enable_Tl())
            {
              T =  node_data->T();
              T.setADValue(node_Tl_offset, 1.0);
            }

            // electron temperature if required
            if(semi_region->get_advanced_model()->enable_Tn())
            {
              AutoDScalar nTn = node_data->n()*node_data->Tn();
              nTn.setADValue(node_Tn_offset, 1.0);
              Tn = nTn/n;
            }

            // hole temperature if required
            if(semi_region->get_advanced_model()->enable_Tp())
            {
              AutoDScalar pTp = node_data->p()*node_data->Tp();
              pTp.setADValue(node_Tp_offset, 1.0);
              Tp = pTp/p;
            }

            // the electrode potential in current iteration
            genius_assert( local_offset()!=invalid_uint );
            AutoDScalar Ve = this->ext_circuit()->Vac();             Ve.setADValue(n_node_var, 1.0);

            // the insert position
            std::vector<PetscInt> row_re, col_re;
            std::vector<PetscInt> row_im, col_im;
            for(unsigned int nv=0; nv<n_node_var; ++nv)
            {
              row_re.push_back(fvm_nodes[i]->global_offset() + nv);
              row_im.push_back(fvm_nodes[i]->global_offset() + n_node_var + nv);
            }

            col_re = row_re;
            col_im = row_im;

            col_re.push_back(bc_global_offset_re);
            col_im.push_back(bc_global_offset_im);

            AutoDScalar nie = semi_region->material()->band->nie(p, n, T);
            AutoDScalar Nc  = semi_region->material()->band->Nc(T);
            AutoDScalar Nv  = semi_region->material()->band->Nv(T);
            AutoDScalar Eg  = semi_region->material()->band->Eg(T);


            //governing equation of pis/n/p for Ohmic contact boundary
            AutoDScalar ff1,ff2,ff3;
            if(semi_region->get_advanced_model()->Fermi) //Fermi
            {
              AutoDScalar Ec =  -(e*V + node_data->affinity() );
              AutoDScalar Ev =  -(e*V + node_data->affinity() + Eg);

              // the quasi-fermi potential equals to electrode Vapp
              AutoDScalar phin = Ve;
              AutoDScalar phip = Ve;

              AutoDScalar etan = (-e*phin-Ec)/kb/T;
              AutoDScalar etap = (Ev+e*phip)/kb/T;

              ff1 =  Nc*fermi_half(etan) - Nv*fermi_half(etap) -node_data->Net_doping();
              ff2 =  n - Nc*fermi_half(etan);
              ff3 =  p - Nv*fermi_half(etap);
            }
            else //Boltzmann
            {

              ff1 = V - kb*T/e*adtl::asinh(node_data->Net_doping()/(2*nie))
                    + Eg/(2*e)
                    + kb*T*log(Nc/Nv)/(2*e)
                    + node_data->affinity()/e
                    - Ve;

              AutoDScalar  electron_density;
              AutoDScalar  hole_density;
              PetscScalar  net_dpoing = node_data->Net_doping();
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
            if(regions[i]->get_advanced_model()->enable_Tl() && (node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion)) )
            {
              // fill A with J
              regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[0], TEMPERATURE, A, b, J, omega, add_value_flag);

              // if this ohmic bc is external boundary, set heat flux here
              PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
              AutoDScalar fT = Heat_Transfer*(T_external()-T)*S;
              MatSetValues(A, 1, &row_re[node_Tl_offset], col_re.size(), &col_re[0], fT.getADValue(),  ADD_VALUES);
              MatSetValues(A, 1, &row_im[node_Tl_offset], col_im.size(), &col_im[0], fT.getADValue(),  ADD_VALUES);
            }

            // electron temperature if required
            if(regions[i]->get_advanced_model()->enable_Tn())
            {
              AutoDScalar fTn = (n*(Tn - T));
              MatSetValues(A, 1, &row_re[node_Tn_offset], col_re.size(), &col_re[0], fTn.getADValue(),  ADD_VALUES);
              MatSetValues(A, 1, &row_im[node_Tn_offset], col_im.size(), &col_im[0], fTn.getADValue(),  ADD_VALUES);
            }

            // hole temperature if required
            if(regions[i]->get_advanced_model()->enable_Tp())
            {
              AutoDScalar fTp = (p*(Tp - T));
              MatSetValues(A, 1, &row_re[node_Tp_offset], col_re.size(), &col_re[0], fTp.getADValue(),  ADD_VALUES);
              MatSetValues(A, 1, &row_im[node_Tp_offset], col_im.size(), &col_im[0], fTp.getADValue(),  ADD_VALUES);
            }

            break;
          }
          // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
          // the boundary for the corner point may be Ohmic.
          case InsulatorRegion:
          // conductor region which has an interface with OhmicContact boundary to semiconductor region
          case MetalRegion:
          case ElectrodeRegion:
          {

            //let psi and T of node in this region equal to psi(T) of node in the insulator region
            regions[i]->DDMAC_Force_equal(fvm_nodes[i], POTENTIAL, A, add_value_flag, regions[0], fvm_nodes[0]);

            // the T of this node is equal to corresponding T of insulator node
            if(regions[i]->get_advanced_model()->enable_Tl())
            {
              regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], TEMPERATURE, A, b, J, omega, add_value_flag, regions[0], fvm_nodes[0]);
              regions[i]->DDMAC_Force_equal(fvm_nodes[i], TEMPERATURE, A, add_value_flag, regions[0], fvm_nodes[0]);
            }


            break;
          }
          case VacuumRegion:
          break;
          default: genius_error(); //we should never reach here
      }
    }

  }


  // the extra equation of ohmic boundary
  // For ac scan
  //
  //          _____  (Z1)          Ve
  //    -----|_____|----/\/\/\/\-------> to ohmic electrode (Ve, I)
  //    |       R          L       |
  //   Vac                      C === (Y2)
  //    |__________________________|
  //           GND
  //

  // we should set matrix entries related with obmic electrode current
  // which is required by external electrode circuit
  {
    const std::complex<PetscScalar> mna_scaling = ext_circuit()->mna_ac_scaling(omega);

    std::vector< std::vector<PetscInt> >    cols_re, cols_im;
    std::vector< std::vector<PetscScalar> > current_jacobian_re, current_jacobian_im;


    BoundaryCondition::const_node_iterator node_it;
    BoundaryCondition::const_node_iterator end_it = nodes_end();
    for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // only process nodes belong to this processor
      if( (*node_it)->processor_id() != Genius::processor_id() ) continue;

      // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
      // but belong to different regions in logic.
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
      for(; rnode_it!=end_rnode_it; ++rnode_it  )
      {
        const SimulationRegion * region = (*rnode_it).second.first;
        const FVM_Node *  fvm_node = (*rnode_it).second.second;

        if ( region->type() == SemiconductorRegion)
        {
          const FVM_NodeData * node_data = fvm_node->node_data();

          unsigned int n_node_var    = region->ebm_n_variables();
          unsigned int node_n_offset = region->ebm_variable_offset(ELECTRON);
          unsigned int node_p_offset = region->ebm_variable_offset(HOLE);

          std::vector<PetscScalar> A1(n_node_var), A2(n_node_var);
          std::vector<PetscScalar> BJ1(n_node_var), BJ2(n_node_var), BJ3(n_node_var), BJ4(n_node_var);
          std::vector<PetscInt>    row_re(n_node_var), row_im(n_node_var);

          for(unsigned int nv=0; nv<n_node_var; ++nv)
          {
            row_re[nv] = fvm_node->global_offset() + nv;
            row_im[nv] = fvm_node->global_offset() + n_node_var + nv;
          }

          // get derivative item of obmic electrode current from Jacobian matrx
          MatGetValues(J, 1, &row_re[node_n_offset], n_node_var, &row_re[0], &A1[0]);
          MatGetValues(J, 1, &row_re[node_p_offset], n_node_var, &row_re[0], &A2[0]);

          // build AC matrix items
          for(unsigned int nv=0; nv<n_node_var; ++nv)
          {
            BJ1[nv] =  (mna_scaling*(A1[nv]-A2[nv])*current_scale).real(); // (real , real )
            BJ2[nv] = -(mna_scaling*(A1[nv]-A2[nv])*current_scale).imag(); // (real , image)
            BJ3[nv] =  (mna_scaling*(A1[nv]-A2[nv])*current_scale).imag(); // (image, real )
            BJ4[nv] =  (mna_scaling*(A1[nv]-A2[nv])*current_scale).real(); // (image, image)
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

          // get the derivative of electrode current to neighbors of ohmic node

          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
          FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
          for(; nb_it != nb_it_end; ++nb_it)
          {
            const FVM_Node *  fvm_nb_node = (*nb_it).first;
            std::vector<PetscInt>    col_re(n_node_var), col_im(n_node_var);
            for(unsigned int i=0; i<n_node_var; ++i)
            {
              col_re[i] = fvm_nb_node->global_offset() + i;
              col_im[i] = fvm_nb_node->global_offset() + n_node_var + i;
            }

            MatGetValues(J, 1, &row_re[node_n_offset], n_node_var, &col_re[0], &A1[0]);
            MatGetValues(J, 1, &row_re[node_p_offset], n_node_var, &col_re[0], &A2[0]);

            for(unsigned int nv=0; nv<n_node_var; ++nv)
            {
              BJ1[nv] =  (mna_scaling*(A1[nv]-A2[nv])*current_scale).real(); // (real , real )
              BJ2[nv] = -(mna_scaling*(A1[nv]-A2[nv])*current_scale).imag(); // (real , image)
              BJ3[nv] =  (mna_scaling*(A1[nv]-A2[nv])*current_scale).imag(); // (image, real )
              BJ4[nv] =  (mna_scaling*(A1[nv]-A2[nv])*current_scale).real(); // (image, image)
            }

            // the real equation location of ohmic bc in matrix entry
            cols_re.push_back(col_re);
            current_jacobian_re.push_back(BJ1);
            cols_re.push_back(col_im);
            current_jacobian_re.push_back(BJ2);

            // the image equation location of ohmic bc in matrix entry
            cols_im.push_back(col_re);
            current_jacobian_im.push_back(BJ3);
            cols_im.push_back(col_im);
            current_jacobian_im.push_back(BJ4);

            // process displacement current. we can't neglect displacement current when frequency is high.
            adtl::AutoDScalar::numdir =2;

            // the psi of this node
            AutoDScalar  V = node_data->psi(); V.setADValue(0, 1.0);
            // the psi of neighbor node
            AutoDScalar V_nb = fvm_nb_node->node_data()->psi(); V_nb.setADValue(1, 1.0);

            // distance from nb node to this node
            PetscScalar distance = fvm_node->distance(fvm_nb_node);

            // area of out surface of control volume related with neighbor node
            PetscScalar cv_boundary = fvm_node->cv_surface_area(fvm_nb_node);
            AutoDScalar D = node_data->eps()*(V-V_nb)/distance;

            // the 1/dt is replaced by j*omega.
            std::complex <PetscScalar> dJdisp_dV  = mna_scaling*cv_boundary*D.getADValue(0)*j*omega*current_scale;
            std::complex <PetscScalar> dJdisp_dVn = mna_scaling*cv_boundary*D.getADValue(1)*j*omega*current_scale;

            // V
            MatSetValue(A, bc_global_offset_re, fvm_node->global_offset(), dJdisp_dV.real(), ADD_VALUES);
            MatSetValue(A, bc_global_offset_re, fvm_node->global_offset()+region->ebm_n_variables(), -dJdisp_dV.imag(), ADD_VALUES);

            MatSetValue(A, bc_global_offset_im, fvm_node->global_offset(), dJdisp_dV.imag(), ADD_VALUES);
            MatSetValue(A, bc_global_offset_im, fvm_node->global_offset()+region->ebm_n_variables(),  dJdisp_dV.real(), ADD_VALUES);

            // V_nb
            MatSetValue(A, bc_global_offset_re, fvm_nb_node->global_offset(), dJdisp_dVn.real(), ADD_VALUES);
            MatSetValue(A, bc_global_offset_re, fvm_nb_node->global_offset()+region->ebm_n_variables(), -dJdisp_dVn.imag(), ADD_VALUES);

            MatSetValue(A, bc_global_offset_im, fvm_nb_node->global_offset(), dJdisp_dVn.imag(), ADD_VALUES);
            MatSetValue(A, bc_global_offset_im, fvm_nb_node->global_offset()+region->ebm_n_variables(),  dJdisp_dVn.real(), ADD_VALUES);
          }
        }
      }
    }
    // d(current)/d(independent variables of bd node and its neighbors)
    for(unsigned int n=0; n<cols_re.size(); ++n)
      MatSetValues(A, 1, &bc_global_offset_re, cols_re[n].size(), &(cols_re[n])[0], &(current_jacobian_re[n])[0], ADD_VALUES);

    for(unsigned int n=0; n<cols_im.size(); ++n)
      MatSetValues(A, 1, &bc_global_offset_im, cols_im[n].size(), &(cols_im[n])[0], &(current_jacobian_im[n])[0], ADD_VALUES);
  }


  if(Genius::is_last_processor())
  {
    // here we process the external circuit, we do not use AD here

    std::complex <PetscScalar> K = this->ext_circuit()->mna_ac_jacobian(omega);

    MatSetValue(A, bc_global_offset_re, bc_global_offset_re,  K.real(), ADD_VALUES);
    MatSetValue(A, bc_global_offset_re, bc_global_offset_im, -K.imag(), ADD_VALUES);
    MatSetValue(A, bc_global_offset_im, bc_global_offset_re,  K.imag(), ADD_VALUES);
    MatSetValue(A, bc_global_offset_im, bc_global_offset_im,  K.real(), ADD_VALUES);

    VecSetValue(b, bc_global_offset_re, this->ext_circuit()->Vac(), ADD_VALUES);
    VecSetValue(b, bc_global_offset_im, 0.0, ADD_VALUES);
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


void OhmicContactBC::DDMAC_Update_Solution(const PetscScalar * lxx, const Mat J, const PetscScalar omega)
{

  std::complex<PetscScalar> Iac(0.0, 0.0);
  std::complex<PetscScalar> j(0.0, 1.0);


  // for 2D mesh, system().z_width() is the device dimension in Z direction; for 3D mesh, system().z_width() is 1.0
  std::complex<PetscScalar> current_scale = system().z_width();

  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(; rnode_it!=end_rnode_it; ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node *  fvm_node = (*rnode_it).second.second;

      if ( region->type() == SemiconductorRegion)
      {
        const FVM_NodeData * node_data = fvm_node->node_data();

        unsigned int n_node_var    = region->ebm_n_variables();
        unsigned int node_n_offset = region->ebm_variable_offset(ELECTRON);
        unsigned int node_p_offset = region->ebm_variable_offset(HOLE);

        std::vector<PetscInt> row;
        for(unsigned int i=0; i<n_node_var; ++i)
          row.push_back(fvm_node->global_offset() + i);

        std::vector<PetscScalar> Vec1(row.size()), Vec2(row.size()), M(row.size()*row.size());

        PetscScalar V_re = lxx[fvm_node->local_offset() + region->ebm_variable_offset(POTENTIAL)];
        PetscScalar V_im = lxx[fvm_node->local_offset() + region->ebm_n_variables() + region->ebm_variable_offset(POTENTIAL)];
        std::complex<PetscScalar> V = std::complex<PetscScalar>(V_re, V_im);

        // the neighbor node to ohmic bd node in semiconductor region
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
        for(; nb_it != nb_it_end; ++nb_it)
        {
          const FVM_Node *  fvm_nb_node = (*nb_it).first;
          std::vector<PetscInt> col;

          for(unsigned int i=0; i<n_node_var; ++i)
            col.push_back(fvm_nb_node->global_offset() + i);

          MatGetValues(J, row.size(), &row[0], col.size(), &col[0], &M[0]);

          for(unsigned int i=0; i<n_node_var; ++i)
          {
            Vec1[i] = lxx[fvm_nb_node->local_offset() + i];
            Vec2[i] = lxx[fvm_nb_node->local_offset() + n_node_var + i];
          }

          TNT::Matrix<PetscScalar> A(row.size(), row.size(), &M[0]);
          TNT::Vector<PetscScalar> v1(col.size(), &Vec1[0]);
          TNT::Vector<PetscScalar> v2(col.size(), &Vec2[0]);

          TNT::Vector<PetscScalar> b1 = A*v1;
          TNT::Vector<PetscScalar> b2 = A*v2;

          // conductance current
          Iac += std::complex<PetscScalar>(b1[node_n_offset]-b1[node_p_offset], b2[node_n_offset]-b2[node_p_offset])*current_scale;

          // displacement current

          // distance from nb node to this node
          PetscScalar distance = (*(fvm_node->root_node()) - *(fvm_nb_node->root_node())).size();

          PetscScalar Vn_re = lxx[fvm_nb_node->local_offset() + region->ebm_variable_offset(POTENTIAL)];
          PetscScalar Vn_im = lxx[fvm_nb_node->local_offset() + region->ebm_n_variables() + region->ebm_variable_offset(POTENTIAL)];
          std::complex<PetscScalar> Vn = std::complex<PetscScalar>(Vn_re, Vn_im);

          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = fvm_node->cv_surface_area(fvm_nb_node);

          Iac += node_data->eps()*(V-Vn)/distance*j*omega*cv_boundary*current_scale;
        }

        MatGetValues(J, row.size(), &row[0], row.size(), &row[0], &M[0]);
        for(unsigned int i=0; i<n_node_var; ++i)
        {
          Vec1[i] = lxx[fvm_node->local_offset() + i];
          Vec2[i] = lxx[fvm_node->local_offset() + n_node_var + i];
        }

        TNT::Matrix<PetscScalar> A(row.size(), row.size(), &M[0]);
        TNT::Vector<PetscScalar> v1(row.size(), &Vec1[0]);
        TNT::Vector<PetscScalar> v2(row.size(), &Vec2[0]);

        TNT::Vector<PetscScalar> b1 = A*v1;
        TNT::Vector<PetscScalar> b2 = A*v2;
        // conductance current
        Iac += std::complex<PetscScalar>(b1[node_n_offset]-b1[node_p_offset], b2[node_n_offset]-b2[node_p_offset])*current_scale;
      }
    }
  }

  Parallel::sum(Iac);

  this->ext_circuit()->current_ac() = Iac;
  this->ext_circuit()->potential_ac() = std::complex<PetscScalar>(lxx[this->local_offset()], lxx[this->local_offset()+1]);
}

