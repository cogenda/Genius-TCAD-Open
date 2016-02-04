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

// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_schottky.h"
#include "parallel.h"
#include "petsc_utils.h"

using PhysicalUnit::e;

void SchottkyContactBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const PetscScalar omega, InsertMode & add_value_flag )
{


  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);
  }

  PetscInt bc_global_offset_re = this->global_offset();
  PetscInt bc_global_offset_im = this->global_offset()+1;

  const PetscScalar Work_Function = this->scalar("workfunction");
  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");

  // impedance at frequency omega
  std::complex <PetscScalar> j(0.0, 1.0);

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();


  // loop again
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
          // Semiconductor Region of course owns Schottky Contact BC
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

            const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

            // mapping this node to material library
            semi_region->material()->mapping(fvm_nodes[i]->root_node(), node_data, SolverSpecify::clock);

            //the indepedent variable number, we only need n_node_var+1 here.
            adtl::AutoDScalar::numdir=n_node_var+1;
            //synchronize with material database
            semi_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

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

            //Schotty Barrier Lowerring
            PetscScalar deltaVB = semi_region->material()->band->SchottyBarrierLowerring(node_data->eps(), node_data->E().size());

            //Schottky current
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            AutoDScalar Fn = semi_region->material()->band->SchottyJsn(n, T, Work_Function - node_data->affinity() - deltaVB) * S;
            AutoDScalar Fp = semi_region->material()->band->SchottyJsp(p, T, Work_Function - node_data->affinity() + deltaVB) * S;

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


            // schottky boundary condition of poisson's equation
            AutoDScalar f_phi = V + Work_Function - deltaVB - Ve;
            MatSetValues(A, 1, &row_re[node_psi_offset], col_re.size(), &col_re[0], f_phi.getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &row_im[node_psi_offset], col_im.size(), &col_im[0], f_phi.getADValue(), ADD_VALUES);

            // process the Jacobian of Schottky current
            regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], ELECTRON, A, b, J, omega, add_value_flag);
            MatSetValues(A, 1, &row_re[node_n_offset],   col_re.size(), &col_re[0], Fn.getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &row_im[node_n_offset],   col_im.size(), &col_im[0], Fn.getADValue(), ADD_VALUES);

            regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], HOLE, A, b, J, omega, add_value_flag);
            MatSetValues(A, 1, &row_re[node_p_offset],   col_re.size(), &col_re[0], (-Fp).getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &row_im[node_p_offset],   col_im.size(), &col_im[0], (-Fp).getADValue(), ADD_VALUES);

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


            // process the Jacobian of current flow out of schottky electrode

            // compute the schottky thermal emit current
            {
              std::complex <PetscScalar> mna_scaling = ext_circuit()->mna_ac_scaling(omega);
              AutoDScalar current_emit = -(Fn+Fp)*current_scale;

              std::vector< PetscScalar >  A1(n_node_var+1), A2(n_node_var+1), A3(n_node_var+1), A4(n_node_var+1);
              for(unsigned int nv=0; nv<n_node_var+1; ++nv)
              {
                std::complex<PetscScalar> C = mna_scaling*current_emit.getADValue(nv);
                A1[nv] =  C.real();
                A2[nv] = -C.imag();
                A3[nv] =  C.imag();
                A4[nv] =  C.real();
              }

              MatSetValues(A, 1, &bc_global_offset_re, col_re.size(), &col_re[0], &A1[0], ADD_VALUES);
              MatSetValues(A, 1, &bc_global_offset_re, col_im.size(), &col_im[0], &A2[0], ADD_VALUES);

              MatSetValues(A, 1, &bc_global_offset_im, col_re.size(), &col_re[0], &A3[0], ADD_VALUES);
              MatSetValues(A, 1, &bc_global_offset_im, col_im.size(), &col_im[0], &A4[0], ADD_VALUES);
            }

            // compute displacement current

            FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
            for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
            {
              const FVM_Node *nb_node = (*nb_it).first;

              // the psi of this node
              AutoDScalar  V = node_data->psi(); V.setADValue(0, 1.0);
              // the psi of neighbor node
              AutoDScalar V_nb = nb_node->node_data()->psi(); V_nb.setADValue(1, 1.0);

              // distance from nb node to this node
              PetscScalar distance = fvm_nodes[i]->distance(nb_node);

              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node);
              AutoDScalar dE = (V-V_nb)/distance;

              // the 1/dt is replaced by j*omega.
              std::complex <PetscScalar> mna_scaling = ext_circuit()->mna_ac_scaling(omega);
              std::complex <PetscScalar> dJdisp_dV  = mna_scaling*cv_boundary*node_data->eps()*dE.getADValue(0)*j*omega*current_scale;
              std::complex <PetscScalar> dJdisp_dVn = mna_scaling*cv_boundary*node_data->eps()*dE.getADValue(1)*j*omega*current_scale;

              // V
              MatSetValue(A, bc_global_offset_re, fvm_nodes[i]->global_offset(), dJdisp_dV.real(), ADD_VALUES);
              MatSetValue(A, bc_global_offset_re, fvm_nodes[i]->global_offset()+regions[i]->ebm_n_variables(), -dJdisp_dV.imag(), ADD_VALUES);

              MatSetValue(A, bc_global_offset_im, fvm_nodes[i]->global_offset(), dJdisp_dV.imag(), ADD_VALUES);
              MatSetValue(A, bc_global_offset_im, fvm_nodes[i]->global_offset()+regions[i]->ebm_n_variables(),  dJdisp_dV.real(), ADD_VALUES);

              // V_nb
              MatSetValue(A, bc_global_offset_re, nb_node->global_offset(), dJdisp_dVn.real(), ADD_VALUES);
              MatSetValue(A, bc_global_offset_re, nb_node->global_offset()+regions[i]->ebm_n_variables(), -dJdisp_dVn.imag(), ADD_VALUES);

              MatSetValue(A, bc_global_offset_im, nb_node->global_offset(), dJdisp_dVn.imag(), ADD_VALUES);
              MatSetValue(A, bc_global_offset_im, nb_node->global_offset()+regions[i]->ebm_n_variables(),  dJdisp_dVn.real(), ADD_VALUES);
            }

            break;
          }
          // conductor region which has an interface with Schottky Contact boundary to semiconductor region
          case ElectrodeRegion:
          // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
          // the boundary for the corner point may be Schottky .
          case InsulatorRegion:
          {

            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            unsigned int n_variables     = regions[i]->ebm_n_variables();
            unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
            unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

            unsigned int semiregion_n_variables     = regions[0]->ebm_n_variables();
            unsigned int semiregion_node_psi_offset = regions[0]->ebm_variable_offset(POTENTIAL);
            unsigned int semiregion_node_Tl_offset  = regions[0]->ebm_variable_offset(TEMPERATURE);

            {
              // psi of this node
              AutoDScalar  V = fvm_nodes[i]->node_data()->psi();  V.setADValue(0,1.0);
              // psi for corresponding semiconductor region
              AutoDScalar  V_semi = fvm_nodes[0]->node_data()->psi(); V_semi.setADValue(1,1.0);
              // the psi of this node is equal to corresponding psi of semiconductor node
              AutoDScalar  ff1 = V - V_semi;

              PetscInt row_re = fvm_nodes[i]->global_offset() + node_psi_offset;
              PetscInt row_im = fvm_nodes[i]->global_offset() + n_variables + node_psi_offset;

              PetscInt col_re[2] = {fvm_nodes[i]->global_offset() + node_psi_offset, fvm_nodes[0]->global_offset() + semiregion_node_psi_offset};
              PetscInt col_im[2] = {fvm_nodes[i]->global_offset() + n_variables + node_psi_offset,
                                    fvm_nodes[0]->global_offset() + semiregion_n_variables + semiregion_node_psi_offset};

              MatSetValues(A, 1, &row_re, 2, &col_re[0], ff1.getADValue(), ADD_VALUES);
              MatSetValues(A, 1, &row_im, 2, &col_im[0], ff1.getADValue(), ADD_VALUES);

            }

            if(regions[i]->get_advanced_model()->enable_Tl())
            {
              regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], TEMPERATURE, A, b, J, omega, add_value_flag, regions[0], fvm_nodes[0]);

              AutoDScalar  T = fvm_nodes[i]->node_data()->T();      T.setADValue(0,1.0); // lattice temperature of this node
              AutoDScalar  T_semi = fvm_nodes[0]->node_data()->T(); T_semi.setADValue(1,1.0);
              // the T of this node is equal to corresponding T of semiconductor node
              AutoDScalar  ff2 = T - T_semi;

              PetscInt row_re = fvm_nodes[i]->global_offset() + node_Tl_offset;
              PetscInt row_im = fvm_nodes[i]->global_offset() + n_variables + node_Tl_offset;

              PetscInt col_re[2] = {fvm_nodes[i]->global_offset() + node_Tl_offset, fvm_nodes[0]->global_offset() + semiregion_node_Tl_offset};
              PetscInt col_im[2] = {fvm_nodes[i]->global_offset() + n_variables + node_Tl_offset,
                                    fvm_nodes[0]->global_offset() + semiregion_n_variables + semiregion_node_Tl_offset};

              MatSetValues(A, 1, &row_re, 2, &col_re[0], ff2.getADValue(), ADD_VALUES);
              MatSetValues(A, 1, &row_im, 2, &col_im[0], ff2.getADValue(), ADD_VALUES);
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



/*---------------------------------------------------------------------
 * update electrode potential
 */
void SchottkyContactBC::DDMAC_Update_Solution(const PetscScalar * lxx, const Mat, const PetscScalar omega)
{
  PetscScalar IacRe = 0;
  PetscScalar IacIm = 0;

  // for 2D mesh, system().z_width() is the device dimension in Z direction; for 3D mesh, system().z_width() is 1.0
  PetscScalar current_scale = system().z_width();

  const PetscScalar Work_Function = this->scalar("workfunction");
  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");

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
        const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
        const FVM_NodeData * node_data = fvm_node->node_data();

        // mapping this node to material library
        semi_region->material()->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

        //the indepedent variable number, we only need 2/3 here.
        adtl::AutoDScalar::numdir=2;

        //synchronize with material database
        semi_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

        AutoDScalar n = node_data->n();    n.setADValue(0, 1.0);  // electron density
        AutoDScalar p = node_data->p();    p.setADValue(1, 1.0);  // hole density
        AutoDScalar T = node_data->T();                           // lattice temperature of this node

        if(region->get_advanced_model()->enable_Tl())
        {
          adtl::AutoDScalar::numdir=3;
          T.setADValue(2, 1.0);
        }

        //Schotty Barrier Lowerring
        PetscScalar deltaVB = semi_region->material()->band->SchottyBarrierLowerring(node_data->eps(), node_data->E().size());

        //Schottky current
        PetscScalar S  = fvm_node->outside_boundary_surface_area();
        AutoDScalar Fn = semi_region->material()->band->SchottyJsn(n, T, Work_Function - node_data->affinity() - deltaVB) * S;
        AutoDScalar Fp = semi_region->material()->band->SchottyJsp(p, T, Work_Function - node_data->affinity() + deltaVB) * S;

        std::vector <PetscScalar> x_real, x_imag;
        {
          x_real.push_back(lxx[fvm_node->local_offset() + region->ebm_variable_offset(ELECTRON)]);
          x_real.push_back(lxx[fvm_node->local_offset() + region->ebm_variable_offset(HOLE)]);

          x_imag.push_back(lxx[fvm_node->local_offset() + region->ebm_n_variables() + region->ebm_variable_offset(ELECTRON)]);
          x_imag.push_back(lxx[fvm_node->local_offset() + region->ebm_n_variables() + region->ebm_variable_offset(HOLE)]);

          if(region->get_advanced_model()->enable_Tl())
          {
            x_real.push_back(lxx[fvm_node->local_offset() + region->ebm_variable_offset(TEMPERATURE)]);
            x_imag.push_back(lxx[fvm_node->local_offset() + region->ebm_n_variables() + region->ebm_variable_offset(TEMPERATURE)]);
          }
        }

        for(unsigned n=0; n<adtl::AutoDScalar::numdir; ++n)
        {
          IacRe += -(Fn.getADValue(n)+Fp.getADValue(n))*x_real[n]*current_scale;
          IacIm += -(Fn.getADValue(n)+Fp.getADValue(n))*x_imag[n]*current_scale;
        }

        // displacement current

        PetscScalar V_re = lxx[fvm_node->local_offset() + region->ebm_variable_offset(POTENTIAL)];
        PetscScalar V_im = lxx[fvm_node->local_offset() + region->ebm_n_variables() + region->ebm_variable_offset(POTENTIAL)];

        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
        {
          const FVM_Node *nb_node = (*nb_it).first;

          // distance from nb node to this node
          PetscScalar distance = (*(fvm_node->root_node()) - *(nb_node->root_node())).size();

          PetscScalar Vn_re = lxx[nb_node->local_offset() + region->ebm_variable_offset(POTENTIAL)];
          PetscScalar Vn_im = lxx[nb_node->local_offset() + region->ebm_n_variables() + region->ebm_variable_offset(POTENTIAL)];

          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = fvm_node->cv_surface_area(nb_node);

          IacRe += -omega*node_data->eps()*(V_im-Vn_im)/distance*cv_boundary*current_scale;
          IacIm +=  omega*node_data->eps()*(V_re-Vn_re)/distance*cv_boundary*current_scale;
        }

      }
    }
  }

  Parallel::sum(IacRe);
  Parallel::sum(IacIm);

  this->ext_circuit()->current_ac() = std::complex<PetscScalar>(IacRe, IacIm);


  PetscScalar VacRe = lxx[this->local_offset()];
  PetscScalar VacIm = lxx[this->local_offset()+1];

  this->ext_circuit()->potential_ac() = std::complex<PetscScalar>(VacRe, VacIm);

}
