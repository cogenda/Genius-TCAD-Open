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
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_gate.h"
#include "parallel.h"
#include "petsc_utils.h"



void GateContactBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const PetscScalar omega, InsertMode & add_value_flag )
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

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
      case InsulatorRegion:
        {
          unsigned int n_variables     = regions[i]->ebm_n_variables();
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          {
            // psi of this node
            AutoDScalar  V = node_data->psi();               V.setADValue(0, 1.0);

            // the electrode potential in current iteration
            genius_assert( local_offset()!=invalid_uint );
            AutoDScalar Ve = this->ext_circuit()->Vac();     Ve.setADValue(1, 1.0);

            // the governing equation of potential
            AutoDScalar ff = V + Work_Function - Ve;

            // the insert position
            PetscInt row_re  = fvm_nodes[i]->global_offset() + node_psi_offset;
            PetscInt cols_re[2] = {row_re, bc_global_offset_re};
            PetscInt row_im  = fvm_nodes[i]->global_offset() + n_variables + node_psi_offset;
            PetscInt cols_im[2] = {row_im, bc_global_offset_im};
            // process the Jacobian of governing equation of potential
            MatSetValues(A, 1, &row_re, 2, cols_re, ff.getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &row_im, 2, cols_im, ff.getADValue(), ADD_VALUES);
          }

          // process the Jacobian of equation of T
          // if this gate bc is external boundary, set heat flux here
          if( regions[i]->get_advanced_model()->enable_Tl() && (node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion)) )
          {
            regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], TEMPERATURE, A, b, J, omega, add_value_flag);

            AutoDScalar T = node_data->T();         T.setADValue(0, 1.0); // T of this node
            PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
            AutoDScalar fT = Heat_Transfer*(T_external()-T)*S;

            PetscInt index_re = fvm_nodes[i]->global_offset() + node_Tl_offset;
            PetscInt col_re   = index_re;
            MatSetValue(A, index_re, col_re, fT.getADValue(0), ADD_VALUES);

            PetscInt index_im = fvm_nodes[i]->global_offset() + n_variables + node_Tl_offset;
            PetscInt col_im   = index_im;
            MatSetValue(A, index_im, col_im, fT.getADValue(0), ADD_VALUES);
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
            PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();

            // area of out surface of control volume related with neighbor node
            PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node);
            AutoDScalar D = node_data->eps()*(V-V_nb)/distance;

            // the 1/dt is replaced by j*omega.
            std::complex <PetscScalar> mna_scaling = ext_circuit()->mna_ac_scaling(omega);
            std::complex <PetscScalar> dJdisp_dV  = mna_scaling*cv_boundary*D.getADValue(0)*j*omega*current_scale;
            std::complex <PetscScalar> dJdisp_dVn = mna_scaling*cv_boundary*D.getADValue(1)*j*omega*current_scale;

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

          //FIXME tunneling current should be considered here

          break;
        }
        // conductor region (gate) which has an interface with insulator region
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



  // the extra equation of gate boundary
  // For ac scan
  //
  //          _____  (Z1)          Ve
  //    -----|_____|----/\/\/\/\-------> to gate electrode (Ve, I)
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
 * update electrode IV
 */
void GateContactBC::DDMAC_Update_Solution(const PetscScalar * lxx, const Mat, const PetscScalar omega)
{

  std::complex<PetscScalar> Vac(lxx[this->local_offset()], lxx[this->local_offset()+1]);
  std::complex<PetscScalar> Iac(0.0, 0.0);
  std::complex<PetscScalar> j(0.0, 1.0);

  // for 2D mesh, system().z_width() is the device dimension in Z direction; for 3D mesh, system().z_width() is 1.0
  PetscScalar current_scale = system().z_width();

  // we use AD again. no matter it is overkill here.
  //the indepedent variable number, we only need 2 here.
  adtl::AutoDScalar::numdir=2;

  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion *region = (*rnode_it).second.first;

      switch ( region->type() )
      {
      case InsulatorRegion:
        {
          const FVM_Node * fvm_node = (*rnode_it).second.second;
          const FVM_NodeData * node_data = fvm_node->node_data();

          // displacement current

          PetscScalar V_re = lxx[fvm_node->local_offset() + region->ebm_variable_offset(POTENTIAL)];
          PetscScalar V_im = lxx[fvm_node->local_offset() + region->ebm_n_variables() + region->ebm_variable_offset(POTENTIAL)];
          std::complex<PetscScalar> V = std::complex<PetscScalar>(V_re, V_im);

          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
          for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
          {
            const FVM_Node *nb_node = (*nb_it).first;

            // distance from nb node to this node
            PetscScalar distance = (*(fvm_node->root_node()) - *(nb_node->root_node())).size();

            PetscScalar Vn_re = lxx[nb_node->local_offset() + region->ebm_variable_offset(POTENTIAL)];
            PetscScalar Vn_im = lxx[nb_node->local_offset() + region->ebm_n_variables() + region->ebm_variable_offset(POTENTIAL)];
            std::complex<PetscScalar> Vn = std::complex<PetscScalar>(Vn_re, Vn_im);

            // area of out surface of control volume related with neighbor node
            PetscScalar cv_boundary = fvm_node->cv_surface_area(nb_node);

            Iac += node_data->eps()*(V-Vn)/distance*j*omega*cv_boundary*current_scale;
          }

          break;
        }
      case MetalRegion:
      case ElectrodeRegion:
        break;

      default: genius_error(); //we should never reach here
      }
    }
  }

  Parallel::sum(Iac);

  this->ext_circuit()->current_ac()   = Iac;
  this->ext_circuit()->potential_ac() = Vac;

}
