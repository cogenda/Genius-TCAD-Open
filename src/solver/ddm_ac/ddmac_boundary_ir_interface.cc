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


#include "simulation_system.h"
#include "resistance_region.h"
#include "insulator_region.h"
#include "boundary_condition_ir.h"
#include "petsc_utils.h"


void ResistanceInsulatorBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag )
{
  std::complex <PetscScalar> j(0.0, 1.0);

  // after that, set new Jacobian entrance to source rows
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const SimulationRegion * resistance_region = get_fvm_node_region((*node_it), MetalRegion);
    const FVM_Node * resistance_fvm_node = get_region_fvm_node((*node_it), MetalRegion);
    unsigned int resistance_fvm_node_re = resistance_fvm_node->global_offset();
    unsigned int resistance_fvm_node_im = resistance_region->ebm_n_variables() + resistance_fvm_node->global_offset();

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
        // Electrode-Insulator interface at Insulator side
      case InsulatorRegion:
        {
          genius_assert(i==0);

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

          // compute displacement current
          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
          for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
          {
            const FVM_Node *nb_node = (*nb_it).second;

            // the psi of this node
            AutoDScalar  V = node_data->psi(); V.setADValue(0, 1.0);
            // the psi of neighbor node
            AutoDScalar V_nb = nb_node->node_data()->psi(); V_nb.setADValue(1, 1.0);

            // distance from nb node to this node
            PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();

            // area of out surface of control volume related with neighbor node
            PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node->root_node());
            AutoDScalar D = node_data->eps()*(V-V_nb)/distance;

            // the 1/dt is replaced by j*omega.
            std::complex <PetscScalar> dJdisp_dV  = -cv_boundary*D.getADValue(0)*j*omega;
            std::complex <PetscScalar> dJdisp_dVn = -cv_boundary*D.getADValue(1)*j*omega;

            // V
            MatSetValue(A, resistance_fvm_node_re, fvm_nodes[i]->global_offset(), dJdisp_dV.real(), ADD_VALUES);
            MatSetValue(A, resistance_fvm_node_re, fvm_nodes[i]->global_offset()+regions[i]->ebm_n_variables(), -dJdisp_dV.imag(), ADD_VALUES);

            MatSetValue(A, resistance_fvm_node_im, fvm_nodes[i]->global_offset(), dJdisp_dV.imag(), ADD_VALUES);
            MatSetValue(A, resistance_fvm_node_im, fvm_nodes[i]->global_offset()+regions[i]->ebm_n_variables(),  dJdisp_dV.real(), ADD_VALUES);

            // V_nb
            MatSetValue(A, resistance_fvm_node_re, nb_node->global_offset(), dJdisp_dVn.real(), ADD_VALUES);
            MatSetValue(A, resistance_fvm_node_re, nb_node->global_offset()+regions[i]->ebm_n_variables(), -dJdisp_dVn.imag(), ADD_VALUES);

            MatSetValue(A, resistance_fvm_node_im, nb_node->global_offset(), dJdisp_dVn.imag(), ADD_VALUES);
            MatSetValue(A, resistance_fvm_node_im, nb_node->global_offset()+regions[i]->ebm_n_variables(),  dJdisp_dVn.real(), ADD_VALUES);
          }

          break;
        }
        // Electrode-Insulator interface at Conductor side
       case MetalRegion:
        {
          regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], POTENTIAL, A, b, J, omega, add_value_flag);
          // the psi of insulator region equals to metal region
          regions[0]->DDMAC_Force_equal(fvm_nodes[0], POTENTIAL, A, add_value_flag, regions[i], fvm_nodes[i]);

          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            // load Jacobian entry of this node from J and fill into AC matrix A
            regions[0]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[0], TEMPERATURE, A, b, J, omega, add_value_flag, regions[i], fvm_nodes[i]);
            regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], TEMPERATURE, A, b, J, omega, add_value_flag);
            regions[0]->DDMAC_Force_equal(fvm_nodes[0], TEMPERATURE, A, add_value_flag, regions[i], fvm_nodes[i]);
          }

          break;
        }
      case VacuumRegion:
        break;
      default: genius_error(); //we should never reach here
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
