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

    // we should only have one resistance region
    const SimulationRegion * resistance_region = get_fvm_node_region((*node_it), MetalRegion);
    const FVM_Node * resistance_fvm_node = get_region_fvm_node((*node_it), MetalRegion);
    unsigned int resistance_fvm_node_re = resistance_fvm_node->global_offset();
    unsigned int resistance_fvm_node_im = resistance_region->ebm_n_variables() + resistance_fvm_node->global_offset();

    resistance_region->DDMAC_Fill_Nodal_Matrix_Vector(resistance_fvm_node, POTENTIAL, A, b, J, omega, add_value_flag);
    if(resistance_region->get_advanced_model()->enable_Tl())
      resistance_region->DDMAC_Fill_Nodal_Matrix_Vector(resistance_fvm_node, TEMPERATURE, A, b, J, omega, add_value_flag);

    // however, this resistance region may contact several insulator region

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      if( region->type() == InsulatorRegion )
      {
        //the indepedent variable number, we need 2 here.
        adtl::AutoDScalar::numdir=2;

        const FVM_Node * insulator_fvm_node = (*rnode_it).second.second;
        const FVM_NodeData * insulator_fvm_node_data = insulator_fvm_node->node_data();

        // compute displacement current
        FVM_Node::fvm_neighbor_node_iterator nb_it = insulator_fvm_node->neighbor_node_begin();
        for(; nb_it != insulator_fvm_node->neighbor_node_end(); ++nb_it)
        {
          const FVM_Node *nb_node = (*nb_it).first;

          // the psi of this node
          AutoDScalar  V = insulator_fvm_node_data->psi(); V.setADValue(0, 1.0);
          // the psi of neighbor node
          AutoDScalar V_nb = nb_node->node_data()->psi(); V_nb.setADValue(1, 1.0);

          // distance from nb node to this node
          PetscScalar distance = insulator_fvm_node->distance(nb_node);

          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = insulator_fvm_node->cv_surface_area(nb_node);
          AutoDScalar D = insulator_fvm_node_data->eps()*(V-V_nb)/distance;

          // the 1/dt is replaced by j*omega.
          std::complex <PetscScalar> dJdisp_dV  = -cv_boundary*D.getADValue(0)*j*omega;
          std::complex <PetscScalar> dJdisp_dVn = -cv_boundary*D.getADValue(1)*j*omega;

          // V
          MatSetValue(A, resistance_fvm_node_re, insulator_fvm_node->global_offset(), dJdisp_dV.real(), ADD_VALUES);
          MatSetValue(A, resistance_fvm_node_re, insulator_fvm_node->global_offset()+region->ebm_n_variables(), -dJdisp_dV.imag(), ADD_VALUES);

          MatSetValue(A, resistance_fvm_node_im, insulator_fvm_node->global_offset(), dJdisp_dV.imag(), ADD_VALUES);
          MatSetValue(A, resistance_fvm_node_im, insulator_fvm_node->global_offset()+region->ebm_n_variables(),  dJdisp_dV.real(), ADD_VALUES);

          // V_nb
          MatSetValue(A, resistance_fvm_node_re, nb_node->global_offset(), dJdisp_dVn.real(), ADD_VALUES);
          MatSetValue(A, resistance_fvm_node_re, nb_node->global_offset()+region->ebm_n_variables(), -dJdisp_dVn.imag(), ADD_VALUES);

          MatSetValue(A, resistance_fvm_node_im, nb_node->global_offset(), dJdisp_dVn.imag(), ADD_VALUES);
          MatSetValue(A, resistance_fvm_node_im, nb_node->global_offset()+region->ebm_n_variables(),  dJdisp_dVn.real(), ADD_VALUES);
        }

        // the psi of insulator region equals to metal region
        region->DDMAC_Force_equal(insulator_fvm_node, POTENTIAL, A, add_value_flag, resistance_region, resistance_fvm_node);

        if(region->get_advanced_model()->enable_Tl())
        {
          // load Jacobian entry of this node from J and fill into AC matrix A
          region->DDMAC_Fill_Nodal_Matrix_Vector(insulator_fvm_node, TEMPERATURE, A, b, J, omega, add_value_flag, resistance_region, resistance_fvm_node);
          region->DDMAC_Force_equal(insulator_fvm_node, TEMPERATURE, A, add_value_flag, resistance_region, resistance_fvm_node);
        }
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
