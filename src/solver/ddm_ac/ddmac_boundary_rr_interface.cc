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
#include "boundary_condition_rr.h"
#include "petsc_utils.h"



void ResistanceResistanceBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag )
{
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * resistance_fvm_node_1  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * resistance_fvm_node_2  = get_region_fvm_node ( ( *node_it ), _r2 );

    _r1->DDMAC_Fill_Nodal_Matrix_Vector(resistance_fvm_node_1, A, b, J, omega, add_value_flag);
    _r2->DDMAC_Fill_Nodal_Matrix_Vector(resistance_fvm_node_2, A, b, J, omega, add_value_flag, _r1, resistance_fvm_node_1);
    _r2->DDMAC_Force_equal(resistance_fvm_node_2, A, add_value_flag, _r1, resistance_fvm_node_1);

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      if ( region->type() == InsulatorRegion )
      {
        const FVM_Node * insulator_fvm_node = (*rnode_it).second.second;
        region->DDMAC_Fill_Nodal_Matrix_Vector(insulator_fvm_node, A, b, J, omega, add_value_flag, _r1, resistance_fvm_node_1);
        // phi_insulator = 0.5*(phi_r1 + phi_r2)
        {
          PetscInt real_row = insulator_fvm_node->global_offset() + region->ebm_variable_offset(POTENTIAL);
          PetscInt imag_row = insulator_fvm_node->global_offset() + region->ebm_n_variables() + region->ebm_variable_offset(POTENTIAL);
          PetscInt real_col[2]={real_row, resistance_fvm_node_1->global_offset() + _r1->ebm_variable_offset(POTENTIAL)};
          PetscInt imag_col[2]={imag_row, resistance_fvm_node_1->global_offset() + _r1->ebm_n_variables() + _r1->ebm_variable_offset(POTENTIAL)};
          PetscScalar diag[2] = {1.0, -0.5};
          MatSetValues(A, 1, &real_row, 2, real_col, diag, ADD_VALUES);
          MatSetValues(A, 1, &imag_row, 2, imag_col, diag, ADD_VALUES);
        }
        {
          PetscInt real_row = insulator_fvm_node->global_offset() + region->ebm_variable_offset(POTENTIAL);
          PetscInt imag_row = insulator_fvm_node->global_offset() + region->ebm_n_variables() + region->ebm_variable_offset(POTENTIAL);
          PetscInt real_col[2]={real_row, resistance_fvm_node_2->global_offset() + _r2->ebm_variable_offset(POTENTIAL)};
          PetscInt imag_col[2]={imag_row, resistance_fvm_node_2->global_offset() + _r2->ebm_n_variables() + _r2->ebm_variable_offset(POTENTIAL)};
          PetscScalar diag[2] = {1.0, -0.5};
          MatSetValues(A, 1, &real_row, 2, real_col, diag, ADD_VALUES);
          MatSetValues(A, 1, &imag_row, 2, imag_col, diag, ADD_VALUES);
        }
      }
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
