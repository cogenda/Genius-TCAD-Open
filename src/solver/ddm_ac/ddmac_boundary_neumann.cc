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


// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_neumann.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;



void NeumannBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag )
{

  // Neumann boundary condition is processed here

  // note, we will use ADD_VALUES to set values of matrix A
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( add_value_flag != ADD_VALUES && add_value_flag != NOT_SET_VALUES)
  {
    MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);
  }

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const SimulationRegion * region = (*region_node_begin(*node_it)).second.first;
    const FVM_Node * fvm_node = (*region_node_begin(*node_it)).second.second;

    region->DDMAC_Fill_Nodal_Matrix_Vector(fvm_node, A, b, J, omega, add_value_flag);

    // process heat flow out of neumann boundary
    if(region->get_advanced_model()->enable_Tl())
    {
      //the indepedent variable number, we only need 1 here.
      adtl::AutoDScalar::numdir=1;

      // process governing equation of T, which should consider heat exchange to entironment
      AutoDScalar T = fvm_node->node_data()->T();  T.setADValue(0, 1.0); // lattice temperature

      // add heat flux out of Neumann boundary to lattice temperature equatiuon
      PetscScalar S  = fvm_node->outside_boundary_surface_area();
      PetscScalar h = this->Heat_Transfer();
      AutoDScalar fT = h*(T_external()-T)*S;
      PetscInt real_row = fvm_node->global_offset() + region->ebm_variable_offset(TEMPERATURE);
      PetscInt imag_row = fvm_node->global_offset() + region->ebm_n_variables() + region->ebm_variable_offset(TEMPERATURE) ;

      MatSetValue(A, real_row, real_row, fT.getADValue(0),  ADD_VALUES);
      MatSetValue(A, imag_row, imag_row, fT.getADValue(0),  ADD_VALUES);
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}
