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


#include <numeric>

#include "fvm_flex_pde_solver.h"



void FVM_FlexPDESolver::set_serial_dof_map()
{
  // set all the global/local offset to invalid_uint
  for(unsigned int n=0; n<_system.n_regions(); ++n)
  {
    SimulationRegion * region = _system.region(n);

    SimulationRegion::local_node_iterator it = region->on_local_nodes_begin();
    SimulationRegion::local_node_iterator it_end = region->on_local_nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = (*it);
      fvm_node->set_local_offset(invalid_uint);
      fvm_node->set_global_offset(invalid_uint);
    }
  }

  // the local index of dof
  n_local_dofs = 0;

  //search for all the regions to build the index of nodal dof
  for(unsigned int n=0; n<_system.n_regions(); ++n)
  {
    SimulationRegion * region = _system.region(n);
    const unsigned int region_node_dofs = this->node_dofs( region );

    SimulationRegion::local_node_iterator it = region->on_local_nodes_begin();
    SimulationRegion::local_node_iterator it_end = region->on_local_nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = (*it);

      fvm_node->set_local_offset(n_local_dofs);
      fvm_node->set_global_offset(n_local_dofs);
      n_local_dofs += region_node_dofs;
    }
  }

  // the total node's dof number
  n_global_node_dofs = n_local_dofs;


  // build the global and local index arrays, they are the same for serial
  local_index_array.clear();
  global_index_array.clear();
  for(unsigned int i=0; i<n_local_dofs; ++i )
  {
    local_index_array.push_back(i);
    global_index_array.push_back(i);
  }


  // we compute the dofs of boundary condition here.
  // these dofs will be added at the end of global_node_dofs
  // so it will not affect previous result.
  // as a result, only the last processor will hold boundary dofs
  n_global_bc_dofs = 0;
  if(_system.get_bcs()!=NULL)
  {
    for(unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n )
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(n);
      unsigned int bc_dofs = this->bc_dofs( bc );

      if( bc_dofs >0 )
      {
        // set the global/local offset of extra bc governing equation
        bc->set_global_offset( n_global_node_dofs + n_global_bc_dofs );
        bc->set_local_offset ( n_global_node_dofs + n_global_bc_dofs );
        bc->set_array_offset ( n_global_node_dofs + n_global_bc_dofs );

        // the extra bc variable should in scatter list
        for(unsigned int i=0; i<bc_dofs; ++i)
        {
          global_index_array.push_back (n_global_node_dofs + n_global_bc_dofs + i);
          local_index_array.push_back  (n_global_node_dofs + n_global_bc_dofs + i);
        }

        n_global_bc_dofs +=  bc_dofs;
      }
      // no extra equation for this bc
      else
      {
        bc->set_global_offset( invalid_uint );
        bc->set_local_offset( invalid_uint );
      }
    }
  }

  unsigned int n_extra_dofs = this->extra_dofs();
  // all the processor should know this value
  n_global_dofs = n_global_node_dofs + n_global_bc_dofs + n_extra_dofs;
  n_local_dofs  = n_global_dofs;
  for(unsigned int i=0; i<n_extra_dofs; ++i )
  {
    local_index_array.push_back(n_global_dofs - n_extra_dofs + i);
    global_index_array.push_back(n_global_dofs - n_extra_dofs +i);
  }

  // for extra dofs
  this->set_extra_matrix_nonzero_pattern();
}
