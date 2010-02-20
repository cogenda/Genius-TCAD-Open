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

#include "boundary_info.h"
#include "fvm_pde_solver.h"



void FVM_PDESolver::set_serial_dof_map()
{

  // set all the global/local offset to invalid_uint
  for(unsigned int n=0; n<_system.n_regions(); ++n)
  {
    SimulationRegion * region = _system.region(n);

    SimulationRegion::node_iterator it = region->nodes_begin();
    SimulationRegion::node_iterator it_end = region->nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = (*it).second;
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

    SimulationRegion::node_iterator it = region->nodes_begin();
    SimulationRegion::node_iterator it_end = region->nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * node = (*it).second;
      assert(node!=NULL);

      //if this node belongs to this processor, set the local offset
      // then we can make sure that each partition has a continuous block
      assert ( node->root_node()->processor_id() == Genius::processor_id() );
      {
        node->set_local_offset(n_local_dofs);
        node->set_global_offset(n_local_dofs);
        n_local_dofs += this->node_dofs( region );
      }
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


  // compute the nonzero pattern of matrix
  // search for all the regions...
  n_nz.resize(n_local_dofs, 0);

  for(unsigned int n=0; n<_system.n_regions(); ++n)
  {
    const SimulationRegion * region = _system.region(n);

    SimulationRegion::const_node_iterator it = region->nodes_begin();
    SimulationRegion::const_node_iterator it_end = region->nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = (*it).second;

      // skip node not belongs to this processor
      if( fvm_node->root_node()->processor_id() != Genius::processor_id() ) continue;

      unsigned int local_offset = fvm_node->local_offset();
      unsigned int local_node_dofs = this->node_dofs( region );
      assert(local_offset!=invalid_uint);

      std::vector<std::pair<unsigned int, unsigned int> > v_region_nodes;
      std::vector<std::pair<unsigned int, unsigned int> >::iterator itn;
      unsigned int node_dofs=0;

      //only one processor? should have no off_processor_dof
      unsigned int off_processor_node_dofs=0;

      // all the nodes involved
      fvm_node->PDE_node_pattern(v_region_nodes, this->all_neighbor_elements_involved(region));
      for(itn=v_region_nodes.begin(); itn!=v_region_nodes.end(); ++itn)
      {
        const SimulationRegion * _region = _system.region((*itn).first);
        unsigned int dof = this->node_dofs( _region );
        unsigned int node_num = (*itn).second;
        node_dofs += node_num*dof;
      }


      // set the nonzero pattern
      for(unsigned int i=0; i<local_node_dofs; ++i)
      {
        n_nz[local_offset + i] = node_dofs-off_processor_node_dofs;
      }

      // not a boundary fvm_node? that's all
      if( fvm_node->boundary_id()==BoundaryInfo::invalid_id ) continue;

      // for boundary node, we need to consider extra dofs contributed by equ of boundary condition
      unsigned int bc_index = _system.get_bcs()->get_bc_index_by_bd_id(fvm_node->boundary_id());
      const BoundaryCondition * bc = _system.get_bcs()->get_bc(bc_index);
      // the dof of this boundary condition
      unsigned int bc_dofs = this->bc_dofs( bc );

      // reserve for bc_dofs
      for(unsigned int i=0; i<local_node_dofs; ++i)
        n_nz[local_offset + i] += bc_dofs;
    }
  }


  //  set n_nz and n_oz for boundary extra equation
  if(_system.get_bcs()!=NULL)
  {
    unsigned int id = 0;

    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); ++b )
    {
      // get the boundary condition
      const BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      // get the nodes belongs to this boundary condition
      const std::vector<const Node *> & bc_nodes = bc->nodes();

      // the dofs of this boundary condition
      unsigned int bc_dofs = this->bc_dofs( bc );

      // the bandwidth of this boundary condition
      unsigned int bc_bandwidth = this->bc_bandwidth( bc );

      // statistic the matrix bandwidth contributed by boundary node
      unsigned int on_processor_dofs = 0;

      // here has some over kill
      for(unsigned int n=0; n<bc_nodes.size(); ++n  )
      {
        unsigned int neighbors = bc->n_node_neighbors(bc_nodes[n]);
        on_processor_dofs += (neighbors+1)*this->bc_node_dofs( bc ) ;
      }

      // prevent overflow, this may be happened for very small problems.
      if ( on_processor_dofs + bc_bandwidth > n_local_dofs )
      { on_processor_dofs = n_local_dofs - bc_bandwidth; }

      // assign to n_nz
      for(unsigned int i=0; i<bc_dofs; ++i, ++id)
      {
        //n_local_dofs - n_global_bc_dofs is the beginning offset of boundary dofs
        n_nz[n_local_dofs - n_global_bc_dofs + id] = on_processor_dofs + bc_bandwidth;
      }
    }
  }

  // set n_nz and n_oz for extra dofs
  if(this->extra_dofs())
    this->set_extra_matrix_nonzero_pattern();
}
