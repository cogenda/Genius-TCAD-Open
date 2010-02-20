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



#include "genius_common.h"

#include <numeric>

#include "mesh.h"
#include "boundary_info.h"
#include "fem_pde_solver.h"
#include "parallel.h"




void FEM_PDESolver::set_parallel_dof_map()
{

  Mesh & mesh = _system.mesh();

  // the local index of dof
  n_local_dofs = 0;

  // the map<Node *, dof of current node> for ghost node dofs
  std::map<Node *, unsigned int> ghost_map;

  //search for all the nodes to build the local index of nodal dof
  MeshBase::node_iterator nd = mesh.active_nodes_begin();
  MeshBase::node_iterator nd_end = mesh.active_nodes_end();
  for ( ; nd!=nd_end; ++nd )
  {
    Node * node = (*nd);

    //if this node belongs to this processor, set the local dof
    // then we can make sure that each partition has a continuous block
    if( node->processor_id() == Genius::processor_id() )
    {
      node->local_dof_id() = n_local_dofs;
      n_local_dofs += this->node_dofs();
    }
    // node is not on processor but on_local() is true, consider it as a ghost node
    else
    {
      if( node->on_local() )
        ghost_map.insert(std::make_pair(node, this->node_dofs()));
    }
  }


  // after the local index are set, we should know the block size on each processor,
  // then the global offset of each local block is known.
  std::vector<int> block_size;
  block_size.push_back(n_local_dofs);
  Parallel::allgather(block_size);
  genius_assert( block_size.size() == Genius::n_processors() );

  // the total node's dof number
  n_global_dofs = std::accumulate(block_size.begin(), block_size.end(), 0 );

  // the offset of local block at global dof array
  global_offset = std::accumulate(block_size.begin(), block_size.begin()+Genius::processor_id(), 0 );

  // build the global and local index arrays,
  // now only contains dof index belongs to this partition
  local_index_array.clear();
  global_index_array.clear();
  for(unsigned int i=0; i<n_local_dofs; ++i )
  {
    local_index_array.push_back(i);
    global_index_array.push_back(global_offset + i);
  }

  // since we konw the block_size, we can assign global offset to each region node
  // search for all the regions again...
  std::vector<unsigned int> offset_in_partition;
  offset_in_partition.resize(Genius::n_processors(), 0);
  for (nd = mesh.active_nodes_begin(); nd!=nd_end; ++nd )
  {
    Node * node = (*nd);
    genius_assert(node!=NULL);

    unsigned int processor_id = node->processor_id();
    // the global offset of this partition
    unsigned int partition_offset = std::accumulate(block_size.begin(), block_size.begin()+processor_id, 0 );

    node->global_dof_id() = partition_offset + offset_in_partition[processor_id];

    offset_in_partition[processor_id] += this->node_dofs();
  }


  // insert ghost node at the end of local_index_array and global_index_array
  std::map<Node *, unsigned int>::iterator map_it = ghost_map.begin();
  for(; map_it!=ghost_map.end(); ++map_it)
  {
    unsigned int local_offset = local_index_array.size();
    unsigned int global_offset = (*map_it).first->global_dof_id();
    (*map_it).first->local_dof_id() = local_offset;

    for(unsigned int i=0; i<(*map_it).second; ++i)
    {
      local_index_array.push_back(local_offset + i);
      global_index_array.push_back(global_offset + i);
    }
  }


  // compute the nonzero pattern of matrix
  {
    n_nz.resize(n_local_dofs, 0);
    n_oz.resize(n_local_dofs, 0);

    std::map<Node *, std::set<Node *> > node_connect_map;

    MeshBase::const_element_iterator       elem_it     = mesh.local_elements_begin();
    const MeshBase::const_element_iterator elem_it_end = mesh.local_elements_end();
    for ( ; elem_it != elem_it_end; ++elem_it)
    {
      const Elem *elem  = (*elem_it);
      // search for all the nodes
      for(unsigned int n=0; n<elem->n_nodes(); ++n)
      {
        Node * node = elem->get_node(n);
        if( node->processor_id() == Genius::processor_id() )
        {
          for (unsigned int m=0; m<elem->n_nodes(); ++m)
            if ( n!=m && elem->node_node_connect(n,m) )
            {
              genius_assert(elem->get_node(m)->on_local());
              node_connect_map[node].insert(elem->get_node(m));
            }
        }
      }
    }

    std::map<Node *, std::set<Node *> >::iterator it = node_connect_map.begin();
    std::map<Node *, std::set<Node *> >::iterator it_end = node_connect_map.end();
    for(; it!=it_end; ++it)
    {
      Node * node = (*it).first;
      //n_nz for this node
      for(unsigned int i=0; i<this->node_dofs(); ++i )
        n_nz[node->local_dof_id()+i] += this->node_dofs();

      //n_nz/n_oz for neighbor nodes
      std::set<Node *> & neighbors = (*it).second;

      std::set<Node *>::iterator neighbor_it = neighbors.begin();
      for(; neighbor_it!=neighbors.end(); ++neighbor_it )
      {
        Node * neighbor_node = (*neighbor_it);
        if( neighbor_node->processor_id() == Genius::processor_id() )
        {
          for(unsigned int i=0; i<this->node_dofs(); ++i )
            n_nz[node->local_dof_id()+i] += this->node_dofs();
        }
        else
        {
          for(unsigned int i=0; i<this->node_dofs(); ++i )
            n_oz[node->local_dof_id()+i] += this->node_dofs();
        }
      }
    }
  }

}



void FEM_PDESolver::set_serial_dof_map()
{

  Mesh & mesh = _system.mesh();

  // the local index of dof
  n_local_dofs = 0;

  //search for all the nodes to build the local index of nodal dof
  MeshBase::node_iterator nd = mesh.active_nodes_begin();
  MeshBase::node_iterator nd_end = mesh.active_nodes_end();
  for ( ; nd!=nd_end; ++nd )
  {
    Node * node = (*nd);
    genius_assert( node->processor_id() == Genius::processor_id() );

    node->local_dof_id()  = n_local_dofs;
    node->global_dof_id() = node->local_dof_id();
    n_local_dofs += this->node_dofs();
  }

  global_offset = 0;
  n_global_dofs = n_local_dofs;


  // build the global and local index arrays,
  // now only contains dof index belongs to this partition
  local_index_array.clear();
  global_index_array.clear();
  for(unsigned int i=0; i<n_local_dofs; ++i )
  {
    local_index_array.push_back(i);
    global_index_array.push_back(global_offset + i);
  }


  // compute the nonzero pattern of matrix
  {
    n_nz.resize(n_local_dofs, 0);
    n_oz.resize(n_local_dofs, 0);

    std::map<Node *, std::set<Node *> > node_connect_map;

    MeshBase::const_element_iterator       elem_it     = mesh.local_elements_begin();
    const MeshBase::const_element_iterator elem_it_end = mesh.local_elements_end();
    for ( ; elem_it != elem_it_end; ++elem_it)
    {
      const Elem *elem  = (*elem_it);
      // search for all the nodes
      for(unsigned int n=0; n<elem->n_nodes(); ++n)
      {
        Node * node = elem->get_node(n);
        for (unsigned int m=0; m<elem->n_nodes(); ++m)
          if ( n!=m && elem->node_node_connect(n,m) )
            node_connect_map[node].insert(elem->get_node(m));
      }
    }

    std::map<Node *, std::set<Node *> >::iterator it = node_connect_map.begin();
    std::map<Node *, std::set<Node *> >::iterator it_end = node_connect_map.end();
    for(; it!=it_end; ++it)
    {
      Node * node = (*it).first;
      unsigned int n_neighbors = (*it).second.size();

      for(unsigned int i=0; i<this->node_dofs(); ++i )
        n_nz[node->local_dof_id()+i] = (n_neighbors+1) * this->node_dofs();
    }
  }

}


void FEM_PDESolver::build_dof_map()
{
#ifdef COGENDA_COMMERCIAL_PRODUCT
  // for commercial version
  set_parallel_dof_map();
#else
  // for open source version
  set_serial_dof_map();
#endif
}


void FEM_PDESolver::build_dof_indices (const Elem* const elem, std::vector<PetscInt>& di) const
{
  const unsigned int n_nodes = elem->n_nodes();
  const unsigned int n_vars  = this->node_dofs();

  // Clear the DOF indices vector
  di.clear();

  for (unsigned int n=0; n<n_nodes; n++)
  {
    const Node* node = elem->get_node(n);
    // Get the dof numbers
    for (unsigned int v=0; v<n_vars; v++)
    {
      // Get the node-based DOF numbers
      di.push_back(static_cast<PetscInt>(node->global_dof_id()+v));
    }
  }
}
