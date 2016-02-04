// $Id: partitioner.cc,v 1.4 2008/05/22 04:53:44 gdiso Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// C++ Includes   -----------------------------------


// Local Includes -----------------------------------
#include "mesh_base.h"
#include "partitioner.h"
#include "mesh_base.h"
#include "perf_log.h"
#include "elem.h"


#if defined(HAVE_TR1_UNORDERED_SET)
#include <tr1/unordered_set>
#elif defined(HAVE_TR1_UNORDERED_SET_WITH_STD_HEADER) || defined(HAVE_UNORDERED_SET)
#include <unordered_set>
#endif


Partitioner::Partitioner()
{}


Partitioner::~Partitioner()
{
  this->_clear_cluster();
}


// ------------------------------------------------------------
// Partitioner implementation


void Partitioner::partition (MeshBase& mesh, const std::vector<std::vector<unsigned int> > *cluster, const unsigned int n)
{
  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n;

  // clear cluster if data exist
  this->_clear_cluster();

  this->_build_cluster(mesh, cluster);

  // Call the partitioning function
  this->_do_partition(mesh,n);

  // clear cluster for saving memory
  this->_clear_cluster();

  // Set the node's processor ids
  this->_set_node_processor_ids(mesh);
}





void Partitioner::repartition(MeshBase& mesh, const std::vector<std::vector<unsigned int> > *cluster, const unsigned int n)
{
  // Set the number of partitions in the mesh
  mesh.set_n_partitions()=n;

  // clear cluster if data exist
  this->_clear_cluster();

  this->_build_cluster(mesh, cluster);

  // Call the partitioning function
  this->_do_repartition(mesh,n);

  // clear cluster for saving memory
  this->_clear_cluster();

  // Set the node's processor ids
  this->_set_node_processor_ids(mesh);
}





void Partitioner::single_partition (MeshBase& mesh)
{
  START_LOG("partition()", "Single Partitioner");

  // Loop over all the elements and assign them to processor 0.
  MeshBase::element_iterator       elem_it  = mesh.elements_begin();
  const MeshBase::element_iterator elem_end = mesh.elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
  {
    (*elem_it)->processor_id() = 0;
    (*elem_it)->on_local() = true;
  }

  // For a single partition, all the nodes are on processor 0
  MeshBase::node_iterator       node_it  = mesh.nodes_begin();
  const MeshBase::node_iterator node_end = mesh.nodes_end();

  for ( ; node_it != node_end; ++node_it)
  {
    (*node_it)->processor_id() = 0;
    (*node_it)->on_local() = true;
  }

  STOP_LOG("partition()", "Single Partitioner");
}




void Partitioner::_set_node_processor_ids(MeshBase& mesh)
{
  START_LOG("set_node_processor_ids()", "Partitioner");

  // Unset any previously-set node processor ids
  // (maybe from previous partitionings).
  MeshBase::node_iterator       node_it  = mesh.nodes_begin();
  const MeshBase::node_iterator node_end = mesh.nodes_end();

  for ( ; node_it != node_end; ++node_it)
  {
    (*node_it)->invalidate_processor_id();
    (*node_it)->on_local() = false;
  }

  // Loop over all the elements
  MeshBase::element_iterator       elem_it  = mesh.elements_begin();
  const MeshBase::element_iterator elem_end = mesh.elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
  {
    Elem* elem = *elem_it;
    elem->on_local() = false;

    // For each node, set the processor ID to the min of
    // its current value and this Element's processor id.
    for (unsigned int n=0; n<elem->n_nodes(); ++n)
      elem->get_node(n)->processor_id() = std::min(elem->get_node(n)->processor_id(),
                                          elem->processor_id());
  }

  // elem has on processor node
#if defined(HAVE_UNORDERED_SET)
  std::unordered_set<const Elem *> elem_has_on_process_node;
#elif defined(HAVE_TR1_UNORDERED_SET) || defined(HAVE_TR1_UNORDERED_SET_WITH_STD_HEADER)
  std::tr1::unordered_set<const Elem *> elem_has_on_process_node;
#else
  std::set<const Elem *> elem_has_on_process_node;
#endif
  for (elem_it  = mesh.elements_begin() ; elem_it != elem_end; ++elem_it)
  {
    const Elem* elem = *elem_it;
    for( unsigned int n=0; n<elem->n_nodes(); n++ )
      if( elem->get_node(n)->processor_id() == Genius::processor_id() )
      { elem_has_on_process_node.insert(elem); break; }
  }


  // now we determine local elements and nodes

  for (elem_it  = mesh.elements_begin() ; elem_it != elem_end; ++elem_it)
  {
    Elem* elem = *elem_it;

    // the element is belongs to this processor
    if(elem->processor_id() == Genius::processor_id())
      goto on_local;

    // otherwise, it has at lease one node belongs to this processor
    // or my neighbor on this processor
    else
    {
      if( elem_has_on_process_node.find(elem) != elem_has_on_process_node.end() )
        goto on_local;
    }

    // also, if elem has a on_local neighbor
    for( unsigned int s=0; s<elem->n_sides(); s++ )
    {
      Elem* neighbor_elem = elem->neighbor(s);
      if(!neighbor_elem) continue;

      if(neighbor_elem->processor_id() == Genius::processor_id())
        goto on_local;
      else
      {
        if( elem_has_on_process_node.find(neighbor_elem) != elem_has_on_process_node.end() )
          goto on_local;
      }
    }

    continue;

  on_local:
    {
      elem->on_local() = true;
      // if we find it, make all its node on_local
      for( unsigned int n=0; n<elem->n_nodes(); n++ )
        elem->get_node(n)->on_local() = true;
    }
  }


#if 0
  MeshBase::element_iterator       el  = mesh.local_elements_begin();
  const MeshBase::element_iterator end = mesh.local_elements_end();
  for (; el != end; ++el)
  {
    Elem * elem = *el;
    genius_assert(elem->on_local());

    for (unsigned int n=0; n<elem->n_nodes(); n++)
    {
      Node * node = elem->get_node(n);
      genius_assert( node->on_local() );
    }
  }
#endif

  STOP_LOG("set_node_processor_ids()", "Partitioner");
}



void Partitioner::_build_cluster(MeshBase& mesh, const std::vector<std::vector<unsigned int> > *cluster_elems)
{
  if(!cluster_elems)
  {
    Cluster * cluster = new Cluster;
    MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
    const MeshBase::element_iterator elem_end = mesh.active_elements_end();
    for (unsigned int n=0; elem_it != elem_end; ++elem_it, ++n)
    {
      const Elem* elem = *elem_it;
      cluster->elems.push_back(elem);
      cluster->elem_id_map.insert(std::make_pair(elem, n));
    }
    _clusters.push_back(cluster);
  }
  else
  {
    for(unsigned int c=0; c<cluster_elems->size(); ++c)
    {
       Cluster * cluster = new Cluster;
       const std::vector<unsigned int> elems = (*cluster_elems)[c];
       for(unsigned int n=0; n<elems.size(); ++n)
       {
          const Elem * elem = mesh.elem(elems[n]);
          cluster->elems.push_back(elem);
          cluster->elem_id_map.insert(std::make_pair(elem, n));
       }
       _clusters.push_back(cluster);
    }
  }
}



void Partitioner::_clear_cluster()
{
  for(unsigned int n=0; n<_clusters.size(); ++n)
    delete _clusters[n];
  _clusters.clear();
}





