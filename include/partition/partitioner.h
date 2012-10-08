// $Id: partitioner.h,v 1.2 2008/05/15 02:36:46 gdiso Exp $

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



#ifndef __partitioner_h__
#define __partitioner_h__

// C++ Includes   -----------------------------------
#include <vector>
#include <map>

// Local Includes -----------------------------------
#include "genius_common.h"
#include "genius_env.h"


// Forward Declarations
class Elem;
class MeshBase;



/**
 * The \p Partitioner class provides a uniform interface for
 * partitioning algorithms.  It takes a reference to a \p MeshBase
 * object as input, which it will partition into a number of
 * subdomains.
 */

// ------------------------------------------------------------
// Partitioner class definition
class Partitioner
{
 public:

  /**
   * Constructor.
   */
   Partitioner ();

  /**
   * Destructor. Virtual so that we can derive from this class.
   */
   virtual ~Partitioner();


  /**
   * Partition the \p MeshBase into \p n parts.  If the
   * user does not specify a number of pieces into which the
   * mesh should be partitioned, then the default behavior
   * of the partitioner is to partition according to the number
   * of processors defined in libMesh::n_processors().
   * The partitioner currently does not modify the subdomain_id
   * of each element.  This number is reserved for things like
   * material properties, etc.
   */
   void partition (MeshBase& mesh, const std::vector<std::vector<unsigned int> > * =NULL, const unsigned int =Genius::n_processors());

  /**
   * Repartitions the \p MeshBase into \p n parts.  This
   * is required since some partitoning algorithms can repartition
   * more efficiently than computing a new partitioning from scratch.
   * The default behavior is to simply call this->partition(n)
   */
   void repartition (MeshBase& mesh, const std::vector<std::vector<unsigned int> > * =NULL, const unsigned int =Genius::n_processors());


protected:

  /**
   * Trivially "partitions" the mesh for one processor.
   * Simply loops through the elements and assigns all of them
   * to processor 0.  It is provided as a separate function
   * so that derived classes may use it without reimplementing it.
   */
  void single_partition (MeshBase& mesh);

  /**
   * This is the actual partitioning method which must be overloaded
   * in derived classes.  It is called via the public partition()
   * method above by the user.
   */
  virtual void _do_partition(MeshBase& mesh, const unsigned int n) = 0;

  /**
   * This is the actual re-partitioning method which can be overloaded
   * in derived classes.  Note that the default behavior is to simply
   * call the partition function.
   */
  virtual void _do_repartition (MeshBase& mesh, const unsigned int n) { this->_do_partition (mesh, n); }

  /**
   * This function is called after partitioning to set the processor IDs
   * for the nodes.  By definition, a Node's processor ID is the minimum
   * processor ID for all of the elements which share the node.
   */
  void _set_node_processor_ids(MeshBase& mesh);

  /**
   * cluster of mesh elements, the elems belongs to the same cluster
   * will always be partitioned into the same block
   */
  struct Cluster
  {
    unsigned int id;
    std::vector<const Elem *> elems;
  };

  /**
   * all the clusters
   */
  std::vector<Cluster *> _clusters;

  /**
   * map elem to cluster _elem_cluster_map[elem->id]
   */
  std::vector<const Cluster *> _elem_cluster_map;

  /**
   * build the cluster, can be rewrite by derived class
   */
  virtual void _build_flat_cluster(MeshBase& mesh);

  /**
   * merge elems into cluster
   */
  virtual void _merge_elem_to_cluster(MeshBase& mesh, const std::vector<std::vector<unsigned int> > *);

  /**
   * get neighbor elems of a cluster
   */
  std::vector<const Elem *> cluster_neighbor_elem(const Cluster *) const;

  /**
   * clear the cluster
   */
  void _clear_cluster();
};




#endif
