// $Id: metis_partitioner.cc,v 1.5 2008/05/22 04:53:44 gdiso Exp $

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
#include "elem.h"
#include "mesh_base.h"
#include "metis_partitioner.h"
#include "perf_log.h"
#include "parallel.h"


#if (defined(PETSC_HAVE_PARMETIS) || defined(PETSC_HAVE_METIS))

namespace Metis
{
  extern "C"
  {
#     include "metis.h"
  }
}

#endif

#include "linear_partitioner.h"


// ------------------------------------------------------------
// MetisPartitioner implementation
void MetisPartitioner::_do_partition (MeshBase& mesh,
                                      const unsigned int n_pieces)
{
  assert (n_pieces > 0);

  // Check for an easy return
  if (n_pieces == 1)
  {
    this->single_partition (mesh);
    return;
  }


#if (defined(PETSC_HAVE_PARMETIS) || defined(PETSC_HAVE_METIS))

  START_LOG("partition()", "MetisPartitioner");

  for(unsigned int c=0; c<_clusters.size(); ++c)
  {
    Cluster * cluster = _clusters[c];

    const unsigned int n_elem        = cluster->elems.size();
    std::vector<int> part(n_elem);  // here stores the partition vector of the graph
    int metis_error=0;

    // only the first process do the partition
    if(_serial_partition || Genius::is_first_processor())
    {
      // build the graph
      std::vector<int> xadj;          // the adjacency structure of the graph
      std::vector<int> adjncy;        // the adjacency structure of the graph
      std::vector<int> options(8);
      std::vector<int> vwgt(n_elem);  // the weights of the vertices

      xadj.reserve(n_elem+1);

      int n = static_cast<int>(n_elem);  // number of "nodes" (elements) in the graph
      int ncon    = 1;                          // The number of balancing constraints. It should be at least 1.
      int wgtflag = 2;                          // weights on vertices only, none on edges
      int numflag = 0;                          // C-style 0-based numbering
      int nparts  = static_cast<int>(n_pieces); // number of subdomains to create
      int edgecut = 0;                          // the numbers of edges cut by the resulting partition

      // Set the options
      options[0] = 0; // use default options

      // build the graph in CSR format.  Note that
      // the edges in the graph will correspond to
      // face neighbors
      {

        adjncy.reserve (n_elem*6);

        for (unsigned int n=0; n<n_elem; ++n)
        {
          const Elem * elem = cluster->elems[n];

          // The weight is used to define what a balanced graph is
          vwgt[n] = mesh.subdomain_weight( elem->subdomain_id () ) * elem->n_nodes();

          // The beginning of the adjacency array for this elem
          xadj.push_back(adjncy.size());

          // Loop over the element's neighbors.  An element
          // adjacency corresponds to a face neighbor
          for (unsigned int ms=0; ms<elem->n_neighbors(); ms++)
          {
            const Elem * neighbor = elem->neighbor(ms);
            if(neighbor && neighbor->subdomain_id() == elem->subdomain_id())
              adjncy.push_back (cluster->elem_id_map.find(neighbor)->second);
          }
        }

        // The end of the adjacency array for the last elem
        xadj.push_back(adjncy.size());
      } // done building the graph

      genius_assert (xadj.size() == n_elem+1);

      if (adjncy.empty())
        adjncy.push_back(0);



#if PETSC_VERSION_GE(3,3,0)
      // METIS-5 interface
      if (n_pieces <= 8)
        metis_error = Metis::METIS_PartGraphRecursive(&n, &ncon, &xadj[0], &adjncy[0], &vwgt[0], NULL/*vsize*/, NULL/*adjwgt*/,
                                             &nparts, NULL, NULL, NULL, &edgecut, &part[0]);
      else
        metis_error = Metis::METIS_PartGraphKway     (&n, &ncon, &xadj[0], &adjncy[0], &vwgt[0], NULL/*vsize*/, NULL/*adjwgt*/,
                                             &nparts, NULL, NULL, NULL, &edgecut, &part[0]);

      if(metis_error == Metis::METIS_OK) metis_error = 0;
#else
      // old METIS-4 interface
      Metis::METIS_PartGraphKway(&n, &xadj[0], &adjncy[0], &vwgt[0], NULL, &wgtflag, &numflag,
                                             &nparts, &options[0], &edgecut, &part[0]);
#endif
    }

    // broadcast partition info to all the processores
    if(!_serial_partition)
    {
      Parallel::broadcast(metis_error, 0);
    }
    else
    {
      Parallel::sum(metis_error);
    }

    // Assign the returned processor ids.  The part array contains
    // the processor id for each active element, but in terms of
    // the contiguous indexing we defined above

    if( !metis_error )
    {
      if(!_serial_partition)
        Parallel::broadcast(part, 0);

      for (unsigned int n=0; n<n_elem; ++n)
      {
        Elem * elem = const_cast<Elem *>(cluster->elems[n]);
        short int processor_id = static_cast<short int>(part[n]);
        elem->processor_id() = processor_id;
      }
    }
    else // linear partition
    {
      const unsigned int n_elem     = cluster->elems.size();
      const unsigned int blksize    = n_elem/n_pieces;

      for (unsigned int n=0; n<n_elem; ++n)
      {
        short int processor_id;
        if ((n/blksize) < n_pieces)
          processor_id = n/blksize;
        else
          processor_id = 0;

        Elem * elem = const_cast<Elem *>(cluster->elems[n]);
        elem->processor_id() = processor_id;
      }
    }

  }

  STOP_LOG("partition()", "MetisPartitioner");

#else

  LinearPartitioner *lp = new LinearPartitioner;
  lp->partition(mesh, 0, n_pieces);
  delete lp;

#endif
}


