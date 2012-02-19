// $Id: mesh_communication.cc,v 1.12 2008/06/15 04:49:30 gdiso Exp $

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
#include <cstring>

// Local Includes -----------------------------------
#include "genius_env.h"
#include "genius_common.h"
#include "log.h"
#include "perf_log.h"
#include "mesh_base.h"
#include "mesh_tools.h"
#include "boundary_info.h"
#include "mesh_communication.h"
#include "parallel.h"
#include "elem.h"
#include "sphere.h"



namespace
{

#ifdef ENABLE_AMR
  const unsigned int packed_elem_header_size = 10;
#else
  const unsigned int packed_elem_header_size = 4;
#endif

}


// ------------------------------------------------------------
// MeshCommunication class members
void MeshCommunication::clear ()
{
  _neighboring_processors.clear();
}


#ifdef HAVE_MPI
void MeshCommunication::find_neighboring_processors (const MeshBase& mesh)
{
  // Don't need to do anything if there is
  // only one processor.
  if (Genius::n_processors() == 1)
    return;

  _neighboring_processors.clear();

  // Get the bounding sphere for the local processor
  Sphere bounding_sphere =
    MeshTools::processor_bounding_sphere (mesh, Genius::processor_id());

  // Just to be sure, increase its radius by 10%.  Sure would suck to
  // miss a neighboring processor!
  bounding_sphere.radius() *= 1.1;

  // Collect the bounding spheres from all processors, test for intersection
  {
    std::vector<float>
    send (4,                         0),
    recv (4*Genius::n_processors(),  0);

    send[0] = bounding_sphere.center()(0);
    send[1] = bounding_sphere.center()(1);
    send[2] = bounding_sphere.center()(2);
    send[3] = bounding_sphere.radius();

    MPI_Allgather (&send[0], send.size(), MPI_FLOAT,
                   &recv[0], send.size(), MPI_FLOAT,
                   PETSC_COMM_WORLD);


    for (unsigned int proc=0; proc<Genius::n_processors(); proc++)
    {
      const Point center (recv[4*proc+0],
                          recv[4*proc+1],
                          recv[4*proc+2]);

      const Real radius = recv[4*proc+3];

      const Sphere proc_sphere (center, radius);

      if (bounding_sphere.intersects(proc_sphere))
        _neighboring_processors.push_back(proc);
    }

    // Print out the _neighboring_processors list
    std::cout << "Processor " << Genius::processor_id()
    << " intersects:" << std::endl;
    for (unsigned int p=0; p<_neighboring_processors.size(); p++)
      std::cout << " " << _neighboring_processors[p] << std::endl;
  }
}
#else
void MeshCommunication::find_neighboring_processors (const MeshBase&)
{}
#endif



void MeshCommunication::broadcast (MeshBase& mesh) const
{
  // Don't need to do anything if there is
  // only one processor.
  if (Genius::n_processors() == 1)
    return;

  MESSAGE<<"Synchronize mesh with all the processors..."<<std::endl;  RECORD();

  this->broadcast_mesh (mesh);
  this->broadcast_bcs  (mesh, *(mesh.boundary_info));

  MESSAGE<<"Mesh synchronization finished.\n"<<std::endl;  RECORD();
}



#ifdef HAVE_MPI
void MeshCommunication::broadcast_mesh (MeshBase& mesh) const
#else // avoid spurious gcc warnings
void MeshCommunication::broadcast_mesh (MeshBase&) const
#endif
{
  // Don't need to do anything if there is
  // only one processor.
  if (Genius::n_processors() == 1)
    return;

#ifdef HAVE_MPI

  START_LOG("broadcast_mesh()","MeshCommunication");

  // Explicitly clear the mesh on all but processor 0.
  if (Genius::processor_id() != 0)
    mesh.clear();

  // broadcast magic number
  Parallel::broadcast (mesh.magic_num());

  // Get important sizes
  unsigned int n_nodes      = mesh.n_nodes();
  unsigned int n_elem       = mesh.n_elem();
  unsigned int n_levels     = MeshTools::n_levels(mesh);
  unsigned int total_weight = MeshTools::total_weight(mesh);
  unsigned int n_subdomains = mesh.n_subdomains ();

  // Broadcast the sizes
  {
    std::vector<unsigned int> buf (4);

    buf[0] = n_nodes;
    buf[1] = n_elem;
    buf[2] = total_weight;
    buf[3] = n_subdomains;

    // Broadcast
    Parallel::broadcast (buf);

    if (Genius::processor_id() != 0)
    {
      n_nodes      = buf[0];
      n_elem       = buf[1];
      total_weight = buf[2];
      mesh.set_n_subdomains () = buf[3];
    }
  }



  // First build up the pts vector which contains
  // the spatial locations of all the nodes
  {
    std::vector<Real> pts;

    // If we are processor 0, we must populate this vector and
    // broadcast it to the other processors.
    if (Genius::processor_id() == 0)
    {
      pts.reserve (3*n_nodes);

      MeshBase::node_iterator       it     = mesh.nodes_begin();
      const MeshBase::node_iterator it_end = mesh.nodes_end();

      for (; it != it_end; ++it)
      {
        assert (*it != NULL);
        assert ((*it)->id()*3 == pts.size());

        const Point& p = **it;

        pts.push_back ( p(0) ); // x
        pts.push_back ( p(1) ); // y
        pts.push_back ( p(2) ); // z

      }

    }
    else
      pts.resize (3*n_nodes);

    // Sanity check for all processors
    assert (pts.size() == (3*n_nodes));

    // Broadcast the pts vector
    Parallel::broadcast (pts);

    // Add the nodes we just received if we are not
    // processor 0.
    if (Genius::processor_id() != 0)
    {
      assert (mesh.n_nodes() == 0);

      for (unsigned int i=0; i<pts.size(); i += 3)
      {
        mesh.add_point (Point(pts[i+0],
                              pts[i+1],
                              pts[i+2]),
                        i/3);

      }
    }

    assert (mesh.n_nodes() == n_nodes);
  } // Done distributing the nodes


  // Now build up the elements vector which
  // contains the element types and connectivity
  {
    // The conn array contains the information needed to construct each element.
    // Pack all this information into one communication to avoid two latency hits
    // For each element it is of the form
    // [ level p_level r_flag p_flag etype subdomain_id
    //   self_ID parent_ID which_child node_0 node_1 ... node_n]
    // We cannot use unsigned int because parent_ID can be negative
    std::vector<int> conn;

    // If we are processor 0, we must populate this vector and
    // broadcast it to the other processors.
    if (Genius::processor_id() == 0)
    {
      conn.reserve (packed_elem_header_size*n_elem + total_weight);

      // We start from level 0. This is a bit simpler than in xdr_io.C
      // because we do not have to worry about economizing by group elements
      // of the same type. Element type is simply specified as an
      // entry in the connectivity vector, "conn".
      // By filling conn in order of levels, parents should exist before children
      // are built when we reconstruct the elements on the other processors.

      for (unsigned int level=0; level<=n_levels; ++level)
      {
        MeshBase::element_iterator it = mesh.level_elements_begin(level);
        const MeshBase::element_iterator it_end = mesh.level_elements_end(level);

        for (; it != it_end; ++it)
        {
          assert (*it);
          const Elem* elem = *it;
          pack_element (conn, elem);
        }
      }
    }
    else
      conn.resize (packed_elem_header_size*n_elem + total_weight);

    // Sanity check for all processors
    assert (conn.size() == (packed_elem_header_size*n_elem + total_weight));

    // Broadcast the element connectivity
    Parallel::broadcast (conn);

    // Build the elements we just received if we are not
    // processor 0.
    if (Genius::processor_id() != 0)
    {
      assert (mesh.n_elem() == 0);

      unsigned int cnt = 0;

      // This map keeps track of elements we've previously added to the mesh
      // to avoid O(n) lookup times for parent pointers.
      std::map<unsigned int, Elem*> parents;

      while (cnt < conn.size())
      {
        // Declare the element that we will add
        Elem* elem = NULL;

        // Unpack the element header
#ifdef ENABLE_AMR
        const int level             = conn[cnt++];
        const int p_level           = conn[cnt++];
        const Elem::RefinementState refinement_flag =
          static_cast<Elem::RefinementState>(conn[cnt++]);
        const Elem::RefinementState p_refinement_flag =
          static_cast<Elem::RefinementState>(conn[cnt++]);
#endif
        const ElemType elem_type    = static_cast<ElemType>(conn[cnt++]);
        const unsigned int elem_PID = conn[cnt++];
        const int subdomain_ID      = conn[cnt++];
        const int self_ID           = conn[cnt++];
#ifdef ENABLE_AMR
        const int parent_ID         = conn[cnt++];
        const int which_child       = conn[cnt++];

        if (parent_ID != -1) // Do a log(n) search for the parent
        {
          Elem* my_parent = parents.count(parent_ID) ? parents[parent_ID] : NULL;

          // If the parent was not previously added, we cannot continue.
          if (my_parent == NULL)
          {
            std::cerr << "Parent element with ID " << parent_ID
            << " not found." << std::endl;
            genius_error();
          }

          assert (my_parent->refinement_flag() == Elem::INACTIVE);

          elem = Elem::build(elem_type, my_parent).release();
          my_parent->add_child(elem);

          assert (my_parent->fvm_compatible_type(my_parent->type())==elem->fvm_compatible_type(elem->type()));
          assert (my_parent->child(which_child) == elem);
        }

        else // level 0 element has no parent
        {
          assert (level == 0);
#endif

          // should be able to just use the integer elem_type
          elem = Elem::build(elem_type).release();
#ifdef ENABLE_AMR

        }

        // Assign the IDs
        assert (elem->level() == static_cast<unsigned int>(level));
        elem->set_refinement_flag(refinement_flag);
        elem->set_p_refinement_flag(p_refinement_flag);
        elem->set_p_level(p_level);
#endif
        elem->processor_id() = elem_PID;
        elem->subdomain_id() = subdomain_ID;
        elem->set_id() = self_ID;

        // Add elem to the map of parents, since it may have
        // children to be added later
        parents.insert(std::make_pair(self_ID,elem));

        // Assign the connectivity
        for (unsigned int n=0; n<elem->n_nodes(); n++)
        {
          assert (cnt < conn.size());

          elem->set_node(n) = mesh.node_ptr (conn[cnt++]);
        }
        elem->prepare_for_fvm();
      } // end while cnt < conn.size

      // Iterate in ascending elem ID order
      for (std::map<unsigned int, Elem *>::iterator i =
             parents.begin();
           i != parents.end(); ++i)
      {
        Elem *elem = i->second;
        if (elem)
          mesh.add_elem(elem);
        else
          // We can probably handle this, but we don't expect it
          genius_error();
      }

    } // end if iam != cpu 0


    assert (mesh.n_elem() == n_elem);
  } // Done distributing the elements

  // distribut subdomain information
  {

    std::vector<std::string> labels;
    std::vector<std::string> materials;

    for(unsigned int n_sub = 0; n_sub < mesh.n_subdomains (); n_sub++)
    {
      if (Genius::processor_id() == 0)
      {
        labels.push_back(mesh.subdomain_label_by_id(n_sub));
        materials.push_back(mesh.subdomain_material(n_sub));
      }
    }
    Parallel::broadcast (labels);
    Parallel::broadcast (materials);

    for(unsigned int n_sub = 0; n_sub < mesh.n_subdomains (); n_sub++)
    {
      if (Genius::processor_id() != 0)
      {
        mesh.set_subdomain_label(n_sub, labels[n_sub]);
        mesh.set_subdomain_material(n_sub, materials[n_sub]);
      }
    }
  } // Done distribut subdomain information

  // now we have serial mesh
  mesh.set_serial(true);

  STOP_LOG("broadcast_mesh()","MeshCommunication");

#else

  // no MPI but multiple processors? Huh??
  genius_error();

#endif
}



#ifdef HAVE_MPI
void MeshCommunication::broadcast_bcs (MeshBase& mesh,
                                       BoundaryInfo& boundary_info) const
#else // avoid spurious gcc warnings
void MeshCommunication::broadcast_bcs (MeshBase&,
                                       BoundaryInfo&) const
#endif
{
  // Don't need to do anything if there is
  // only one processor.
  if (Genius::n_processors() == 1)
    return;


#ifdef HAVE_MPI

  START_LOG("broadcast_bcs()","MeshCommunication");

  // Explicitly clear the boundary conditions on all
  // but processor 0.
  if (Genius::processor_id() != 0)
    boundary_info.clear();

  // Build up the list of elements with boundary conditions
  {
    std::vector<unsigned int>       el_id;
    std::vector<unsigned short int> side_id;
    std::vector<short int>          bc_id;

    if (Genius::processor_id() == 0)
      boundary_info.build_side_list (el_id, side_id, bc_id);

    assert (el_id.size() == side_id.size());
    assert (el_id.size() == bc_id.size());

    unsigned int n_bcs = el_id.size();

    // Broadcast the number of bcs to expect from processor 0.
    Parallel::broadcast (n_bcs);

    // Only continue if we have element BCs
    if (n_bcs > 0)
    {
      // Allocate space. On CPU 0, these vectors should already have size n_bcs.
      el_id.resize   (n_bcs);
      side_id.resize (n_bcs);
      bc_id.resize   (n_bcs);

      // Broadcast the element identities
      Parallel::broadcast (el_id);

      // Broadcast the side ids for those elements
      Parallel::broadcast (side_id);

      // Broadcast the bc ids for each side
      Parallel::broadcast (bc_id);

      // Build the boundary_info structure if we aren't processor 0
      if (Genius::processor_id() != 0)
        for (unsigned int e=0; e<n_bcs; e++)
        {
          assert (el_id[e] < mesh.n_elem());

          const Elem* elem = mesh.elem(el_id[e]);

          assert (elem != NULL);

          // sanity: be sure that the element returned by mesh.elem() really has id()==el_id[e]
          genius_assert(elem->id() == el_id[e]);

          assert (side_id[e] < elem->n_sides());

          boundary_info.add_side (elem, side_id[e], bc_id[e]);
        }
    }
  }


  // distribute boundary ids
  std::set<short int> & boundary_ids = boundary_info.get_boundary_ids();
  Parallel::broadcast (boundary_ids);

  // distribute boundary labels
  {
    std::vector<std::string> labels;
    std::vector<std::string> descriptions;
    std::vector<bool> user_defined;

    std::set<short int>::iterator it=boundary_ids.begin();
    for(; it!=boundary_ids.end(); ++it)
    {
      if (Genius::processor_id() == 0)
      {
        labels.push_back(boundary_info.get_label_by_id(*it));
        descriptions.push_back(boundary_info.get_description_by_id(*it));
        user_defined.push_back(boundary_info.boundary_id_has_user_defined_label(*it));
      }
    }

    Parallel::broadcast (labels);
    Parallel::broadcast (descriptions);
    Parallel::broadcast (user_defined);

    it=boundary_ids.begin();
    for(unsigned int n=0; it!=boundary_ids.end(); ++n, ++it)
    {
      if (Genius::processor_id() != 0)
      {
        boundary_info.set_label_to_id(*it, labels[n], user_defined[n]);
        boundary_info.set_description_to_id(*it, descriptions[n]);
      }
    }
  }


  // distribute extra boundary descriptions
  {
    std::vector<std::string> & extra_descriptions = boundary_info.extra_descriptions();
    Parallel::broadcast(extra_descriptions);
  }

  // Build up the list of nodes with boundary conditions
  {
    std::vector<unsigned int> node_id;
    std::vector<short int>    bc_id;

    if (Genius::processor_id() == 0)
      boundary_info.build_node_list (node_id, bc_id);


    assert (node_id.size() == bc_id.size());

    unsigned int n_bcs = node_id.size();

    // Broadcast the number of bcs to expect from processor 0.
    Parallel::broadcast (n_bcs);

    // Only continue if we have nodal BCs
    if (n_bcs > 0)
    {
      // Allocate space, again on CPU 0 this should be a no-op.
      node_id.resize (n_bcs);
      bc_id.resize   (n_bcs);

      // Broadcast the node ids
      Parallel::broadcast (node_id);

      // Broadcast the bc ids for each side
      Parallel::broadcast (bc_id);

      // Build the boundary_info structure if we aren't processor 0
      if (Genius::processor_id() != 0)
        for (unsigned int n=0; n<n_bcs; n++)
        {
          assert (node_id[n] < mesh.n_nodes());

          const Node* node = mesh.node_ptr (node_id[n]);

          assert (node != NULL);

          // sanity: be sure that the node returned by mesh.node_ptr() really has id()==node_id[n]
          genius_assert(node->id() == node_id[n]);

          boundary_info.add_node (node, bc_id[n]);
        }
    }
  }

  STOP_LOG("broadcast_bcs()","MeshCommunication");

#else

  // no MPI but multiple processors? Huh??
  genius_error();

#endif
}


// Pack all this information into one communication to avoid two latency hits
// For each element it is of the form
// [ level p_level r_flag p_flag etype subdomain_id
//   self_ID parent_ID which_child node_0 node_1 ... node_n]
// We cannot use unsigned int because parent_ID can be negative
void MeshCommunication::pack_element (std::vector<int> &conn, const Elem* elem) const
{
  assert (elem != NULL);

#ifdef ENABLE_AMR
  conn.push_back (static_cast<int>(elem->level()));
  conn.push_back (static_cast<int>(elem->p_level()));
  conn.push_back (static_cast<int>(elem->refinement_flag()));
  conn.push_back (static_cast<int>(elem->p_refinement_flag()));
#endif
  conn.push_back (static_cast<int>(elem->type()));
  conn.push_back (static_cast<int>(elem->processor_id()));
  conn.push_back (static_cast<int>(elem->subdomain_id()));
  conn.push_back (elem->id());

#ifdef ENABLE_AMR
  // use parent_ID of -1 to indicate a level 0 element
  if (elem->level() == 0)
  {
    conn.push_back(-1);
    conn.push_back(-1);
  }
  else
  {
    conn.push_back(elem->parent()->id());
    conn.push_back(elem->parent()->which_child_am_i(elem));
  }
#endif

  for (unsigned int n=0; n<elem->n_nodes(); n++)
    conn.push_back (elem->node(n));
}
