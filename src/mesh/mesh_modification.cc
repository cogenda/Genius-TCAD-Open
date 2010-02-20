// $Id: mesh_modification.cc,v 1.2 2008/04/17 10:09:16 gdiso Exp $

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



// C++ includes
#include <cmath> // for std::acos()
#include <algorithm>
#include <map>

// Local includes
#include "mesh_tools.h"
#include "mesh_modification.h"
#include "unstructured_mesh.h"
#include "perf_log.h"
#include "boundary_info.h"



// ------------------------------------------------------------
// MeshTools::Modification functions for mesh modification
void MeshTools::Modification::distort (MeshBase& mesh,
                                       const Real factor,
                                       const bool perturb_boundary)
{
  assert (mesh.mesh_dimension() != 1);
  assert (mesh.n_nodes());
  assert (mesh.n_elem());
  assert ((factor >= 0.) && (factor <= 1.));

  START_LOG("distort()", "MeshTools::Modification");



  // First find nodes on the boundary and flag them
  // so that we don't move them
  // on_boundary holds false (not on boundary) and true (on boundary)
  std::vector<bool> on_boundary (mesh.n_nodes(), false);

  if (!perturb_boundary) MeshTools::find_boundary_nodes (mesh, on_boundary);

  // Now calculate the minimum distance to
  // neighboring nodes for each node.
  // hmin holds these distances.
  std::vector<float> hmin (mesh.n_nodes(), 1.e20);

  MeshBase::element_iterator       el  = mesh.active_elements_begin();
  const MeshBase::element_iterator end = mesh.active_elements_end();

  for (; el!=end; ++el)
    for (unsigned int n=0; n<(*el)->n_nodes(); n++)
      hmin[(*el)->node(n)] = std::min(hmin[(*el)->node(n)],
                                      static_cast<float>((*el)->hmin()));


  // Now actually move the nodes
  {
    const unsigned int seed = 123456;

    // seed the random number generator
    std::srand(seed);

    // If the node is on the boundary or
    // the node is not used by any element (hmin[n]<1.e20)
    // then we should not move it.
    // [Note: Testing for (in)equality might be wrong
    // (different types, namely float and double)]
    for (unsigned int n=0; n<mesh.n_nodes(); n++)
      if (!on_boundary[n] && (hmin[n] < 1.e20) )
      {
        // the direction, random but unit normalized

        Point dir( static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX),
                   static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX),
                   ((mesh.mesh_dimension() == 3) ?
                    static_cast<Real>(std::rand())/static_cast<Real>(RAND_MAX) :
                    0.)
                 );

        dir(0) = (dir(0)-.5)*2.;
        dir(1) = (dir(1)-.5)*2.;
        if (mesh.mesh_dimension() == 3)
          dir(2) = (dir(2)-.5)*2.;

        dir = dir.unit();

        mesh.node(n)(0) += dir(0)*factor*hmin[n];
        mesh.node(n)(1) += dir(1)*factor*hmin[n];
        if (mesh.mesh_dimension() == 3)
          mesh.node(n)(2) += dir(2)*factor*hmin[n];
      }
  }


  // All done
  STOP_LOG("distort()", "MeshTools::Modification");
}



void MeshTools::Modification::translate (MeshBase& mesh,
    const Real xt,
    const Real yt,
    const Real zt)
{
  const Point p(xt, yt, zt);

  for (unsigned int n=0; n<mesh.n_nodes(); n++)
    mesh.node(n) += p;
}


// void MeshTools::Modification::rotate2D (MeshBase& mesh,
//                                         const Real alpha)
// {
//   assert (mesh.mesh_dimension() != 1);

//   const Real pi = std::acos(-1);
//   const Real  a = alpha/180.*pi;
//   for (unsigned int n=0; n<mesh.n_nodes(); n++)
//     {
//       const Point p = mesh.node(n);
//       const Real  x = p(0);
//       const Real  y = p(1);
//       const Real  z = p(2);
//       mesh.node(n) = Point(std::cos(a)*x - std::sin(a)*y,
//                            std::sin(a)*x + std::cos(a)*y,
//                            z);
//     }

// }



void MeshTools::Modification::rotate (MeshBase& mesh,
                                      const Real phi,
                                      const Real theta,
                                      const Real psi)
{
  assert (mesh.mesh_dimension() != 1);

  const Real pi = std::acos(-1.);
  const Real  p = -phi/180.*pi;
  const Real  t = -theta/180.*pi;
  const Real  s = -psi/180.*pi;
  const Real sp = std::sin(p), cp = std::cos(p);
  const Real st = std::sin(t), ct = std::cos(t);
  const Real ss = std::sin(s), cs = std::cos(s);

  for (unsigned int n=0; n<mesh.n_nodes(); n++)
  {
    const Point p = mesh.node(n);
    const Real  x = p(0);
    const Real  y = p(1);
    const Real  z = p(2);
    mesh.node(n) = Point(( cp*cs-sp*ct*ss)*x + ( sp*cs+cp*ct*ss)*y + (st*ss)*z,
                         (-cp*ss-sp*ct*cs)*x + (-sp*st+cp*ct*cs)*y + (st*cs)*z,
                         ( sp*st)*x          + (-cp*st)*y          + (ct)*z   );
  }
}


void MeshTools::Modification::scale (MeshBase& mesh,
                                     const Real xs,
                                     const Real ys,
                                     const Real zs)
{
  const Real x_scale = xs;
  Real y_scale       = ys;
  Real z_scale       = zs;

  if (ys == 0.)
  {
    assert (zs == 0.);

    y_scale = z_scale = x_scale;
  }

  // Scale the x coordinate in all dimensions
  for (unsigned int n=0; n<mesh.n_nodes(); n++)
    mesh.node(n)(0) *= x_scale;


  // Only scale the y coordinate in 2 and 3D
  if (mesh.spatial_dimension() > 1)
  {

    for (unsigned int n=0; n<mesh.n_nodes(); n++)
      mesh.node(n)(1) *= y_scale;

    // Only scale the z coordinate in 3D
    if (mesh.spatial_dimension() == 3)
    {
      for (unsigned int n=0; n<mesh.n_nodes(); n++)
        mesh.node(n)(2) *= z_scale;
    }
  }
}



void MeshTools::Modification::smooth (MeshBase& mesh,
                                      const unsigned int n_iterations,
                                      const Real power)
{
  /**
   * This implementation assumes every element "side" has only 2 nodes.
   */
  assert (mesh.mesh_dimension() == 2);

  /*
   * find the boundary nodes
   */
  std::vector<bool>  on_boundary;
  MeshTools::find_boundary_nodes(mesh, on_boundary);

  for (unsigned int iter=0; iter<n_iterations; iter++)

  {
    /*
     * loop over the mesh refinement level
     */
    unsigned int n_levels = MeshTools::n_levels(mesh);
    for (unsigned int refinement_level=0; refinement_level != n_levels;
         refinement_level++)
    {
      // initialize the storage (have to do it on every level to get empty vectors
      std::vector<Point> new_positions;
      std::vector<Real>   weight;
      new_positions.resize(mesh.n_nodes());
      weight.resize(mesh.n_nodes());

      {
        /*
         * Loop over the elements to calculate new node positions
         */
        MeshBase::element_iterator       el  = mesh.level_elements_begin(refinement_level);
        const MeshBase::element_iterator end = mesh.level_elements_end(refinement_level);

        for (; el != end; ++el)
        {
          /*
           * Constant handle for the element
           */
          const Elem* elem = *el;

          /*
           * We relax all nodes on level 0 first
           * If the element is refined (level > 0), we interpolate the
           * parents nodes with help of the embedding matrix
           */
          if (refinement_level == 0)
          {
            for (unsigned int s=0; s<elem->n_neighbors(); s++)
            {
              /*
               * Only operate on sides which are on the
               * boundary or for which the current element's
               * id is greater than its neighbor's.
               * Sides get only built once.
               */
              if ((elem->neighbor(s) != NULL) &&
                  (elem->id() > elem->neighbor(s)->id()) )
              {
                AutoPtr<Elem> side(elem->build_side(s));

                Node* node0 = side->get_node(0);
                Node* node1 = side->get_node(1);

                Real node_weight = 1.;
                // calculate the weight of the nodes
                if (power > 0)
                {
                  Point diff = (*node0)-(*node1);
                  node_weight = std::pow( diff.size(), power );
                }

                const unsigned int id0 = node0->id(), id1 = node1->id();
                new_positions[id0].add_scaled( *node1, node_weight );
                new_positions[id1].add_scaled( *node0, node_weight );
                weight[id0] += node_weight;
                weight[id1] += node_weight;
              }
            } // element neighbor loop
          }
#ifdef ENABLE_AMR
          else   // refinement_level > 0
          {
            /*
             * Find the positions of the hanging nodes of refined elements.
             * We do this by calculating their position based on the parent
             * (one level less refined) element, and the embedding matrix
             */

            const Elem* parent = elem->parent();

            /*
             * find out which child I am
             */
            for (unsigned int c=0; c < parent->n_children(); c++)
            {
              if (parent->child(c) == elem)
              {
                /*
                 *loop over the childs (that is, the current elements) nodes
                 */
                for (unsigned int nc=0; nc < elem->n_nodes(); nc++)
                {
                  /*
                   * the new position of the node
                   */
                  Point point;
                  for (unsigned int n=0; n<parent->n_nodes(); n++)
                  {
                    /*
                     * The value from the embedding matrix
                     */
                    const float em_val = parent->embedding_matrix(c,nc,n);

                    if (em_val != 0.)
                      point.add_scaled (parent->point(n), em_val);
                  }

                  const unsigned int id = elem->get_node(nc)->id();
                  new_positions[id] = point;
                  weight[id] = 1.;
                }

              } // if parent->child == elem
            } // for parent->n_children
          } // if element refinement_level
#endif // #ifdef ENABLE_AMR

        } // element loop

        /*
         * finally reposition the vertex nodes
         */
        for (unsigned int nid=0; nid<mesh.n_nodes(); ++nid)
          if (!on_boundary[nid] && weight[nid] > 0.)
            mesh.node(nid) = new_positions[nid]/weight[nid];
      }

      {
        /*
         * Now handle the additional second_order nodes by calculating
         * their position based on the vertex postitions
         * we do a second loop over the level elements
         */
        MeshBase::element_iterator       el  = mesh.level_elements_begin(refinement_level);
        const MeshBase::element_iterator end = mesh.level_elements_end(refinement_level);

        for (; el != end; ++el)
        {
          /*
           * Constant handle for the element
           */
          const Elem* elem = *el;
          const unsigned int son_begin = elem->n_vertices();
          const unsigned int son_end   = elem->n_nodes();
          for (unsigned int n=son_begin; n<son_end; n++)
          {
            const unsigned int n_adjacent_vertices =
              elem->n_second_order_adjacent_vertices(n);

            Point point;
            for (unsigned int v=0; v<n_adjacent_vertices; v++)
              point.add(elem->point( elem->second_order_adjacent_vertex(n,v) ));

            const unsigned int id = elem->get_node(n)->id();
            mesh.node(id) = point/n_adjacent_vertices;
          }
        }
      }

    } // refinement_level loop

  } // end iteration
}



#ifdef ENABLE_AMR
void MeshTools::Modification::flatten(MeshBase& mesh)
{
  // Algorithm:
  // .) For each active element in the mesh: construct a
  //    copy which is the same in every way *except* it is
  //    a level 0 element.  Store the pointers to these in
  //    a separate vector. Save any boundary information as well.
  //    Delete the active element from the mesh.
  // .) Loop over all (remaining) elements in the mesh, delete them.
  // .) Add the level-0 copies back to the mesh

  // Temporary storage for old and new element pointers
  std::vector<Elem*> old_elements, new_elements;

  // BoundaryInfo Storage for element ids, sides, and BC ids
  std::vector<Elem*>              saved_boundary_elements;
  std::vector<short int>          saved_bc_ids;
  std::vector<unsigned short int> saved_bc_sides;

  // Reserve a reasonable amt. of space for each
  new_elements.reserve(mesh.n_active_elem());
  saved_boundary_elements.reserve(mesh.boundary_info->n_boundary_conds());
  saved_bc_ids.reserve(mesh.boundary_info->n_boundary_conds());
  saved_bc_sides.reserve(mesh.boundary_info->n_boundary_conds());

  {
    MeshBase::element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::element_iterator end = mesh.active_elements_end();

    for (; it != end; ++it)
    {
      Elem* elem = *it;

      // Make a new element of the same type
      Elem* copy = Elem::build(elem->type()).release();

      // Set node pointers (they point to nodes in the original mesh)
      // The copy's element ID will be set upon adding it to the mesh
      for(unsigned int n=0; n<elem->n_nodes(); n++)
        copy->set_node(n) = elem->get_node(n);

      // set the subdomain id
      copy->subdomain_id () = elem->subdomain_id ();

      // This element could have boundary info as well.  We need
      // to save the (elem, side, bc_id) triples
      if( elem->on_boundary () )
        for (unsigned int s=0; s<elem->n_sides(); s++)
        {
          short int bc_id = mesh.boundary_info->boundary_id (elem,s);
          // make sure it is a boundary face
          if (bc_id != BoundaryInfo::invalid_id && elem->neighbor(s) == NULL )
          {
	    saved_boundary_elements.push_back(copy);
            saved_bc_ids.push_back(bc_id);
            saved_bc_sides.push_back(s);
          }
        }

      if( elem->on_interface () )
        for (unsigned int s=0; s<elem->n_sides(); s++)
        {
          short int bc_id = mesh.boundary_info->boundary_id (elem,s);
          // make sure it is an interface face
          if (bc_id != BoundaryInfo::invalid_id )
            if(elem->neighbor(s) && elem->neighbor(s)->subdomain_id () != elem->subdomain_id () )
	    {
	      saved_boundary_elements.push_back(copy);
              saved_bc_ids.push_back(bc_id);
              saved_bc_sides.push_back(s);
            }
        }

      // We're done with this element
      old_elements.push_back(elem);

      // But save the copy
      new_elements.push_back(copy);
    }

    // delete old elements
    {
      for (std::vector<Elem*>::iterator it = old_elements.begin();
           it != old_elements.end();
           ++it)
        mesh.delete_elem(*it);
    }


    // Make sure we saved the same number of boundary conditions
    // in each vector.
    assert (saved_boundary_elements.size() == saved_bc_ids.size());
    assert (saved_bc_ids.size()            == saved_bc_sides.size());
  }


  // Loop again, delete any remaining elements
  {
    MeshBase::element_iterator       it  = mesh.elements_begin();
    const MeshBase::element_iterator end = mesh.elements_end();

    for (; it != end; ++it)
      mesh.delete_elem( *it );
  }


  // Add the copied (now level-0) elements back to the mesh
  {
    for (std::vector<Elem*>::iterator it = new_elements.begin();
         it != new_elements.end();
         ++it)
      mesh.add_elem(*it);
  }

  // Finally, also add back the saved boundary information
  for (unsigned int e=0; e<saved_boundary_elements.size(); ++e)
    mesh.boundary_info->add_side(saved_boundary_elements[e],
                                 saved_bc_sides[e],
                                 saved_bc_ids[e]);

  // Trim unused nodes and elements
  //mesh.prepare_for_use();
  mesh.renumber_nodes_and_elements ();
}
#endif // #ifdef ENABLE_AMR

