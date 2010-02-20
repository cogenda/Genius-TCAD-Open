// $Id: unstructured_mesh.cc,v 1.14 2008/07/10 04:03:46 gdiso Exp $

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
#include <fstream>


// Local includes
#include "boundary_info.h"
#include "unstructured_mesh.h"
#include "mesh_communication.h"
#include "mesh_tools.h" // For n_levels
#include "perf_log.h"
#include "elem.h"

#if   defined(HAVE_HASH_MAP)
# include <hash_map>
#elif defined(HAVE_EXT_HASH_MAP)
# include <ext/hash_map>
#endif


// ------------------------------------------------------------
// UnstructuredMesh class member functions
UnstructuredMesh::UnstructuredMesh (unsigned int d) :
    MeshBase (d)
{
}



void UnstructuredMesh::copy_nodes_and_elements
(const UnstructuredMesh& other_mesh)
{
  // We're assuming our subclass data needs no copy
  genius_assert(_n_sbd == other_mesh._n_sbd);
  genius_assert(_n_parts == other_mesh._n_parts);
  genius_assert(_dim == other_mesh._dim);
  genius_assert(_is_prepared == other_mesh._is_prepared);

  //Copy in Nodes
  {
    //Preallocate Memory if necessary
    this->reserve_nodes(other_mesh.n_nodes());

    const_node_iterator it = other_mesh.nodes_begin();
    const_node_iterator end = other_mesh.nodes_end();

    for (; it != end; ++it)
      this->add_point(*(*it)); //Add new nodes in old node Point locations
  }

  //Copy in Elements
  {
    //Preallocate Memory if necessary
    this->reserve_elem(other_mesh.n_elem());

    // Loop over the elements
    MeshBase::const_element_iterator it = other_mesh.elements_begin();
    const MeshBase::const_element_iterator end = other_mesh.elements_end();

    // FIXME: Where do we set element IDs??
    for (; it != end; ++it)
    {
      //Look at the old element
      Elem *old = *it;
      //Build a new element
      Elem *newparent = old->parent() ?
                        this->elem(old->parent()->id()) : NULL;
      AutoPtr<Elem> ap = Elem::build(old->type(), newparent);
      Elem * elem = ap.release();

#ifdef ENABLE_AMR
      //Create the parent's child pointers if necessary
      if (newparent)
      {
        // Make sure we have space for those child pointers
        newparent->add_child(elem);

        // We'd better be adding these in the correct order
        assert (newparent->which_child_am_i(elem) ==
                old->parent()->which_child_am_i(old));
      }

      // Copy the refinement flags
      elem->set_refinement_flag(old->refinement_flag());
      elem->set_p_refinement_flag(old->p_refinement_flag());
#endif // #ifdef ENABLE_AMR

      //Assign all the nodes
      for(unsigned int i=0;i<elem->n_nodes();i++)
        elem->set_node(i) = &this->node(old->node(i));

      //Hold onto it
      this->add_elem(elem);
    }
  }

  //Finally prepare the Mesh for use
  this->prepare_for_use();
}



UnstructuredMesh::~UnstructuredMesh ()
{
  //  this->clear ();  // Nothing to clear at this level


}





void UnstructuredMesh::find_neighbors()
{
  genius_assert(this->n_nodes() != 0);
  genius_assert(this->n_elem()  != 0);

  START_LOG("find_neighbors()", "Mesh");


  //TODO:[BSK] This should be removed later?!
  const element_iterator el_end = this->elements_end();
  for (element_iterator el = this->elements_begin(); el != el_end; ++el)
  {
    Elem* elem = *el;
    for (unsigned int s=0; s<elem->n_neighbors(); s++)
      elem->set_neighbor(s,NULL);
  }

  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors
  {
    // data structures -- Use the hash_multimap if available
    typedef unsigned int                    key_type;
    typedef std::pair<Elem*, unsigned char> val_type;
    typedef std::pair<key_type, val_type>   key_val_pair;

#if   defined(HAVE_HASH_MAP)
    typedef std::hash_multimap<key_type, val_type> map_type;
#elif defined(HAVE_EXT_HASH_MAP)
# if    (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
    typedef std::hash_multimap<key_type, val_type> map_type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
    typedef __gnu_cxx::hash_multimap<key_type, val_type> map_type;
# else
    // XLC and who knows what other compilers get here.
    // Try the most standard thing we can:
    typedef std::multimap<key_type, val_type>  map_type;
# endif
#else
    typedef std::multimap<key_type, val_type>  map_type;
#endif

    // A map from side keys to corresponding elements & side numbers
    map_type side_to_elem_map;



    for (element_iterator el = this->elements_begin(); el != el_end; ++el)
    {
      Elem* element = *el;

      for (unsigned int ms=0; ms<element->n_neighbors(); ms++)
      {
      next_side:

        if (element->neighbor(ms) == NULL)
        {
          // Get the key for the side of this element
          const unsigned int key = element->key(ms);

          // Look for elements that have an identical side key
          std::pair <map_type::iterator, map_type::iterator>
          bounds = side_to_elem_map.equal_range(key);

          // May be multiple keys, check all the possible
          // elements which _might_ be neighbors.
          if (bounds.first != bounds.second)
          {
            // Get the side for this element
            const AutoPtr<DofObject> my_side(element->side(ms));

            // Look at all the entries with an equivalent key
            while (bounds.first != bounds.second)
            {
              // Get the potential element
              Elem* neighbor = bounds.first->second.first;

              // Get the side for the neighboring element
              const unsigned int ns = bounds.first->second.second;
              const AutoPtr<DofObject> their_side(neighbor->side(ns));
              //assert (my_side.get() != NULL);
              //assert (their_side.get() != NULL);

              // If found a match wbounds.firsth my side
              //
              // We need a special case here for 1D, since parents
              // and children have an equal side (i.e. a node),
              // so need to check ns != ms in 1D as well
              if( (*my_side == *their_side) && ((_dim != 1) || (ns != ms)) )
              {
                // So share a side.  Is this a mixed pair
                // of subactive and active/ancestor
                // elements?
                // If not, then we're neighbors.
                // If so, then the subactive's neighbor is

                if (element->subactive() == neighbor->subactive())
                {
                  // an element is only subactive if it has
                  // been coarsened but not deleted
                  element->set_neighbor (ms,neighbor);
                  neighbor->set_neighbor(ns,element);
                }
                else if (element->subactive())
                {
                  element->set_neighbor(ms,neighbor);
                }
                else if (neighbor->subactive())
                {
                  neighbor->set_neighbor(ns,element);
                }
                side_to_elem_map.erase (bounds.first);

                // get out of this nested crap
                goto next_side;
              }

              ++bounds.first;
            }
          }

          // didn't find a match...
          // Build the map entry for this element
          key_val_pair kvp;

          kvp.first         = key;
          kvp.second.first  = element;
          kvp.second.second = ms;

          // use the lower bound as a hint for
          // where to put it.
#if defined(HAVE_HASH_MAP) || defined(HAVE_EXT_HASH_MAP)
          side_to_elem_map.insert (kvp);
#else
          side_to_elem_map.insert (bounds.first,kvp);
#endif

        }
      }
    }
  }



#ifdef ENABLE_AMR

  /**
   * Here we look at all of the child elements.
   * If a child element has a NULL neighbor it is
   * either because it is on the boundary or because
   * its neighbor is at a different level.  In the
   * latter case we must get the neighbor from the
   * parent.
   *
   * Furthermore, that neighbor better be active,
   * otherwise we missed a child somewhere.
   */
  element_iterator end = this->not_level_elements_end(0);
  for (element_iterator el = this->not_level_elements_begin(0);
       el != end; ++el)
  {
    Elem* elem = *el;

    assert (elem->parent() != NULL);
    for (unsigned int s=0; s < elem->n_neighbors(); s++)
      if (elem->neighbor(s) == NULL)
      {
        elem->set_neighbor(s, elem->parent()->neighbor(s));
      }
  }

#endif // AMR

  STOP_LOG("find_neighbors()", "Mesh");
}




// ------------------------------------------------------------
// UnstructuredMesh class member functions for mesh modification
void UnstructuredMesh::all_first_order ()
{
  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
    this->renumber_nodes_and_elements ();

  START_LOG("all_first_order()", "Mesh");

  /**
   * Loop over the high-ordered elements.
   * First make sure they _are_ indeed high-order, and then replace
   * them with an equivalent first-order element.
   */
  const_element_iterator endit = elements_end();
  for (const_element_iterator it = elements_begin();
       it != endit; ++it)
  {
    Elem* so_elem = *it;

    assert (so_elem != NULL);

    if ( Elem::first_order_equivalent_type(so_elem->type()) == so_elem->type() ) continue;

    /*
     * build the first-order equivalent, add to
     * the new_elements list.
     */
    Elem *newparent = so_elem->parent();
    Elem* lo_elem = Elem::build
                    (Elem::first_order_equivalent_type
                     (so_elem->type()), newparent).release();

#ifdef ENABLE_AMR
    /*
     * Add this element to it's parent if it has one
     */
    if (newparent)
      newparent->add_child(lo_elem);


    /*
     * Reset the parent links of any child elements
     */
    if (so_elem->has_children())
    {
      for (unsigned int c=0; c != so_elem->n_children(); ++c)
        so_elem->child(c)->set_parent(lo_elem);
    }

    /*
     * Copy as much data to the new element as makes sense
     */
    lo_elem->set_p_level(so_elem->p_level());
    lo_elem->set_refinement_flag(so_elem->refinement_flag());
    lo_elem->set_p_refinement_flag(so_elem->p_refinement_flag());
#endif

    assert (lo_elem->n_vertices() == so_elem->n_vertices());

    /*
     * By definition the vertices of the linear and
     * second order element are identically numbered.
     * transfer these.
     */
    for (unsigned int v=0; v < so_elem->n_vertices(); v++)
      lo_elem->set_node(v) = so_elem->get_node(v);

    /*
     * set the subdomain id
     */
    lo_elem->subdomain_id () = so_elem->subdomain_id ();


    /**
     * If the second order element had any boundary conditions they
     * should be transfered to the first-order element.  The old
     * boundary conditions will be removed from the BoundaryInfo
     * data structure by insert_elem.
     */
    assert (lo_elem->n_sides() == so_elem->n_sides());

    for (unsigned int s=0; s<so_elem->n_sides(); s++)
    {
      const short int boundary_id =
        this->boundary_info->boundary_id (so_elem, s);

      if (boundary_id != this->boundary_info->invalid_id)
        this->boundary_info->add_side (lo_elem, s, boundary_id);
    }

    /*
     * The new first-order element is ready.
     * Inserting it into the mesh will replace and delete
     * the second-order element.
     */
    lo_elem->set_id(so_elem->id());
    this->insert_elem(lo_elem);
  }

  STOP_LOG("all_first_order()", "Mesh");

  // delete or renumber nodes, etc
  this->prepare_for_use();
}



void UnstructuredMesh::all_second_order (const bool full_ordered)
{
  /*
   * when the mesh is not prepared,
   * at least renumber the nodes and
   * elements, so that the node ids
   * are correct
   */
  if (!this->_is_prepared)
    this->renumber_nodes_and_elements ();

  /*
   * If the mesh is empty or already second order
   * then we have nothing to do
   */
  if (!this->n_elem() ||
      (*(this->elements_begin()))->default_order() != FIRST)
    return;

  // does this work also in parallel?
  // assert (this->n_processors() == 1);

  START_LOG("all_second_order()", "Mesh");

  /*
   * this map helps in identifying second order
   * nodes.  Namely, a second-order node:
   * - edge node
   * - face node
   * - bubble node
   * is uniquely defined through a set of adjacent
   * vertices.  This set of adjacent vertices is
   * used to identify already added higher-order
   * nodes.  We are safe to use node id's since we
   * make sure that these are correctly numbered.
   */
  std::map<std::vector<unsigned int>, Node*> adj_vertices_to_so_nodes;

  /*
   * for speed-up of the \p add_point() method, we
   * can reserve memory.  Guess the number of additional
   * nodes for different dimensions
   */
  switch (this->mesh_dimension())
  {
  case 1:
    /*
     * in 1D, there can only be order-increase from Edge2
     * to Edge3.  Something like 1/2 of n_nodes() have
     * to be added
     */
    this->reserve_nodes(static_cast<unsigned int>(1.5*this->n_nodes()));
    break;

  case 2:
    /*
     * in 2D, either refine from Tri3 to Tri6 (double the nodes)
     * or from Quad4 to Quad8 (again, double) or Quad9 (2.25 that much)
     */
    this->reserve_nodes(static_cast<unsigned int>(2*this->n_nodes()));
    break;


  case 3:
    /*
     * in 3D, either refine from Tet4 to Tet10 (factor = 2.5) up to
     * Hex8 to Hex27 (something  > 3).  Since in 3D there _are_ already
     * quite some nodes, and since we do not want to overburden the memory by
     * a too conservative guess, use the lower bound
     */
    this->reserve_nodes(static_cast<unsigned int>(2.5*this->n_nodes()));
    break;

  default:
    // Hm?
    genius_error();
  }



  /*
   * form a vector that will hold the node id's of
   * the vertices that are adjacent to the son-th
   * second-order node.  Pull this outside of the
   * loop so that silly compilers don't repeatedly
   * create and destroy the vector.
   */
  std::vector<unsigned int> adjacent_vertices_ids;


  /**
   * Loop over the low-ordered elements in the _elements vector.
   * First make sure they _are_ indeed low-order, and then replace
   * them with an equivalent second-order element.  Don't
   * forget to delete the low-order element, or else it will leak!
   */
  const_element_iterator endit = elements_end();
  for (const_element_iterator it = elements_begin();
       it != endit; ++it)
  {
    // the linear-order element
    Elem* lo_elem = *it;

    assert (lo_elem != NULL);

    // make sure it is linear order
    if (lo_elem->default_order() != FIRST)
    {
      std::cerr << "ERROR: This is not a linear element: type="
      << lo_elem->type() << std::endl;
      genius_error();
    }

    // this does _not_ work for refined elements
    assert (lo_elem->level () == 0);

    /*
     * build the second-order equivalent, add to
     * the new_elements list.  Note that this here
     * is the only point where \p full_ordered
     * is necessary.  The remaining code works well
     * for either type of seconrd-order equivalent, e.g.
     * Hex20 or Hex27, as equivalents for Hex8
     */
    Elem* so_elem =
      Elem::build (Elem::second_order_equivalent_type(lo_elem->type(),
                   full_ordered) ).release();

    assert (lo_elem->n_vertices() == so_elem->n_vertices());


    /*
     * By definition the vertices of the linear and
     * second order element are identically numbered.
     * transfer these.
     */
    for (unsigned int v=0; v < lo_elem->n_vertices(); v++)
      so_elem->set_node(v) = lo_elem->get_node(v);

    /*
     * set the subdomain id
     */
    so_elem->subdomain_id () = lo_elem->subdomain_id ();

    /*
     * Now handle the additional mid-side nodes.  This
     * is simply handled through a map that remembers
     * the already-added nodes.  This map maps the global
     * ids of the vertices (that uniquely define this
     * higher-order node) to the new node.
     * Notation: son = second-order node
     */
    const unsigned int son_begin = so_elem->n_vertices();
    const unsigned int son_end   = so_elem->n_nodes();


    for (unsigned int son=son_begin; son<son_end; son++)
    {
      const unsigned int n_adjacent_vertices =
        so_elem->n_second_order_adjacent_vertices(son);

      adjacent_vertices_ids.resize(n_adjacent_vertices);

      for (unsigned int v=0; v<n_adjacent_vertices; v++)
        adjacent_vertices_ids[v] =
          so_elem->node( so_elem->second_order_adjacent_vertex(son,v) );

      /*
       * \p adjacent_vertices_ids is now in order of the current
       * side.  sort it, so that comparisons  with the
       * \p adjacent_vertices_ids created through other elements'
       * sides can match
       */
      std::sort(adjacent_vertices_ids.begin(),
                adjacent_vertices_ids.end());


      // does this set of vertices already has a mid-node added?
      std::pair<std::map<std::vector<unsigned int>, Node*>::iterator,
      std::map<std::vector<unsigned int>, Node*>::iterator>
      pos = adj_vertices_to_so_nodes.equal_range (adjacent_vertices_ids);

      // no, not added yet
      if (pos.first == pos.second)
      {
        /*
         * for this set of vertices, there is no
         * second_order node yet.  Add it.
         *
         * compute the location of the new node as
         * the average over the adjacent vertices.
         */
        Point new_location = this->point(adjacent_vertices_ids[0]);
        for (unsigned int v=1; v<n_adjacent_vertices; v++)
          new_location += this->point(adjacent_vertices_ids[v]);

        new_location /= static_cast<Real>(n_adjacent_vertices);

        // add the new point to the mesh
        Node* so_node = this->add_point (new_location);

        /*
         * insert the new node with its defining vertex
         * set into the map, and relocate pos to this
         * new entry, so that the so_elem can use
         * \p pos for inserting the node
         */
        adj_vertices_to_so_nodes.insert(pos.first,
                                        std::make_pair(adjacent_vertices_ids,
                                                       so_node));

        so_elem->set_node(son) = so_node;
      }
      // yes, already added.
      else
      {
        assert (pos.first->second != NULL);

        so_elem->set_node(son) = pos.first->second;
      }
    }


    /**
     * If the linear element had any boundary conditions they
     * should be transfered to the second-order element.  The old
     * boundary conditions will be removed from the BoundaryInfo
     * data structure by insert_elem.
     */
    assert (lo_elem->n_sides() == so_elem->n_sides());

    for (unsigned int s=0; s<lo_elem->n_sides(); s++)
    {
      const short int boundary_id =
        this->boundary_info->boundary_id (lo_elem, s);

      if (boundary_id != this->boundary_info->invalid_id)
        this->boundary_info->add_side (so_elem, s, boundary_id);
    }

    /*
     * The new second-order element is ready.
     * Inserting it into the mesh will replace and delete
     * the first-order element.
     */
    so_elem->set_id(lo_elem->id());
    this->insert_elem(so_elem);
  }


  // we can clear the map
  adj_vertices_to_so_nodes.clear();


  STOP_LOG("all_second_order()", "Mesh");

  // renumber nodes, elements etc
  this->prepare_for_use();
}



bool UnstructuredMesh::all_fvm_elem ()
{

  // make sure all the element is first order (FEM element)
  // note, mesh is prepared here!

  this->all_first_order ();


  // here we convert all the active FEM element to FVM element, maybe only element belongs to local
  // procesor needs to be converted.
  const_element_iterator endit = active_elements_end();
  for (const_element_iterator it = active_elements_begin();  it != endit; ++it )
  {
    Elem* fem_elem = *it;

    assert (fem_elem != NULL);

    // skip elements not on_local
    if( !fem_elem->on_local() ) continue;

    // can this element be used in FVM?
    if ( fem_elem->fvm_compatible_test() == false ) return false;

    /*
     * build the FVM compatible element, add to
     * the new_elements list.
     */
    Elem *newparent = fem_elem->parent();
    Elem *fvm_elem = Elem::build (Elem::fvm_compatible_type(fem_elem->type()), newparent).release();

#ifdef ENABLE_AMR

    /*
     * replace FEM element with FVM element to it's parent if it has one
     */
    if (newparent)
    {
      newparent->delete_child(fem_elem);
      newparent->add_child(fvm_elem);
    }

    /*
     * Reset the parent links of any child elements
     */

    if (fem_elem->has_children())
    {
      for (unsigned int c=0; c != fem_elem->n_children(); ++c)
        fem_elem->child(c)->set_parent(fvm_elem);
    }

    /*
     * Copy as much data to the new element as makes sense
     */
    fvm_elem->set_p_level(fem_elem->p_level());
    fvm_elem->set_refinement_flag(fem_elem->refinement_flag());
    fvm_elem->set_p_refinement_flag(fem_elem->p_refinement_flag());
#endif
    assert (fvm_elem->n_vertices() == fem_elem->n_vertices());

    /*
     * The FVM element and first order FEM elem has the same node
     */
    for (unsigned int v=0; v < fem_elem->n_vertices(); v++)
      fvm_elem->set_node(v) = fem_elem->get_node(v);

    /*
     * build cell's geometry information for FVM usage
     */
    fvm_elem->prepare_for_fvm();

    /*
     * set the subdomain id
     */
    fvm_elem->subdomain_id () = fem_elem->subdomain_id ();

    /*
     * set processor_id and on_local information
     */
    fvm_elem->processor_id () = fem_elem->processor_id ();
    fvm_elem->on_local ()     = fem_elem->on_local ();

    /**
     * If the FEM element had any boundary conditions they
     * should be transfered to the FVM element.  The old
     * boundary conditions will be removed from the BoundaryInfo
     * data structure by later this->insert_elem(fvm_elem) call.
     */
    assert (fvm_elem->n_sides() == fem_elem->n_sides());


    for (unsigned int s=0; s<fem_elem->n_sides(); s++)
    {
      // only search fem_elem itself in the boundary_info structure,
      // send false to function boundary_id() avoid search the top parent of fem_elem in the boundary_info structure!
      const short int boundary_id = this->boundary_info->boundary_id (fem_elem, s, false);

      if (boundary_id != BoundaryInfo::invalid_id)
      {
        this->boundary_info->add_side ( fvm_elem, s, boundary_id );
      }
    }


    /*
     * The new FVM element is ready.
     * Inserting it into the mesh will replace and delete
     * the first-order FEM element.
     */
    fvm_elem->set_id(fem_elem->id());

    // delete old fem element,
    this->insert_elem(fvm_elem);

  }

  this->boundary_info->rebuild_ids();

  // rebuild neighbor information
  this->find_neighbors();

  return true;
}



void UnstructuredMesh::create_pid_mesh(UnstructuredMesh& pid_mesh,
                                       const unsigned int pid) const
{

  // Issue a warning if the number the number of processors
  // currently available is less that that requested for
  // partitioning.  This is not necessarily an error since
  // you may run on one processor and still partition the
  // mesh into several partitions.
#ifdef DEBUG
  if (this->n_processors() < pid)
  {
    std::cout << "WARNING:  You are creating a "
    << "mesh for a processor id (="
    << pid
    << ") greater than "
    << "the number of processors available for "
    << "the calculation. (="
    << Genius::n_processors()
    << ")."
    << std::endl;
  }
#endif

  // Create iterators to loop over the list of elements
  //   const_active_pid_elem_iterator       it(this->elements_begin(),   pid);
  //   const const_active_pid_elem_iterator it_end(this->elements_end(), pid);

  const_element_iterator       it     = this->active_pid_elements_begin(pid);
  const const_element_iterator it_end = this->active_pid_elements_end(pid);

  this->create_submesh (pid_mesh, it, it_end);
}







void UnstructuredMesh::create_submesh (UnstructuredMesh& new_mesh,
                                       const_element_iterator& it,
                                       const const_element_iterator& it_end) const
{
  // Just in case the subdomain_mesh already has some information
  // in it, get rid of it.
  new_mesh.clear();

  // Fail if (*this == new_mesh), we cannot create a submesh inside ourself!
  // This may happen if the user accidently passes the original mesh into
  // this function!  We will check this by making sure we did not just
  // clear ourself.
  assert (this->n_nodes() != 0);
  assert (this->n_elem()  != 0);

  // How the nodes on this mesh will be renumbered to nodes
  // on the new_mesh.
  std::vector<unsigned int> new_node_numbers (this->n_nodes());

  std::fill (new_node_numbers.begin(),
             new_node_numbers.end(),
             invalid_uint);



  // the number of nodes on the new mesh, will be incremented
  unsigned int n_new_nodes = 0;
  unsigned int n_new_elem  = 0;

  for (; it != it_end; ++it)
  {
    // increment the new element counter
    n_new_elem++;

    const Elem* old_elem = *it;

    // Add an equivalent element type to the new_mesh
    Elem* new_elem =
      new_mesh.add_elem (Elem::build(old_elem->type()).release());

    assert (new_elem != NULL);

    // Loop over the nodes on this element.
    for (unsigned int n=0; n<old_elem->n_nodes(); n++)
    {
      assert (old_elem->node(n) < new_node_numbers.size());

      if (new_node_numbers[old_elem->node(n)] == invalid_uint)
      {
        new_node_numbers[old_elem->node(n)] = n_new_nodes;

        // Add this node to the new mesh
        new_mesh.add_point (old_elem->point(n));

        // Increment the new node counter
        n_new_nodes++;
      }

      // Define this element's connectivity on the new mesh
      assert (new_node_numbers[old_elem->node(n)] < new_mesh.n_nodes());

      new_elem->set_node(n) = new_mesh.node_ptr (new_node_numbers[old_elem->node(n)]);
    }

    // Maybe add boundary conditions for this element
    for (unsigned int s=0; s<old_elem->n_sides(); s++)
      if (old_elem->neighbor(s) == NULL)
        if (this->boundary_info->boundary_id (old_elem, s) !=
            this->boundary_info->invalid_id)
          new_mesh.boundary_info->add_side (new_elem,
                                            s,
                                            this->boundary_info->boundary_id (old_elem, s));
  } // end loop over elements


  // Prepare the new_mesh for use
  new_mesh.prepare_for_use();

}



#ifdef ENABLE_AMR
bool UnstructuredMesh::contract ()
{
  START_LOG ("contract()", "Mesh");

  // Flag indicating if this call actually changes the mesh
  bool mesh_changed = false;

  element_iterator in        = elements_begin();
  element_iterator out       = elements_begin();
  const element_iterator end = elements_end();

#ifdef DEBUG
  for ( ; in != end; ++in)
    if (*in != NULL)
    {
      Elem* elem = *in;
      genius_assert(elem->active() || elem->subactive() || elem->ancestor());
    }
  in = elements_begin();
#endif

  // Loop over the elements.
  for ( ; in != end; ++in)
    if (*in != NULL)
    {
      Elem* elem = *in;

      // Delete all the subactive ones
      if (elem->subactive())
      {
        // Huh?  no level-0 element should be subactive
        assert (elem->level() != 0);

        // Make sure we dealt with parents first
        if (elem->parent()->has_children())
        {
          std::cerr << "Element being deleted is still a child." << std::endl;
        }

        // Delete the element
        // This just sets a pointer to NULL, and doesn't
        // invalidate any iterators
        this->delete_elem(elem);

        // the mesh has certainly changed
        mesh_changed = true;
      }
      else
      {
        // Compress all the active ones
        if (elem->active())
          elem->contract();
        else
          assert (elem->ancestor());
      }
    }

  // Strip any newly-created NULL voids out of the element array
  this->renumber_nodes_and_elements();

  STOP_LOG ("contract()", "Mesh");

  return mesh_changed;
}
#endif // #ifdef ENABLE_AMR



