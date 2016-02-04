// $Id: serial_mesh.h,v 1.4 2008/05/17 07:25:13 gdiso Exp $

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



#ifndef __serial_mesh_h__
#define __serial_mesh_h__

// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "unstructured_mesh.h"

/**
 * The \p SerialMesh class is derived from the \p MeshBase class,
 * and currently represents the default Mesh implementation.
 * Most methods for this class are found in MeshBase, and most
 * implementation details are found in UnstructuredMesh.
*/

// ------------------------------------------------------------
// Mesh class definition
class SerialMesh : public UnstructuredMesh
{
 public:

  /**
   * Constructor.  Requires the dimension and optionally
   * a processor id.  Note that \p proc_id should always
   * be provided for multiprocessor applications.
   */
  SerialMesh (unsigned int d);

  /**
   * Copy-constructor.  This should be able to take a
   * serial or parallel mesh.
   */
  SerialMesh (const UnstructuredMesh& other_mesh);

  /**
   * Copy-constructor, possibly specialized for a
   * serial mesh.
   */
  SerialMesh (const SerialMesh& other_mesh);

  /**
   * Virtual copy-constructor, creates a copy of this mesh
   */
  virtual AutoPtr<MeshBase> clone () const
    { return AutoPtr<MeshBase>(new SerialMesh(*this)); }

  /**
   * Destructor.
   */
  virtual ~SerialMesh();

  /**
   * Clear all internal data.
   */
  virtual void clear();

  /**
   * @returns \p true if all elements and nodes of the mesh
   * exist on the current processor, \p false otherwise
   */
  virtual bool is_serial () const
  { return _is_serial; }

  /**
   * set the state of mesh
   */
  virtual void set_serial (bool flag)
  {_is_serial = flag; }

  /**
   * broadcast all elements and nodes of the mesh onto
   * other processor
   */
  virtual void broadcast (unsigned int root_id=0);

  /**
   * Gathers all elements and nodes of the mesh onto
   * root processor. mesh can be totally distributed
   */
  virtual void gather(unsigned int root_id=0);

  /**
   * Gathers all elements and nodes of the mesh onto
   * every processor. mesh can be totally distributed
   */
  virtual void allgather();

  /**
   * When supported, deletes all nonlocal elements of the mesh
   * except for "ghosts" which touch a local element, and deletes
   * all nodes which are not part of a local or ghost element
   */
  virtual void delete_remote_elements (bool volume_elem=true, bool surface_elem=true);

  /**
   * pack all the mesh node location (x, y, z) one by one into an real array
   * with size 3*n_nodes(), should be executed in parallel
   */
  virtual void pack_nodes(std::vector<Real> &) const;

  /**
   * pack all the mesh edge info, should be executed in parallel
   */
  virtual void pack_egeds(std::vector< std::pair<unsigned int, unsigned int> > &) const;

  /**
   * pack all the mesh elem info, should be executed in parallel
   */
  virtual void pack_elems(std::vector<int> &) const;

  /**
   * pack all the mesh boundary face info, should be executed in parallel
   */
  virtual void pack_boundary_faces (std::vector<unsigned int> &,
                                    std::vector<unsigned short int> &,
                                    std::vector<short int> &) const;

  /**
   * pack all the mesh boundary edge info, should be executed in parallel
   */
  virtual void pack_boundary_egeds(std::vector< std::pair<unsigned int, unsigned int> > &) const;

  /**
   * pack all the mesh boundary node info, should be executed in parallel
   */
  virtual void pack_boundary_nodes (std::vector<unsigned int> &,
                                    std::vector<short int> &) const;

  /**
   * @Return the clone of \f$ i^{th} \f$ node, even it is not local
   * must be executed in parallel
   * use AutoPtr to prevent memory lost.
   */
  virtual AutoPtr<Node> node_clone (const unsigned int i) const;

  /**
   * @Return the clone of \f$ i^{th} \f$ element, even it is not local
   * NOTE this function also clone elem nodes
   * must be executed in parallel
   * use AutoPtr to prevent memory lost.
   */
  virtual AutoPtr<Elem> elem_clone (const unsigned int i) const;

  /**
   * Remove NULL elements from arrays
   */
  virtual void renumber_nodes_and_elements ();


  virtual unsigned int n_nodes () const { return _nodes.size(); }
  virtual unsigned int max_node_id () const { return _nodes.size(); }
  virtual void reserve_nodes (const unsigned int nn) { _nodes.reserve (nn); }
  virtual unsigned int n_elem ()  const { return _elements.size(); }
  virtual unsigned int max_elem_id ()  const { return _elements.size(); }
  virtual void reserve_elem (const unsigned int ne) { _elements.reserve (ne); }

  /**
   * For meshes that don't store points/elems, these functions may be an issue!
   */
  virtual const Point& point (const unsigned int i) const ;
  virtual const Node&  node  (const unsigned int i) const ;
  virtual Node& node (const unsigned int i) ;
  virtual const Node* node_ptr (const unsigned int i) const ;
  virtual Node* & node_ptr (const unsigned int i) ;
  virtual Elem* elem (const unsigned int i) const ;

  /**
   * functions for adding /deleting nodes elements.
   */
  virtual Node* add_point (const Point& p,
               const unsigned int id =
                 DofObject::invalid_id,
               const unsigned int proc_id =
                 DofObject::invalid_processor_id);
  virtual Node* add_node (Node* n) ;
  virtual void delete_node (Node* n) ;
  virtual Elem* add_elem (Elem* e) ;
  virtual Elem* add_elem (Elem* e, unsigned int id) ;
  virtual Elem* insert_elem (Elem* e) ;
  virtual void delete_elem (Elem* e) ;


  /**
   * functions for reordering elems
   */
  virtual bool reorder_elems(std::string &err);

  /**
   * functions for reordering nodes
   */
  virtual bool reorder_nodes (std::string &err);

  /**
   * generate all boundary elem-side pair with given boundary id
   */
  virtual void generate_boundary_info(short int id);

  /**
   * the subdomain interconnect graph in CSR format
   */
  virtual void subdomain_graph(std::vector<std::vector<unsigned int> >&) const;


  /**
   * approx memory usage
   */
  virtual size_t memory_size() const;

public:
  /**
   * Elem iterator accessor functions.
   */
  element_iterator elements_begin ();
  element_iterator elements_end   ();

  element_iterator active_elements_begin ();
  element_iterator active_elements_end   ();

  element_iterator subactive_elements_begin ();
  element_iterator subactive_elements_end   ();

  element_iterator not_active_elements_begin ();
  element_iterator not_active_elements_end   ();

  element_iterator not_subactive_elements_begin ();
  element_iterator not_subactive_elements_end   ();

  element_iterator this_pid_elements_begin ();
  element_iterator this_pid_elements_end   ();

  element_iterator not_this_pid_elements_begin ();
  element_iterator not_this_pid_elements_end   ();

  element_iterator active_this_pid_elements_begin ();
  element_iterator active_this_pid_elements_end   ();

  element_iterator active_not_this_pid_elements_begin ();
  element_iterator active_not_this_pid_elements_end   ();

  element_iterator local_elements_begin             ();
  element_iterator local_elements_end               ();

  element_iterator ghost_elements_begin             ();
  element_iterator ghost_elements_end               ();

  element_iterator level_elements_begin (const unsigned int level);
  element_iterator level_elements_end   (const unsigned int level);

  element_iterator not_level_elements_begin (const unsigned int level);
  element_iterator not_level_elements_end   (const unsigned int level);

  element_iterator this_pid_level_elements_begin (const unsigned int level);
  element_iterator this_pid_level_elements_end   (const unsigned int level);

  element_iterator this_pid_not_level_elements_begin (const unsigned int level);
  element_iterator this_pid_not_level_elements_end   (const unsigned int level);

  element_iterator pid_elements_begin (const unsigned int proc_id);
  element_iterator pid_elements_end   (const unsigned int proc_id);

  element_iterator subdomain_elements_begin (const unsigned int subdomain_id);
  element_iterator subdomain_elements_end   (const unsigned int subdomain_id);

  element_iterator type_elements_begin (const ElemType type);
  element_iterator type_elements_end   (const ElemType type);

  element_iterator active_type_elements_begin (const ElemType type);
  element_iterator active_type_elements_end   (const ElemType type);

  element_iterator active_pid_elements_begin (const unsigned int proc_id);
  element_iterator active_pid_elements_end   (const unsigned int proc_id);

  element_iterator unpartitioned_elements_begin ();
  element_iterator unpartitioned_elements_end ();


  /**
   * const Elem iterator accessor functions.
   */
  const_element_iterator elements_begin() const;
  const_element_iterator elements_end()   const;

  const_element_iterator active_elements_begin() const;
  const_element_iterator active_elements_end()   const;

  const_element_iterator subactive_elements_begin() const;
  const_element_iterator subactive_elements_end()   const;

  const_element_iterator not_active_elements_begin() const;
  const_element_iterator not_active_elements_end()   const;

  const_element_iterator not_subactive_elements_begin() const;
  const_element_iterator not_subactive_elements_end()   const;

  const_element_iterator this_pid_elements_begin () const;
  const_element_iterator this_pid_elements_end   () const;

  const_element_iterator not_this_pid_elements_begin () const;
  const_element_iterator not_this_pid_elements_end   () const;

  const_element_iterator active_this_pid_elements_begin () const;
  const_element_iterator active_this_pid_elements_end   () const;

  const_element_iterator active_not_this_pid_elements_begin () const;
  const_element_iterator active_not_this_pid_elements_end   () const;

  const_element_iterator local_elements_begin             () const;
  const_element_iterator local_elements_end               () const;

  const_element_iterator ghost_elements_begin             () const;
  const_element_iterator ghost_elements_end               () const;

  const_element_iterator level_elements_begin (const unsigned int level) const;
  const_element_iterator level_elements_end   (const unsigned int level) const;

  const_element_iterator not_level_elements_begin (const unsigned int level) const;
  const_element_iterator not_level_elements_end   (const unsigned int level) const;

  const_element_iterator this_pid_level_elements_begin (const unsigned int level) const;
  const_element_iterator this_pid_level_elements_end   (const unsigned int level) const;

  const_element_iterator this_pid_not_level_elements_begin (const unsigned int level) const;
  const_element_iterator this_pid_not_level_elements_end   (const unsigned int level) const;

  const_element_iterator pid_elements_begin (const unsigned int proc_id) const;
  const_element_iterator pid_elements_end   (const unsigned int proc_id) const;

  const_element_iterator subdomain_elements_begin (const unsigned int subdomain_id) const;
  const_element_iterator subdomain_elements_end   (const unsigned int subdomain_id) const;

  const_element_iterator type_elements_begin (const ElemType type) const;
  const_element_iterator type_elements_end   (const ElemType type) const;

  const_element_iterator active_type_elements_begin (const ElemType type) const;
  const_element_iterator active_type_elements_end   (const ElemType type) const;

  const_element_iterator active_pid_elements_begin (const unsigned int proc_id) const;
  const_element_iterator active_pid_elements_end   (const unsigned int proc_id) const;

  const_element_iterator unpartitioned_elements_begin () const;
  const_element_iterator unpartitioned_elements_end () const;







  /**
   * non-const Node iterator accessor functions.
   */
  node_iterator nodes_begin();
  node_iterator nodes_end();

  node_iterator active_nodes_begin();
  node_iterator active_nodes_end();

  node_iterator this_pid_nodes_begin  ();
  node_iterator this_pid_nodes_end    ();

  node_iterator pid_nodes_begin (const unsigned int proc_id);
  node_iterator pid_nodes_end   (const unsigned int proc_id);

  node_iterator local_nodes_begin             () ;
  node_iterator local_nodes_end               () ;

  /**
   * const Node iterator accessor functions.
   */
  const_node_iterator nodes_begin() const;
  const_node_iterator nodes_end()   const;

  const_node_iterator active_nodes_begin() const;
  const_node_iterator active_nodes_end()   const;

  const_node_iterator this_pid_nodes_begin  () const;
  const_node_iterator this_pid_nodes_end    () const;

  const_node_iterator pid_nodes_begin (const unsigned int proc_id) const;
  const_node_iterator pid_nodes_end   (const unsigned int proc_id) const;

  const_node_iterator local_nodes_begin             () const;
  const_node_iterator local_nodes_end               () const;

protected:
  /**
   * The verices (spatial coordinates) of the mesh.
   */
  std::vector<Node*> _nodes;

  /**
   * The elements in the mesh.
   */
  std::vector<Elem*> _elements;

  /**
   * A boolean remembering whether we're serialized or not
   */
  bool _is_serial;

private:

  /**
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Elem*>.
   */
  typedef std::vector<Elem*>::iterator             elem_iterator_imp;
  typedef std::vector<Elem*>::const_iterator const_elem_iterator_imp;

  /**
   * Typedefs for the container implementation.  In this case,
   * it's just a std::vector<Node*>.
   */
  typedef std::vector<Node*>::iterator             node_iterator_imp;
  typedef std::vector<Node*>::const_iterator const_node_iterator_imp;

private:

  void _unpack_mesh (const std::vector<Real> &pts, const std::vector<int> &conn);
  void _unpack_bc_faces (const std::vector<unsigned int> &el_id,
                         const std::vector<unsigned short int> &side_id,
                         const std::vector<short int> &bc_id);
  void _unpack_bc_nodes (const std::vector<unsigned int> &node_id,
                         const std::vector<short int> &bc_id);

};





#endif
