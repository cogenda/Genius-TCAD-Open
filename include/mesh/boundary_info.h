// $Id: boundary_info.h,v 1.10 2008/06/20 05:59:01 gdiso Exp $

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



#ifndef __boundary_info_h__
#define __boundary_info_h__

// C++ includes
#include <vector>
#include <map>
#include <set>

// Local includes
#include "genius_common.h"


// Forward declarations
class Elem;
class Node;
class MeshBase;
class BoundaryMesh;



/**
 * The \p BoundaryInfo class contains information relevant
 * to boundary conditions: it does not hold actual boundary
 * condition data , but can mark
 * element faces and nodes with ids useful for identifying the
 * type of boundary condtion.  It can also build a mesh that
 * just includes boundary elements/faces.
 *
 * TODO[JWP]: Generalize this to work with MeshBase again.
 */

//------------------------------------------------------
// BoundaryInfo class definition
class BoundaryInfo
{
protected:
  friend class MeshBase;

  /**
   * Constructor.  Takes a reference to the mesh.
   * The BoundaryInfo class is only used internally
   * by the Mesh class.  A user should never instantiate
   * this class.  Therefore the constructor is protected.
   */
  BoundaryInfo (const MeshBase& m);

public:

  /**
   * Destructor.  Not much to do.
   */
  ~BoundaryInfo ();

  /**
   * Clears the underlying data structures.
   * Returns the object to a pristine state
   * with no data stored.
   */
  void clear ();

  /**
   * Close the data structures and prepare for use.
   * Synchronizes the \p boundary_mesh
   * data structures with the \p mesh data structures.
   * Allows the \p boundary_mesh to be used like any other mesh.
   * Before this is called the \p boundary_mesh data structure is
   * empty.
   */
  void sync (BoundaryMesh& boundary_mesh);


  /**
   * Create boundary_mesh with special boundary id.
   * Allows the \p boundary_mesh to be used like any other mesh.
   * Before this is called the \p boundary_mesh data structure is
   * empty.
   */
  //void create_boundary_mesh (BoundaryMesh& boundary_mesh, const short int id);


  /**
   * Add \p Node \p node with boundary id \p id to the boundary
   * information data structures.
   */
  void add_node (const Node* node,
         const short int id);

  /**
   * Add node number \p node with boundary id \p id to the boundary
   * information data structures.
   */
  void add_node (const unsigned int node,
         const short int id);

  /**
   * Add side \p side of element number \p elem with boundary id \p id
   * to the boundary information data structure.
   */
  void add_side (const unsigned int elem,
         const unsigned short int side,
         const short int id);

  /**
   * Add side \p side of element \p elem with boundary id \p id
   * to the boundary information data structure.
   */
  void add_side (const Elem* elem,
         const unsigned short int side,
         const short int id);


  /**
   * build _boundary_node_id map with a priority order of boundary id
   * from boundary side. i.e. one node lies on the interface of two sides with
   * different boundary id, the id of this node will be determined by an priority order.
   */
  void build_node_ids_from_priority_order(const std::map<short int, unsigned int> & order);


  /**
   * Removes the boundary conditions associated with node \p node,
   * if any exist.
   */
  void remove (const Node* node);

  /**
   * Removes the boundary conditions associated with element \p elem,
   * if any exist.
   */
  void remove (const Elem* elem);

  /**
   * Removes the boundary conditions associated with side of element \p elem,
   * if any exist.
   */
  void remove (const Elem* elem, unsigned short int side);

  /**
   * Returns the number of user-specified boundary ids.
   */
  unsigned int n_boundary_ids () const { return _boundary_ids.size(); }

  /**
   * Returns the boundary id associated with \p Node \p node.
   * Returns \p invalid_id if the node is not found, so \p invalid_id
   * can be thought of as a "default" boundary id.
   */
  short int boundary_id (const Node* node) const;

  /**
   * Returns the boundary id associated with the \p side side of
   * element \p elem.  Note that only one id per side is allowed,
   * however multiple sides per element are allowed.  Returns \p invalid_id
   * if the \p side does not have an associated boundary id, hence
   * \p invalid_id can be used as the default boundary id.
   */
  short int boundary_id (const Elem* const elem,
             const unsigned short int side, bool to_top_parent=true) const;

  /**
   * Returns a side of element \p elem whose associated boundary id is
   * \p boundary_id if such a side exists.
   * If multiple sides of \p elem have the same id, only the lowest numbered
   * such side is returned.
   *
   * Returns \p invalid_uint if no side has the requested boundary id.
   */
  unsigned int side_with_boundary_id(const Elem* const elem,
                     short int boundary_id) const;

  /**
   * @return true when given elem is a boundary element
   */
  bool is_boundary_elem(const Elem *e) const
  { return _boundary_side_id.find(e) != _boundary_side_id.end(); }


  /**
   * @returns the number of element-based boundary conditions.
   */
  unsigned int n_boundary_conds () const
  { return _boundary_side_id.size(); }


  /**
   * @return a list of nodes on boundary side, nodes are sorted by their id, node on the interface of two boundary will exist on both side
   */
  void boundary_side_nodes_with_id ( std::map<short int, std::set<const Node *> > & boundary_side_nodes_id_map) const;


  // NOTE: node are created by build_node_ids_from_priority_order for the following functions.

  /**
   * Creates a list of nodes and ids for those nodes.
   */
  void build_node_list (std::vector<unsigned int>& nl, std::vector<short int>&    il) const;

  /**
   * Creates a list of on processor nodes and ids for those nodes.
   */
  void build_on_processor_node_list (std::vector<unsigned int>& nl, std::vector<short int>&    il) const;

  /**
   * @return a list of nodes with the same ids, nodes are sorted by their id
   */
  void nodes_with_boundary_id (std::vector<unsigned int>& nl, short int boundary_id) const;

  /**
   * @return a list of node pointer with the same ids, nodes are sorted by their id
   */
  void nodes_with_boundary_id (std::vector<const Node *>& nl, short int boundary_id) const;

  /**
   * @return a list of nodes with the same ids, nodes are sorted by their id
   */
  void nodes_with_boundary_id ( std::map<short int, std::vector<const Node *> > & node_boundary_id_map) const;

  /**
   * @return the regions node belongs to
   */
  void node_in_regions ( std::map<const Node *, std::set<unsigned int > > & node_region_map) const;

  /**
   * @return active elem and its side with given boundary_id
   */
  void active_elem_with_boundary_id (std::vector<const Elem *>& el, std::vector<unsigned int>& sl, short int boundary_id) const;

  /**
   * @return active elem and its side with given boundary_id
   */
  void active_elem_with_boundary_id ( std::map<short int, std::vector< std::pair<const Elem *, unsigned int> > > & boundary_elem_side_map ) const;

  /**
   * Creates a list of element numbers, sides, and ids for those sides.
   */
  void build_side_list (std::vector<unsigned int>&       el,
            std::vector<unsigned short int>& sl,
            std::vector<short int>&          il) const;

  /**
   * Creates a list of active element numbers, sides, and ids for those sides.
   */
  void build_active_side_list (std::vector<unsigned int>&       el,
                   std::vector<unsigned short int>& sl,
                   std::vector<short int>&          il) const;


  /**
   * Creates a list of on processor element numbers, sides, and ids for those sides.
   */
  void build_on_processor_side_list (std::vector<unsigned int>&       el,
                                     std::vector<unsigned short int>& sl,
                                     std::vector<short int>&          il) const;

  /**
   * @returns the user-specified boundary ids.
   */
  const std::set<short int>& get_boundary_ids () const
  { return _boundary_ids; }

  /**
   * @returns the user-specified boundary ids.
   */
  std::set<short int>& get_boundary_ids ()
  { return _boundary_ids; }



  /**
   * @return the subdomain id of this boundary side lies on.
   * if the side on INTERFACE, sub_id1 and sub_id2 are meanful
   * if the side is on the BOUNDARY, sub_id2 is set to invalid_id
   */
  void get_subdomains_bd_on( std::map<short int,  std::pair<unsigned int, unsigned int > > & boundary_subdomain_map ) const;

  /**
   * @return all the boundary ids belongs to some subdomain.
   */
  void get_boundary_ids_by_subdomain (unsigned int subdomain, std::set<short int> & bd_ids) const;


  /**
   * @return the neighbor subdomains surround me (have common boundary with me)
   */
  void get_subdomain_neighbors( unsigned int subdomain, std::vector<unsigned int> & neighbor_subdomains) const;

  /**
   * Print the boundary information data structure.
   */
  void print_info () const;

  /**
   * Print the boundary label<->id information.
   */
  void print_boundary_label () const;

  /**
   * set the id with label
   */
  void set_label_to_id(short int id, const std::string & label, bool user_defined=true);

  /**
   * get the id by label
   */
  short int get_id_by_label(const std::string & label) const;

  /**
   * get the label by id
   */
  std::string get_label_by_id(short int id) const;

  /**
   * @return ture when boundary id has user defind (not assigned by genius) label
   */
  bool boundary_id_has_user_defined_label(short int id) const;

  /**
   * @return all the boundary ids on interface of two subdomains
   */
  std::set<short int> get_ids_on_region_interface(unsigned int sub1, unsigned int sub2) const;

  /**
   * set the id with description
   */
  void set_description_to_id(short int id, const std::string & description);

  /**
   * get the boundary description by id
   */
  std::string get_description_by_id(short int id) const;

  /**
   * record extra boundary description
   */
  void add_extra_description(const std::string &description)
  { _extra_descriptions.push_back(description); }


  /**
   * get extra boundary descriptions
   */
  std::vector<std::string> & extra_descriptions()
  { return _extra_descriptions; }

  /**
   * get extra boundary descriptions
   */
  const std::vector<std::string> & extra_descriptions() const
  { return _extra_descriptions; }

  /**
   * this function build neighbor information of interface elements
   */
  void find_neighbors();

  /**
   * rebuild _boundary_ids after node or elem remove.
   * prevent empty boundary id (no node or elem belongs to it)
   */
  void rebuild_ids();


  /**
   * Number used for internal use. This is the return value
   * if a boundary condition is not specified.
   */
  static const short int invalid_id;


 private:


  /**
   * The Mesh this boundary info pertains to.
   */
  const MeshBase& _mesh;

  /**
   * Data structure that maps nodes in the mesh
   * to boundary ids.
   */
  std::map<const Node*, short int> _boundary_node_id;

  /**
   * Data structure that maps sides of elements
   * to boundary ids.
   */
  std::multimap<const Elem*,
                std::pair<unsigned short int, short int> >
                                             _boundary_side_id;

  /**
   * A collection of user-specified boundary ids.
   */
  std::set<short int> _boundary_ids;


  /**
   * data structure for convert label to index
   */
  std::map<std::string, const short int> _boundary_labels_to_ids;

  /**
   * data structure for convert index to label
   */
  std::map<short int, const std::string> _boundary_ids_to_labels;

  /**
   * sometimes, not all the boundary_id has user defined label.
   * genius will try to assign an uinque label to each boundary_id,
   * this flag indicate that the boundary_id has user defined label.
   */
  std::set<short int> _boundary_id_has_user_defined_label;

  /**
   * data structure for boundary description
   */
  std::map<short int, const std::string> _boundary_ids_to_descriptions;

  /**
   * extra boundary descriptions
   */
  std::vector<std::string> _extra_descriptions;



//   /**
//    * Functor class for printing a single node's info
//    * To be used with "for_each".
//    */
//   class PrintNodeInfo
//   {
//   public:
//     inline
//     void operator() (const std::pair<const Node*, short int>& np) const
//     {
//       std::cout << "  (" << np.first->id()
//      << ", "  << np.second
//      << ")"  << std::endl;
//     }
//   };


//   /**
//    * Functor class for printing a single side's info.
//    * To be used with "for_each".
//    */
//   class PrintSideInfo
//   {
//   public:
//     PrintSideInfo() {}
//     inline
//     void operator() (const std::pair<const Elem*, std::pair<unsigned short int,short int> >& sp) const
//     {
//       std::cout << "  (" << sp.first->id()
//      << ", "  << sp.second.first
//      << ", "  << sp.second.second
//      << ")"   << std::endl;
//     }
//   };



  /**
   * Functor class for initializing a map.
   * The entries being added to the map
   * increase by exactly one each time.
   * The desctructor also inserts the
   * invalid_id entry.
   */
  class Fill
  {
  public:
    Fill(std::map<short int, unsigned int>& im) : id_map(im), cnt(0) {}

    ~Fill()
    {
      id_map[invalid_id] = cnt;
    }

    inline
    void operator() (const short int & pos)
    {
      id_map[pos] = cnt++;
    }

  private:
    std::map<short int, unsigned int>& id_map;
    unsigned int cnt;
  };

};


#endif
