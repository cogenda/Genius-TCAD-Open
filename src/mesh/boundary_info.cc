// $Id: boundary_info.cc,v 1.13 2008/06/22 02:58:03 gdiso Exp $

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


// Local includes
#include "boundary_info.h"
#include "boundary_mesh.h"
#include "elem.h"
#include "log.h"


#if defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
#elif defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER) || defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
#endif



//------------------------------------------------------
// BoundaryInfo static member initializations
const short int BoundaryInfo::invalid_id = -1234;



//------------------------------------------------------
// BoundaryInfo functions
BoundaryInfo::BoundaryInfo(const MeshBase& m) :
    _mesh (m)
{}



BoundaryInfo::~BoundaryInfo()
{
  this->clear();
}



void BoundaryInfo::clear()
{
  _boundary_node_id.clear();
  _boundary_side_id.clear();
  _boundary_ids.clear();
  _boundary_labels_to_ids.clear();
  _boundary_ids_to_labels.clear();
  _boundary_ids_to_descriptions.clear();
  _extra_descriptions.clear();
}



void BoundaryInfo::sync(BoundaryMesh& boundary_mesh)
{
  boundary_mesh.clear();

  /**
   * Re-create the boundary mesh.
   */

  std::map<unsigned int, unsigned int> new_node_numbers;

  boundary_mesh.set_n_subdomains() = this->n_boundary_ids();


  // Add sides to the structure.
  std::map<short int, unsigned int> id_map;

  // Original Code
  //     unsigned int cnt = 0;
  //     for (std::set<short int>::iterator pos = boundary_ids.begin();
  //     pos != boundary_ids.end(); ++pos)
  //       id_map[*pos] = cnt++;

  //     id_map[invalid_id] = cnt;


  // New code
  // Here we need to use iota() once it is in the
  // Utility namespace.
  std::for_each(_boundary_ids.begin(),
                _boundary_ids.end(),
                Fill(id_map));



  boundary_mesh.set_n_subdomains() = id_map.size();


  // Make individual copies of all the nodes in the current mesh
  // and add them to the boundary mesh.  Yes, this is overkill because
  // all of the current mesh nodes will not end up in the the boundary
  // mesh.  These nodes can be trimmed later via a call to prepare_for_use().
  {
    assert (boundary_mesh.n_nodes() == 0);
    boundary_mesh.reserve_nodes(_mesh.n_nodes());

    MeshBase::const_node_iterator it  = _mesh.nodes_begin();
    MeshBase::const_node_iterator end = _mesh.nodes_end();

    for(; it != end; ++it)
    {
      const Node* node = *it;
      boundary_mesh.add_point(*node); // calls Node::build(Point, id)
    }
  }

  // Add additional sides that aren't flagged with boundary conditions
  MeshBase::const_element_iterator       el     = _mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = _mesh.active_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    for (unsigned int s=0; s<elem->n_sides(); s++)
      if (elem->neighbor(s) == NULL) // on the boundary
      {

        // Build the side - do not use a "proxy" element here:
        // This will be going into the BoundaryMesh and needs to
        // stand on its own.
        AutoPtr<Elem> side (elem->build_side(s, false));

        // Get the top-level parent for this element
        const Elem* top_parent = elem->top_parent();

        // A convenient typedef
        typedef
        std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator
        Iter;

        // Find the right id number for that side
        std::pair<Iter, Iter> pos = _boundary_side_id.equal_range(top_parent);

        while (pos.first != pos.second)
        {
          if (pos.first->second.first == s) // already flagged with a boundary condition
          {
            side->subdomain_id() =
              id_map[pos.first->second.second];

            side->processor_id() =
              side->subdomain_id();
            break;
          }

          ++pos.first;
        }

        // either the element wasn't found or side s
        // doesn't have a boundary condition
        if (pos.first == pos.second)
        {
          side->subdomain_id() = id_map[invalid_id];
        }

        // Add the side
        Elem* new_elem = boundary_mesh.add_elem(side.release());

        // This side's Node pointers still point to the nodes of the  original mesh.
        // We need to re-point them to the boundary mesh's nodes!  Since we copied *ALL* of
        // the original mesh's nodes over, we should be guaranteed to have the same ordering.
        for (unsigned int nn=0; nn<new_elem->n_nodes(); ++nn)
        {
          // Get the correct node pointer, based on the id()
          Node* new_node = boundary_mesh.node_ptr(new_elem->node(nn));

          // sanity check: be sure that the new Nodes global id really matches
          assert (new_node->id() == new_elem->node(nn));

          // Assign the new node pointer
          new_elem->set_node(nn) = new_node;
        }
      }
  } // end loop over active elements



  // Trim any un-used nodes from the Mesh
  boundary_mesh.prepare_for_use();
}



void BoundaryInfo::add_node(const unsigned int node,  const short int id)
{
  this->add_node (_mesh.node_ptr(node), id);
}



void BoundaryInfo::add_node(const Node* node, const short int id)
{
  if (id == invalid_id)
  {
    std::cerr << "ERROR: You may not set a boundary ID of "
    << invalid_id << std::endl
    << " That is reserved for internal use.\n"
    << std::endl;

    genius_error();
  }

  _boundary_node_id[node] = id;
  _boundary_ids.insert(id);
}



void BoundaryInfo::add_side(const unsigned int e,
                            const unsigned short int side,
                            const short int id)
{
  this->add_side (_mesh.elem(e), side, id);
}



void BoundaryInfo::add_side(const Elem* elem,
                            const unsigned short int side,
                            const short int id)
{
  assert (elem != NULL);

  // Only add BCs for level-0 elements.
  assert (elem->level() == 0);

  if (id == invalid_id)
  {
    std::cerr << "ERROR: You may not set a boundary ID of "
    << invalid_id << std::endl
    << " That is reserved for internal use.\n"
    << std::endl;

    genius_error();
  }

  // if the elem/side pair has already been set, skip it
  {
    typedef std::multimap<const Elem*, std::pair<unsigned short int, short int> >::iterator It;
    std::pair<It, It> bound = _boundary_side_id.equal_range(elem);
    while(bound.first!=bound.second)
    {
      if(bound.first->second.first==side && bound.first->second.second==id) return;
      ++bound.first;
    }
  }

  std::pair<unsigned short int, short int> p(side,id);
  std::pair<const Elem*, std::pair<unsigned short int, short int> >  kv (elem, p);

  _boundary_side_id.insert(kv);
  _boundary_ids.insert(id);

  // Possilby add the nodes of the side,
  // no matter they are already there.
  /*
  {
    assert (side < elem->n_sides());

    AutoPtr<Elem> side_elem(elem->build_side(side));

    for (unsigned int n=0; n<side_elem->n_nodes(); n++)
      //if (this->boundary_id(side_elem->get_node(n)) == invalid_id)
      this->add_node(side_elem->get_node(n), id);
  }*/
}

void BoundaryInfo::boundary_side_nodes_with_id ( std::map<short int, std::set<const Node *> > & boundary_side_nodes_id_map) const
{
  boundary_side_nodes_id_map.clear();

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;
  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
  {
    const Elem* elem = (*pos).first;
    unsigned short int side = (*pos).second.first;
    short int bd_id = (*pos).second.second;

    std::vector<const Elem*> family;
    elem->active_family_tree_by_side(family, side);

    for(unsigned int f=0; f<family.size(); ++f)
    {
      const Elem* family_elem = family[f];
      AutoPtr<Elem> side_elem = family_elem->build_side(side);

      for (unsigned int n=0; n<side_elem->n_nodes(); ++n)
      {
        const Node * node = side_elem->get_node(n);
        boundary_side_nodes_id_map[bd_id].insert(node);
      }
    }
  }
}



void BoundaryInfo::build_node_ids_from_priority_order(const std::map<short int, unsigned int> & order)
{
  _boundary_node_id.clear();

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;

  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
  {
    const Elem* elem = (*pos).first;
    unsigned short int side = (*pos).second.first;
    short int bd_id = (*pos).second.second;

    std::vector<const Elem*> family;
    elem->active_family_tree_by_side(family, side);

    for(unsigned int f=0; f<family.size(); ++f)
    {
      const Elem* family_elem = family[f];
      AutoPtr<Elem> side_elem = family_elem->build_side(side);

      for (unsigned int n=0; n<side_elem->n_nodes(); ++n)
      {

        const Node * node = side_elem->get_node(n);

        std::map<const Node*, short int>::iterator it = _boundary_node_id.find(node);

        // the node is already exist
        if( it != _boundary_node_id.end()  )
        {
          // the bd id of existing node
          short int node_bd_id = (*it).second;

          //if current bd_id has a higher priority order, replace existing bd_id
          unsigned int o1 = (*order.find(bd_id)).second;
          unsigned int o2 = (*order.find(node_bd_id)).second;
          if(  o1 > o2  )
            _boundary_node_id[node] = bd_id;

          // or two node has the same priority order, set bd_id to little one
          if( o1 == o2 && bd_id < node_bd_id)
            _boundary_node_id[node] = bd_id;
        }
        // not exist yet? insert it
        else
          _boundary_node_id[node] = bd_id;

      }
    }
  }

}


void BoundaryInfo::remove (const Node* node)
{
  assert (node != NULL);

  // Erase everything associated with node
  _boundary_node_id.erase (node);

  // for efficency reason, we don't do it here.
  // please call rebuild_ids() after all the remove operator

  // after that, rebuild _boundary_ids
  //rebuild_ids();

}


void BoundaryInfo::remove (const Elem* elem)
{
  assert (elem != NULL);

  // Erase everything associated with elem
  _boundary_side_id.erase (elem);


  // for efficency reason, we don't do it here.
  // please call rebuild_ids() after all the remove operator

  // after that, rebuild _boundary_ids
  //rebuild_ids();
}



void BoundaryInfo::remove (const Elem* elem, unsigned short int side)
{
  assert (elem != NULL);

  // Erase element with the given side
  typedef std::multimap<const Elem*, std::pair<unsigned short int, short int> >::iterator It;
  std::pair<It,It> pos = _boundary_side_id.equal_range(elem);

  std::vector<It> its;

  //record all the element should be erased.
  while (pos.first != pos.second)
  {
    if( (*pos.first).second.first == side )
      its.push_back(pos.first);
    ++pos.first;
  }

  // erase here
  for(size_t n=0; n<its.size(); n++)
    _boundary_side_id.erase (its[n]);


  // for efficency reason, we don't do it here.
  // please call rebuild_ids() after all the remove operator

  // after that, rebuild _boundary_ids
  //rebuild_ids();
}



void BoundaryInfo::rebuild_ids()
{
  _boundary_ids.clear();

  std::map<const Node*, short int>::iterator node_it = _boundary_node_id.begin();
  for(; node_it != _boundary_node_id.end(); ++node_it)
    _boundary_ids.insert((*node_it).second);

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::iterator elem_it = _boundary_side_id.begin();
  for(; elem_it != _boundary_side_id.end(); ++elem_it)
    _boundary_ids.insert((*elem_it).second.second);
}




short int BoundaryInfo::boundary_id(const Node* node) const
{
  std::map<const Node*, short int>::const_iterator
  n = _boundary_node_id.find(node);

  // node not in the data structure
  if (n == _boundary_node_id.end())
    return invalid_id;

  return n->second;
}



short int BoundaryInfo::boundary_id(const Elem* const elem,  const unsigned short int side, bool to_top_parent) const
{
  assert (elem != NULL);


  const Elem*  searched_elem = elem;

  // Only level-0 elements store BCs.  If this is not a level-0
  // element get its level-0 parent and infer the BCs when to_top_parent is true.
  if ( to_top_parent && elem->level() != 0)
    searched_elem = elem->top_parent ();

  std::pair<std::multimap<const Elem*,
  std::pair<unsigned short int, short int> >::const_iterator,
  std::multimap<const Elem*,
  std::pair<unsigned short int, short int> >::const_iterator >
  e = _boundary_side_id.equal_range(searched_elem);

  // elem not in the data structure
  if (e.first == e.second)
    return invalid_id;


  // elem is there, maybe multiple occurances
  while (e.first != e.second)
  {
    // if this is true we found the requested side
    // of the element and want to return the id

    if (e.first->second.first == side)
      return e.first->second.second;

    ++e.first;
  }

  // if we get here, we found elem in the data structure but not
  // the requested side, so return the default value
  return invalid_id;
}


unsigned int BoundaryInfo::side_with_boundary_id(const Elem* const elem, short int boundary_id) const
{
  const Elem* searched_elem = elem;
  if (elem->level() != 0)
    searched_elem = elem->top_parent();

  typedef  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator CIter;
  std::pair<CIter,CIter> e = _boundary_side_id.equal_range(searched_elem);

  // elem not in the data structure
  if (e.first == e.second)
    return invalid_uint;

  // elem is there, maybe multiple occurances
  while (e.first != e.second)
  {
    // if this is true we found the requested boundary_id
    // of the element and want to return the side
    if (e.first->second.second == boundary_id)
      return e.first->second.first;

    ++e.first;
  }

  // if we get here, we found elem in the data structure but not
  // the requested boundary id, so return the default value
  return invalid_uint;
}



void BoundaryInfo::build_node_list (std::vector<unsigned int>& nl,
                                    std::vector<short int>&    il) const
{
  // Reserve the size, then use push_back
  nl.reserve (_boundary_node_id.size());
  il.reserve (_boundary_node_id.size());

  std::map<const Node*, short int>::const_iterator pos
  = _boundary_node_id.begin();

  for (; pos != _boundary_node_id.end(); ++pos)
  {
    nl.push_back (pos->first->id());
    il.push_back (pos->second);
  }
}


void BoundaryInfo::build_on_processor_node_list (std::vector<unsigned int>& nl, std::vector<short int>& il) const
{
  nl.clear();
  il.clear();

  std::map<const Node*, short int>::const_iterator pos
      = _boundary_node_id.begin();

  for (; pos != _boundary_node_id.end(); ++pos)
  {
    if(pos->first->processor_id() != Genius::processor_id()) continue;

    nl.push_back (pos->first->id());
    il.push_back (pos->second);
  }
}


void BoundaryInfo::build_side_list (std::vector<unsigned int>&       el,
                                    std::vector<unsigned short int>& sl,
                                    std::vector<short int>&          il) const
{
  el.resize (_boundary_side_id.size());
  sl.resize (_boundary_side_id.size());
  il.resize (_boundary_side_id.size());

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos = _boundary_side_id.begin();

  for (unsigned int n=0; pos != _boundary_side_id.end(); ++pos, ++n)
  {
    el[n] = pos->first->id();
    sl[n] = pos->second.first;
    il[n] = pos->second.second;
  }
}



void BoundaryInfo::build_active_side_list (std::vector<unsigned int>&       el,
                                           std::vector<unsigned short int>& sl,
                                           std::vector<short int>&          il) const
{
  el.clear();
  sl.clear();
  il.clear();

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;

  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end();
       ++pos)
  {
    const Elem * elem = pos->first;
    if (elem->active() )
    {
      el.push_back (pos->first->id());
      sl.push_back (pos->second.first);
      il.push_back (pos->second.second);
    }
    // this element has child
    else
    {
      unsigned short int side = pos->second.first;
      std::vector<const Elem*> family;

      elem->active_family_tree_by_side(family, side);

      for(unsigned int n=0; n<family.size(); ++n)
      {
        el.push_back (family[n]->id());
        sl.push_back (side);
        il.push_back (pos->second.second);
      }

    }
  }

}


void BoundaryInfo::build_on_processor_side_list (std::vector<unsigned int>&       el,
                                                 std::vector<unsigned short int>& sl,
                                                 std::vector<short int>&          il) const
{
  el.clear();
  sl.clear();
  il.clear();

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos = _boundary_side_id.begin();

  for (unsigned int n=0; pos != _boundary_side_id.end(); ++pos, ++n)
  {
    if(pos->first->processor_id() != Genius::processor_id()) continue;

    el.push_back(pos->first->id());
    sl.push_back(pos->second.first);
    il.push_back(pos->second.second);
  }

}



void BoundaryInfo::nodes_with_boundary_id (std::vector<unsigned int>& nl, short int boundary_id) const
{
  nl.clear();
  // use set to reorder the nodes by their id
  std::set<unsigned int> bd_node_set;

  std::map<const Node*, short int>::const_iterator pos = _boundary_node_id.begin();
  for (; pos != _boundary_node_id.end(); ++pos)
  {
    if( pos->second==boundary_id )
      bd_node_set.insert(pos->first->id());
  }

  for(std::set<unsigned int>::iterator it=bd_node_set.begin(); it!=bd_node_set.end(); ++it)
      nl.push_back(*it);
}



void BoundaryInfo::nodes_with_boundary_id (std::vector<const Node *>& nl, short int boundary_id) const
{
  nl.clear();
  // use map to reorder the nodes by their id
  std::map<unsigned int, const Node *> bd_node_map;

  std::map<const Node*, short int>::const_iterator pos = _boundary_node_id.begin();
  for (; pos != _boundary_node_id.end(); ++pos)
  {
    if( pos->second==boundary_id )
      bd_node_map[pos->first->id()] = pos->first;
  }

  for(std::map<unsigned int, const Node *>::iterator it=bd_node_map.begin(); it!=bd_node_map.end(); ++it)
    nl.push_back((*it).second);
}

void BoundaryInfo::nodes_with_boundary_id ( std::map<short int, std::vector<const Node *> > & node_boundary_id_map) const
{
  node_boundary_id_map.clear();

  std::map<const Node*, short int>::const_iterator pos = _boundary_node_id.begin();
  for (; pos != _boundary_node_id.end(); ++pos)
  {
    node_boundary_id_map[pos->second].push_back(pos->first);
  }
}



void BoundaryInfo::node_in_regions ( std::map<const Node *, std::set<unsigned int > > & node_region_map) const
{
  node_region_map.clear();

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;
  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
  {
    const Elem* elem = (*pos).first;
    unsigned short int side = (*pos).second.first;
    short int bd_id = (*pos).second.second;

    std::vector<const Elem*> family;
    elem->active_family_tree_by_side(family, side);

    for(unsigned int f=0; f<family.size(); ++f)
    {
      const Elem* family_elem = family[f];
      AutoPtr<Elem> side_elem = family_elem->build_side(side);

      for (unsigned int n=0; n<side_elem->n_nodes(); ++n)
      {
        const Node * node = side_elem->get_node(n);
        node_region_map[node].insert(elem->subdomain_id());
      }
    }
  }
  assert(node_region_map.size() == _boundary_node_id.size());
}



void BoundaryInfo::active_elem_with_boundary_id (std::vector<const Elem *>& el, std::vector<unsigned int>& sl, short int boundary_id) const
{
  el.clear();
  sl.clear();

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;

  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end();
       ++pos)
  {
    const Elem * elem = pos->first;
    if (elem->active() )
    {
      if( pos->second.second == boundary_id )
      {
        el.push_back (pos->first);
        sl.push_back (pos->second.first);
      }
    }
    // this element has child
    else
    {
      unsigned short int side = pos->second.first;
      std::vector<const Elem*> family;

      elem->active_family_tree_by_side(family, side);

      for(unsigned int n=0; n<family.size(); ++n)
      {
        if( pos->second.second == boundary_id )
        {
          el.push_back (family[n]);
          sl.push_back (side);
        }
      }

    }
  }

}

void BoundaryInfo::active_elem_with_boundary_id (
  std::map<short int, std::vector< std::pair<const Elem *, unsigned int> > > & boundary_elem_side_map ) const
{
  boundary_elem_side_map.clear();

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;
  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
  {
    const Elem * elem = pos->first;
    if (elem->active() )
    {
      boundary_elem_side_map[pos->second.second].push_back(std::make_pair(pos->first, pos->second.first));
    }
    // this element has child
    else
    {
      unsigned short int side = pos->second.first;
      std::vector<const Elem*> family;

      elem->active_family_tree_by_side(family, side);

      for(unsigned int n=0; n<family.size(); ++n)
      {
        boundary_elem_side_map[pos->second.second].push_back(std::make_pair(pos->first, pos->second.first));
      }
    }
  }
}



void BoundaryInfo::get_subdomains_bd_on(short int boundary_id, unsigned int & sub_id1, unsigned int & sub_id2) const
{

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;

  std::set< unsigned int > sub_ids;


  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
  {
    // if not match required  boundary_id, continue
    if( pos->second.second != boundary_id ) continue;

    // only process element with this boundary_id
    const Elem * elem = pos->first;
    sub_ids.insert( elem->subdomain_id() );

    // get the side
    unsigned int side = pos->second.first;

    const Elem * neighbor_elem = elem->neighbor(side);

    if( neighbor_elem == NULL )
      sub_ids.insert( invalid_uint );
    else
      sub_ids.insert( neighbor_elem->subdomain_id() );

  }

  // fatal error
  if( sub_ids.size() != 2 )
  {
    MESSAGE<<std::endl;
    MESSAGE<<"Fatal Error: Genius assert that all the boundary faces (edges) with same"<<std::endl;
    MESSAGE<<"boundary label should either be external boundary of a special region or"<<std::endl;
    MESSAGE<<"interface between two regions."<<std::endl;
    MESSAGE<<"However, the mesh boundary "<<get_label_by_id(boundary_id)<<" breaks this assert."<<std::endl;
    MESSAGE<<"This boundary involves "<<sub_ids.size()<<" regions."<<std::endl;
    MESSAGE<<"Please check the device mesh and try again."<<std::endl;
    RECORD();
    genius_error();
  }

  std::set< unsigned int >::iterator it = sub_ids.begin();

  it = sub_ids.begin();
  sub_id1 = *it;
  ++it;
  sub_id2 = *it;

}


void BoundaryInfo::get_subdomains_bd_on( std::map<short int,  std::pair<unsigned int, unsigned int > > & boundary_subdomain_map ) const
{
  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;

  std::map<short int, std::set< unsigned int > > boundary_sub_ids;
  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
  {
    short int boundary_id = pos->second.second;

    const Elem * elem = pos->first;
    boundary_sub_ids[boundary_id].insert( elem->subdomain_id() );

    // get the side
    unsigned int side = pos->second.first;
    const Elem * neighbor_elem = elem->neighbor(side);
    if( neighbor_elem == NULL )
      boundary_sub_ids[boundary_id].insert( invalid_uint );
    else
      boundary_sub_ids[boundary_id].insert( neighbor_elem->subdomain_id() );
  }

  // fatal error
  std::map<short int, std::set< unsigned int > >::const_iterator boundary_sub_ids_it = boundary_sub_ids.begin();
  std::map<short int, std::set< unsigned int > >::const_iterator boundary_sub_ids_it_end = boundary_sub_ids.end();
  for(; boundary_sub_ids_it != boundary_sub_ids_it_end; ++boundary_sub_ids_it )
  {
    short int boundary_id = boundary_sub_ids_it->first;
    const std::set< unsigned int > & sub_ids = boundary_sub_ids_it->second;
    if( sub_ids.size() != 2 )
    {
      MESSAGE<<std::endl;
      MESSAGE<<"Fatal Error: Genius assert that all the boundary faces (edges) with same"<<std::endl;
      MESSAGE<<"boundary label should either be external boundary of a special region or"<<std::endl;
      MESSAGE<<"interface between two regions."<<std::endl;
      MESSAGE<<"However, the mesh boundary "<<get_label_by_id(boundary_id)<<" breaks this assert."<<std::endl;
      MESSAGE<<"This boundary involves "<<sub_ids.size()<<" regions."<<std::endl;
      MESSAGE<<"Please check the device mesh and try again."<<std::endl;
      RECORD();
      genius_error();
    }
    std::set< unsigned int >::const_iterator it = sub_ids.begin();
    unsigned int sub_id1 = *it;
    ++it;
    unsigned int sub_id2 = *it;
    boundary_subdomain_map[boundary_id] = std::make_pair(sub_id1, sub_id2);
  }
}



void BoundaryInfo::get_boundary_ids_by_subdomain (unsigned int subdomain, std::set<short int> & bd_ids) const
  {
    bd_ids.clear();

    std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;

    for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
    {

      const Elem * elem = pos->first;

      if( elem->subdomain_id() != subdomain) continue;

      // only process element at this subdomain
      bd_ids.insert(pos->second.second);
    }

  }



void BoundaryInfo::get_subdomain_neighbors( unsigned int subdomain, std::vector<unsigned int> & neighbor_subdomains) const
{
  std::set<unsigned int> neighbor_subdomains_set;

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;
  for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
  {

    const Elem * elem = pos->first;

    if( elem->subdomain_id() != subdomain) continue;

    // only process element at this subdomain
    const Elem * neighbor_elem = elem->neighbor(pos->second.first);
    if( neighbor_elem )
      neighbor_subdomains_set.insert(neighbor_elem->subdomain_id());
  }

  // save to vector
  neighbor_subdomains.clear();
  std::set<unsigned int>::iterator it =  neighbor_subdomains_set.begin();
  for(; it!= neighbor_subdomains_set.end(); ++it)
  {
    neighbor_subdomains.push_back(*it);
  }

}


std::set<short int> BoundaryInfo::get_ids_on_region_interface(unsigned int sub1, unsigned int sub2) const
  {
    std::set<short int> bds;
    std::multimap<const Elem*, std::pair<unsigned short int, short int> >::const_iterator pos;
    for (pos=_boundary_side_id.begin(); pos != _boundary_side_id.end(); ++pos)
    {
      const Elem * elem = pos->first;
      if( elem->subdomain_id() != sub1) continue; // only consider elem on subdomain 1

      // if my neighbor on subdomain 2, record it
      const Elem * beighbor_elem = elem->neighbor(pos->second.first);
      if(!beighbor_elem) continue;
      if( beighbor_elem->subdomain_id() == sub2) bds.insert(pos->second.second);
    }
    return bds;
  }


void BoundaryInfo::find_neighbors()
{

  // Find neighboring elements by first finding elements
  // with identical side keys and then check to see if they
  // are neighbors

  // data structures -- Use the unordered_multimap if available
  typedef unsigned int                    key_type;
  typedef std::pair<Elem*, unsigned char> val_type;
  typedef std::pair<key_type, val_type>   key_val_pair;

#if defined(HAVE_UNORDERED_MAP)
  typedef std::unordered_multimap<key_type, val_type> map_type;
#elif defined(HAVE_TR1_UNORDERED_MAP) || defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER)
  typedef std::tr1::unordered_multimap<key_type, val_type> map_type;
#else
  typedef std::multimap<key_type, val_type>  map_type;
#endif

  // A map from side keys to corresponding elements & side numbers
  map_type side_to_elem_map;

  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::iterator el = _boundary_side_id.begin();
  std::multimap<const Elem*, std::pair<unsigned short int, short int> >::iterator el_end = _boundary_side_id.end();
  for (; el != el_end; ++el)
  {
   next_side:

    Elem * element    = const_cast<Elem *>(el->first);
    unsigned int side = el->second.first;
    if (element->neighbor(side) == NULL)
    {
      // Get the key for the side of this element
      const unsigned int key = element->key(side);

      // Look for elements that have an identical side key
      std::pair <map_type::iterator, map_type::iterator>
      bounds = side_to_elem_map.equal_range(key);

      // May be multiple keys, check all the possible
      // elements which _might_ be neighbors.
      if (bounds.first != bounds.second)
      {
        // Get the side for this element
        const AutoPtr<DofObject> my_side(element->side(side));

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
          if( *my_side == *their_side )
          {
            element->set_neighbor (side, neighbor);
            neighbor->set_neighbor(ns, element);
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
      kvp.second.second = side;

      // use the lower bound as a hint for
      // where to put it.
      side_to_elem_map.insert (bounds.first,kvp);
    }

  }
}


void BoundaryInfo::print_info() const
{
  // Print out the nodal BCs
  if (!_boundary_node_id.empty())
  {
    std::cout << "Nodal Boundary conditions:" << std::endl
    << "--------------------------" << std::endl
    << "  (Node No., ID)               " << std::endl;

    //       std::for_each(_boundary_node_id.begin(),
    //          _boundary_node_id.end(),
    //          PrintNodeInfo());

    std::map<const Node*, short int>::const_iterator it        = _boundary_node_id.begin();
    const std::map<const Node*, short int>::const_iterator end = _boundary_node_id.end();

    for (; it != end; ++it)
      std::cout << "  (" << (*it).first->id()
      << ", "  << (*it).second
      << ")"  << std::endl;
  }

  // Print out the element BCs
  if (!_boundary_side_id.empty())
  {
    std::cout << std::endl
    << "Side Boundary conditions:" << std::endl
    << "-------------------------" << std::endl
    << "  (Elem No., Side No., ID)      " << std::endl;

    //       std::for_each(_boundary_side_id.begin(),
    //          _boundary_side_id.end(),
    //              PrintSideInfo());

    std::multimap<const Elem*,
    std::pair<unsigned short int, short int> >::const_iterator it = _boundary_side_id.begin();
    const std::multimap<const Elem*,
    std::pair<unsigned short int, short int> >::const_iterator end = _boundary_side_id.end();

    for (; it != end; ++it)
      std::cout << "  (" << (*it).first->id()
      << ", "  << (*it).second.first
      << ", "  << (*it).second.second
      << ")"   << std::endl;
  }


}


void BoundaryInfo::print_boundary_label() const
{
  std::cout<<"Boundary label in processor " << Genius::processor_id()<<std::endl;
  std::map<const std::string, const short int>::const_iterator it     = _boundary_labels_to_ids.begin();
  std::map<const std::string, const short int>::const_iterator it_end = _boundary_labels_to_ids.end();
  for(; it!=it_end; ++it)
  {
    std::cout<< " Label: " << it->first <<"  ID: " <<it->second<<std::endl;
  }

}


void BoundaryInfo::set_label_to_id(short int id, const std::string & label)
{

  _boundary_labels_to_ids.insert(std::pair<const std::string, const short int>(label,id));
  _boundary_ids_to_labels.insert(std::pair<const short int, const std::string>(id,label));

}



short int BoundaryInfo::get_id_by_label(const std::string & label) const
{
  if ( _boundary_labels_to_ids.find(label) != _boundary_labels_to_ids.end() )
    return (*_boundary_labels_to_ids.find(label)).second;

  return   invalid_id;
}



std::string BoundaryInfo::get_label_by_id(short int id) const
{

  std::string result("invalid_bd");

  if ( _boundary_ids_to_labels.find(id) != _boundary_ids_to_labels.end() )
    return (*_boundary_ids_to_labels.find(id)).second;

  return result  ;
}



void BoundaryInfo::set_description_to_id(short int id, const std::string & description)
{
  _boundary_ids_to_descriptions.insert(std::pair<const short int, const std::string>(id, description));
}


std::string BoundaryInfo::get_description_by_id(short int id) const
{
  if(_boundary_ids_to_descriptions.find(id)!=_boundary_ids_to_descriptions.end())
    return _boundary_ids_to_descriptions.find(id)->second;
  return std::string();
}



