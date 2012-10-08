// $Id: serial_mesh.cc,v 1.5 2008/05/18 03:18:27 gdiso Exp $

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
#include <queue>
#include <algorithm>

// Local includes
#include "boundary_info.h"
#include "elem.h"
#include "perf_log.h"
#include "serial_mesh.h"
#include "mesh_tools.h"
#include "parallel.h"

// ------------------------------------------------------------
// SerialMesh class member functions
SerialMesh::SerialMesh (unsigned int d) :
    UnstructuredMesh (d),_is_serial(true)
{}


SerialMesh::~SerialMesh ()
{
  this->clear();  // Free nodes and elements
}


// This might be specialized later, but right now it's just here to
// make sure the compiler doesn't give us a default (non-deep) copy
// constructor instead.
SerialMesh::SerialMesh (const SerialMesh &other_mesh) :
    UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
}


SerialMesh::SerialMesh (const UnstructuredMesh &other_mesh) :
    UnstructuredMesh (other_mesh)
{
  this->copy_nodes_and_elements(other_mesh);
}


const Point& SerialMesh::point (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() == i); // This will change soon

  return (*_nodes[i]);
}





const Node& SerialMesh::node (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL );
  assert (_nodes[i]->id() == i); // This will change soon

  return (*_nodes[i]);
}





Node& SerialMesh::node (const unsigned int i)
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL );
  assert (_nodes[i]->id() == i); // This will change soon

  return (*_nodes[i]);
}



const Node* SerialMesh::node_ptr (const unsigned int i) const
{
  assert (_nodes[i] == NULL || _nodes[i]->id() == i);

  return _nodes[i];
}




Node* & SerialMesh::node_ptr (const unsigned int i)
{
  assert (_nodes[i] == NULL || _nodes[i]->id() == i);

  return _nodes[i];
}




Elem* SerialMesh::elem (const unsigned int i) const
{
  assert (_elements[i] == NULL || _elements[i]->id() == i);

  return _elements[i];
}




Elem* SerialMesh::add_elem (Elem* e)
{
  if (e != NULL)
    e->set_id (_elements.size());

  _elements.push_back(e);

  return e;
}



Elem* SerialMesh::add_elem (Elem* e, unsigned int id)
{
  genius_assert((e != NULL));
  e->set_id (id);

  genius_assert(id < _elements.size());
  _elements[id] = e;

  return e;
}


Elem* SerialMesh::insert_elem (Elem* e)
{
  unsigned int eid = e->id();
  genius_assert(eid < _elements.size());
  Elem *oldelem = _elements[eid];

  if (oldelem)
  {
    genius_assert(oldelem->id() == eid);
    this->delete_elem(oldelem);
  }

  _elements[e->id()] = e;

  return e;
}



void SerialMesh::delete_elem(Elem* e)
{
  assert (e != NULL);

  // Initialize an iterator to eventually point to the element we want to delete
  std::vector<Elem*>::iterator pos = _elements.end();

  // In many cases, e->id() gives us a clue as to where e
  // is located in the _elements vector.  Try that first
  // before trying the O(n_elem) search.
  assert (e->id() < _elements.size());

  if (_elements[e->id()] == e)
  {
    // We found it!
    pos = _elements.begin();
    std::advance(pos, e->id());
  }

  else
  {
    // This search is O(n_elem)
    pos = std::find (_elements.begin(),
                     _elements.end(),
                     e);
  }

  // Huh? Element not in the vector?
  assert (pos != _elements.end());

  // Remove the element from the BoundaryInfo object
  this->boundary_info->remove(e);

  // delete the element
  delete e;

  // explicitly NULL the pointer
  *pos = NULL;
}



Node* SerialMesh::add_point (const Point& p,
                             const unsigned int id,
                             const unsigned int proc_id)
{
  //   // We only append points with SerialMesh
  //   genius_assert(id == DofObject::invalid_id || id == _nodes.size());
  //   Node *n = Node::build(p, _nodes.size()).release();
  //   n->processor_id() = proc_id;
  //   _nodes.push_back (n);

  Node *n = NULL;

  // If the user requests a valid id, either
  // provide the existing node or resize the container
  // to fit the new node.
  if (id != DofObject::invalid_id)
    if (id < _nodes.size())
      n = _nodes[id];
    else
      _nodes.resize(id+1);
  else
    _nodes.push_back (static_cast<Node*>(NULL));

  // if the node already exists, then assign new (x,y,z) values
  if (n)
    *n = p;
  // otherwise build a new node, put it in the right spot, and return
  // a valid pointer.
  else
  {
    n = Node::build(p, (id == DofObject::invalid_id) ? _nodes.size()-1 : id).release();
    n->processor_id() = proc_id;

    if (id == DofObject::invalid_id)
      _nodes.back() = n;
    else
      _nodes[id] = n;
  }

  // better not pass back a NULL pointer.
  assert (n);

  return n;
}



Node* SerialMesh::add_node (Node* n)
{
  genius_assert(n);
  // We only append points with SerialMesh
  genius_assert(!n->valid_id() || n->id() == _nodes.size());

  n->set_id (_nodes.size());

  _nodes.push_back(n);

  return n;
}



void SerialMesh::delete_node(Node* n)
{
  assert (n != NULL);
  assert (n->id() < _nodes.size());

  // Initialize an iterator to eventually point to the element we want
  // to delete
  std::vector<Node*>::iterator pos;

  // In many cases, e->id() gives us a clue as to where e
  // is located in the _elements vector.  Try that first
  // before trying the O(n_elem) search.
  if (_nodes[n->id()] == n)
  {
    pos = _nodes.begin();
    std::advance(pos, n->id());
  }
  else
  {
    pos = std::find (_nodes.begin(),
                     _nodes.end(),
                     n);
  }

  // Huh? Node not in the vector?
  assert (pos != _nodes.end());

  // Delete the node from the BoundaryInfo object
  this->boundary_info->remove(n);

  // delete the node
  delete n;

  // explicitly NULL the pointer
  *pos = NULL;
}



void SerialMesh::reorder_elems()
{
  START_LOG("reorder_elems()", "Mesh");

  // do it only on serial mesh
  assert(_is_serial);

  {
    const Elem * elem_begin = _elements[0];
    for(unsigned int n=0; n<_elements.size(); ++n)
      if( _elements[n]->centroid().all_less(elem_begin->centroid()) )
        elem_begin = _elements[n];

    // the visited flag array
    std::vector<bool> visit_flag(n_elem(), false);
    unsigned int new_index = 0;
    // a queue for Breadth-First Search
    std::queue<const Elem *> Q;
    std::vector<unsigned int> new_order(n_elem(), invalid_uint);

    // do Breadth-First Search
    // begin at this node.
    visit_flag[elem_begin->id()] = true;
    Q.push( elem_begin );

    while(!Q.empty())
    {
      const Elem * current = Q.front();
      Q.pop();
      new_order[current->id()] = new_index++;

      // here just use a simple neighbor search method.
      // the neighbor is in arbitrary ordered.
      // however, one may use other nodal reordering
      // i.e. sort neighbor nodes by Degree of the neighboring nodes,
      // Physical distance of the neighboring nodes to current node
      // or even the solution based sort
      for(unsigned int e=0; e<current->n_neighbors(); ++e)
      {
        const Elem * neighbor = current->neighbor(e);
        if(neighbor && !visit_flag[neighbor->id()] )
        {
          Q.push(neighbor);
          visit_flag[neighbor->id()] = true;
        }
      }
    }

    // ok, assign ordered index to each elem
    for(unsigned int n=0; n<_elements.size(); ++n)
      _elements[n]->set_id () = new_order[_elements[n]->id()];

    // sort the elems by new ID
    DofObject::Less less;
    std::sort( _elements.begin(), _elements.end(), less );
  }


  // also sort nodes
  {
    std::vector<bool> visit_flag(n_nodes(), false);
    std::vector<unsigned int> new_order(n_nodes(), invalid_uint);
    unsigned int new_index = 0;
    for(unsigned int n=0; n<_elements.size(); ++n)
    {
      const Elem * elem = _elements[n];
      for( unsigned int v=0; v<elem->n_nodes(); ++v)
      {
        const Node * node = elem->get_node(v);

        if( !visit_flag[node->id()] )
        {
          new_order[node->id()] = new_index++;
          visit_flag[node->id()] = true;
        }
      }
    }

    // ok, assign ordered index to each node
    for (unsigned int n=0; n<_nodes.size(); ++n)
      _nodes[n]->set_id() = new_order[_nodes[n]->id()];

    // sort the nodes by new ID
    DofObject::Less less;
    std::sort( _nodes.begin(), _nodes.end(), less );
  }

  STOP_LOG("reorder_elems()", "Mesh");
}



void SerialMesh::reorder_nodes()
{
  START_LOG("reorder_nodes()", "Mesh");

  // do it only on serial mesh
  assert(_is_serial);

  // reorder the node index by Reverse Cuthill-McKee Algorithm
  // which can reduce filling in LU (ILU) Factorization
  node_iterator       node_it  = nodes_begin();
  const node_iterator node_end = nodes_end();

  // some prepare work for finding node neighbors
  std::vector<std::vector<const Elem*> > nodes_to_elem_map;
  MeshTools::build_nodes_to_elem_map (*this, nodes_to_elem_map);

  //find the node which has the minimal [x y z] coordinates.
  node_iterator start_node = nodes_begin();

  for (node_it  = nodes_begin() ; node_it != node_end; ++node_it)
    if( (*(*node_it)).all_less(*(*start_node)))
      start_node = node_it;


  // the visited flag array
  std::vector<bool> visit_flag(n_nodes(), false);
  unsigned int new_index = 0;
  // a queue for Breadth-First Search
  std::queue<const Node *> Q;
  std::vector<unsigned int> new_order(n_nodes(), invalid_uint);

  // do Breadth-First Search
  // begin at this node.
  visit_flag[(*start_node)->id()] = true;
  Q.push( *start_node );

  while(!Q.empty())
  {
    const Node * current = Q.front();
    Q.pop();
    new_order[current->id()] = new_index++;

    std::vector<const Node*> neighbors;
    std::vector<const Node*>::iterator it;
    MeshTools::find_nodal_neighbors(*this, *current, nodes_to_elem_map, neighbors, false);

    // here just use a simple neighbor search method.
    // the neighbor is in arbitrary ordered.
    // however, one may use other nodal reordering
    // i.e. sort neighbor nodes by Degree of the neighboring nodes,
    // Physical distance of the neighboring nodes to current node
    // or even the solution based sort
    for(it=neighbors.begin(); it!=neighbors.end(); it++)
      if( !visit_flag[(*it)->id()] )
      {
        Q.push(*it);
        visit_flag[(*it)->id()] = true;
      }
  }

  // ok, assign ordered index to each node
  for (node_it  = nodes_begin() ; node_it != nodes_end(); ++node_it)
    (*node_it)->set_id () = new_order[(*node_it)->id()];

  // sort the nodes by new ID
  DofObject::Less less;
  std::sort( _nodes.begin(), _nodes.end(), less );

  STOP_LOG("reorder_nodes()", "Mesh");
}



void SerialMesh::clear ()
{
  // Call parent clear function
  MeshBase::clear();


  // Clear our elements and nodes
  {
    std::vector<Elem*>::iterator       it  = _elements.begin();
    const std::vector<Elem*>::iterator end = _elements.end();

    // There is no need to remove the elements from
    // the BoundaryInfo data structure since we
    // already cleared it.
    for (; it != end; ++it)
      delete *it;

    _elements.clear();
  }

  // clear the nodes data structure
  {
    std::vector<Node*>::iterator       it  = _nodes.begin();
    const std::vector<Node*>::iterator end = _nodes.end();

    // There is no need to remove the nodes from
    // the BoundaryInfo data structure since we
    // already cleared it.
    for (; it != end; ++it)
      delete *it;

    _nodes.clear();
  }

}


AutoPtr<Elem> SerialMesh::elem_clone (const unsigned int i) const
{
  const Elem * elem = _elements[i];

  std::vector<Real> pts; // elem node location
  std::vector<int> conn; // elem info
  unsigned int elem_processor_id=invalid_uint; // elem on which processor

  if( elem && elem->processor_id()==Genius::processor_id())
  {
    for(unsigned int n=0; n<elem->n_nodes(); ++n)
    {
      const Point & p = elem->point(n);
      pts.push_back( p[0] );
      pts.push_back( p[1] );
      pts.push_back( p[2] );
    }
    elem->pack_element(conn);
    elem_processor_id = elem->processor_id();
  }

  Parallel::min(elem_processor_id); // sync the processor_id which elem on
  Parallel::broadcast(pts,  elem_processor_id);
  Parallel::broadcast(conn,  elem_processor_id);

  unsigned int cnt = 0;

  const ElemType elem_type    = static_cast<ElemType>(conn[cnt++]);
  AutoPtr<Elem> clone = Elem::build_clone(elem_type);//build ElemClone
  //Elem * clone = new ElemClone<Elem>();
  clone->processor_id() = conn[cnt++];
  clone->subdomain_id() = conn[cnt++];
  clone->set_id() = conn[cnt++];
#ifdef ENABLE_AMR
  const int level             = conn[cnt++];
  const int p_level           = conn[cnt++];
  const Elem::RefinementState refinement_flag =  static_cast<Elem::RefinementState>(conn[cnt++]);
  const Elem::RefinementState p_refinement_flag = static_cast<Elem::RefinementState>(conn[cnt++]);
  const int parent_ID         = conn[cnt++];
  const int which_child       = conn[cnt++];
  clone->set_refinement_flag(refinement_flag);
  clone->set_p_refinement_flag(p_refinement_flag);
  clone->set_p_level(p_level);
#endif
  // Assign the nodes, these nodes will be deleted by ElemClone
  for (unsigned int n=0; n<clone->n_nodes(); n++)
  {
    clone->set_node(n) = new Node(pts[3*n+0], pts[3*n+1], pts[3*n+2], conn[cnt++]);
  }

  return clone;
}


void SerialMesh::renumber_nodes_and_elements ()
{

  START_LOG("renumber_nodes_and_elem()", "Mesh");

  assert(_is_serial);

  // node and element id counters
  unsigned int next_free_elem = 0;
  unsigned int next_free_node = 0;

  // Loop over the elements.  Note that there may
  // be NULLs in the _elements vector from the coarsening
  // process.  Pack the elements in to a contiguous array
  // and then trim any excess.
  {
    std::vector<Elem*>::iterator in        = _elements.begin();
    std::vector<Elem*>::iterator out       = _elements.begin();
    const std::vector<Elem*>::iterator end = _elements.end();

    for (; in != end; ++in)
      if (*in != NULL)
      {
        Elem* elem = *in;

        *out = *in;
        ++out;

        // Increment the element counter
        elem->set_id (next_free_elem++);

        // Loop over this element's nodes.  Number them,
        // if they have not been numbered already.  Also,
        // position them in the _nodes vector so that they
        // are packed contiguously from the beginning.
        for (unsigned int n=0; n<elem->n_nodes(); n++)
          if (elem->node(n) == next_free_node)     // don't need to process
            next_free_node++;                      // [(src == dst) below]

          else if (elem->node(n) > next_free_node) // need to process
          {
            // The source and destination indices
            // for this node
            const unsigned int src_idx = elem->node(n);
            const unsigned int dst_idx = next_free_node++;

            // ensure we want to swap valid nodes
            assert (_nodes[src_idx] != NULL);
            assert (_nodes[dst_idx] != NULL);

            // Swap the source and destination nodes
            std::swap(_nodes[src_idx],
                      _nodes[dst_idx] );

            // Set proper indices
            _nodes[src_idx]->set_id (src_idx);
            _nodes[dst_idx]->set_id (dst_idx);
          }
      }

    // Erase any additional storage. These elements have been
    // copied into NULL voids by the procedure above, and are
    // thus repeated and unnecessary.
    _elements.erase (out, end);
  }

  // Any nodes in the vector >= _nodes[next_free_node]
  // are not connected to any elements and may be deleted
  // if desired.

  // (This code block will erase the unused nodes)
  // Now, delete the unused nodes
  {
    std::vector<Node*>::iterator nd        = _nodes.begin();
    const std::vector<Node*>::iterator end = _nodes.end();

    std::advance (nd, next_free_node);

    for (std::vector<Node*>::iterator it=nd;
         it != end; ++it)
    {
      assert (*it != NULL);

      // remove any boundary information associated with
      // this node
      this->boundary_info->remove (*it);

      // delete the node
      delete *it;
      *it = NULL;
    }

    _nodes.erase (nd, end);
  }


  assert (next_free_elem == _elements.size());
  assert (next_free_node == _nodes.size());

  // c++0x shrink_to_fit
  std::vector< Elem * >(_elements).swap(_elements);
  std::vector< Node * >(_nodes).swap(_nodes);

  STOP_LOG("renumber_nodes_and_elem()", "Mesh");
}


void SerialMesh::delete_remote_elements(bool volume_elem, bool surface_elem)
{
  if(Genius::n_processors() == 1) return;

  std::set<const Elem *> deleted_elems;
  for(unsigned int n=0; n<_elements.size(); ++n)
  {
    Elem * elem = _elements[n];
    if(!elem || elem->on_local()) continue;

    bool be = boundary_info->is_boundary_elem(elem);
    if( (volume_elem && !be) || (surface_elem && be) )
    {
      deleted_elems.insert(elem);
      //will also delete elem in boundary_info
      delete_elem(elem);
    }
  }

  //reset neighbor info, set remote neighbor pointer to NULL
  for(unsigned int n=0; n<_elements.size(); ++n)
  {
    Elem * elem = _elements[n];
    if( elem )
    {
      for (unsigned int ms=0; ms<elem->n_neighbors(); ms++)
      {
        if (elem->neighbor(ms) && deleted_elems.find(elem->neighbor(ms))!=deleted_elems.end())
        {
          elem->set_neighbor(ms, NULL);
        }
      }
    }
  }

  deleted_elems.clear();

  //statistic nodes
  std::set<const Node *> remaining_nodes;
  for(unsigned int n=0; n<_elements.size(); ++n)
  {
    Elem * elem = _elements[n];
    if( elem )
    {
      for (unsigned int md=0; md<elem->n_nodes(); md++)
      {
        remaining_nodes.insert(elem->get_node(md));
      }
    }
  }

  for(unsigned int n=0; n<_nodes.size(); ++n)
  {
    Node * node = _nodes[n];
    if( node && remaining_nodes.find(node) == remaining_nodes.end() )
    {
      // will also delete node in boundary_info
      delete_node(node);
    }
  }

  remaining_nodes.clear();

  _is_serial = false;

}


void SerialMesh::pack_nodes(std::vector<Real> & pts) const
{
  parallel_only();

  // collect local node information
  std::vector<Real> node_locations;
  std::vector<unsigned int> node_ids;
  node_locations.reserve(3*_nodes.size());
  node_ids.reserve(_nodes.size());
  for (unsigned int n=0; n<_nodes.size(); ++n)
  {
    const Node * node = _nodes[n];
    if( node && node->processor_id() == Genius::processor_id())
    {
      node_locations.push_back ( node->x() ); // x
      node_locations.push_back ( node->y() ); // y
      node_locations.push_back ( node->z() ); // z
      node_ids.push_back(n);
    }
  }

  // parallel allgather
  Parallel::allgather(node_ids);
  Parallel::allgather(node_locations);
  assert(node_ids.size() == _nodes.size());

  // write to packed vector
  pts.resize( node_locations.size() );
  for(unsigned int n=0; n<node_ids.size(); ++n)
  {
    pts[ 3*node_ids[n]+0 ] = node_locations[3*n+0];
    pts[ 3*node_ids[n]+1 ] = node_locations[3*n+1];
    pts[ 3*node_ids[n]+2 ] = node_locations[3*n+2];
  }

}


void SerialMesh::pack_egeds(std::vector< std::pair<unsigned int, unsigned int> > &edges) const
{
  parallel_only();

  std::set< std::pair<unsigned int, unsigned int> >  edges_set;
  for (unsigned int n=0; n<_elements.size(); ++n)
  {
    const Elem * elem = _elements[n];
    if(elem && elem->processor_id() == Genius::processor_id())
    {
      for(unsigned int e=0; e<elem->n_edges(); ++e)
      {
        std::pair<unsigned int, unsigned int> edge_nodes;
        elem->nodes_on_edge(e, edge_nodes);
        unsigned int node1_id = elem->get_node(edge_nodes.first)->id();
        unsigned int node2_id = elem->get_node(edge_nodes.second)->id();
        if(node1_id >  node2_id ) std::swap(node1_id, node2_id);
        edges_set.insert( std::make_pair(node1_id, node2_id) );
      }
    }
  }

  std::vector<unsigned int> buffer;
  std::set< std::pair<unsigned int, unsigned int> >::iterator it=  edges_set.begin();
  for(; it!=edges_set.end(); ++it)
  {
    buffer.push_back(it->first);
    buffer.push_back(it->second);
  }
  Parallel::allgather(buffer);

  for(unsigned int n=0; n<buffer.size();)
    edges_set.insert( std::make_pair(buffer[n++], buffer[n++]) );

  for(it=edges_set.begin(); it!=edges_set.end(); ++it)
  {
    edges.push_back(*it);
  }

}



void SerialMesh::pack_elems(std::vector<int> &conn) const
{
  parallel_only();

  std::vector<int> cell_info;
  std::vector<unsigned int> cell_ids;
  cell_info.reserve(24*_elements.size());// a little over kill
  cell_ids.reserve(_elements.size());
  for (unsigned int n=0; n<_elements.size(); ++n)
  {
    const Elem * elem = _elements[n];
    if(elem && elem->processor_id() == Genius::processor_id())
    {
      elem->pack_element(cell_info);
      cell_ids.push_back(n);
    }
  }

  // parallel allgather
  Parallel::allgather(cell_ids);
  Parallel::allgather(cell_info);
  assert(cell_ids.size() == _elements.size());

  // find the right order
  std::vector< std::pair<unsigned int, unsigned int> > cell_offset(cell_ids.size());
  unsigned int cnt = 0;
  for(unsigned int n=0; n<cell_ids.size(); ++n)
  {
    const ElemType elem_type    = static_cast<ElemType>(cell_info[cnt]);
    unsigned int pack_size = Elem::pack_size(elem_type);
    cell_offset[cell_ids[n]] = std::make_pair( cnt,  cnt+pack_size);
    cnt += pack_size;
    assert(cnt<=cell_info.size());
  }

  // write to packed vector
  conn.clear();
  conn.reserve( cell_info.size() );
  for(unsigned int n=0; n<cell_offset.size(); ++n)
  {
    for(unsigned int m=cell_offset[n].first; m < cell_offset[n].second; ++m)
      conn.push_back(cell_info[m]);
  }
}


void SerialMesh::pack_boundary_faces (std::vector<unsigned int> & el_id,
                                      std::vector<unsigned short int> & side_id,
                                      std::vector<short int> & bc_id) const
{
  this->boundary_info->build_on_processor_side_list (el_id, side_id, bc_id);

  Parallel::allgather(el_id);
  Parallel::allgather(side_id);
  Parallel::allgather(bc_id);
}


void SerialMesh::pack_boundary_egeds(std::vector< std::pair<unsigned int, unsigned int> > &edges) const
{
  parallel_only();

  std::vector<unsigned int> el_id;
  std::vector<unsigned short int> side_id;
  std::vector<short int> bc_id;
  this->boundary_info->build_on_processor_side_list (el_id, side_id, bc_id);

  std::set< std::pair<unsigned int, unsigned int> >  edges_set;
  for (unsigned int n=0; n<el_id.size(); ++n)
  {
    const Elem * elem = _elements[el_id[n]];
    if(elem && elem->processor_id() == Genius::processor_id())
    {
      AutoPtr<Elem> face = elem->build_side(side_id[n], false);

      // is an edge
      if(face->dim() == 1)
      {
        unsigned int node1_id = face->get_node(0)->id();
        unsigned int node2_id = face->get_node(1)->id();
        if(node1_id >  node2_id ) std::swap(node1_id, node2_id);
        edges_set.insert( std::make_pair(node1_id, node2_id) );
      }
      else
      {
        for(unsigned int e=0; e<face->n_edges(); ++e)
        {
          std::pair<unsigned int, unsigned int> edge_nodes;
          face->nodes_on_edge(e, edge_nodes);
          unsigned int node1_id = face->get_node(edge_nodes.first)->id();
          unsigned int node2_id = face->get_node(edge_nodes.second)->id();
          if(node1_id >  node2_id ) std::swap(node1_id, node2_id);
          edges_set.insert( std::make_pair(node1_id, node2_id) );
        }
      }
    }
  }

  std::vector<unsigned int> buffer;
  std::set< std::pair<unsigned int, unsigned int> >::iterator it=  edges_set.begin();
  for(; it!=edges_set.end(); ++it)
  {
    buffer.push_back(it->first);
    buffer.push_back(it->second);
  }
  Parallel::allgather(buffer);

  for(unsigned int n=0; n<buffer.size();)
    edges_set.insert( std::make_pair(buffer[n++], buffer[n++]) );

  for(it=edges_set.begin(); it!=edges_set.end(); ++it)
  {
    edges.push_back(*it);
  }
}

void SerialMesh::pack_boundary_nodes (std::vector<unsigned int> & node_id,
                                      std::vector<short int> & bc_id) const
{
  this->boundary_info->build_on_processor_node_list (node_id, bc_id);

  Parallel::allgather(node_id);
  Parallel::allgather(bc_id);
}


void SerialMesh::gather(unsigned int root_id)
{
  if(Genius::n_processors() == 1) return;

  assert(_is_prepared);

  {
    std::vector<Real> pts;
    pack_nodes(pts);

    std::vector<int> conn;
    pack_elems(conn);

    if(Genius::processor_id() == root_id)
      _unpack_mesh(pts, conn);
  }

  {
    std::vector<unsigned int>       el_id;
    std::vector<unsigned short int> side_id;
    std::vector<short int>          bc_id;
    pack_boundary_faces(el_id, side_id, bc_id);

    if(Genius::processor_id() == root_id)
      _unpack_bc_faces(el_id, side_id, bc_id);
  }

  {
    std::vector<unsigned int> node_id;
    std::vector<short int>    bc_id;
    pack_boundary_nodes(node_id, bc_id);

    if(Genius::processor_id() == root_id)
      _unpack_bc_nodes(node_id, bc_id);
  }

  if(Genius::processor_id() == root_id)
    _is_serial = true;
}


void SerialMesh::allgather()
{
  if(Genius::n_processors() == 1) return;

  assert(_is_prepared);

  {
    std::vector<Real> pts;
    pack_nodes(pts);

    std::vector<int> conn;
    pack_elems(conn);

    _unpack_mesh(pts, conn);
  }

  {
    std::vector<unsigned int>       el_id;
    std::vector<unsigned short int> side_id;
    std::vector<short int>          bc_id;
    pack_boundary_faces(el_id, side_id, bc_id);
    _unpack_bc_faces(el_id, side_id, bc_id);
  }

  {
    std::vector<unsigned int> node_id;
    std::vector<short int>    bc_id;
    pack_boundary_nodes(node_id, bc_id);
    _unpack_bc_nodes(node_id, bc_id);
  }

  _is_serial = true;
}


void SerialMesh::broadcast (unsigned int root_id)
{
  if(Genius::n_processors() == 1) return;

  assert(_is_prepared);
  assert(root_id < Genius::n_processors());

  if( Genius::processor_id() == root_id )
    assert( _is_serial );

  // mesh is serial in all the processor
  int serial = _is_serial;
  Parallel::min(serial);
  if(serial) return;

  // broadcast the pts vector
  std::vector<Real> pts;
  if (Genius::processor_id() == root_id)
  {
    pts.reserve (3*_nodes.size());

    for (unsigned int n=0; n<_nodes.size(); ++n)
    {
      assert(_nodes[n]);
      const Point& p = *_nodes[n];
      pts.push_back ( p(0) ); // x
      pts.push_back ( p(1) ); // y
      pts.push_back ( p(2) ); // z
    }
  }
  Parallel::broadcast (pts, root_id);

  // broadcast cells
  std::vector<int> conn;
  if (Genius::processor_id() == root_id)
  {
    for (unsigned int n=0; n<_elements.size(); ++n)
    {
      assert(_elements[n]);
      const Elem * elem = _elements[n];
      elem->pack_element(conn);
    }
  }
  Parallel::broadcast (conn, root_id);

  _unpack_mesh(pts, conn);

  // broadcast bcs
  {
    std::vector<unsigned int>       el_id;
    std::vector<unsigned short int> side_id;
    std::vector<short int>          bc_id;

    if (Genius::processor_id() == root_id)
      this->boundary_info->build_side_list (el_id, side_id, bc_id);

    unsigned int n_bcs = el_id.size();

    // Broadcast the number of bcs to expect from processor 0.
    Parallel::broadcast (n_bcs, root_id);

    // Only continue if we have element BCs
    if (n_bcs > 0)
    {
      // Broadcast the element identities
      Parallel::broadcast (el_id, root_id);

      // Broadcast the side ids for those elements
      Parallel::broadcast (side_id, root_id);

      // Broadcast the bc ids for each side
      Parallel::broadcast (bc_id, root_id);

      _unpack_bc_faces(el_id, side_id, bc_id);
    }
  }

  {
    std::vector<unsigned int> node_id;
    std::vector<short int>    bc_id;

    if (Genius::processor_id() == root_id)
      this->boundary_info->build_node_list (node_id, bc_id);

    unsigned int n_bcs = node_id.size();

    // Broadcast the number of bcs to expect from processor 0.
    Parallel::broadcast (n_bcs, root_id);

    // Only continue if we have nodal BCs
    if (n_bcs > 0)
    {
      // Allocate space, again on CPU 0 this should be a no-op.
      node_id.resize (n_bcs);
      bc_id.resize   (n_bcs);

      // Broadcast the node ids
      Parallel::broadcast (node_id, root_id);

      // Broadcast the bc ids for each side
      Parallel::broadcast (bc_id, root_id);

      _unpack_bc_nodes(node_id, bc_id);
    }
  }
  _is_serial = true;

}



void SerialMesh::_unpack_mesh (const std::vector<Real> &pts, const std::vector<int> &conn)
{

  // unpack nodes
  {
    for (unsigned int n=0; n<_nodes.size(); ++n)
    {
      if(_nodes[n] == NULL)
      {
        Point p( pts[3*n+0], pts[3*n+1], pts[3*n+2] );
        _nodes[n] = new Node(p, n);
      }
    }
  }

  // unpack elems
  {
    unsigned int cnt = 0;

    // This map keeps track of elements we've previously added to the mesh
    // to avoid O(n) lookup times for parent pointers.
    std::map<unsigned int, Elem*> top_elems;
    // This map saved elem neighbor information
    std::map<unsigned int, std::vector<int> > elem_neighbors;

    while (cnt < conn.size())
    {
      // Declare the element that we will add
      Elem* elem = NULL;

      // Unpack the element header
      const ElemType elem_type    = static_cast<ElemType>(conn[cnt++]);
      const unsigned int elem_PID = conn[cnt++];
      const int subdomain_ID      = conn[cnt++];
      const int self_ID           = conn[cnt++];

#ifdef ENABLE_AMR
      const int level             = conn[cnt++];
      const int p_level           = conn[cnt++];
      const Elem::RefinementState refinement_flag =  static_cast<Elem::RefinementState>(conn[cnt++]);
      const Elem::RefinementState p_refinement_flag = static_cast<Elem::RefinementState>(conn[cnt++]);
      const int parent_ID         = conn[cnt++];
      const int which_child       = conn[cnt++];

      if (parent_ID != -1) // Do a log(n) search for the parent
      {
        Elem* my_parent = top_elems.count(parent_ID) ? top_elems[parent_ID] : NULL;

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
      top_elems.insert(std::make_pair(self_ID,elem));

      // Assign the node
      for (unsigned int n=0; n<elem->n_nodes(); n++)
      {
        assert (cnt < conn.size());
        elem->set_node(n) = this->node_ptr (conn[cnt++]);
      }

      // save neighbor information
      std::vector<int> neighbors;
      for (unsigned int n=0; n<elem->n_sides(); n++)
      {
        assert (cnt < conn.size());
        neighbors.push_back(conn[cnt++]);
      }
      elem_neighbors.insert( std::make_pair(elem->id(), neighbors) );

      elem->prepare_for_fvm();
    } // end while cnt < conn.size

    // assign elems to _elements array
    for (unsigned int n=0; n<_elements.size(); ++n)
    {
      if(_elements[n] == NULL)
      {
        _elements[n] = top_elems[n];
        top_elems.erase(n);
      }
      assert(_elements[n]->id() == n);
    }

    // set neighbors
    for (unsigned int n=0; n<_elements.size(); ++n)
    {
      Elem * elem = _elements[n];
      assert(elem_neighbors.find(elem->id()) != elem_neighbors.end() );

      const std::vector<int> & neighbors = elem_neighbors.find(elem->id())->second;
      for (unsigned int s=0; s<elem->n_sides(); s++)
        elem->set_neighbor( s, neighbors[s] != -1 ?  _elements[neighbors[s]] : NULL);
    }

    // delete extra elems
    std::map<unsigned int, Elem*>::iterator it = top_elems.begin();
    for(; it != top_elems.end(); ++it)
      delete it->second;

  }

}


void SerialMesh::_unpack_bc_faces (const std::vector<unsigned int> &el_id,
                                   const std::vector<unsigned short int> &side_id,
                                   const std::vector<short int> &bc_id)
{
  for (unsigned int e=0; e<bc_id.size(); e++)
  {
    const Elem* elem = this->elem(el_id[e]);

    assert (elem != NULL);
    assert (side_id[e] < elem->n_sides());

    if( this->boundary_info->boundary_id(elem, side_id[e]) != bc_id[e])
      this->boundary_info->add_side (elem, side_id[e], bc_id[e]);
  }
}



void SerialMesh::_unpack_bc_nodes (const std::vector<unsigned int> &node_id,
                                   const std::vector<short int> &bc_id)
{
  for (unsigned int n=0; n<bc_id.size(); n++)
  {
    const Node* node = this->node_ptr (node_id[n]);

    assert (node != NULL);

    if(this->boundary_info->boundary_id(node) != bc_id[n])
      this->boundary_info->add_node (node, bc_id[n]);
  }
}

void SerialMesh::subdomain_graph(std::vector<std::vector<unsigned int> >& adjncy) const
{
  std::vector< std::set<unsigned int > > subdomain_neighbors(this->n_subdomains());

  for (unsigned int n=0; n<_elements.size(); ++n)
  {
    const Elem * elem = _elements[n];
    unsigned int subdomain_id = elem->subdomain_id();

    for(unsigned int m=0; m<elem->n_neighbors(); ++m)
    {
      const Elem * neighbor = elem->neighbor(m);
      if(!neighbor) continue;

      unsigned int neighbor_subdomain_id = neighbor->subdomain_id();
      if( neighbor_subdomain_id != subdomain_id )
      {
        subdomain_neighbors[subdomain_id].insert(neighbor_subdomain_id);
      }
    }
  }

  adjncy.clear();
  adjncy.resize(this->n_subdomains());
  for(unsigned int n=0; n<this->n_subdomains(); ++n)
  {
    adjncy[n].insert(adjncy[n].end(), subdomain_neighbors[n].begin(), subdomain_neighbors[n].end());
  }
}



size_t SerialMesh::memory_size() const
{
  size_t counter = sizeof(*this);

  counter += _elements.capacity()*sizeof(Elem *);
  for(unsigned int n=0; n<_elements.size(); ++n)
  {
    Elem * elem = _elements[n];
    if(!elem) continue;

    counter += Elem::memory_size(elem->type());
  }

  counter += _nodes.capacity()*sizeof(Node *);
  for(unsigned int n=0; n<_nodes.size(); ++n)
  {
    Node * node = _nodes[n];
    if(!node) continue;

    counter += sizeof(Node);
  }

  return counter;
}

