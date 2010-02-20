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

// ------------------------------------------------------------
// SerialMesh class member functions
SerialMesh::SerialMesh (unsigned int d) :
  UnstructuredMesh (d)
{
}


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
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() == i); // This will change soon

  return (*_nodes[i]);
}





Node& SerialMesh::node (const unsigned int i)
{
  if (i >= this->n_nodes())
    {
      std::cout << " i=" << i
		<< ", n_nodes()=" << this->n_nodes()
		<< std::endl;
      genius_error();
    }

  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() == i); // This will change soon

  return (*_nodes[i]);
}



const Node* SerialMesh::node_ptr (const unsigned int i) const
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




Node* & SerialMesh::node_ptr (const unsigned int i)
{
  assert (i < this->n_nodes());
  assert (_nodes[i] != NULL);
  assert (_nodes[i]->id() == i); // This will change soon

  return _nodes[i];
}




Elem* SerialMesh::elem (const unsigned int i) const
{
  assert (i < this->n_elem());
  assert (_elements[i] != NULL);
  assert (_elements[i]->id() == i); // This will change soon

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




static bool less_than( const Node * n1, const Node * n2 )
{
  return n1->id() < n2->id();
}



void SerialMesh::reorder_nodes()
{

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
    if( (*(*node_it))(0) < (*(*start_node))(0) &&
        (*(*node_it))(1) < (*(*start_node))(1) &&
	(*(*node_it))(2) < (*(*start_node))(2)
      )
      start_node = node_it;


  // the visited flag array
  std::vector<bool> visit_flag(n_nodes(), false);
  unsigned int new_index = 0;
  // a queue for Breadth-First Search
  std::queue<const Node *> Q;
  std::map<const Node *, unsigned int> new_order;

  // do Breadth-First Search
  // begin at this node.
  visit_flag[(*start_node)->id()] = true;
  Q.push( *start_node );

  while(!Q.empty())
  {
    const Node * current = Q.front();
    Q.pop();
    new_order.insert(std::pair<const Node *, unsigned int>(current, new_index++) );

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
  genius_assert( new_order.size()==n_nodes() );

  // ok, assign ordered index to each node
  for (node_it  = nodes_begin() ; node_it != nodes_end(); ++node_it)
      (*node_it)->set_id () = (*new_order.find(*node_it)).second;

  // if only reset the id, may cause many problems
  // as a result, we must sort the nodes by new ID
  std::sort( _nodes.begin(), _nodes.end(), less_than );

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



void SerialMesh::renumber_nodes_and_elements ()
{

  START_LOG("renumber_nodes_and_elem()", "Mesh");

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

  STOP_LOG("renumber_nodes_and_elem()", "Mesh");
}
