// $Id: tree.cc,v 1.1 2008/05/22 14:13:24 gdiso Exp $

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
#include "tree.h"
#include "mesh_base.h"
#include "mesh_tools.h"
#include "boundary_info.h"


// ------------------------------------------------------------
// Tree class method

// constructor
template <unsigned int N>
Tree<N>::Tree (const MeshBase& m,
               const unsigned int tbs,
               const unsigned int mlevel,
               const Trees::BuildType bt) :
    TreeBase(m),
    root(m,tbs, mlevel),
    build_type(bt)
{
  // Set the root node bounding box equal to the bounding
  // box for the entire domain.
  root.set_bounding_box (mesh.bounding_box());


  if (build_type == Trees::NODES)
  {
    // Add all the nodes to the root node.  It will
    // automagically build the tree for us.
    MeshBase::const_node_iterator       it  = mesh.nodes_begin();
    const MeshBase::const_node_iterator end = mesh.nodes_end();

    for (; it != end; ++it)
      root.insert (*it);

    // Now the tree contains the nodes.
    // However, we want element pointers, so here we
    // convert between the two.
    std::vector<std::vector<const Elem*> > nodes_to_elem;

    MeshTools::build_nodes_to_elem_map (mesh, nodes_to_elem);
    root.transform_nodes_to_elements (nodes_to_elem);
  }

  else if (build_type == Trees::ELEMENTS)
  {
    // Add all active elements to the root node.  It will
    // automatically build the tree for us.
    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    for (; it != end; ++it)
    {
        root.insert (*it);
    }
  }

  else if (build_type == Trees::ELEMENTS_ON_BOUNDARY)
  {
    // Add all active elements on BOUNDARY to the root node.  It will
    // automatically build the tree for us.
    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();

    for (; it != end; ++it)
    {
      if((*it)->on_boundary())
        root.insert (*it);
    }
  }

  else if (build_type == Trees::SURFACE_ELEMENTS)
  {
    std::vector<unsigned int>       elems;
    std::vector<unsigned short int> sides;
    std::vector<short int>          bds;
    mesh.boundary_info->build_active_side_list (elems, sides, bds);

    for(unsigned int n=0; n<elems.size(); ++n)
    {
      const Elem * surface_elem = mesh.elem(elems[n])->build_side(sides[n], false).release();
      root.insert (surface_elem);
      surface_elems.push_back(surface_elem);
    }
  }
}


template <unsigned int N>
Tree<N>::~Tree()
{
  for(unsigned int n=0; n<surface_elems.size(); ++n)
    delete surface_elems[n];
  surface_elems.clear();
}


template <unsigned int N>
const Elem* Tree<N>::find_element(const Point& p) const
{
  return root.find_element(p);
}



template <unsigned int N>
const Elem * Tree<N>::hit_element(const Point & p, const Point & dir) const
{
  if(!hit_boundbox(p, dir)) return NULL;

  const Elem * elem = root.hit_element(p, dir);

#if defined(HAVE_FENV_H) && defined(DEBUG)
  feclearexcept(FE_ALL_EXCEPT);
#endif

  return elem;
}



template <unsigned int N>
bool Tree<N>::hit_boundbox(const Point & p, const Point & dir) const
{
  std::pair<double, double> t;
  bool result = root.hit_boundbox(p, dir, t);

#if defined(HAVE_FENV_H) && defined(DEBUG)
  feclearexcept(FE_ALL_EXCEPT);
#endif

  return result;
}


// ------------------------------------------------------------
// Explicit Instantiations
template class Tree<4>;
template class Tree<8>;









