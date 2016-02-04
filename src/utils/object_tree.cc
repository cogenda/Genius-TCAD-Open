/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/


#include "elem.h"
#include "mesh_base.h"
#include "mesh_tools.h"
#include "tree.h"
#include "object_tree.h"

ObjectTree::ObjectTree(const MeshBase& mesh, Trees::BuildType type)
{
  if (mesh.mesh_dimension() == 2)
    _tree = new Trees::QuadTree (mesh, 16, 10, type);
  else
    _tree = new Trees::OctTree  (mesh, 36, 10, type);
}


ObjectTree::~ObjectTree()
{
  delete this->_tree;
}



const Elem * ObjectTree::hit(const Point & p, const Point & d) const
{
  return this->_tree->hit_element(p, d);
}

bool ObjectTree::hit_boundbox(const Point & p, const Point & d) const
{
  return this->_tree->hit_boundbox(p, d);
}

bool ObjectTree::hit_domain(const Point & p, const Point & d, std::pair<double, double> &t) const
{
  return this->_tree->hit_domain(p, d, t);
}





