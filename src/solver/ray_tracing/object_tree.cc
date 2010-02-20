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
#include "ray_tracing/light_thread.h"
#include "ray_tracing/object_tree.h"

ObjectTree::ObjectTree(const MeshBase& mesh)
{
  // A 1D/2D mesh in 3D space needs special consideration.
  // If the mesh is planar XY, we want to build a QuadTree
  // to search efficiently.  If the mesh is truly a manifold,
  // then we need an octree
  bool is_planar_xy = false;

  // Build the bounding box for the mesh.  If the delta-z bound is
  // negligibly small then we can use a quadtree.
  {
    MeshTools::BoundingBox bbox = MeshTools::bounding_box(mesh);

    const Real
    Dx = bbox.second(0) - bbox.first(0),
         Dz = bbox.second(2) - bbox.first(2);

    if (std::abs(Dz/(Dx + 1.e-20)) < 1e-10)
      is_planar_xy = true;
  }

  if (is_planar_xy)
    _tree = new Trees::QuadTree (mesh, 12, Trees::BOUNDARY_ELEMENTS);
  else
    _tree = new Trees::OctTree  (mesh, 27, Trees::BOUNDARY_ELEMENTS);
}


ObjectTree::~ObjectTree()
{
  delete this->_tree;
}


const Elem * ObjectTree::hit(const LightThread *light) const
{
  return this->_tree->hit_element(light->start_point(), light->dir());
}


const Elem * ObjectTree::hit(const Point & p, const Point & d) const
{
  return this->_tree->hit_element(p, d);
}

bool ObjectTree::hit_boundbox(const Point & p, const Point & d) const
{
  return this->_tree->hit_boundbox(p, d);
}


