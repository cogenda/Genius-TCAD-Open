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


#ifndef __object_tree_h__
#define __object_tree_h__

#include "tree_base.h"

// Forward Declarations
class MeshBase;
class Elem;
class LightThread;

/**
 * an octree/quad tree for fast ray-elem intersection test
 */
class ObjectTree
{
public:

  ObjectTree (const MeshBase& mesh);

  ~ObjectTree();

  bool is_octree() const
  { return _tree->n_leaf() == 8; }

  bool is_quadtree() const
  { return _tree->n_leaf() == 4; }

  /**
   * @return the first elem the light thread hit
   */
  const Elem * hit(const LightThread *) const;

  /**
   * @return the first elem the ray(p,d) hit
   */
  const Elem * hit(const Point & p, const Point & d) const;

  /**
   * @return true if the ray(p,d) hit the bounding box of the mesh
   */
  bool hit_boundbox(const Point & p, const Point & d) const;

private:

  /**
   * Pointer to our tree.  The tree is built at run-time
   * through \p init().  For servant PointLocators (not master),
   * this simply points to the tree of the master.
   */
  TreeBase* _tree;
};

#endif

