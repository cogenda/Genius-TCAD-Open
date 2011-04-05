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


#ifndef __nn_locator_h__
#define __nn_locator_h__

#include <set>
#include <vector>
#include <functional>
#include <kdtree.hpp>

class Node;
class MeshBase;

/**
 * find the nearest nodes in given subdomain by point p
 */
// ------------------------------------------------------------
// NearestNodeLocator class definition
class NearestNodeLocator
{
public:
  /**
   * constructor, build kdtree for each subdomain
   */
  NearestNodeLocator(const MeshBase& mesh);

  ~NearestNodeLocator();

  /**
   * @return the mininal distance to nearest node in given region
   */
  Real distance_to_nearest_node(const Point &p, unsigned int subdomain) const;

  /**
   * @return the nearest nodes of point p inside given radius and in specified region
   */
  std::vector<const Node * > nearest_nodes(const Point &p, Real radius, unsigned int subdomain) const;

  /**
   * @return the nearest nodes of segment p1-p2 inside given radius and in specified region
   */
  std::vector<const Node * > nearest_nodes(const Point &p1, const Point &p2, Real radius, unsigned int subdomain) const;

  /**
   * @return the nearest nodes of segment p1-p2 inside given radius and in specified region
   */
  void nearest_nodes(const Point &p1, const Point &p2, Real radius, unsigned int subdomain, std::set<const Node * >& nns) const;

private:

  const MeshBase& _mesh;

  typedef KDTree::KDTree<3, const Node *, std::pointer_to_binary_function<const Node *, unsigned int, Real> > kdtree_type;

  /**
   * kdtree for each subdomain
   */
  std::vector<kdtree_type *> _kdtrees;
};


#endif

