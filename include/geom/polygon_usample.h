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

#ifndef __polygon_uniform_sample_h__
#define __polygon_uniform_sample_h__

#include <vector>

#include "polygon.h"


typedef std::pair<unsigned int, unsigned int> Idx;


/**
 * build structured mesh which covers the given polygon
 */
class PolygonUSample : public Polygon
{
public:

  /**
   * build the polygon with given contour points
   */
  PolygonUSample( const std::vector<Point> & pts, const Real mesh_size);


  /**
   * @return the max index
   */
  Idx max_index() const;

  /**
   * @return true when p (coplanr and ) inside polygon
   */
  bool has_point(const Point &p) const;

  /**
   * get the index range of intersection area with given area (defined by p and d)
   * @return true when intersection happens
   */
  bool index_range(const Point &p, Real d, std::pair<Idx, Idx> &) const;

  /**
   * @return the in/out flag of given index pair
   */
  bool inside(const Idx & idx) const;

  /**
   * @return the location of given index pair
   */
  Point location(const Idx & idx) const;

private:


  Real _mesh_size;

  Point _project_norm;

  /**
   * cood used for projection
   */
  std::vector<unsigned int> _project_cood;

  /**
   * projected contour points
   */
  std::vector<Point> _projected_contour_points;

  bool _has_projected_point(const Point &projected_p) const;

  /// bounding box
  Point _bl_point, _tr_point;

  unsigned int _max_index_x;

  unsigned int _max_index_y;


  std::vector<bool> _inside;

};


#endif
