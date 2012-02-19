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

#ifndef __polygon_h__
#define __polygon_h__

#include <vector>

#include "point.h"
#include "plane.h"

/**
 * A \p Polygon defines a 3D space plane in figure that is bounded by a closed path
 * only support simple convex polygon
 */
class Polygon
{
public:

  /**
   * Constructor, build empty polygon
   */
  Polygon() {}

  /**
   * build the polygon with given contour points
   */
  Polygon( const std::vector<Point> & pts)
      : _points(pts)
  {}

  virtual ~Polygon() {}

  /**
   * set the contour points of polygon
   */
  void set_contour_points( const std::vector<Point> & pts )
  { _points = pts; }

  /**
   * add point to polygon contour points linst
   */
  void add_contour_point( const Point & pt )
  { _points.push_back(pt); }

  /**
   * @return true when polygon valid
   */
  bool valid() const;

  /**
   * @return the vertex number
   */
  unsigned int n_vertex() const { return _points.size(); }

  /**
   * the norm of the polygon
   */
  Point norm() const;

  /**
   * return the plane of the polygon
   */
  Plane plane() const;

  /**
   * @return the signed area of polygon
   */
  Real signed_area() const;

  /**
   * @return the area of polygon
   */
  Real area() const { return std::abs(signed_area()); }

  /**
   * @return true when p (coplanr and ) inside polygon
   */
  //bool has_point(const Point &p) const;

  /**
   * clip the polygon by given plane, any part below the plane will be droped
   */
  void clip( const Plane &);

  /**
   * debug output
   */
  void print(std::ostream& os) const;

  /**
   * debug output
   */
  friend std::ostream& operator << (std::ostream& os, const Polygon & poly)
  {
    poly.print(os);
    return os;
  }

protected:

  /**
   * polygon contour points, not closed
   */
  std::vector<Point> _points;


};


#endif
