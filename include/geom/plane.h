// $Id: plane.h,v 1.3 2008/07/10 09:39:38 gdiso Exp $

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



#ifndef __plane_h__
#define __plane_h__

// C++ includes

// Local includes
#include "surface.h"


/**
 * This class defines a plane.
 *
 * @author Benjamin S. Kirk, 2002
 */

// ------------------------------------------------------------
// Plane class definition
class Plane : public Surface
{
public:

  /**
   * Dummy Constructor.
   */
  Plane ();

  /**
   * Constructs a plane containing point p with normal n.
   */
  Plane (const Point& p, const Point& n);

  /**
   * Constructs a plane containing the three points.  The
   * normal is determined in a counter-clockwise sense.  See
   * the create_from_three_points method for more details.
   */
  Plane (const Point& p0, const Point& p1, const Point& p2);

  /**
   * Copy-constructor.
   */
  Plane (const Plane& other_plane);

  /**
   * Destructor.  Does nothing at the moment.
   */
  ~Plane ();

  /**
   * Defines a plane containing point p with normal n.
   */
  void create_from_point_normal (const Point& p, const Point& n);

  /**
   * Defines a plane intersecting the three points
   * p0, p1, and p2.  The normal is constructed in a
   * counter-clockwise sense, i.e. (p1-p0)x(p2-p0);
   */
  void create_from_three_points (const Point& p0,
				 const Point& p1,
				 const Point& p2 );

  /**
   * Returns the inversed plane (with inversed norm)
   */
  inline Plane inverse() const;

  /**
   * Returns the point for the plane.
   */
  inline const Point& point () const;


  /**
   * Returns the normal for the plane.
   */
  inline const Point& normal () const;


  /**
   * Returns plane parameter as Ax+By+Cz+D=0
   * normalize as A^2+B^2+C^2=1
   */
  void plane_parameter(Real &A, Real &B, Real &C, Real &D) const;


  /**
   * Creates an XY plane located at z=zpos,
   */
  void xy_plane (const Real zpos=0.);

  /**
   * Creates an XZ plane located at y=ypos,
   */
  void xz_plane (const Real ypos=0.);

  /**
   * Creates an YZ plane located at x=xpos,
   */
  void yz_plane (const Real xpos=0.);

  /**
   * @return true for plane XY
   */
  bool is_xy_plane() const;

  /**
   * @return true for plane XZ
   */
  bool is_xz_plane() const;

  /**
   * @return true for plane YZ
   */
  bool is_yz_plane() const;

  /**
   * @return the signed distance of p to the plane
   */
  Real signed_distance (const Point& p) const;

  /**
   * @returns true if the point p is above the surface,
   * false otherwise.
   */
  bool above_surface (const Point& p) const;

  /**
   * @returns true if the point p is below the surface,
   * false otherwise.
   */
  bool below_surface (const Point& p) const;

  /**
   * @returns true if the point p is on the surface,
   * false otherwise.  Note that the definition of on
   * the surface really means "very close" to account
   * for roundoff error.
   */
  bool on_surface (const Point& p) const;

  /**
   * @returns the closest point on the surface to point p.
   */
  Point closest_point (const Point& p) const;

  /**
   * @returns a unit vector normal to the surface at
   * point p.
   */
  Point unit_normal (const Point& p=Point(0,0,0)) const;

  /**
   * @return true when given vector parallel to this plane
   */
  bool parallel_to( const Point& p ) const;

  /**
   * @return the plane - line segment intersection point
   * line segment should not parallel/on the plane
   */
  bool intersect_point(const Point& v1, const Point& v2, Point *result ) const;

  /**
   * @return the plane - ray intersection point by t
   * ray should not parallel/on the plane
   */
  bool intersect_point(const Point& p, const Point& dir, Real &t ) const;

private:


  /**
   *  The plane is defined by a point and a normal.
   */
  Point _point;

  /**
   *  The plane normal.
   */
  Point _normal;

};



// ------------------------------------------------------------
// Plane class inline members

Plane Plane::inverse() const
{
  return Plane(_point, -_normal);
}


const Point & Plane::point () const
{
  return _point;
}


const Point & Plane::normal () const
{
  return _normal;
}

#endif
