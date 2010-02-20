// $Id: point.h,v 1.2 2008/03/20 02:20:14 gdiso Exp $

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



#ifndef __point_h__
#define __point_h__

// C++ includes
#include <cmath>

// Local includes
#include "type_vector.h"




/**
 * A \p Point defines a location in DIM dimensional Real space.  Points
 * are always real-valued, even if the library is configured with
 * \p --enable-complex.
 *
 * \author Benjamin S. Kirk, 2003.
 */

class Point : public TypeVector<Real>
{
public:

  /**
   * Constructor.  By default sets all entries to 0.  Gives the point 0 in
   * \p DIM dimensions.
   */
  Point  (const Real x=0.,
          const Real y=0.,
          const Real z=0.);

  /**
   * Copy-constructor.
   */
  Point (const Point& p);

  /**
   * Copy-constructor.
   */
  Point (const TypeVector<Real>& p);

  /**
   * Empty.
   */
  virtual ~Point() {}

  //   /**
  //    * @returns a key associated with this point.  Useful for sorting.
  //    */
  //   unsigned int key() const;

  /**
   * get x location
   */
  Real x() const
  { return (*this)(0); }

  /**
   * get reference of x location
   */
  Real & x()
  { return (*this)(0); }


  /**
   * get y location
   */
  Real y() const
  { return (*this)(1); }

  /**
   * get reference of y location
   */
  Real & y()
  { return (*this)(1); }


  /**
   * get z location
   */
  Real z() const
  { return (*this)(2); }

  /**
   * get reference of z location
   */
  Real & z()
  { return (*this)(2); }

protected:


  /**
   * Make the derived class a friend
   */
  friend class Node;
};



//------------------------------------------------------
// Inline functions
inline
Point::Point (const Real x,
              const Real y,
              const Real z) :
    TypeVector<Real> (x,y,z)
{}



inline
Point::Point (const Point& p) :
    TypeVector<Real> (p)
{}



inline
Point::Point (const TypeVector<Real>& p) :
    TypeVector<Real> (p)
{}


/**
 *
 *  Return a positive value if the point pd lies below the
 *  plane passing through pa, pb, and pc; "below" is defined so
 *  that pa, pb, and pc appear in counterclockwise order when
 *  viewed from above the plane.  Returns a negative value if
 *  pd lies above the plane.  Returns zero if the points are
 *  coplanar.  The result is also a rough approximation of six
 *  times the signed volume of the tetrahedron defined by the
 *  four points.
 *
 */
inline Real orient3dfast(const Point &pa, const Point &pb,
                         const Point &pc, const Point &pd)
{

  Point pad = pa-pd;
  Point pbd = pb-pd;
  Point pcd = pc-pd;
  return pad*(pbd.cross(pcd));
}


/**
 *  Return a positive value if the point pe lies inside the
 *  sphere passing through pa, pb, pc, and pd; a negative value
 *  if it lies outside; and zero if the five points are
 *  cospherical.
 *
 *  Original edition:
 *      The points pa, pb, pc, and pd must be ordered
 *      so that they have a positive orientation (as defined by
 *      orient3d()), or the sign of the result will be reversed.
 *
 *  Now the return value is corrected by orient3d() function.
 */
inline Real inspherefast(const Point &pa, const Point &pb,
                         const Point &pc, const Point &pd,
                         const Point &pe)
{
  Real aex, bex, cex, dex;
  Real aey, bey, cey, dey;
  Real aez, bez, cez, dez;
  Real alift, blift, clift, dlift;
  Real ab, bc, cd, da, ac, bd;
  Real abc, bcd, cda, dab;

  Real orient3d = orient3dfast(pa, pb, pc, pd);
  genius_assert(std::abs(orient3d)>TOLERANCE);

  // we must compute 4th order det here...
  // the code should be cleaned later.
  aex = pa[0] - pe[0];
  bex = pb[0] - pe[0];
  cex = pc[0] - pe[0];
  dex = pd[0] - pe[0];
  aey = pa[1] - pe[1];
  bey = pb[1] - pe[1];
  cey = pc[1] - pe[1];
  dey = pd[1] - pe[1];
  aez = pa[2] - pe[2];
  bez = pb[2] - pe[2];
  cez = pc[2] - pe[2];
  dez = pd[2] - pe[2];

  ab = aex * bey - bex * aey;
  bc = bex * cey - cex * bey;
  cd = cex * dey - dex * cey;
  da = dex * aey - aex * dey;

  ac = aex * cey - cex * aey;
  bd = bex * dey - dex * bey;

  abc = aez * bc - bez * ac + cez * ab;
  bcd = bez * cd - cez * bd + dez * bc;
  cda = cez * da + dez * ac + aez * cd;
  dab = dez * ab + aez * bd + bez * da;

  alift = aex * aex + aey * aey + aez * aez;
  blift = bex * bex + bey * bey + bez * bez;
  clift = cex * cex + cey * cey + cez * cez;
  dlift = dex * dex + dey * dey + dez * dez;

  return orient3d*((dlift * abc - clift * dab) + (blift * cda - alift * bcd));
}

#endif // #define __point_h__
