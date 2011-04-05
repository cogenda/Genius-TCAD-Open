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


#include "polygon.h"


bool Polygon::valid() const
{
  if(_points.size() < 3) return false;
  // or we have just an degenerated polygon
  return true;
}


Point Polygon::norm() const
{
  if( !valid() ) return Point();

  const Point e0 = _points[1] - _points[0];
  const Point e1 = _points[2] - _points[0];
  const Point n  = e0.cross(e1);

  return  n.unit(true);
}


Real Polygon::signed_area() const
{
  if( !valid() ) return 0.0;

  Point cross;
  for(unsigned int n=0; n<_points.size(); ++n)
  {
    unsigned int v1 = n;
    unsigned int v2 = (n+1)%_points.size();
    cross += _points[v1].cross(_points[v2]);
  }

  return 0.5*norm()*cross;
}


void Polygon::clip( const Plane &p )
{
  if( !valid() ) return;

  unsigned int above = 0;
  unsigned int below = 0;

  // flag of point locator. 1 above, -1 below, 0 on
  std::vector<int> location(_points.size());

  // test for point location
  for(unsigned int n=0; n<_points.size(); ++n)
  {
    Real d = p.signed_distance(_points[n]);
    if( d > 1e-10 )
    {
      location[n] = 1;
      above++;
    }
    else
    {
      if( d < -1e-10 )
      {
        location[n] = -1;
        below++;
      }
      else
      {
        location[n] = 0;
      }
    }
  }

  // all the points are above/on plane
  if( below == 0 ) { return; }
  // all the points are below/on plane
  if( above == 0 ) { _points.clear(); return; }

  std::vector<Point> result;

  unsigned int previous = _points.size() - 1;
  for (unsigned int index = 0; index < _points.size(); index++)
  {
    int loc = location[index];
    if (loc == -1) // this point below the plane
    {
      if (location[previous] == 1) // previous point above the plane
      {
        const Point & v1 = _points[previous];
        const Point & v2 = _points[index];
        Point v;
        p.intersect_point(v1, v2, &v);
        result.push_back(v);
      }
    }
    else // above or on the plane
    {
      const Point& v1 = _points[index];
      if ((loc == 1) && (location[previous] == -1))
      {
        const Point & v2 = _points[previous];
        Point v;
        p.intersect_point(v1, v2, &v);
        result.push_back(v);
      }

      result.push_back(v1);
    }

    previous = index;
  }

  _points = result;

}


void Polygon::print(std::ostream& os) const
{
  os << "Polygon " << std::endl;
  for(unsigned int n=0; n<_points.size(); ++n)
    os << "  " <<_points[n];
  os << "  area = " << area() << std::endl;
  os << std::endl;
}
