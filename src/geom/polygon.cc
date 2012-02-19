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

#include <map>

#include "polygon.h"


bool Polygon::valid() const
{
  if(_points.size() < 3) return false;
  // or we have just an degenerated polygon

  const Point n  = (_points[1] - _points[0]).cross(_points[2] - _points[0]);
  if( n.size() == 0.0 ) return false;

  Plane plane(_points[0], n.unit() );
  for(unsigned int n=0; n<_points.size(); ++n)
  {
    const Point & p = _points[n];
    if( std::abs(plane.signed_distance(p)) > 1e-8 )
      return false;
  }

  return true;
}


Point Polygon::norm() const
{
  const Point n  = (_points[1] - _points[0]).cross(_points[2] - _points[0]);
  if( n.size() == 0.0 ) return Point(0,0,0);

  return  n.unit();
}


Plane Polygon::plane() const
{
  return Plane(_points[0],  norm() );
}


Real Polygon::signed_area() const
{
  if(_points.size() < 3) return 0.0;

  Point cross;
  for(unsigned int n=0; n<_points.size(); ++n)
  {
    unsigned int v1 = n;
    unsigned int v2 = (n+1)%_points.size();
    cross += _points[v1].cross(_points[v2]);
  }

  return 0.5*norm()*cross;
}

/*
bool Polygon::has_point(const Point &p) const
{
  Point norm = this->norm();

  // project to 2D
  std::vector<unsigned int> project_cood;
  std::multimap<Real, unsigned int> sorter;
  sorter.insert( std::make_pair(norm[0], 0) );
  sorter.insert( std::make_pair(norm[1], 1) );
  sorter.insert( std::make_pair(norm[2], 2) );

  unsigned int min_cood = sorter.begin()->second;
  for(unsigned int i=0; i<=2; i++)
    if( i!= min_cood ) project_cood.push_back(i);
  assert(project_cood.size() ==2);

  std::vector<Point> projected_contour_points;
  Point projected_p(p[project_cood[0]], p[project_cood[1]], 0.0);
  for(unsigned int n=0; n<_points.size(); ++n)
  {
    projected_contour_points.push_back( Point(_points[n][project_cood[0]], _points[n][project_cood[1]], 0.0) );
  }

  // test
  unsigned int i, j;
  unsigned int nvert = projected_contour_points.size();
  Real testx = projected_p.x();
  Real testy = projected_p.y();
  bool c = false;
  for (i = 0, j = nvert-1; i < nvert; j = i++)
  {
    Real vi_x = projected_contour_points[i].x();
    Real vi_y = projected_contour_points[i].y();
    Real vj_x = projected_contour_points[j].x();
    Real vj_y = projected_contour_points[j].y();
    if ( ((vi_y>testy) != (vj_y>testy)) && (testx < (vj_x-vi_x) * (testy-vi_y) / (vj_y-vi_y) + vi_x) )
      c = !c;
  }
  return c;
}
*/


void Polygon::clip( const Plane &p )
{
  if(_points.size() < 3) return;

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
