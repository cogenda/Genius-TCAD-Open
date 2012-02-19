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

#include "polygon_usample.h"



PolygonUSample::PolygonUSample(const std::vector<Point> & pts, const Real mesh_size)
  :Polygon(pts)
{
  Point norm = this->norm();

  // project to 2D
  std::multimap<Real, unsigned int> sorter;
  sorter.insert( std::make_pair(std::abs(norm[0]), 0) );
  sorter.insert( std::make_pair(std::abs(norm[1]), 1) );
  sorter.insert( std::make_pair(std::abs(norm[2]), 2) );

  unsigned int min_cood = sorter.rbegin()->second;
  for(unsigned int i=0; i<=2; i++)
    if( i!= min_cood ) _project_cood.push_back(i);
  assert(_project_cood.size() == 2);
  _project_cood.push_back(min_cood);
  _project_norm[min_cood] = norm[min_cood];

  for(unsigned int n=0; n<_points.size(); ++n)
  {
    Point p(_points[n][_project_cood[0]], _points[n][_project_cood[1]], _points[n][_project_cood[2]]);
    _projected_contour_points.push_back( p );
  }

  _mesh_size = mesh_size * _project_norm.size() / norm.size();

  // build bounding box in XY
  _bl_point = _projected_contour_points[0];
  _tr_point = _projected_contour_points[0];

  for(unsigned int n=1; n<_projected_contour_points.size(); ++n)
  {
    Point point = _projected_contour_points[n];

    _bl_point.x() = std::min(_bl_point.x() , point.x() );
    _bl_point.y() = std::min(_bl_point.y() , point.y() );

    _tr_point.x() = std::max(_tr_point.x() , point.x() );
    _tr_point.y() = std::max(_tr_point.y() , point.y() );
  }


  _max_index_x = static_cast<unsigned int>(ceil((_tr_point.x() - _bl_point.x())/_mesh_size));
  _max_index_y = static_cast<unsigned int>(ceil((_tr_point.y() - _bl_point.y())/_mesh_size));

  _tr_point = _bl_point + Point(_max_index_x*_mesh_size, _max_index_y*_mesh_size);

  _inside.resize( (_max_index_x+1) * (_max_index_y+1) );
  for(unsigned int ix=0; ix<=_max_index_x; ++ix)
    for(unsigned int iy=0; iy<=_max_index_y; ++iy)
    {
      Point loc = _bl_point + Point( ix*_mesh_size, iy*_mesh_size  );
      _inside[ix*(_max_index_y+1)+iy] = _has_projected_point(loc);
    }
}



bool PolygonUSample::has_point(const Point &p) const
{
  Point projected_p(p[_project_cood[0]], p[_project_cood[1]], p[_project_cood[2]]);
  return _has_projected_point(projected_p);
}



bool PolygonUSample::index_range(const Point &p, Real d, std::pair<Idx, Idx> & ids) const
{
  Point projected_p(p[_project_cood[0]], p[_project_cood[1]], p[_project_cood[2]]);
  Point b_low(projected_p[0]-d, projected_p[1]-d, projected_p[2]);
  Point b_high(projected_p[0]+d, projected_p[1]+d, projected_p[2]);

  Point bl(std::max(b_low.x(), _bl_point.x()), std::max(b_low.y(), _bl_point.y()) );
  Point tr(std::min(b_high.x(), _tr_point.x()), std::min(b_high.y(), _tr_point.y()) );

  if( bl.x() >  _tr_point.x() || bl.y() >  _tr_point.y() ) return false;
  if( tr.x() <  _bl_point.x() || tr.y() <  _bl_point.y() ) return false;

  ids.first.first   = static_cast<unsigned int>(floor((bl.x() - _bl_point.x())/_mesh_size));
  ids.first.second  = static_cast<unsigned int>(floor((bl.y() - _bl_point.y())/_mesh_size));
  ids.second.first  = _max_index_x - static_cast<unsigned int>(floor((_tr_point.x() - tr.x())/_mesh_size));
  ids.second.second = _max_index_y - static_cast<unsigned int>(floor((_tr_point.y() - tr.y())/_mesh_size));

  assert(ids.second.first >= ids.first.first);
  assert(ids.second.second >= ids.first.second);
  return true;
}



Idx PolygonUSample::max_index() const
{
  return std::make_pair(_max_index_x, _max_index_y);
}

bool PolygonUSample::inside(const Idx & idx) const
{
  return _inside[ idx.first*(_max_index_y+1) + idx.second ];
}

Point PolygonUSample::location(const Idx & idx) const
{
  Point projected_point = _bl_point + Point( idx.first*_mesh_size, idx.second*_mesh_size  );
  Point p;
  p[_project_cood[0]] = projected_point[0];
  p[_project_cood[1]] = projected_point[1];
  p[_project_cood[2]] = projected_point[2];

  Real t;
  Plane plane = this->plane();
  plane.intersect_point(p, _project_norm, t);

  //std::cout<<projected_point;
  //std::cout<< p + _project_norm*t;
  return p + _project_norm*t;
}


bool PolygonUSample::_has_projected_point(const Point &projected_p) const
{
  // test
  unsigned int i, j;
  unsigned int nvert = _projected_contour_points.size();
  Real testx = projected_p.x();
  Real testy = projected_p.y();
  bool c = false;
  for (i = 0, j = nvert-1; i < nvert; j = i++)
  {
    Real vi_x = _projected_contour_points[i].x();
    Real vi_y = _projected_contour_points[i].y();
    Real vj_x = _projected_contour_points[j].x();
    Real vj_y = _projected_contour_points[j].y();
    if ( ((vi_y>testy) != (vj_y>testy)) && (testx < (vj_x-vi_x) * (testy-vi_y) / (vj_y-vi_y) + vi_x) )
      c = !c;
  }
  return c;
}

