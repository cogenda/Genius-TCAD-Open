#include <cassert>
#include <cmath>

#include "genius_common.h"
#include "asinh.hpp"
#include "interpolation_2d_nn.h"
#include "parallel.h"

#include "log.h"

Interpolation2D_NN::Interpolation2D_NN()
{}


Interpolation2D_NN::~Interpolation2D_NN()
{
  this->clear();
}


void Interpolation2D_NN::clear()
{
  std::map<int, DATA>::iterator it=field.begin();
  for(; it!=field.end(); ++it)
  {
    NN::delaunay_destroy( it->second.d );
    NN::lpi_destroy(it->second.li);
  }
  field.clear();
}



void Interpolation2D_NN::add_scatter_data(const Point & point, int group, double value)
{
  field[group].x.push_back(point.x());
  field[group].y.push_back(point.y());

  InterpolationType type = _interpolation_type[group];
  switch(type)
  {
  case Linear    : break;
  case SignedLog : value = (value>0 ? 1.0 : -1.0)*std::log(1.0+std::abs(value)); break;
  case Asinh     : value = boost::math::asinh(value); break;
  }

  field[group].f.push_back(value);
  field[group].n = field[group].f.size();
}



void Interpolation2D_NN::setup(int group)
{
  DATA & data = field[group];
  std::vector<NN::point> points;
  for( unsigned int i=0; i<data.n; ++i )
  {
    NN::point  p = { data.x[i], data.y[i], data.f[i] };
    points.push_back(p);
  }
  data.d = NN::delaunay_build(data.n, &points[0], 0, 0, 0, 0 );
  data.li = NN::lpi_build(data.d);
}


void Interpolation2D_NN::broadcast(unsigned int root)
{
  std::vector<int> groups;
  std::map<int, DATA>::iterator it=field.begin();
  for(; it!=field.end(); ++it)
    groups.push_back(it->first);
  Parallel::broadcast(groups, root);

  for(unsigned int n=0; n<groups.size(); ++n)
  {
    DATA & data = field[groups[n]];
    Parallel::broadcast(data.x, root);
    Parallel::broadcast(data.y, root);
    Parallel::broadcast(data.f, root);
  }
}


double Interpolation2D_NN::get_interpolated_value(const Point & point, int group) const
{
  const DATA & data = field.find(group)->second;
  NN::point  p = { point.x(), point.y(), 0 };
  NN::lpi_interpolate_point(data.li, &p);
  return p.z;
}

