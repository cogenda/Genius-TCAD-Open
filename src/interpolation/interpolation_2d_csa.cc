#include <cassert>
#include <cmath>

#include "config.h"
#include "asinh.hpp"
#include "interpolation_2d_csa.h"
#include "parallel.h"

Interpolation2D_CSA::Interpolation2D_CSA()
{}


Interpolation2D_CSA::~Interpolation2D_CSA()
{
  clear();
}


void Interpolation2D_CSA::clear()
{
  std::map<int, CSA::csa *>::iterator it = field_map.begin();
  for(; it != field_map.end(); ++it)
    CSA::csa_destroy(it->second);
  csa_points.clear();
}



void Interpolation2D_CSA::add_scatter_data(const Point & point, int group, double value)
{
  CSA::point p;
  p.x = point.x();
  p.y = point.y();

  if(field_limit.find(group)==field_limit.end())
    field_limit[group] = std::make_pair(value, value);
  else
  {
    field_limit[group].first  = std::min(field_limit[group].first, value);
    field_limit[group].second = std::max(field_limit[group].second, value);
  }

  InterpolationType type = _interpolation_type[group];

  switch(type)
  {
  case Linear    : p.z = value; break;
  case SignedLog : p.z = (value>0 ? 1.0 : -1.0)*std::log(1.0+std::abs(value)); break;
  case Asinh     : p.z = boost::math::asinh(value); break;
  }

  csa_points[group].push_back(p);
}



void Interpolation2D_CSA::setup(int group)
{
  CSA::csa * field=CSA::csa_create();
  field_map[group] = field;
  CSA::csa_addpoints(field, csa_points[group].size(), &(csa_points[group][0]));
  CSA::csa_calculatespline(field);
#ifdef HAVE_FENV_H
  feclearexcept(FE_INVALID);
#endif
}


void Interpolation2D_CSA::broadcast(unsigned int root)
{
  std::vector<int> groups;
  std::map<int, std::vector<CSA::point> >::iterator it=csa_points.begin();
  for(; it!=csa_points.end(); ++it)
    groups.push_back(it->first);
  Parallel::broadcast(groups, root);

  for(unsigned int n=0; n<groups.size(); ++n)
  {
     std::vector<CSA::point> & csa_point = csa_points[groups[n]];

     std::vector<double> point_location_x;
     std::vector<double> point_location_y;
     std::vector<double> point_location_z;
     for(unsigned int i=0; i<csa_point.size(); ++i)
     {
       point_location_x.push_back(csa_point[i].x);
       point_location_y.push_back(csa_point[i].y);
       point_location_z.push_back(csa_point[i].z);
     }
     Parallel::broadcast(point_location_x,  root);
     Parallel::broadcast(point_location_y,  root);
     Parallel::broadcast(point_location_z,  root);

     csa_point.clear();

     for(unsigned int i=0; i<point_location_x.size(); ++i)
     {
       CSA::point p;
       p.x = point_location_x[i];
       p.y = point_location_y[i];
       p.z = point_location_z[i];
       csa_point.push_back(p);
     }

     std::pair<double, double> & limits = field_limit[groups[n]];
     Parallel::broadcast(limits.first,  root);
     Parallel::broadcast(limits.second, root);
  }
}


double Interpolation2D_CSA::get_interpolated_value(const Point & point, int group) const
{
  CSA::csa * field = field_map.find(group)->second;

  CSA::point pout;
  pout.x = point.x();
  pout.y = point.y();
  CSA::csa_approximatepoints(field, 1, &pout);

  InterpolationType type = _interpolation_type.find(group)->second;

  switch(type)
  {
  case Linear    : break;
  case SignedLog : pout.z = pout.z>0 ? exp(pout.z)-1 : 1-exp(-pout.z) ; break;
  case Asinh     : pout.z = sinh(pout.z); break;
  }

  double vmin = field_limit.find(group)->second.first;
  double vmax = field_limit.find(group)->second.second;
  if(pout.z<vmin) pout.z = vmin;
  if(pout.z>vmax) pout.z = vmax;

#ifdef HAVE_FENV_H
  feclearexcept(FE_INVALID);
#endif

  return pout.z;
}

