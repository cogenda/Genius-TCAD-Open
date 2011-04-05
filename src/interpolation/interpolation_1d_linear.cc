#include <cmath>

#include "genius_env.h"
#include "interpolation_1d_linear.h"
#include "asinh.hpp"
#include "parallel.h"


void Interpolation1D_Linear::clear()
{
  _coordinate.clear();
  _values.clear();
  _pre_computed_value.clear();
}


void Interpolation1D_Linear::add_scatter_data(const Point & point, int group, double value)
{
  if(_coordinate.size())
    genius_assert( point.x() > _coordinate[_coordinate.size()-1] );
  _coordinate.push_back(point.x());
  _values[group].push_back(value);
}


void Interpolation1D_Linear::setup(int group)
{
  const std::vector<double> & value = _values[group];
  std::vector<double> & y = _pre_computed_value[group];

  InterpolationType type = _interpolation_type[group];
  for(unsigned int u=0; u<value.size(); ++u)
  {
    switch(type)
    {
    case Linear    : y.push_back(value[u]); break;
    case SignedLog : y.push_back((value[u]>0 ? 1.0 : -1.0)*std::log(1.0+std::abs(value[u]))); break;
    case Asinh     : y.push_back(boost::math::asinh(value[u])); break;
    }
  }
}


void Interpolation1D_Linear::broadcast(unsigned int root)
{
  Parallel::broadcast(_coordinate, root);

  std::vector<int> groups;
  std::map<int, std::vector<double> >::iterator it=_values.begin();
  for(; it!=_values.end(); ++it)
    groups.push_back(it->first);
  Parallel::broadcast(groups, root);

  for(unsigned int n=0; n<groups.size(); ++n)
  {
     std::vector<double> & value = _values[groups[n]];
     Parallel::broadcast(value, root);
  }
}


double Interpolation1D_Linear::get_interpolated_value(const Point & point, int group) const
{
  int n = _coordinate.size();
  const std::vector<double> & t = _coordinate;
  const std::vector<double> & y = _pre_computed_value.find(group)->second;

  double tval = point.x();
  double yval = 0.0;

  {
    int begin=0, end=n-1;
    while(end-begin>2)
    {
      int mid = (begin+end)/2;
      if( tval < t[mid] )
        end = mid;
      else
        begin = mid;
    }

    int ival = n - 2;
    for (int i = begin; i < end; i++ )
    {
      if ( tval < t[i+1] )
      {
        ival = i;
        break;
      }
    }
    genius_assert(ival>=0 && ival <n);

    double dt = tval - t[ival];
    double h = t[ival+1] - t[ival];
    yval = y[ival]+ dt * ( y[ival+1] - y[ival] ) / h;
  }

  InterpolationType type = _interpolation_type.find(group)->second;
  switch(type)
  {
  case Linear    : break;
  case SignedLog : yval = yval>0 ? exp(yval)-1 : 1-exp(-yval) ; break;
  case Asinh     : yval = sinh(yval); break;
  }


  return  yval;
}


