#include <cmath>

#include "genius_env.h"
#include "interpolation_1d.h"
#include "asinh.hpp"
#include "parallel.h"


void Interpolation1D::clear()
{
  _coordinate.clear();
  _values.clear();
  _pre_computed_value.clear();
  _pre_computed_value_pp.clear();
}


void Interpolation1D::add_scatter_data(const Point & point, int group, double value)
{
  if(_coordinate.size())
    genius_assert( point.x() > _coordinate[_coordinate.size()-1] );
  _coordinate.push_back(point.x());
  _values[group].push_back(value);
}


void Interpolation1D::setup(int group)
{
  int n = _coordinate.size();

  std::vector<double> a(3*n);
  std::vector<double>      b(n);

  const std::vector<double> & t = _coordinate;
  const std::vector<double> & value = _values[group];
  std::vector<double> & y = _pre_computed_value[group];
  std::vector<double> & ypp = _pre_computed_value_pp[group];

  InterpolationType type = _interpolation_type[group];
  for(unsigned int u=0; u<value.size(); ++u)
  {
    switch(type)
    {
    case Linear    : y.push_back(value[u]); break;
    case SignedLog : y.push_back((value[u]>0 ? 1.0 : -1.0)*std::log(1.0+value[u])); break;
    case Asinh     : y.push_back(boost::math::asinh(value[u])); break;
    }
  }

  //
  //  Set up the first equation.
  //
  {
    b[0] = 0.0;
    a[1+0*3] = 1.0;
    a[0+1*3] = 0.0;
  }

  //
  //  Set up the intermediate equations.
  //
  for (int i = 1; i < n-1; i++ )
  {
    b[i] = ( y[i+1] - y[i] ) / ( t[i+1] - t[i] ) - ( y[i] - y[i-1] ) / ( t[i] - t[i-1] );
    a[2+(i-1)*3] = ( t[i] - t[i-1] ) / 6.0;
    a[1+ i   *3] = ( t[i+1] - t[i-1] ) / 3.0;
    a[0+(i+1)*3] = ( t[i+1] - t[i] ) / 6.0;
  }
  //
  //  Set up the last equation.
  //
  {
    b[n-1] = 0.0;
    a[2+(n-2)*3] = 0.0;
    a[1+(n-1)*3] = 1.0;
  }

  //
  //  Solve the linear system.
  //
  if ( n == 2 )
  {
    ypp.resize(n);

    ypp[0] = 0.0;
    ypp[1] = 0.0;
  }
  else
  {
    ypp.resize(n);

    for (int i = 0; i < n; i++ )
    {
      ypp[i] = b[i];
    }

    for (int  i = 1; i < n; i++ )
    {
      double xmult = a[2+(i-1)*3] / a[1+(i-1)*3];
      a[1+i*3] = a[1+i*3] - xmult * a[0+i*3];
      ypp[i] = ypp[i] - xmult * ypp[i-1];
    }

    ypp[n-1] = ypp[n-1] / a[1+(n-1)*3];
    for (int  i = n-2; 0 <= i; i-- )
    {
      ypp[i] = ( ypp[i] - a[0+(i+1)*3] * ypp[i+1] ) / a[1+i*3];
    }
  }
}


void Interpolation1D::broadcast(unsigned int root)
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


double Interpolation1D::get_interpolated_value(const Point & point, int group) const
{
  int n = _coordinate.size();
  const std::vector<double> & t = _coordinate;
  const std::vector<double> & y = _pre_computed_value.find(group)->second;
  const std::vector<double> & ypp = _pre_computed_value_pp.find(group)->second;

  double tval = point.x();

  //
  //  Determine the interval [ double(I), double(I+1) ] that contains TVAL.
  //  Values below double[0] or above double[N-1] use extrapolation.
  //
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

  //
  //  In the interval I, the polynomial is in terms of a normalized
  //  coordinate between 0 and 1.
  //
  double dt = tval - t[ival];
  double h = t[ival+1] - t[ival];

  double yval = y[ival]+ dt * ( ( y[ival+1] - y[ival] ) / h
                                - ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * h
                                + dt * ( 0.5 * ypp[ival] + dt * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0 * h ) ) ) );
  InterpolationType type = _interpolation_type.find(group)->second;

  switch(type)
  {
  case Linear    : break;
  case SignedLog : yval = yval>0 ? exp(yval)-1 : 1-exp(-yval) ; break;
  case Asinh     : yval = sinh(yval); break;
  }


  return  yval;
}


