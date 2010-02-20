#ifndef __interpolation_base_h__
#define __interpolation_base_h__

#include <cassert>
#include <map>
#include <string>

#include "point.h"

/**
 * this class offers an unique interface to various interpolation method
 * in 1D/2D/3D.
 */
class InterpolationBase
{
public:
  InterpolationBase() {}

  virtual ~InterpolationBase()
  {
    _interpolation_type.clear();
    _variable_group_map.clear();
  }

  /**
   * clear internal interpolation data
   */
  virtual void clear()=0;

  /**
   * give the variable string a group code
   */
  int set_group_code(const std::string & var)
  {
    if(_variable_group_map.find(var)!=_variable_group_map.end())
      return _variable_group_map.find(var)->second;
    int size = _variable_group_map.size();
    _variable_group_map.insert(std::make_pair(var, size));
    return size;
  }

  /**
   * get the group code by variable string
   */
  int group_code(const std::string & var) const
  { return _variable_group_map.find(var)->second; }


  /**
   * add the data with GROUP_ID group in (N+1)D for interpolation
   */
  virtual void add_scatter_data(const Point & point, int group, double value)=0;

  /**
   * build internal data structure
   */
  virtual void setup(int /* group */)=0;

  /**
   * virtual function to set the options of Interpolation.
   * each dirived class can override it
   */
  virtual void set_option(const std::string & /* option */)  {/* do nothing */ }

  /**
   * get interpolated value with GROUP_ID group in location point
   */
  virtual double get_interpolated_value(const Point & point, int group)const=0;

  /**
   * InterpolationType, should support linear (for potential, etc) and asinh (doping concentration and carrier density)
   */
  enum InterpolationType{Linear, SignedLog, Asinh};

  /**
   * broadcast data to all the processor
   */
  virtual void broadcast(unsigned int=0)=0;

  /**
   * how data interpolated to given point, supported type can be linear, signedlog and asinh
   */
  void set_interpolation_type(int group, InterpolationType type)
  { _interpolation_type[group] = type; }

protected:

  std::map<int, InterpolationType> _interpolation_type;

  std::map<std::string, int> _variable_group_map;

  //how to store the point and their value?

  inline double scaleValue(InterpolationType type, const double value) const
  {
    switch(type)
    {
      case Linear    : return value;
      case SignedLog : return (value>0 ? 1.0 : -1.0)*std::log(1.0+std::abs(value));
      case Asinh     : return boost::math::asinh(value);
    }
    return 0.0; //prevent warning
  }

  inline double unscaleValue(InterpolationType type, const double value) const
  {
    switch(type)
    {
      case Linear    : return value;
      case SignedLog : return value>0 ? exp(value)-1 : 1-exp(-value) ;
      case Asinh     : return sinh(value);
    }
    return 0.0;//prevent warning
  }

};


#endif

