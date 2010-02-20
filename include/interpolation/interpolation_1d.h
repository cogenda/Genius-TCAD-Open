#ifndef __interpolation_1d_h__
#define __interpolation_1d_h__

#include <vector>

#include "interpolation_base.h"

class Interpolation1D   : public InterpolationBase
{
public:
  Interpolation1D() {}

  ~Interpolation1D()
  { clear(); }

  /**
   * clear internal interpolation data
   */
  void clear();

  /**
   * build internal data structure
   */
  void setup(int group);

  /**
   * broadcast data to all the processor
   */
  virtual void broadcast(unsigned int root=0);

  /**
   * add the data with GROUP_ID group in 2D for interpolation
   * NOTE: only the first coordinate of point is used.
   */
  void add_scatter_data(const Point & point, int group, double value);

  /**
   * get interpolated value with GROUP_ID group in location point
   */
  double get_interpolated_value(const Point & point, int group) const;

private:

  /**
   * the coordinate should be sorted: coordinate[n+1]>coordinate[n]
   */
  std::vector<double> _coordinate;

  /**
   * the values
   */
  std::map<int, std::vector<double> > _values;

  /**
   * precompute value
   */
  std::map<int, std::vector<double> > _pre_computed_value;

  /**
   * precompute second order derivative of value
   */
  std::map<int, std::vector<double> > _pre_computed_value_pp;
};


#endif

