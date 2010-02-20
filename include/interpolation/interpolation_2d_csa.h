#ifndef __interpolation_2d_csa_h__
#define __interpolation_2d_csa_h__

#include <vector>

#include "interpolation_base.h"
#include "csa.h"

/**
 * 2D Interpolation by CSA
 */
class Interpolation2D_CSA   : public InterpolationBase
{
public:
  Interpolation2D_CSA();

  ~Interpolation2D_CSA();

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
  virtual void broadcast(unsigned int  root=0);

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

  std::map<int, CSA::csa *>  field_map;  //

  std::map<int, std::vector<CSA::point> > csa_points;

 /*
  * record the min/max value of the field. the interpolated value should be limited by this value
  */
  std::map<int, std::pair<double, double> > field_limit;
};


#endif

