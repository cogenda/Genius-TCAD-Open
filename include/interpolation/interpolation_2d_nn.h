#ifndef __interpolation_2d_nn_h__
#define __interpolation_2d_nn_h__

#include <vector>

#include "interpolation_base.h"
#include "nn.h"



class Interpolation2D_NN : public InterpolationBase
{
public:
  Interpolation2D_NN ();

  virtual ~Interpolation2D_NN();

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


  struct DATA
  {
    /// the number of nodes and associated data values
    unsigned int n;
    /// x location
    std::vector<double> x;
    /// y location
    std::vector<double> y;
    /// f(x,y,z)
    std::vector<double> f;
    /// delaunay diagram for Natural Neighbours
    NN::delaunay * d;
    /// linear interpolator
    NN::lpi      *li;
  };

  std::map<int, DATA> field;


};




#endif
