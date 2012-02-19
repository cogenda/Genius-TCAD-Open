#ifndef __interpolation_3d_qshep_h__
#define __interpolation_3d_qshep_h__

#if 0

#include <vector>

#include "interpolation_base.h"


/**
 * 3D Interpolation by qshep3d.f90
 */
class Interpolation3D_qshep   : public InterpolationBase
{
public:
  Interpolation3D_qshep();

  ~Interpolation3D_qshep();

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

  struct QSHEP
  {
    /// the number of nodes and associated data values
    int n;
    /// x location
    std::vector<double> x;
    /// y location
    std::vector<double> y;
    /// z location
    std::vector<double> z;
    /// f(x,y,z)
    std::vector<double> f;
    /// the number of data points to be used in the least squares fit
    int nq;
    /// the number of rows, columns, and planes in the cell grid defined
    int nr;
    /// the number of nodes within (and defining) the radii of influence R(K)
    int nw;
    /// nodal indices associated with cells
    std::vector<int> lcell;
    /// next-node indices
    std::vector<int> lnext;
    /// minimum nodal coordinates
    Point xyzmin;
    /// cell dimensions
    Point xyzdel;
    /// square root of the largest element in RSQ
    double rmax;
    /// the squares of the radii r(k) which enter into the weights W(K).
    std::vector<double> rsq;
    /// the coefficients for quadratic nodal function
    std::vector<double> A;
    /// error code
    int ier;
  };

  std::map<int, QSHEP> field;

 /*
  * record the min/max value of the field. the interpolated value should be limited by this value
  */
  std::map<int, std::pair<double, double> > field_limit;
};


#endif


#endif

