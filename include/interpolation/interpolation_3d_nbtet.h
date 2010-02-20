#ifndef __interpolation_3d_nb_h__
#define __interpolation_3d_nb_h__

#include <vector>
#include "ANN/ANN.h"

#include "interpolation_base.h"

typedef std::vector<std::pair<int, double> > ANNResultType; // array of (index, dist) pairs

class ANNSession
{
public:
  ANNSession(int dim);
  ~ANNSession();

  void clear();

  void broadcast(unsigned int root);

  void setup();

  bool addPoint(const Point &pt);

  ANNResultType search(const Point &pt, const unsigned int k) const;

  Point getPointCoord(unsigned int i) const;

  long size() const;

private:
  int _dim;
  bool _tree_built;

  std::vector<Point> _pts;

  ANNpointArray   _datapts;   // hold pointers to points
  ANNpoint        _databuf;   // hold coordinates of points

  ANNkd_tree*     _kdTree;

};

class Interpolation3D_nbtet : public InterpolationBase
{
public:
  Interpolation3D_nbtet ();

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
  ANNSession _ann;
  std::vector<double> _field;

};




#endif
