/******************************************************************************
 *
 * File:           nn.h
 *
 * Created:        04/08/2000
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        Header file for nn library
 *
 * Description:    None
 *
 * Revisions:      None
 *
 *****************************************************************************/

#if !defined(_NN_H)
#define _NN_H

namespace NN
{

  typedef enum { SIBSON, NON_SIBSONIAN } NN_RULE;

  /* "point" is a basic data structure in this package.
   */
#if !defined(_POINT_STRUCT)
#define _POINT_STRUCT
  typedef struct
  {
    double x;
    double y;
    double z;
  }
  point;
#endif

  /* Constructors for interpolators in this package require Delaunay
   * triangulation of the input data.
   */
#if !defined(_DELAUNAY_STRUCT)
#define _DELAUNAY_STRUCT
  struct delaunay;
  typedef struct delaunay delaunay;
#endif

  /** Builds Delaunay triangulation of the given array of points.
   *
   * @param np Number of points
   * @param points Array of points [np] (input)
   * @param ns Number of forced segments
   * @param segments Array of (forced) segment endpoint indices [2*ns]
   * @param nh Number of holes
   * @param holes Array of hole (x,y) coordinates [2*nh]
   * @return Delaunay triangulation structure with triangulation results
   */
  delaunay* delaunay_build(int np, point points[], int ns, int segments[], int nh, double holes[]);

  /** Destroys Delaunay triangulation.
   *
   * @param d Structure to be destroyed
   */
  void delaunay_destroy(delaunay* d);


  /** `lpi' -- "Linear Point Interpolator" is a structure for linear
   ** interpolation of data on a "point-to-point" basis.
   *
   * `lpi' interpolates linearly within each triangle resulted from the Delaunay
   * triangluation of the input data. `lpi' is much faster than all Natural
   * Neighbours interpolators below.
   */
  struct lpi;
  typedef struct lpi lpi;

  /** Builds linear interpolator.
   *
   * @param d Delaunay triangulation
   * @return Linear interpolator
   */
  lpi* lpi_build(delaunay* d);

  /** Destroys linear interpolator.
   *
   * @param l Structure to be destroyed
   */
  void lpi_destroy(lpi* l);

  /** Finds linearly interpolated value in a point.
   *
   * @param l Linear point interpolator
   * @param p Point to be interpolated (p->x, p->y -- input; p->z -- output)
   */
  void lpi_interpolate_point(lpi* l, point* p);

  /** Linearly interpolates data in an array of points.
   *
   * @param nin Number of input points
   * @param pin Array of input points [pin]
   * @param nout Number of ouput points
   * @param pout Array of output points [nout]
   */
  void lpi_interpolate_points(int nin, point pin[], int nout, point pout[]);


}

#endif                          /* _NN_H */
