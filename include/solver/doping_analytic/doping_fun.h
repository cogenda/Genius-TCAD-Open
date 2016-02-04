/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/

//  $Id: doping_fun.h,v 1.2 2008/07/09 05:58:16 gdiso Exp $

#ifndef __doping_fun_h__
#define __doping_fun_h__

#include "config.h"

#include <vector>

#include "point.h"
#include "plane.h"
#include "enum_direction.h"


enum PROFILE_CARD_ERROR
{
  PROFILE_NO_ERROR,
  PROFILE_UNKNOW_PARAMETER,
  PROFILE_ION_TYPE,
  PROFILE_XLOCATION,
  PROFILE_YLOCATION,
  PROFILE_NEGTIVE_DOSE,
  PROFILE_NEGTIVE_CONCENTR,
  PROFILE_DOPING_CHAR_LENGTH
};





/**
 * the base class for analytic doping function
 */
class DopingFunction
{
public:
  /**
   * constructor, do nothing
   */
  DopingFunction(double ion)
      : _ion(ion)
  {}

  /**
   * destructor, do nothing
   */
  virtual ~DopingFunction() {}

  /**
   * virtual compute function
   */
  virtual double profile(double x, double y, double z)=0;

protected:
  /**
   * impurity ion type N-ion or P-ion
   */
  double _ion;


};



//------------------------------------------------------------------



/**
 * uniform doping function
 */
class UniformDopingFunction : public DopingFunction
{
public:

  /**
   * constructor
   */
  UniformDopingFunction(double ion, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double N)
  : DopingFunction(ion), _peak(N), _xmin(xmin), _xmax(xmax), _ymin(ymin), _ymax(ymax), _zmin(zmin), _zmax(zmax)
  {}

  /**
   * compute doping concentration by point location
   */
  double profile(double x,double y,double z);

private:
  /**
   * the peak value of doping concentration
   */
  double _peak;

  /**
   * bound box of doping region
   */
  double _xmin;    // bound box

  /**
   * bound box of doping region
   */
  double _xmax;

  /**
   * bound box of doping region
   */
  double _ymin;

  /**
   * bound box of doping region
   */
  double _ymax;

  /**
   * bound box of doping region
   */
  double _zmin;

  /**
   * bound box of doping region
   */
  double _zmax;
};


//------------------------------------------------------------------


/**
 * linear doping function
 */
class LinearDopingFunction : public DopingFunction
{
public:

  /**
   * constructor
   */
  LinearDopingFunction(double ion, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                       double doping_begin, double doping_slope, double base_plane, double end_plane, Direction dir)
      : DopingFunction(ion),   _xmin(xmin), _xmax(xmax), _ymin(ymin), _ymax(ymax), _zmin(zmin), _zmax(zmax),
      _doping_begin(doping_begin), _doping_slope(doping_slope), _base_plane(base_plane), _end_plane(end_plane), _dir(dir)
  {}

  /**
   * destructor, do nothing
   */
  virtual ~LinearDopingFunction() {}

  /**
   * compute profile
   */
  double profile(double x,double y,double z);


private:

  
  /**
   * bound box of doping region
   */
  double _xmin;    // bound box

  /**
   * bound box of doping region
   */
  double _xmax;

  /**
   * bound box of doping region
   */
  double _ymin;

  /**
   * bound box of doping region
   */
  double _ymax;

  /**
   * bound box of doping region
   */
  double _zmin;

  /**
   * bound box of doping region
   */
  double _zmax;
  
  /**
   * doping fraction at the _base_plane
   */
  double _doping_begin;

  /**
   * the slope of the doping fraction for graded compounds
   */
  double _doping_slope;

  /**
   * the reference base plane, at which the doping fraction is _doping_begin
   */
  double _base_plane;

  /**
   * the end plane, between which the doping fraction is _doping_begin + x*_doping_slope
   */
  double _end_plane;

  /**
   * the doping fraction grading direction
   */
  Direction _dir;

};


//------------------------------------------------------------------



/**
 * analytic (gauss or erf) doping function
 */
class AnalyticDopingFunction : public DopingFunction
{
public:

  /**
   * constructor
   * @param XERFC true to use erfc distribution in x direction
   * @param ZERFC true to use erfc distribution in z direction
   */
  AnalyticDopingFunction(double ion, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double N,
                         double ax, double ay, double az, bool XERFC=false, bool YERFC=false, bool ZERFC=false)
  : DopingFunction(ion), _xmin(xmin), _xmax(xmax), _ymin(ymin), _ymax(ymax), _zmin(zmin), _zmax(zmax),
      _peak(N), _XCHAR(ax), _YCHAR(ay), _ZCHAR(az), _XERFC(XERFC), _YERFC(YERFC), _ZERFC(ZERFC)
  {}

  /**
   * compute doping concentration by point location
   */
  double profile(double x, double y, double z);

private:
  /**
   * the peak value of doping concentration
   */
  double _peak;

  /**
   * bound box of doping region
   */
  double _xmin;    // bound box

  /**
   * bound box of doping region
   */
  double _xmax;

  /**
   * bound box of doping region
   */
  double _ymin;

  /**
   * bound box of doping region
   */
  double _ymax;

  /**
   * bound box of doping region
   */
  double _zmin;

  /**
   * bound box of doping region
   */
  double _zmax;

  /**
   * characteristic length along x direction
   */
  double _XCHAR;

  /**
   * characteristic length along y direction
   */
  double _YCHAR;

  /**
   * characteristic length along z direction
   */
  double _ZCHAR;

  /**
   * when it is true, erfc distribution for x direction is used instead of gauss distribution
   */
  bool   _XERFC;

  /**
   * when it is true, erfc distribution for y direction is used instead of gauss distribution
   */
  bool   _YERFC;

  /**
   * when it is true, erfc distribution for z direction is used instead of gauss distribution
   */
  bool   _ZERFC;

};



//------------------------------------------------------------------
#include "polygon_usample.h"
/**
 * poly mask based doping function
 */
class PolyMaskDopingFunction : public DopingFunction
{
public:

  /**
   * constructor
   */
  PolyMaskDopingFunction(double ion, const std::vector<Point> &poly, double theta, double phi,
                     double rmin, double rmax, double npeak, double char_depth, double char_lateral,
                     double resolution_factor=4.0);

  /**
   * compute doping concentration by point location
   */
  double profile(double x, double y, double z);

private:

  /**
   * direction of profile line
   */
  Point _dir;

  /**
   * min short range of ion, should be positive
   */
  double _rmin;

  /**
   *  max short range of ion, should be positive
   */
  double _rmax;

  /**
   * the peak value of doping concentration
   */
  double _peak;

  /**
   * characteristic length along ion  direction
   */
  double _char_depth;

  /**
   * lateral characteristic length
   */
  double _char_lateral;

  /**
   * control the density of doping line
   */
  double _resolution_factor;


  PolygonUSample _mask_mesh;


  std::vector<double> prof_func_r(const std::vector<double> &rs) const;


  std::vector<double> prof_func_l(const std::vector<double> &dists) const;


  double prof_func_lnorm() const;


};



/**
 * rectangle mask based doping function
 */
class RecMaskDopingFunction : public DopingFunction
{
public:

  /**
   * constructor
   */
  RecMaskDopingFunction(double ion, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double theta, double phi,
                     double rmin, double rmax, double npeak, double char_depth, double char_lateral, double resolution_factor=4.0)
  : DopingFunction(ion),
    _rmin(rmin), _rmax(rmax), _peak(npeak), _char_depth(char_depth), _char_lateral(char_lateral), _resolution_factor(resolution_factor)
  {
    _init( xmin,  xmax,  ymin,  ymax,  zmin,  zmax,  theta,  phi);
  }

  /**
   * compute doping concentration by point location
   */
  double profile(double x, double y, double z);

private:

  /**
   * mask plane
   */
  Plane _mask;

  /**
   * mask box
   */
  Point _mask_min;

  /**
   * mask box
   */
  Point _mask_max;

  /**
   * direction of profile line
   */
  Point _dir;

  /**
   * min short range of ion, should be positive
   */
  double _rmin;

  /**
   *  max short range of ion, should be positive
   */
  double _rmax;

  /**
   * the peak value of doping concentration
   */
  double _peak;

  /**
   * characteristic length along ion  direction
   */
  double _char_depth;

  /**
   * lateral characteristic length
   */
  double _char_lateral;

  /**
   * doping bounding box
   */
  Point _doping_min;

  /**
   * doping bounding box
   */
  Point _doping_max;

  /**
   * control the density of doping line
   */
  double _resolution_factor;

  void _init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double theta, double phi);

  double _profile_r(double r) const;

};

#endif
