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
//win32 does not have erfc function
#ifdef CYGWIN
  #include "mathfunc.h"
#endif

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
  DopingFunction(double ion, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
      : _ion(ion), _xmin(xmin), _xmax(xmax), _ymin(ymin), _ymax(ymax), _zmin(zmin), _zmax(zmax)
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
      : DopingFunction(ion, xmin, xmax, ymin, ymax, zmin, zmax), _PEAK(N)
  {}

  /**
   * compute doping concentration by point location
   */
  double profile(double x,double y,double z)
  {
    if( x >= _xmin-1e-6 && x <= _xmax+1e-6 &&
        y >= _ymin-1e-6 && y <= _ymax+1e-6 &&
        z >= _zmin-1e-6 && z <= _zmax+1e-6 )
      return _ion*_PEAK;
    else
      return 0.0;
  }

private:
  /**
   * the peak value of doping concentration
   */
  double _PEAK;
};



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
      : DopingFunction(ion, xmin, xmax, ymin, ymax, zmin, zmax),
      _PEAK(N), _XCHAR(ax), _YCHAR(ay), _ZCHAR(az), _XERFC(XERFC), _YERFC(YERFC), _ZERFC(ZERFC)
  {}

  /**
   * compute doping concentration by point location
   */
  double profile(double x, double y, double z)
  {
    double dx, dy, dz;
    if ( _XERFC )
#ifdef CYGWIN
      dx = (Erfc((x-_xmax)/_XCHAR)-Erfc((x-_xmin)/_XCHAR))/2.0;
#else
      dx = (erfc((x-_xmax)/_XCHAR)-erfc((x-_xmin)/_XCHAR))/2.0;
#endif
    else
    {
      if(x<_xmin)
        dx = exp(-(x-_xmin)*(x-_xmin)/(_XCHAR*_XCHAR));
      else if(x>=_xmin&&x<=_xmax)
        dx = 1.0;
      else
        dx = exp(-(x-_xmax)*(x-_xmax)/(_XCHAR*_XCHAR));
    }

    if ( _YERFC )
#ifdef CYGWIN
      dy = (Erfc((y-_ymax)/_YCHAR)-Erfc((y-_ymin)/_YCHAR))/2.0;
#else
      dy = (erfc((y-_ymax)/_YCHAR)-erfc((y-_ymin)/_YCHAR))/2.0;
#endif
    else
    {
      if(y<_ymin)
        dy = exp(-(y-_ymin)*(y-_ymin)/(_YCHAR*_YCHAR));
      else if(y>=_ymin&&y<=_ymax)
        dy = 1.0;
      else
        dy = exp(-(y-_ymax)*(y-_ymax)/(_YCHAR*_YCHAR));
    }

    if ( _ZERFC )
#ifdef CYGWIN
      dz = (Erfc((z-_zmax)/_ZCHAR)-Erfc((z-_zmin)/_ZCHAR))/2.0;
#else
      dz = (erfc((z-_zmax)/_ZCHAR)-erfc((z-_zmin)/_ZCHAR))/2.0;
#endif
    else
    {
      if(z<_zmin)
        dz = exp(-(z-_zmin)*(z-_zmin)/(_ZCHAR*_ZCHAR));
      else if(z>=_zmin&&z<=_zmax)
        dz = 1.0;
      else
        dz = exp(-(z-_zmax)*(z-_zmax)/(_ZCHAR*_ZCHAR));
    }

    return _ion*_PEAK*dx*dy*dz;
  }

private:
  /**
   * the peak value of doping concentration
   */
  double _PEAK;

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


#endif
