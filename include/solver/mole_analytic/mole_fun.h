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

#ifndef __mole_fun_h__
#define __mole_fun_h__

#include "enum_direction.h"

/**
 * the base class for mole function
 */
class MoleFunction
{
public:
  /**
   * constructor, do nothing
   */
  MoleFunction(const std::string & region, double xmin=0, double xmax=0, double ymin=0, double ymax=0, double zmin=0, double zmax=0)
      : _region(region), _xmin(xmin), _xmax(xmax), _ymin(ymin), _ymax(ymax), _zmin(zmin), _zmax(zmax)
  {}

  /**
   * destructor, do nothing
   */
  virtual ~MoleFunction() {}

  /**
   * virtual compute function
   */
  virtual double profile(double x, double y, double z)=0;

  /**
   * @return region name as const reference
   */
  const std::string & region() const
  { return _region;}

  /**
   * @return \p true if this is first mole fraction for signle compound material
   */
  virtual bool mole_x()=0;

  /**
   * @return \p true if this second mole fraction for complex compound material
   */
  virtual bool mole_y()=0;

protected:

  /**
   * the compound semicondcutor region
   */
  std::string _region;

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
 * linear mole function
 */
class LinearMoleFunction : public MoleFunction
{
public:

  /**
   * constructor
   */
  LinearMoleFunction(const std::string & region, double base_plane, Direction dir, double mole_begin, double mole_slope, bool mole_x)
      : MoleFunction(region),  _mole_begin(mole_begin), _mole_slope(mole_slope), _base_plane(base_plane), _dir(dir)
  {
    _mole_x = mole_x;
    _mole_y = !mole_x;
  }

  /**
   * destructor, do nothing
   */
  virtual ~LinearMoleFunction() {}

  /**
   * compute mole fraction by point location
   */
  double profile(double x,double y,double z)
  {

    switch(_dir)
    {
      case X_Direction:
         return  _mole_begin + (x-_base_plane)*_mole_slope;
      case Y_Direction:
         return  _mole_begin + (y-_base_plane)*_mole_slope;
      case Z_Direction:
         return  _mole_begin + (z-_base_plane)*_mole_slope;
      default: genius_error();
    }
    return 0;
  
  }

  /**
   * @return \p true if this is first mole fraction for signle compound material
   */
  virtual bool mole_x()
  { return _mole_x;}

  /**
   * @return \p true if this second mole fraction for complex compound material
   */
  virtual bool mole_y()
  { return _mole_y;}

private:

  /**
   * indicate iff first mole fraction for single compound material
   */
  bool   _mole_x;

  /**
   * indicate iff second mole fraction for complex compound material
   */
  bool   _mole_y;

  /**
   * mole fraction at the _base_plane
   */
  double _mole_begin;

  /**
   * the slope of the mole fraction for graded compounds
   */
  double _mole_slope;

  /**
   * the reference base plane, at which the mole fraction is _mole_begin
   */
  double _base_plane;

  /**
   * the mole fraction grading direction
   */
  Direction _dir;

};





#endif // __mole_fun_h__
