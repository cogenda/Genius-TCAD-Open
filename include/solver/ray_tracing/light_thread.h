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


#ifndef __light_thread_h__
#define __light_thread_h__


//local include
#include "point.h"
#include "elem_intersection.h"

class Elem;
class ARCoatings;

/**
 * class to define a light
 */
class LightThread
{
public:

  LightThread(const Point & p, const Point & dir, const Point & E_dir, double wavelength, double init_power, double power)
      :_p(p), _dir(dir.unit()), _E_dir(E_dir), _wavelength(wavelength),
      _init_power(init_power), _power(power)
  {
    hit_elem = NULL;
  }

  /**
   * light advance to a new point. recalculate the light power.
   * @param  p_end new light point
   * @param  a_band band to band absorption
   * @param  a_tail band tail absorption
   * @param  a_fc free carrier absorption
   * @return power loss due to band/tail/fc
   */
  std::vector<double> advance_to(const Point & p_end, double a_band, double a_tail, double a_fc);

  /**
   * when light energy less than some factor (default 1e-3) of origin energy, it is dead
   */
  bool is_dead() const { return _power < dead_factor*_init_power; }

  /**
   * @return writable reference to start point of the light
   */
  Point & start_point()  { return _p; }

  /**
   * @return the start point of the light
   */
  const Point & start_point() const    { return _p; }

  /**
   * @return the writable reference to direction of the light
   */
  Point & dir()  { return _dir; }

  /**
   * @return the direction of the light
   */
  const Point & dir() const    { return _dir; }

  /**
   * @return writable reference to direction of E field
   */
  Point & E_dir() { return _E_dir; }

  /**
   * @return the direction of E field
   */
  const Point & E_dir() const { return _E_dir; }

  /**
   *@return the wavelength of the light
   */
  double wavelength() const { return _wavelength; }

  /**
   * @return the intial power
   */
  double init_power() const
    { return _init_power; }

  /**
   * @return light power
   */
  double power() const
    { return _power; }

  /**
   * writable reference to power
   */
  double & power()
  { return _power; }

  /**
   * set the light as reflection by mirror
   */
  LightThread * reflection(const Point & in_p, const Point & norm) const;

  /**
   * generate reflection and transmission light at interface
   * @param in_p incident point
   * @param norm the norm of interface, form material 2 to material 1
   * @param n1   refraction index of material 1
   * @param n2   refraction index of material 2
   * @param arc  anti-reflection coating
   * @param inv  the layer order of anti-reflection coating seen by the light
   */
  std::pair<LightThread *, LightThread *> interface_light_gen_linear_polarized
      (const Point & in_p, const Point & norm, double n1, double n2, const ARCoatings *arc=0, bool inv=false) const;

  /**
   * pointer to the elem this light hit
   */
  const Elem * hit_elem;

  /**
   * the intersection result with  elem
   */
  IntersectionResult result;

  /**
   * factor to determine it is dead
   */
  static double dead_factor;

private:

  /**
   * starting point of this thread
   */
  Point _p;

  /**
   * directional vector of this thread
   */
  Point _dir;

  /**
   * direction vector of electrical field
   */
  Point _E_dir;

  /**
   * wave length of the light
   */
  const double _wavelength;

  /**
   * initial power of this thread
   */
  const double _init_power;

  /**
   * current power of this thread
   */
  double _power;

  /// light refraction on simple interface
  std::pair<LightThread *, LightThread *> _interface_light_gen_linear_polarized_simple
      (const Point & in_p, const Point & norm, double n1, double n2 ) const;

  /// light refraction on stacked interface
  std::pair<LightThread *, LightThread *> _interface_light_gen_linear_polarized_stack
      (const Point & in_p, const Point & norm, double n1, double n2, const ARCoatings *arc, bool inv ) const;


};

#endif


