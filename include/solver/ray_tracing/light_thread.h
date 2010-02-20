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

#include <cmath>

//local include
#include "point.h"
#include "plane.h"

class Elem;

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
   * @param  im_index loss parameter
   * @return power loss
   */
  double advance_to(const Point & p_end, double im_index)
  {
    const double pi = 3.14159265358979;
    double loss_rate = 4*pi*im_index/_wavelength;

    double length = (_p - p_end).size();
    _p = p_end;

    double power_start = _power;
    _power *= exp(-loss_rate*length);
    return  power_start - _power;
  }

  /**
   * when light energy less than 1/1000 of origin energy, it is dead
   */
  bool is_dead() const
    { return _power < 1e-3*_init_power; }

  /**
   * @return the start point of the light
   */
  Point & start_point()
  { return _p; }

  /**
   * @return the start point of the light
   */
  const Point & start_point() const
  { return _p; }

  /**
   * @return the direction of the light
   */
  const Point & dir() const
  { return _dir; }

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
   * generate reflection and transmission light at interface
   * @param in_p incident point
   * @param norm the norm of interface, form material 2 to material 1
   * @param n1   refraction index of material 1
   * @param n2   refraction index of material 2
   */
  std::pair<LightThread *, LightThread *> interface_light_gen_linear_polarized(const Point & in_p, const Point & norm, double n1, double n2)
  {
    double n = n2/n1;

    //for incidence angle
    double in_angle = acos(fabs(norm.dot(_dir)));

    //for reflection angle
    Point reflect_dir = (_dir-2*(norm.dot(_dir))*norm).unit();

    //for refraction angle
    double refract_angle;
    Point refract_dir;
    if (1/n*(1-pow(norm.dot(_dir),2))<=1)
    {
      double cos_theta_1 = -norm.dot(_dir);
      double cos_theta_2 = sqrt(1-1/n*(1-pow(norm.dot(_dir),2)));
      refract_angle = acos(fabs(cos_theta_2));
      refract_dir = (1/n*_dir-(cos_theta_2-1/n*cos_theta_1)*norm).unit();
    }
    else
      refract_angle = 1.0e12; //indicate full reflection

    // build reflect_light and refract_light
    if (in_angle>=1e-9) //if not perpendicular incidence
    {
      Plane incident_plane(in_p, norm.cross(_dir));
      //project electrical vector to the plane consisted by norm.cross(_dir)
      // the E_perpendicular direction will not changed in reflect/refract
      Point E_perpendicular = (_E_dir*incident_plane.unit_normal())*incident_plane.unit_normal();
      //however, the E_parallel should be defined as light_dir cross incident_plane.unit_normal

      if (refract_angle < 1.0e5) //if not full reflection
      {
        double _polarization_angle = asin(E_perpendicular.size());

        //for reflection efficiency
        double reflect_parallel = pow(tan(in_angle-refract_angle),2)/pow(tan(in_angle+refract_angle),2);
        double reflect_perpendicular = pow(sin(in_angle-refract_angle),2)/pow(sin(in_angle+refract_angle),2);
        double reflect_eff = reflect_parallel*pow(cos(_polarization_angle),2)+reflect_perpendicular*pow(sin(_polarization_angle),2);

        //for refraction efficiency
        double refract_parallel = sin(2*in_angle)*sin(2*refract_angle)/pow(sin(in_angle+refract_angle)*cos(in_angle-refract_angle),2);
        double refract_perpendicular = sin(2*in_angle)*sin(2*refract_angle)/pow(sin(in_angle+refract_angle),2);
        double refract_eff = refract_parallel*pow(cos(_polarization_angle),2)+refract_perpendicular*pow(sin(_polarization_angle),2);

        Point _E_dir_reflect;
        if(n2>n1)
          _E_dir_reflect = (sqrt(reflect_parallel)*reflect_dir.cross(incident_plane.unit_normal()) - sqrt(reflect_perpendicular)*E_perpendicular).unit();
        else
          _E_dir_reflect = (sqrt(reflect_parallel)*reflect_dir.cross(incident_plane.unit_normal()) + sqrt(reflect_perpendicular)*E_perpendicular).unit();

        Point _E_dir_refract = (sqrt(refract_parallel)*refract_dir.cross(incident_plane.unit_normal()) + sqrt(refract_perpendicular)*E_perpendicular).unit();

        LightThread * reflect_light = new LightThread(in_p, reflect_dir, _E_dir_reflect, _wavelength, _init_power, reflect_eff*_power);
        LightThread * refract_light = new LightThread(in_p, refract_dir, _E_dir_refract, _wavelength, _init_power, refract_eff*_power);

        return std::make_pair(reflect_light, refract_light);
      }
      else //for full reflection
      {
        Point _E_dir_reflect = (reflect_dir.cross(incident_plane.unit_normal()) + E_perpendicular).unit();
        LightThread * reflect_light = new LightThread(in_p, reflect_dir, _E_dir_reflect, _wavelength, _init_power, _power);

        return std::make_pair(reflect_light, (LightThread *)0);
      }
    }
    else //if perpendicular incidence, for this instance, full reflection will never happened
    {
      double reflect_eff = (n-1)*(n-1)/((n+1)*(n+1));
      double refract_eff = 4*n/((n+1)*(n+1));

      LightThread * reflect_light = new LightThread(in_p, -_dir, -_E_dir, _wavelength, _init_power, reflect_eff*_power);
      LightThread * refract_light = new LightThread(in_p,  _dir,  _E_dir, _wavelength, _init_power, refract_eff*_power);

      return std::make_pair(reflect_light, refract_light);
    }

  }

  /**
   * pointer to the elem this light hit
   */
  const Elem * hit_elem;

  /**
   * the intersection result with  elem
   */
  IntersectionResult result;

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

};

#endif


