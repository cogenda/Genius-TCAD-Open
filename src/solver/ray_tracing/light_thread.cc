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


#include <cmath>

#include "plane.h"
#include "light_thread.h"


double LightThread::dead_factor = 1e-3;


std::vector<double> LightThread::advance_to(const Point & p_end, double a_band, double a_tail, double a_fc)
{
  const double length = (_p - p_end).size();
  _p = p_end;

  const double l_band = exp(-a_band*length);
  const double l_tail = exp(-a_tail*length);
  const double l_fc   = exp(-a_fc*length);

  const double power_start = _power;
  _power *= l_band*l_tail*l_fc;
  double dpower = power_start - _power;

  std::vector<double> power_loss(3);
  power_loss[0] = dpower*a_band/(a_band+a_tail+a_fc+1e-30);
  power_loss[1] = dpower*a_tail/(a_band+a_tail+a_fc+1e-30);
  power_loss[2] = dpower*a_fc/(a_band+a_tail+a_fc+1e-30);

  return power_loss;
}


std::pair<LightThread *, LightThread *> LightThread::interface_light_gen_linear_polarized(const Point & in_p, const Point & norm, double n1, double n2)
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

      LightThread * reflect_light = 0;
      LightThread * refract_light = 0;

      if( reflect_eff >=1e-9 )
      {
        Point _E_dir_reflect;
        if(n2>n1)
          _E_dir_reflect = (sqrt(reflect_parallel)*reflect_dir.cross(incident_plane.unit_normal()) - sqrt(reflect_perpendicular)*E_perpendicular).unit();
        else
          _E_dir_reflect = (sqrt(reflect_parallel)*reflect_dir.cross(incident_plane.unit_normal()) + sqrt(reflect_perpendicular)*E_perpendicular).unit();
        reflect_light = new LightThread(in_p, reflect_dir, _E_dir_reflect, _wavelength, _init_power, reflect_eff*_power);
      }

      if( refract_eff >=1e-9 )
      {
        Point _E_dir_refract = (sqrt(refract_parallel)*refract_dir.cross(incident_plane.unit_normal()) + sqrt(refract_perpendicular)*E_perpendicular).unit();
        refract_light = new LightThread(in_p, refract_dir, _E_dir_refract, _wavelength, _init_power, refract_eff*_power);
      }

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

