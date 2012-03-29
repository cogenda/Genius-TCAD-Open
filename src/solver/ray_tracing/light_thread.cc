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
#include <algorithm>

#include "plane.h"
#include "light_thread.h"
#include "anti_reflection_coating.h"
#include "physical_unit.h"

using PhysicalUnit::eps0;
using PhysicalUnit::mu0;


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


LightThread * LightThread::reflection(const Point & in_p, const Point & norm) const
{
  Point reflect_dir = (_dir-2*(norm.dot(_dir))*norm).unit();
  return  new LightThread(in_p, reflect_dir, _E_dir, _wavelength, _init_power, _power);
}

#if 0
std::pair<LightThread *, LightThread *> LightThread::_interface_light_gen_linear_polarized_simple(const Point & in_p, const Point & norm, double n1, double n2 ) const
{
  double n = n2/n1;

  // ref http://en.wikipedia.org/wiki/Snell%27s_law

  //for incidence angle
  double in_angle = acos(fabs(norm.dot(_dir)));

  //for reflection angle
  Point reflect_dir = (_dir-2*(norm.dot(_dir))*norm).unit(true);

  //for refraction angle
  bool full_reflection;
  double refract_angle;
  Point refract_dir;
  if (1/n*(1-pow(norm.dot(_dir),2))<=1)
  {
    full_reflection = false;
    double cos_theta_1 = -norm.dot(_dir);
    double cos_theta_2 = sqrt(1-1/n*(1-pow(norm.dot(_dir),2)));
    refract_angle = acos(fabs(cos_theta_2));
    refract_dir = (1/n*_dir-(cos_theta_2-1/n*cos_theta_1)*norm).unit(true);
  }
  else
  {
    full_reflection = true; //indicate full reflection
  }

  // build reflect_light and refract_light
  if (in_angle>=1e-9) //if not perpendicular incidence
  {
    Plane incident_plane(in_p, norm.cross(_dir));
    // project electrical vector to the plane consisted by norm.cross(_dir)
    // the E_perpendicular direction will not changed in reflect/refract
    Point E_perpendicular = (_E_dir*incident_plane.unit_normal())*incident_plane.unit_normal();
    Point E_parallel = _E_dir - E_perpendicular;//TM
    Point B_perpendicular = _dir.cross(E_parallel).unit(true);//TM,also, B_perpendicular keeps const in the reflect/refract

    if (!full_reflection) //if not full reflection
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

      // the E_perpendicular keeps the same
      //however, the E_parallel direction should be re-defined as light_dir cross incident_plane.unit_normal
      if( reflect_eff >=1e-9 )
      {
        Point _E_dir_reflect;
        Point E_parallel_reflect = E_parallel.size()*B_perpendicular.cross(reflect_dir);
        if(n2>n1)
          _E_dir_reflect = (sqrt(reflect_parallel)*E_parallel_reflect - sqrt(reflect_perpendicular)*E_perpendicular).unit(true);
        else
          _E_dir_reflect = (sqrt(reflect_parallel)*E_parallel_reflect + sqrt(reflect_perpendicular)*E_perpendicular).unit(true);
        reflect_light = new LightThread(in_p, reflect_dir, _E_dir_reflect, _wavelength, _init_power, reflect_eff*_power);
      }

      if( refract_eff >=1e-9 )
      {
        Point E_parallel_refract = E_parallel.size()*B_perpendicular.cross(refract_dir);
        Point _E_dir_refract = (sqrt(refract_parallel)*E_parallel_refract + sqrt(refract_perpendicular)*E_perpendicular).unit(true);
        refract_light = new LightThread(in_p, refract_dir, _E_dir_refract, _wavelength, _init_power, refract_eff*_power);
      }

      return std::make_pair(reflect_light, refract_light);
    }
    else //for full reflection
    {
      Point _E_dir_reflect = (reflect_dir.cross(incident_plane.unit_normal()) + E_perpendicular).unit(true);
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
#endif


#if 1
std::pair<LightThread *, LightThread *> LightThread::_interface_light_gen_linear_polarized_simple(const Point & in_p, const Point & norm, double n1, double n2 ) const
{
  double n = n2/n1;

  // ref http://en.wikipedia.org/wiki/Snell%27s_law

  //for incidence angle
  double cos_in_angle = -norm.dot(_dir);
  double in_angle = acos(cos_in_angle);

  //for reflection angle
  Point reflect_dir = (_dir-2*(norm.dot(_dir))*norm).unit(true);

  //for refraction angle
  bool full_reflection;
  double refract_angle;
  Point refract_dir;
  if (1/n*(1-pow(norm.dot(_dir),2))<=1)
  {
    full_reflection = false;
    double cos_theta_1 = -norm.dot(_dir);
    double cos_theta_2 = sqrt(1-1/n*(1-pow(norm.dot(_dir),2)));
    refract_angle = acos(fabs(cos_theta_2));
    refract_dir = (1/n*_dir-(cos_theta_2-1/n*cos_theta_1)*norm).unit(true);
  }
  else
  {
    full_reflection = true; //indicate full reflection
  }

  // build reflect_light and refract_light
  if (in_angle>=1e-9) //if not perpendicular incidence
  {
    Plane incident_plane(in_p, norm.cross(_dir));
    // project electrical vector to the plane consisted by norm.cross(_dir)
    // the E_perpendicular direction will not changed in reflect/refract
    Point E_perpendicular = (_E_dir*incident_plane.unit_normal())*incident_plane.unit_normal();
    Point E_parallel = _E_dir - E_perpendicular;
    Point B_perpendicular = _dir.cross(E_parallel).unit(true);//TM,also, B_perpendicular keeps const in the reflect/refract

    if (!full_reflection) //if not full reflection
    {
      double _polarization_angle = asin(E_perpendicular.size());

      //for reflection efficiency
      double reflect_parallel = tan(refract_angle-in_angle)/tan(in_angle+refract_angle);//Rp
      double reflect_perpendicular = sin(refract_angle-in_angle)/sin(in_angle+refract_angle);//Rs
      double reflect_eff = pow(reflect_parallel*cos(_polarization_angle),2)+pow(reflect_perpendicular*sin(_polarization_angle),2);

      //for refraction efficiency
      double refract_parallel = 2*sin(refract_angle)*cos(in_angle)/(sin(in_angle+refract_angle)*cos(in_angle-refract_angle));//Tp
      double refract_perpendicular = 2*sin(refract_angle)*cos(in_angle)/sin(in_angle+refract_angle);//Ts
      double refract_eff = pow(refract_parallel*cos(_polarization_angle),2)+pow(refract_perpendicular*sin(_polarization_angle),2);

      LightThread * reflect_light = 0;
      LightThread * refract_light = 0;

      // the E_perpendicular keeps the same
      //however, the E_parallel direction should be re-defined as light_dir cross incident_plane.unit_normal
      if( reflect_eff >=1e-9 )
      {
        Point _E_dir_reflect;
        Point E_parallel_reflect = E_parallel.size()*B_perpendicular.cross(reflect_dir);
        _E_dir_reflect = (reflect_parallel*E_parallel_reflect + reflect_perpendicular*E_perpendicular).unit(true);
        reflect_light = new LightThread(in_p, reflect_dir, _E_dir_reflect, _wavelength, _init_power, reflect_eff*_power);
      }

      if( refract_eff >=1e-9 )
      {
        Point E_parallel_refract = E_parallel.size()*B_perpendicular.cross(refract_dir);
        Point _E_dir_refract = (refract_parallel*E_parallel_refract + refract_perpendicular*E_perpendicular).unit(true);
        refract_light = new LightThread(in_p, refract_dir, _E_dir_refract, _wavelength, _init_power, refract_eff*_power);
      }

      return std::make_pair(reflect_light, refract_light);
    }
    else //for full reflection
    {
      Point E_parallel_reflect = E_parallel.size()*B_perpendicular.cross(reflect_dir);
      Point _E_dir_reflect = (E_parallel_reflect + E_perpendicular).unit(true);
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
#endif


// TNT matrix-vector library
#include <TNT/tnt.h>
std::pair<LightThread *, LightThread *> LightThread::_interface_light_gen_linear_polarized_stack(const Point & in_p, const Point & norm,
    double n1, double n2, const ARCoatings *arc, bool inv ) const
{
  Point incident_plane_norm = norm.cross(_dir);
  if(incident_plane_norm.size() == 0.0)
    incident_plane_norm = _E_dir;
  Plane incident_plane(in_p, incident_plane_norm);

  // project electrical vector to the plane consisted by norm.cross(_dir)
  // the E_perpendicular direction will not changed in reflect/refract
  Point E_perpendicular = (_E_dir*incident_plane.unit_normal())*incident_plane.unit_normal();//TE
  Point E_parallel = _E_dir - E_perpendicular;//TM
  Point B_perpendicular = _dir.cross(E_parallel).unit(true);//TM

  double TE_power = _power*E_perpendicular.size_sq()/(E_perpendicular.size_sq()+E_parallel.size_sq());
  double TM_power = _power*E_parallel.size_sq()/(E_perpendicular.size_sq()+E_parallel.size_sq());

  // stack parameter
  unsigned int layers = arc->layers;
  std::vector< std::complex<double> > layer_refractive_index = arc->layer_refractive_index;
  std::vector<double> layer_thickness = arc->layer_thickness;
  if(inv)
  {
    std::reverse(layer_refractive_index.begin(), layer_refractive_index.end());
    std::reverse(layer_thickness.begin(), layer_thickness.end());
  }

  TNT::Array2D<std::complex<double> > TE_M(2,2);
  TE_M[0][0] = 1.0; TE_M[1][1] = 1.0;
  TNT::Array2D<std::complex<double> > TM_M(2,2);
  TM_M[0][0] = 1.0; TM_M[1][1] = 1.0;

  //for incidence angle
  double n_in  = n1;
  double n_out = n2;
  std::complex<double> cos_in_angle = std::abs(norm.dot(_dir));
  std::complex<double> cos_out_angle;
  for(unsigned int i=0; i<layers; i++)
  {
    double n = layer_refractive_index[i].real();
    std::complex<double> cos_angle = std::sqrt(1.0-n_in*n_in/(n*n)*(1.0-cos_in_angle*cos_in_angle));

    std::complex<double> N_TE = layer_refractive_index[i]*cos_angle;
    std::complex<double> N_TM = layer_refractive_index[i]/cos_angle;

    //optical admittance
    std::complex<double> Y_TE = N_TE;
    std::complex<double> Y_TM = N_TM;

    //std::cout<<layer_thickness[i]*cos_angle/_wavelength << " " << layer_refractive_index[i] << " "<< sqrt(n1*n2)<<std::endl;

    // phase delay of the light passing through layer
    std::complex<double> D = 2*3.14159265359*layer_refractive_index[i]*layer_thickness[i]*cos_angle/_wavelength;

    //TE characteristic matrix
    TNT::Array2D<std::complex<double> > TE_matrix(2,2);
    TE_matrix[0][0] = std::cos(D);  TE_matrix[0][1] = std::complex<double>(0,1)*std::sin(D)/Y_TE;
    TE_matrix[1][0] = std::complex<double>(0,1)*std::sin(D)*Y_TE;  TE_matrix[1][1] = std::cos(D);
    TE_M = TNT::matmult(TE_M, TE_matrix);

    //TM characteristic matrix
    TNT::Array2D<std::complex<double> > TM_matrix(2,2);
    TM_matrix[0][0] = std::cos(D);  TM_matrix[0][1] = std::complex<double>(0,1)*std::sin(D)/Y_TM;
    TM_matrix[1][0] = std::complex<double>(0,1)*std::sin(D)*Y_TM;  TM_matrix[1][1] = std::cos(D);
    TM_M = TNT::matmult(TM_M, TE_matrix);

    n_in = n;
    cos_in_angle = cos_angle;
  }

  n_out = n2;
  cos_out_angle = std::sqrt(1.0-n_in*n_in/(n_out*n_out)*(1.0-cos_in_angle*cos_in_angle));
  n_in  = n1;
  cos_in_angle = std::abs(norm.dot(_dir));


  std::complex<double> Yin_TE  = n_in*cos_in_angle;
  std::complex<double> Yout_TE = n_out*cos_out_angle;
  std::complex<double> r_TE = (Yin_TE*TE_M[0][0]+Yin_TE*Yout_TE*TE_M[0][1]-TE_M[1][0]-Yout_TE*TE_M[1][1])
                             /(Yin_TE*TE_M[0][0]+Yin_TE*Yout_TE*TE_M[0][1]+TE_M[1][0]+Yout_TE*TE_M[1][1]);
  std::complex<double> t_TE = 2.0*Yin_TE/(Yin_TE*TE_M[0][0]+Yin_TE*Yout_TE*TE_M[0][1]+TE_M[1][0]+Yout_TE*TE_M[1][1]);
  double R_TE = std::norm(r_TE);
  double T_TE = Yout_TE.real()/Yin_TE.real()*std::norm(t_TE);

  std::complex<double> Yin_TM  = n_in/cos_in_angle;
  std::complex<double> Yout_TM = n_out/cos_out_angle;
  std::complex<double> r_TM = (Yin_TM*TM_M[0][0]+Yin_TM*Yout_TM*TM_M[0][1]-TM_M[1][0]-Yout_TM*TM_M[1][1])
                             /(Yin_TM*TM_M[0][0]+Yin_TM*Yout_TM*TM_M[0][1]+TM_M[1][0]+Yout_TM*TM_M[1][1]);
  std::complex<double> t_TM = 2.0*Yin_TM/(Yin_TM*TM_M[0][0]+Yin_TM*Yout_TM*TM_M[0][1]+TM_M[1][0]+Yout_TM*TM_M[1][1]);
  double R_TM = std::norm(r_TM);
  double T_TM = Yout_TM.real()/Yin_TM.real()*std::norm(t_TM);

  LightThread * reflect_light = 0;
  LightThread * refract_light = 0;

  //for reflection angle
  Point reflect_dir = (_dir-2*(norm.dot(_dir))*norm).unit(true);
  //for refraction angle
  Point refract_dir =  _dir;
  if(acos(cos_in_angle.real()) > 0)
  {
    double n = sin(acos(cos_out_angle.real()))/sin(acos(cos_in_angle.real()));
    refract_dir = (1/n*_dir-(cos_out_angle.real()-1/n*cos_in_angle.real())*norm).unit(true);
  }

  Point E_dir_reflect = (E_perpendicular*r_TE.real() + E_parallel.size()*B_perpendicular.cross(reflect_dir)*r_TM.real()).unit(true);
  Point E_dir_refract = (E_perpendicular*t_TE.real() + E_parallel.size()*B_perpendicular.cross(refract_dir)*t_TM.real()).unit(true);;

  double reflect_power = TE_power*R_TE + TM_power*R_TM;
  double refract_power = TE_power*T_TE + TM_power*T_TM;

  //std::cout<<"COFF " << R_TE << " " << R_TM << " " << T_TE << " " << T_TM << std::endl;

  if( reflect_power > 1e-10*_power )
    reflect_light = new LightThread(in_p, reflect_dir, E_dir_reflect, _wavelength, _init_power, reflect_power);

  if( refract_power > 1e-10*_power )
    refract_light = new LightThread(in_p, refract_dir, E_dir_refract, _wavelength, _init_power, refract_power);

  return std::make_pair(reflect_light, refract_light);
}





std::pair<LightThread *, LightThread *> LightThread::interface_light_gen_linear_polarized(const Point & in_p, const Point & norm,
    double n1, double n2, const ARCoatings *arc, bool inv) const
{
  if(norm.dot(_dir) == 0.0 ) return std::make_pair((LightThread *)0, (LightThread *)0);

  if(arc)
    return _interface_light_gen_linear_polarized_stack(in_p, norm, n1, n2, arc, inv);
  else
    return _interface_light_gen_linear_polarized_simple(in_p, norm, n1, n2);
}

