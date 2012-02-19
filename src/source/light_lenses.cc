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
/*  Author: xianghua zhang   zhangxih@163.com                                   */
/*                                                                              */
/********************************************************************************/
#include "config.h"

#ifdef WINDOWS
#define _USE_MATH_DEFINES // M_PI for VC
#endif

#include <cmath>

#include "parser.h"
#include "plane.h"
#include "light_thread.h"
#include "light_lenses.h"
#include "physical_unit.h"

using  PhysicalUnit::um;

LightLenses::LightLenses(Parser::InputParser &decks)
{
  for( decks.begin(); !decks.end(); decks.next() )
  {
    Parser::Card c = decks.get_current_card();

    if( c.key() == "LENS" )   // It's a LENS card
    {
      Lens * lens = new Lens;
      lens->id = c.get_string("id", "");
      if( c.is_parameter_exist("center") )
      {
        std::vector<double> center = c.get_array<double>("center");
        lens->center[0] = center[0]*um;
        lens->center[1] = center[1]*um;
        lens->center[2] = center[2]*um;
      }
      else
      {
        lens->center[0] = c.get_real("center.x",0.0)*um;
        lens->center[1] = c.get_real("center.y",0.0)*um;
        lens->center[2] = c.get_real("center.z",0.0)*um;
      }
      if( c.is_parameter_exist("norm") )
      {
        std::vector<double> norm = c.get_array<double>("norm");
        lens->norm[0] = norm[0];
        lens->norm[1] = norm[1];
        lens->norm[2] = norm[2];
      }
      else
      {
        lens->norm[0] = c.get_real("norm.x",0.0);
        lens->norm[1] = c.get_real("norm.y",1.0);
        lens->norm[2] = c.get_real("norm.z",0.0);
      }
      lens->norm.to_unit();
      lens->radius  = c.get_real("radius",1e4)*um;
      if(c.is_parameter_exist("f"))
      {
        double f = c.get_real("f",1.0)*um;
        lens->A = 1.0;
        lens->B = 0.0;
        lens->C = -1.0/f;
        lens->D = 1.0;
      }
      else
      {
        lens->A = c.get_real("a",1.0);
        lens->B = c.get_real("b",1.0)*um;
        lens->C = c.get_real("c",0.0)/um;
        lens->D = c.get_real("d",1.0);
      }

      _lenses.insert( std::make_pair(lens->id, lens) );
    }
  }
}


LightLenses::~LightLenses()
{
  std::map<std::string, Lens *>::const_iterator it= _lenses.begin();
  for(; it != _lenses.end(); ++it)
    delete it->second;
  _lenses.clear();
}


void LightLenses::set_lenses(const std::vector<std::string> & lens)
{
  _effect_lenses = lens;
}


Sphere LightLenses::bounding_sphere() const
{
  assert(!_effect_lenses.empty());
  const Lens * first_lens = _lenses.find(_effect_lenses[0])->second;
  Sphere s1(first_lens->center, first_lens->radius);

  for(unsigned int i=1; i<_effect_lenses.size(); ++i)
  {
    const Lens * lens = _lenses.find(_effect_lenses[i])->second;
    Sphere s2(lens->center, lens->radius);

    Point c1 = s1.center();
    Point c2 = s2.center();
    Point s_end1 = c1 + (c1-c2).unit(true)*s1.radius();
    Point s_end2 = c2 + (c2-c1).unit(true)*s2.radius();
    s1 = Sphere(0.5*(s_end1 + s_end2), 0.5*(s_end1-s_end2).size());
  }

  return s1;
}



LightThread * LightLenses::operator << (LightThread *light) const
{
  const double energy = light->power();
  do
  {
    //find intersection lens
    double t=1e30;
    const Lens * active_lens=0;
    for(unsigned int n=0; n<_effect_lenses.size(); ++n)
    {
      const Lens * lens = _lenses.find(_effect_lenses[n])->second;
      Plane lens_plane(lens->center, lens->norm);
      double dist;
      if(lens_plane.intersect_point(light->start_point(), light->dir(), dist) )
      {
        if( (lens->center - (light->start_point() + dist*light->dir()) ).size() > lens->radius) continue;
        if(dist > 0 && dist<t)
        {
          t = dist;
          active_lens = lens;
        }
      }
    }

    if(!active_lens)
    {
      return light;
    }

    // do ABCD transform
    {

      // E_p           y.axis
      // | E_t         |
      // .____rayin____| aim.point
      //               |\ray_out
      //               | \
      //               |  \   / E_p
      // aim_distance  |   \ /
      //               |    . E_t
      //               |     \
      //            z  .----------------------- lens.norm (x.axis)
      //           lens.center


      // the light incident point
      Point aim_point = light->start_point() + t*light->dir();
      // set y axis as lens center to light incident point
      Point y_axis = (aim_point - active_lens->center).unit(true);
      Point z_axis = active_lens->norm.cross(y_axis);
      // the parallel and perpendicular E
      const double Et = light->E_dir()*z_axis;
      const double Ep = light->E_dir()*(z_axis.cross(light->dir()));
      //
      double aim_distance = (aim_point - active_lens->center).size();
      double aim_angle = active_lens->norm.angle(light->dir());
      if(y_axis*light->dir() < 0.0)
        aim_angle = 2*M_PI-aim_angle;
      // ABCD transform
      double t_distance = active_lens->A*aim_distance + active_lens->B*aim_angle;
      double t_angle = active_lens->C*aim_distance + active_lens->D*aim_angle;
      // recalculate startpoint, dir and E_dir of light
      light->dir() =  cos(t_angle)*active_lens->norm + sin(t_angle)*y_axis;
      light->start_point() = active_lens->center + t_distance*y_axis + 1e-6*active_lens->radius*light->dir();
      light->E_dir() =  Ep*z_axis.cross(light->dir()) + Et*z_axis;
    }
  } while(light->power() > 0.01*energy);

  return 0;
}


