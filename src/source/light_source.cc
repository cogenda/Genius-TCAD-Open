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

#include <fstream>
#include <iomanip>


#include "simulation_system.h"
#include "simulation_region.h"
#include "interpolation_2d_csa.h"
#include "interpolation_2d_nn.h"
//#include "interpolation_3d_qshep.h"
#include "interpolation_3d_nbtet.h"
#include "light_source.h"
#include "mathfunc.h"
#include "log.h"

using PhysicalUnit::s;
using PhysicalUnit::V;
using PhysicalUnit::um;
using PhysicalUnit::cm;
using PhysicalUnit::m;
using PhysicalUnit::J;
using PhysicalUnit::eps0;
using PhysicalUnit::h;

//Light_Source_From_File::Light_Source_From_File(SimulationSystem &system, const Parser::Card &c):Light_Source(system)
//{
//  MESSAGE<<"Loading Light Source"; RECORD();
//
//  genius_assert(c.key() == "LIGHT");
//
//  std::string fname = c.get_string("profile.file", "");
//  _wave_length   = c.get_real("wavelength", 0.532, "lambda") * um;   // wave length
//  _power         = c.get_real("intensity",  0.0) * J/s/cm/cm;      // incident power;
//  _eta           = c.get_real("quan.eff", 1.0);
//  _eta_auto      = !c.is_parameter_exist("quan.eff");
//
//
//  // PROFILE card should have "type" parameter
//  if( !c.is_parameter_exist("profile") )
//  {
//    MESSAGE<<"ERROR at " << c.get_fileline() <<" LIGHT: Should have a 'profile' parameter."<<std::endl; RECORD();
//    genius_error();
//  }
//
//  MESSAGE<<" at lambda="<< std::setprecision(3) << std::setw(10) << _wave_length/um << "um from "<< fname << " ..."; RECORD();
//
//
//  int skip_line      = c.get_int("skipline", 0);
//
//  // set Light profile here.
//  if(c.is_enum_value("profile","fromfile2d"))
//    set_particle_profile_fromfile2d(c, fname, skip_line);
//
//  if(c.is_enum_value("profile","fromfile3d"))
//    set_particle_profile_fromfile3d(c, fname, skip_line);
//
//  MESSAGE<<"ok\n"<<std::endl; RECORD();
//
//}

Light_Source_From_File::Light_Source_From_File(SimulationSystem &system, const Parser::Card &c,
               const std::string &fname_ext,
               const double wave_length, const double power,
               const double eta, const bool eta_auto):Light_Source(system)
{

  genius_assert(c.key() == "LIGHT");

  _wave_length = wave_length;
  _power       = power;
  _eta         = eta;
  _eta_auto    = eta_auto;

  // PROFILE card should have "type" parameter
  if( !c.is_parameter_exist("profile") )
  {
    MESSAGE<<"ERROR at " << c.get_fileline() <<" LIGHT: Should have a 'profile' parameter."<<std::endl; RECORD();
    genius_error();
  }


  _fname = c.get_string("profile.file", "");
  if (fname_ext.length()>0)
  {
    _fname += ".";
    _fname += fname_ext;
  }

  _skip_line      = c.get_int("skipline", 0);

  // set Light profile here.
  if(c.is_enum_value("profile","fromfile2d"))
    _dim = 2;
  if(c.is_enum_value("profile","fromfile3d"))
    _dim = 3;

  if (c.is_enum_value("lunit", "m"))
    _LUnit = 1.0;
  else if (c.is_enum_value("lunit", "cm"))
    _LUnit = 1e-2;
  else if (c.is_enum_value("lunit", "um"))
    _LUnit = 1e-6;
  else if (c.is_enum_value("lunit", "nm"))
    _LUnit = 1e-9;
  else
    _LUnit = 1e-6;

  if (c.is_enum_value("funit", "m"))
    _FUnit = 1.0;
  else if (c.is_enum_value("funit", "cm"))
    _FUnit = 1e2;
  else if (c.is_enum_value("funit", "um"))
    _FUnit = 1e6;
  else if (c.is_enum_value("funit", "nm"))
    _FUnit = 1e9;
  else
    _FUnit = 1.0;

  if( c.is_parameter_exist("translate") )
  {
    std::vector<double> dummy = c.get_array<double>("translate");
    _translate = VectorValue<double>(&dummy[0])*um;
  }
  else
  {
    _translate = VectorValue<double>(c.get_real("translate.x", 0.0)*um,
                                     c.get_real("translate.y", 0.0)*um,
                                     c.get_real("translate.z", 0.0)*um);
  }

  if( c.is_parameter_exist("transform") )
  {
    std::vector<double> dummy = c.get_array<double>("transform");
    _transform = TensorValue<double>(&dummy[0]);
  }
  else
  {
    _transform = TensorValue<double>(c.get_real("transform.xx", 1.0), c.get_real("transform.xy", 0.0), c.get_real("transform.xz", 0.0),
                                     c.get_real("transform.yx", 0.0), c.get_real("transform.yy", 1.0), c.get_real("transform.yz", 0.0),
                                     c.get_real("transform.zx", 0.0), c.get_real("transform.zy", 0.0), c.get_real("transform.zz", 1.0));
  }

}

int Light_Source_From_File::load_light_profile_fromfile(InterpolationBase * interpolator, const std::string &fname, const int skip_line)
{

  if(Genius::processor_id()==0)
  {
    std::ifstream in(fname.c_str());
    if(!in.good())
    {
      MESSAGE << "Error: E-Field: file "<<fname<<" can't be opened."<<std::endl; RECORD();
      genius_error();
    }

    Point p;
    double EField;

    int i;
    for (i=0; i<skip_line && !in.eof(); i++)
    {
      std::string line;
      std::getline(in, line);
    }

    int data_cnt=0;
    while(!in.eof() && in.good())
    {
      in >> p[0];// read x location
      in >> p[1];// read y location
      if (_dim==3)
        in >> p[2];// read z location

      p *= PhysicalUnit::m*_LUnit;     // scale to length unit
      p = _transform*p + _translate; //do transform & translate
      in >> EField;
      EField *= V/PhysicalUnit::m * _FUnit;
      if(!in.fail())
      {
        interpolator->add_scatter_data(p, 0, EField);
        data_cnt ++;
      }
      i++;
    }
    if(!in.eof())
    {
      MESSAGE<<"ERROR : error reading data file " << fname;
      MESSAGE<<" at line " << i << "." << std::endl; RECORD();
      genius_error();
    }
    in.close();
    return data_cnt;
  }
  return 0;

}


void Light_Source_From_File::update_system()
{

  // load light E-field profile from File and fill into the interpolator
  MESSAGE<<"Loading Light Source"; RECORD();
  MESSAGE<<" at lambda="<< std::setprecision(3) << std::setw(10) << _wave_length/um << "um from "<< _fname << " ..."; RECORD();

  InterpolationBase * interpolator;
  if (_dim == 2)
    interpolator = new Interpolation2D_NN;
    //interpolator = new Interpolation2D_CSA;
  else
    interpolator = new Interpolation3D_nbtet;
    //interpolator = new Interpolation3D_qshep;

  interpolator->set_interpolation_type(0, InterpolationBase::Linear);

  if(Genius::processor_id()==0)
  {
    int count=0;
    count = load_light_profile_fromfile(interpolator, _fname, _skip_line);
    MESSAGE<<count<<" points ";
  }

  interpolator->broadcast(0);
  interpolator->setup(0);

  MESSAGE<<"ok"<<std::endl; RECORD();

  double c0 = 299792458*m/s;
  double P0 = 0.5 * eps0 * c0 * (1*V/m) * (1*V/m);
  double scale = _power/P0;

  // interpolate to mesh node
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);

    double eta;
    if (_eta_auto)
      eta = floor(h*c0/_wave_length / region->get_optical_Eg(region->T_external()));
    else
      eta = _eta;

    Complex r = region->get_optical_refraction(_wave_length);
    Complex eps(r.real()*r.real()-r.imag()*r.imag(), -2*r.real()*r.imag()) ;

    SimulationRegion::processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = (*it);
      FVM_NodeData * node_data = fvm_node->node_data();

      double E_field = interpolator->get_interpolated_value(*(fvm_node->root_node()), 0);
      node_data->OptG() += eta*M_PI/h*eps0*(-eps.imag())* E_field *E_field * scale;
    }
  }

  delete interpolator;
}



void Light_Source_Uniform::update_system()
{

  double g = _card.get_real("optical.gen", 0.0)*pow(cm, -3)/s;

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);

    SimulationRegion::processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = (*it);
      FVM_NodeData * node_data = fvm_node->node_data();

      node_data->OptG() += g;
    }
  }

}


//-------------------------------------------------------------------------------------------

#include "ray_tracing/ray_tracing.h"
void Light_Source_RayTracing::update_system()
{
  RayTraceSolver * solver = new RayTraceSolver(_system, _card);
  solver->create_solver();
  solver->solve();
  solver->destroy_solver();
  delete solver;
}


//-------------------------------------------------------------------------------------------

#include "emfem2d/emfem2d.h"
void Light_Source_EMFEM2D::update_system()
{
  EMFEM2DSolver * solver = new EMFEM2DSolver(_system, _card);
  solver->create_solver();
  solver->solve();
  solver->destroy_solver();
  delete solver;
}
