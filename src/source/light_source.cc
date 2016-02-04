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
#include "semiconductor_region.h"
#include "interpolation_2d_csa.h"
#include "interpolation_2d_nn.h"
//#include "interpolation_3d_qshep.h"
#include "interpolation_3d_nbtet.h"
#include "light_source.h"
#include "solver_specify.h"
#include "material.h"
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



void Light_Source::carrier_generation(double t)
{
  double optical_gen_waveform = 1.0;
  if(SolverSpecify::TimeDependent)
  {
    if(_global_waveform)
      optical_gen_waveform = 0.5*(_global_waveform->waveform(t) + _global_waveform->waveform(t-SolverSpecify::dt));
    else if(_waveform)
      optical_gen_waveform = 0.5*(_waveform->waveform(t) + _waveform->waveform(t-SolverSpecify::dt));
  }

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);

    if( region->type() != SemiconductorRegion ) continue;

    SimulationRegion::processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = (*it);
      FVM_NodeData * node_data = fvm_node->node_data();

      node_data->OptG() += _fvm_node_particle_deposit[fvm_node]*optical_gen_waveform;
    }
  }
}

double Light_Source::limit_dt(double time, double dt, double dt_min) const
{
  if(_global_waveform)
  {
    do
    {
      double a1 = _global_waveform->waveform(time);
      double am = _global_waveform->waveform(time+0.5*dt);
      double a2 = _global_waveform->waveform(time+dt);
      if( 0.5*fabs(a1+a2)<1e-3 || fabs(am-0.5*a1-0.5*a2)<0.02*fabs(a1+a2) ) break;
      if( dt < dt_min ) break;
    } while( dt*=0.9 );
    
    _global_waveform->dt_critial_limit(time, dt, dt_min);
  }
  else if(_waveform)
  {
    do
    {
      double a1 = _waveform->waveform(time);
      double am = _waveform->waveform(time+0.5*dt);
      double a2 = _waveform->waveform(time+dt);
      if( 0.5*fabs(a1+a2)<1e-3 || fabs(am-0.5*a1-0.5*a2)<0.02*fabs(a1+a2) ) break;
      if( dt < dt_min ) break;
    } while( dt*=0.9 );
    
    _waveform->dt_critial_limit(time, dt, dt_min);
  }
  return dt;
}

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

  if(c.is_enum_value("profile","efile2d") || c.is_enum_value("profile","efile3d")) _field_type="efield";
  if(c.is_enum_value("profile","pfile2d") || c.is_enum_value("profile","pfile3d")) _field_type="pfield";


  _fname = c.get_string("data.file", "");
  if (fname_ext.length()>0)
  {
    _fname += ".";
    _fname += fname_ext;
  }

  _skip_line      = c.get_int("skipline", 0);

  // set Light profile here.
  if(c.is_enum_value("profile","efile2d") || c.is_enum_value("profile","pfile2d"))
    _dim = 2;
  if(c.is_enum_value("profile","efile3d") || c.is_enum_value("profile","pfile3d"))
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

int Light_Source_From_File::load_light_elec_profile_fromfile(InterpolationBase * interpolator, const std::string &fname, const int skip_line)
{

  if(Genius::processor_id()==0)
  {
    std::ifstream in(fname.c_str());
    if(!in.good())
    {
      MESSAGE << "Light Source Error: file "<<fname<<" can't be opened."<<std::endl; RECORD();
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


int Light_Source_From_File::load_light_pow_profile_fromfile(InterpolationBase * interpolator, const std::string &fname, const int skip_line)
{

  if(Genius::processor_id()==0)
  {
    std::ifstream in(fname.c_str());
    if(!in.good())
    {
      MESSAGE << "Light Source Error: file "<<fname<<" can't be opened."<<std::endl; RECORD();
      genius_error();
    }

    Point p;
    double Power;

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
      in >> Power;
      Power *= J/pow(PhysicalUnit::m*_FUnit, 3.0);
      if(!in.fail())
      {
        interpolator->add_scatter_data(p, 0, Power);
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


void Light_Source_From_File::update_source()
{

  // load lightprofile from File and fill into the interpolator
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
    if(_field_type=="efield")
      count = load_light_elec_profile_fromfile(interpolator, _fname, _skip_line);
    if(_field_type=="pfield")
      count = load_light_pow_profile_fromfile(interpolator, _fname, _skip_line);
    MESSAGE<<count<<" points ";
  }

  interpolator->broadcast(0);
  interpolator->setup(0);

  MESSAGE<<"ok"<<std::endl; RECORD();

  double pi = 3.1415927;
  double c0 = 299792458*m/s; // light speed
  double nu     = c0/_wave_length; // light frequency
  double photon = h*nu; //photon energy

  // Poynting vector E \cross H at E=1V/m(peak value), 0.5 is the factor for RMS value
  double P0 = 0.5 * eps0 * c0 * (1*V/m) * (1*V/m);

  double scale = _power/P0;



  // interpolate to mesh node
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    const SimulationRegion * region = _system.region(n);

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = (*it);

      double eta;
      if (_eta_auto)
        eta = floor(photon / region->get_optical_Eg(fvm_node));
      else
        eta = _eta;

      Complex r = region->get_optical_refraction(fvm_node, _wave_length); // n+ik
      Complex eps = r*r;// n^2-k^2 + 2nki

      if(_field_type=="efield")
      {
        double E_field = interpolator->get_interpolated_value(*(fvm_node->root_node()), 0);
        // optical energy is J \cdot E, here J=\sigma E,
        // and \sigma has the relationship with eps imag part as \eps^{''} = \frac{\sigma}{\omega}, here  \omega = 2\pi\nu
        // so J \cdot E can be expressed as J \cdot E = \sigma E^2 = 2\pi\nu\eps^{''} E^2
        double energy =   2*pi*nu*(eps0*(eps.imag())*0.5*E_field*E_field)*scale; //0.5 is peak value to RMS value converter
        _fvm_node_particle_deposit[fvm_node] = eta*energy/photon;
        //node_data->OptE() += energy;
      }

      if(_field_type=="pfield")
      {
        double Power = interpolator->get_interpolated_value(*(fvm_node->root_node()), 0);
        // Power = eps^{'}eps_0E^2, here E is RMS value
        double energy = 2*pi*nu*Power*eps.imag()/eps.real()*scale;
        _fvm_node_particle_deposit[fvm_node] = eta*energy/photon;
        //node_data->OptE() += energy;
      }
    }
  }

  delete interpolator;
}


//-------------------------------------------------------------------------------------------


void Light_Source_Uniform::update_source()
{
  if(_card.is_parameter_exist("optical.gen"))
  {
    double g = _card.get_real("optical.gen", 0.0)*pow(cm, -3)/s;

    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      const SimulationRegion * region = _system.region(n);

      SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        const FVM_Node * fvm_node = (*it);
        _fvm_node_particle_deposit[fvm_node] = g;
      }
    }
  }

  if(_card.is_parameter_exist("optical.power"))
  {
    double wave_length   = _card.get_real("wavelength", 0.532, "lambda") * um;   // wave length
    double power = _card.get_real("optical.power", 0.0)*J/(cm*cm*cm);
    double c0 = 299792458*m/s;
    double photon = h*c0/wave_length; //photon energy

    // interpolate to mesh node
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      const SimulationRegion * region = _system.region(n);

      SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        const FVM_Node * fvm_node = (*it);

        double eta = floor(photon / region->get_optical_Eg(fvm_node));

        Complex r = region->get_optical_refraction(fvm_node, wave_length); // n+ik
        Complex eps = r*r;// n^2-k^2 + 2nki

        _fvm_node_particle_deposit[fvm_node] = eta*(power*eps.imag()/eps.real())/photon/(1*s);
      }
    }
  }
}


//-------------------------------------------------------------------------------------------

void Light_Source_Xray::update_source()
{
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    const SimulationRegion * region = _system.region(n);
    if( region->type() != SemiconductorRegion ) continue;

    const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);

    Material::MaterialSemiconductor * m = semi_region->material();
    double density = m->basic->Density(region->T_external());
    double q = m->band->ParticleQuantumEffect(region->T_external());

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = (*it);
      _fvm_node_particle_deposit[fvm_node] = _doserate*density/q;
    }
  }
}




//-------------------------------------------------------------------------------------------

#ifdef TCAD_SOLVERS

#include "ray_tracing/ray_tracing.h"
void Light_Source_RayTracing::update_source()
{
  RayTraceSolver * solver = new RayTraceSolver(_system, _card);
  solver->create_solver();
  solver->solve();
  solver->destroy_solver();
  delete solver;


  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    const SimulationRegion * region = _system.region(n);

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = (*it);
      const FVM_NodeData * fvm_node_data = fvm_node->node_data();
      _fvm_node_particle_deposit[fvm_node] = fvm_node_data->OptG();
    }
  }
}


//-------------------------------------------------------------------------------------------

#include "emfem2d/emfem2d.h"
void Light_Source_EMFEM2D::update_source()
{
  EMFEM2DSolver * solver = new EMFEM2DSolver(_system, _card);
  solver->create_solver();
  solver->solve();
  solver->destroy_solver();
  delete solver;

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    const SimulationRegion * region = _system.region(n);

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = (*it);
      const FVM_NodeData * fvm_node_data = fvm_node->node_data();
      _fvm_node_particle_deposit[fvm_node] = fvm_node_data->OptG();
    }
  }
}

#else
void Light_Source_RayTracing::update_source() {}
void Light_Source_EMFEM2D::update_source() {}
#endif
