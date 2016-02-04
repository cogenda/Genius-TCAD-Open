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

#include <limits>
#include <fstream>

#include "parser.h"
#include "mesh_base.h"
#include "particle_source.h"
#include "simulation_system.h"
#include "simulation_region.h"
#include "semiconductor_region.h"
#include "insulator_region.h"
#include "interpolation_2d_csa.h"
//#include "interpolation_3d_qshep.h"
#include "interpolation_3d_nbtet.h"
#include "nearest_node_locator.h"
#include "object_tree.h"
#include "parallel.h"
#include "mathfunc.h"
#include "solver_specify.h"
#include "log.h"

//#define DEBUG

using PhysicalUnit::s;
using PhysicalUnit::um;
using PhysicalUnit::eV;
using PhysicalUnit::cm;
using PhysicalUnit::g;






void Particle_Source::carrier_generation(double t)
{
  double ct = 0.5*(carrier_generation_t(t+0.5*SolverSpecify::dt) + carrier_generation_t(t-0.5*SolverSpecify::dt));

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

      node_data->PatG() += _fvm_node_particle_deposit[fvm_node]*ct;
    }
  }
}


double Particle_Source::limit_dt(double time, double dt, double dt_min) const
{
  const double total = 0.5*sqrt(3.1415926536)*_t_char*(1 - Erf((_t0-_t_max)/_t_char));
  bool before = time+dt < _t_max;
  bool after  = time > _t_max;
  if( (before || after) && carrier_generation_t(time) < 0.001*total && carrier_generation_t(time+dt) < 0.001*total ) return dt;

  do{
    const double interval = 0.5*sqrt(3.1415926536)*_t_char*(Erf((time+dt-_t_max)/_t_char) - Erf((time-_t_max)/_t_char));
    double error = std::abs(0.5*(carrier_generation_t(time) + carrier_generation_t(time+dt))*dt - interval);
    if( error < 0.001*total || error < 0.001*interval || dt < 0.05*_t_char || dt < dt_min) break;
    dt *= 0.9;
  }
  while(1);

  return dt;
}



double Particle_Source::quan_eff(const SimulationRegion * region) const
{
  switch(region->type())
  {
    case SemiconductorRegion:
    {
      const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
      return semi_region->material()->band->ParticleQuantumEffect(_system.T_external());
    }
    case InsulatorRegion:
    {
      const InsulatorSimulationRegion * insu_region = dynamic_cast<const InsulatorSimulationRegion *>(region);
      return insu_region->material()->band->ParticleQuantumEffect(_system.T_external());
    }
    default: return 3.6*eV;
  }

  return 3.6*eV;
}



double Particle_Source::carrier_generation_t(double t) const
{
  if( t>= _t0 && (t-_t_max)*(t-_t_max)/(_t_char*_t_char)<30)
    return exp(-(t-_t_max)*(t-_t_max)/(_t_char*_t_char));
  return 0;
}

//-------------------------------------------------------------------------------------------------------------------------

Particle_Source_DataFile::Particle_Source_DataFile(SimulationSystem &system, const Parser::Card &c):Particle_Source(system)
{
  MESSAGE<<"Setting Radiation Source from data file..."; RECORD();

  genius_assert(c.key() == "PARTICLE");

  // set Particle profile here.
  if(c.is_enum_value("profile","fromfile2d"))
    set_particle_profile_fromfile2d(c);

  if(c.is_enum_value("profile","fromfile3d"))
    set_particle_profile_fromfile3d(c);

  _t0     = c.get_real("t0", 0.0)*s;
  _t_max  = c.get_real("tmax", 0.0)*s;
  _t_char = c.get_real("t.char", 2e-12)*s;


  MESSAGE<<"ok\n"<<std::endl; RECORD();

}



void Particle_Source_DataFile::set_particle_profile_fromfile2d(const Parser::Card &c)
{
  interpolator = AutoPtr<InterpolationBase>(new Interpolation2D_CSA);
  interpolator->set_interpolation_type(0, InterpolationBase::Asinh);

  if(Genius::processor_id()==0)
  {
    VectorValue<double> translate;
    TensorValue<double> transform;
    if( c.is_parameter_exist("translate") )
    {
      std::vector<double> dummy = c.get_array<double>("translate");
      translate = VectorValue<double>(&dummy[0])*um;
    }
    else
    {
      translate = VectorValue<double>(c.get_real("translate.x", 0.0)*um,
                                      c.get_real("translate.y", 0.0)*um,
                                      c.get_real("translate.z", 0.0)*um);
    }

    if( c.is_parameter_exist("transform") )
    {
      std::vector<double> dummy = c.get_array<double>("transform");
      transform = TensorValue<double>(&dummy[0]);
    }
    else
    {
      transform = TensorValue<double>(c.get_real("transform.xx", 1.0), c.get_real("transform.xy", 0.0), c.get_real("transform.xz", 0.0),
                                      c.get_real("transform.yx", 0.0), c.get_real("transform.yy", 1.0), c.get_real("transform.yz", 0.0),
                                      c.get_real("transform.zx", 0.0), c.get_real("transform.zy", 0.0), c.get_real("transform.zz", 1.0));
    }

    std::string file = c.get_string("profile.file", "");
    std::ifstream in(file.c_str());
    if(!in.good())
    {
      MESSAGE<<"ERROR at " << c.get_fileline() <<" PARTICLE: file "<<file<<" can't be opened."<<std::endl; RECORD();
      genius_error();
    }

    Point p;
    double energy;

    while(!in.eof())
    {
      in >> p[0];// read x location
      in >> p[1];// read y location
      p *= um;     // scale to um
      p = transform*p + translate; //do transform & translate
      in >> energy;
      energy *= eV/(pow(um,3));
      if(!in.fail())
        interpolator->add_scatter_data(p, 0, energy);
    }
    if(!in.eof())
    {
      MESSAGE<<"ERROR at " << c.get_fileline() <<": error reading data file " << file << std::endl; RECORD();
      genius_error();
    }
    in.close();
  }

  interpolator->broadcast(0);
  interpolator->setup(0);
}


void Particle_Source_DataFile::set_particle_profile_fromfile3d(const Parser::Card &c)
{
  //interpolator = AutoPtr<InterpolationBase>(new Interpolation3D_qshep);
  interpolator = AutoPtr<InterpolationBase>(new Interpolation3D_nbtet);


  interpolator->set_interpolation_type(0, InterpolationBase::Asinh);


  VectorValue<double> translate;
  TensorValue<double> transform;
  if( c.is_parameter_exist("translate") )
  {
    std::vector<double> dummy = c.get_array<double>("translate");
    translate = VectorValue<double>(&dummy[0])*um;
  }
  else
  {
    translate = VectorValue<double>(c.get_real("translate.x", 0.0)*um,
                                    c.get_real("translate.y", 0.0)*um,
                                    c.get_real("translate.z", 0.0)*um);
  }

  if( c.is_parameter_exist("transform") )
  {
    std::vector<double> dummy = c.get_array<double>("transform");
    transform = TensorValue<double>(&dummy[0]);
  }
  else
  {
    transform = TensorValue<double>(c.get_real("transform.xx", 1.0), c.get_real("transform.xy", 0.0), c.get_real("transform.xz", 0.0),
                                    c.get_real("transform.yx", 0.0), c.get_real("transform.yy", 1.0), c.get_real("transform.yz", 0.0),
                                    c.get_real("transform.zx", 0.0), c.get_real("transform.zy", 0.0), c.get_real("transform.zz", 1.0));
  }

  if(Genius::processor_id()==0)
  {
    std::string file = c.get_string("profile.file", "");
    std::ifstream in(file.c_str());
    if(!in.good())
    {
      MESSAGE<<"ERROR at " << c.get_fileline() <<" PARTICLE: file "<<file<<" can't be opened."<<std::endl; RECORD();
      genius_error();
    }

    Point p;
    double energy;

    while(!in.eof())
    {
      in >> p[0];// read x location
      in >> p[1];// read y location
      in >> p[2];// read z location
      p *= um;     // scale to um
      p = transform*p + translate; //do transform & translate
      in >> energy;
      energy *= eV/(pow(um,3));
      if(!in.fail())
        interpolator->add_scatter_data(p, 0, energy);
    }
    if(!in.eof())
    {
      MESSAGE<<"ERROR at " << c.get_fileline() <<": error reading data file " << file << std::endl; RECORD();
      genius_error();
    }
    in.close();
  }
  interpolator->broadcast(0);
  interpolator->setup(0);
}


void Particle_Source_DataFile::update_source()
{
  // interpolate to mesh node
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    const SimulationRegion * region = _system.region(n);

    if( region->type() != SemiconductorRegion ) continue;
    double _quan_eff = quan_eff(region);

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = (*it);

      double E = interpolator->get_interpolated_value(*(fvm_node->root_node()), 0);
      _fvm_node_particle_deposit[fvm_node] = 2*E/_quan_eff/_t_char/sqrt(3.1415926536)/(1+Erf((_t_max-_t0)/_t_char));
    }
  }
}



//-------------------------------------------------------------------------------------------------------------------------

Particle_Source_Analytic::Particle_Source_Analytic(SimulationSystem &system, const Parser::Card &c):Particle_Source(system)
{
  MESSAGE<<"Setting Radiation Source from analytic expression..."<<std::endl; RECORD();

  genius_assert(c.key() == "PARTICLE");

  _start.x() = c.get_real("start.x", 0.0)*um;
  _start.y() = c.get_real("start.y", 0.0)*um;
  _start.z() = c.get_real("start.z", 0.0)*um;

  if( c.is_parameter_exist("dir.x") || c.is_parameter_exist("dir.y") || c.is_parameter_exist("dir.z") )
  {
    _dir.x() = c.get_real("dir.x", 0.0);
    _dir.y() = c.get_real("dir.y", 0.0);
    _dir.z() = c.get_real("dir.z", 0.0);
  }
  else
  {
    // ray direction
    double phi   = c.get_real("k.phi",   0.0)/180.0*3.14159265358979323846;
    double theta = c.get_real("k.theta", 0.0)/180.0*3.14159265358979323846;
    _dir = Point(sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta));
  }

  _dir.to_unit();

  _t0     = c.get_real("t0", 0.0)*s;
  _t_max  = c.get_real("tmax", 0.0)*s;
  _t_char = c.get_real("t.char", 2e-12)*s;


  _lateral_char = c.get_real("lateral.char", 0.1)*um;
  _length = c.get_real("length", 50.0)*um;

  if( c.is_parameter_exist("dedx") )
    _dEdx = c.get_real("let", 8.0)*1e6*eV/um;
  else if ( c.is_parameter_exist("let") )
    _dEdx = c.get_real("let", 35.0)*1e6*eV*cm*cm/(g/1e3)*(2.32*g/(cm*cm*cm));
  else
    _dEdx = 0.0;
}





void Particle_Source_Analytic::update_source()
{
  const double pi = 3.1415926536;
  bool is_2d = (_system.mesh().mesh_dimension() == 2) && (!_system.cylindrical_mesh());
  double z_width = is_2d ? _system.z_width() : 1.0;

  // make sure the track not exceeds the device domain
  bool hit;
  std::pair<double, double> t;
  if(Genius::processor_id()==0)
  {
    // processor 0 has all the mesh elements, infact, boundary elements should exist on all the processors
    ObjectTree * surface_elem_tree = new ObjectTree(_system.mesh(), Trees::ELEMENTS_ON_BOUNDARY);
    hit = surface_elem_tree->hit_domain(_start, _dir, t);
    delete surface_elem_tree;
  }

  Parallel::broadcast(hit);
  Parallel::broadcast(t.first);
  Parallel::broadcast(t.second);

  MESSAGE<<"Process carrier generation from particle analytic expression..."<<std::endl; RECORD();
  if(hit && t.second>0.0)
  {
    if(t.first > 0.0)
    {
      _start = _start + t.first*_dir;
      _length = std::min(t.second - t.first, _length);

      MESSAGE<<"  Particle dEdx " << _dEdx/(1e6*eV/um) << " MeV/um," << std::endl; RECORD();
      MESSAGE<<"  Particle hit device at " << _start/um << " um." << std::endl; RECORD();
      MESSAGE<<"  Particle effective length limited to " << _length/um << " um due to the size of device." <<std::endl; RECORD();
    }
    else
    {
      _length = std::min(t.second, _length);

      MESSAGE<<"  Particle dEdx " << _dEdx/(1e6*eV/um) << " MeV/um" << std::endl; RECORD();
      MESSAGE<<"  Particle start at " << _start/um << " um." << std::endl; RECORD();
      MESSAGE<<"  Particle effective length limited to " << _length/um << " um due to the size of device." <<std::endl; RECORD();
    }
  }
  else
  {
    MESSAGE<<"  Particle does not hit device, skip." << std::endl; RECORD();
    return;
  }


  const double dEdx = _dEdx;
  const double total_energy = dEdx*_length;

  //used for energy conservation
  double track_energy = 0.0;
  std::map<const FVM_Node *, double> track_energy_density;

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    const SimulationRegion * region = _system.region(n);

    double lateral_norm = pi*_lateral_char*_lateral_char;

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = (*it);

      Point loc = fvm_node->position();
      Point loc_pp = _start + (loc-_start)*_dir*_dir;

      double r = (loc-loc_pp).size();
      double e_r = exp(-r*r/(_lateral_char*_lateral_char));
      double e_z = Erf((loc_pp-_start)*_dir/_lateral_char) - Erf((loc_pp-(_start+_length*_dir))*_dir/_lateral_char);
      double energy_density = _dEdx/(2*lateral_norm)*e_r*e_z;

      track_energy += energy_density*fvm_node->volume()*z_width;
      track_energy_density[fvm_node] += energy_density;
    }
  }

  Parallel::sum(track_energy);


  if(track_energy > 0.0)
  {
    double alpha = total_energy/track_energy; //used for keep energy conservation

    std::map<const FVM_Node *, double>::const_iterator it = track_energy_density.begin();
    for(; it != track_energy_density.end(); ++it)
    {
      const FVM_Node * fvm_node = it->first;
      if(!fvm_node) continue;

      const SimulationRegion * region = _system.region(fvm_node->subdomain_id());
      if( region->type() != SemiconductorRegion) continue;

      double energy_density = it->second;
      double _quan_eff = quan_eff(region);

      _fvm_node_particle_deposit[fvm_node] += alpha*energy_density/_quan_eff/(_t_char/2.0*sqrt(pi)*(1+Erf((_t_max-_t0)/_t_char)));
    }
  }

}



//-------------------------------------------------------------------------------------------------------------------------

Particle_Source_Track::Particle_Source_Track(SimulationSystem &system, const Parser::Card &c):Particle_Source(system)
{
  genius_assert(c.key() == "PARTICLE");
  genius_assert(c.is_enum_value("profile", "track"));

  // get parameters for the particle track
  _t0     = c.get_real("t0", 0.0)*s;
  _t_max  = c.get_real("tmax", 0.0)*s;
  _t_char = c.get_real("t.char", 2e-12)*s;


  std::string track_file = c.get_string("profile.file", "");
  std::string hdf5_file = c.get_string("profile.hdf5", "");
  std::string event_path = c.get_string("hdf5.path", "");

  MESSAGE<<"Setting Radiation Source from particle event file " << track_file << "..."; RECORD();
  if(!track_file.empty())
    _read_particle_profile_track_evt(track_file, c.get_real("lateral.char", 0.1));
#ifdef HAVE_HDF5
  if(!hdf5_file.empty())
    _read_particle_profile_track_hdf5(hdf5_file, event_path, c.get_real("lateral.char", 0.1));
#endif
  MESSAGE<<"ok\n"<<std::endl; RECORD();
}



void Particle_Source_Track::_read_particle_profile_track_evt(const std::string & filename, Real lateral_char)
{
  std::vector<double> meta_data;
  if(Genius::processor_id()==0)
  {
    std::ifstream in(filename.c_str());
    if(!in.good())
    {
      MESSAGE<<"ERROR PARTICLE: file "<<filename<<" can't be opened."<<std::endl; RECORD();
      genius_error();
    }

    while(!in.eof())
    {
      std::string context;
      std::getline(in, context);
      if(context.empty()) continue;
      if(context[0] == '#') continue;

      std::stringstream ss(context);
      std::string particle;
      ss >> particle;
      if(particle.empty()) continue;

      double p1[3], p2[3], energy, sigma=lateral_char;
      ss >>  p1[0] >> p1[1] >> p1[2] >> p2[0] >> p2[1] >> p2[2] >> energy;
      if(!ss.eof()) { ss >> sigma;}

      // skip track with very little energy
      // a very low energy such as 1.38073e-315 may break some floating point compare result
      if(energy < 1e-12) continue;
      // skip very short track
      Point begin(p1[0]*um,p1[1]*um,p1[2]*um);
      Point end(p2[0]*um,p2[1]*um,p2[2]*um);
      if( (begin-end).size() < 1e-6*um ) continue;

      meta_data.push_back(p1[0]*um);
      meta_data.push_back(p1[1]*um);
      meta_data.push_back(p1[2]*um);
      meta_data.push_back(p2[0]*um);
      meta_data.push_back(p2[1]*um);
      meta_data.push_back(p2[2]*um);
      meta_data.push_back(energy*1e6*eV);
      meta_data.push_back(sigma*um);
    }

    in.close();
  }

  Parallel::broadcast(meta_data);
  for(unsigned int n=0; n<meta_data.size()/8; ++n)
  {
    track_t track;
    track.start.x()    = meta_data[8*n+0];
    track.start.y()    = meta_data[8*n+1];
    track.start.z()    = meta_data[8*n+2];
    track.end.x()      = meta_data[8*n+3];
    track.end.y()      = meta_data[8*n+4];
    track.end.z()      = meta_data[8*n+5];
    track.energy       = meta_data[8*n+6];
    track.lateral_char = meta_data[8*n+7];
    _tracks.push_back(track);
  }

  Parallel::verify(_tracks.size());

}

#ifdef HAVE_HDF5

#include "ParticleEvent.h"
void Particle_Source_Track::_read_particle_profile_track_hdf5(const std::string & filename, const std::string & path, Real lateral_char)
{
  std::vector<double> meta_data;
  if(Genius::processor_id()==0)
  {
    hid_t    file_handle;

    file_handle = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_handle<0)
    {
      MESSAGE<<"ERROR PARTICLE: file "<<filename<<" can't be opened."<<std::endl; RECORD();
      genius_error();
    }

    H5G_info_t ginfo;
    H5Gget_info(file_handle, &ginfo);

    CogendaHDF5::LookupTable region_lut;
    std::vector<std::string> material_list;
    {
      hid_t grp;
      grp = H5Gopen(file_handle, "target", H5P_DEFAULT);
      material_list = CogendaHDF5::getAttribute< std::vector<std::string> >(grp, ".", "MaterialList");
      region_lut   = CogendaHDF5::getAttribute< CogendaHDF5::LookupTable >(grp, ".", "RegionList");
    }

    //for(unsigned int n=0; n<material_list.size(); ++n)
    //  std::cout<<material_list[n]<<std::endl;

    //std::string event_path = "./" + path;
    CogendaHDF5::GSeat::Event evt_hdf5(region_lut, material_list);
    bool rc = evt_hdf5.readData(file_handle, path);
    if (!rc)
    {
      MESSAGE<<"ERROR PARTICLE: file "<<filename<<" path " << path << " can not be opened." <<std::endl; RECORD();
      genius_error();
    }

    for (CogendaHDF5::GSeat::Event::TrackIter iTrack = evt_hdf5.track_begin();  iTrack != evt_hdf5.track_end(); ++iTrack)
    {
      const CogendaHDF5::GSeat::Track & track_hdf5 = *iTrack;

      std::vector<double> p = track_hdf5.StartPoint();
      Point StartPoint(p[0], p[1], p[2]);

      CogendaHDF5::GSeat::Event::StepIter iStep = evt_hdf5.step_begin(track_hdf5);
      for (unsigned int c=0; iStep != evt_hdf5.step_end(track_hdf5); ++iStep, ++c)
      {
        const CogendaHDF5::GSeat::Step & step_hdf5 = *iStep;

        std::vector<double> p = step_hdf5.EndPoint();
        Point EndPoint(p[0], p[1], p[2]);
        double energy = step_hdf5.EnergyDeposit();

        meta_data.push_back(StartPoint[0]*um);
        meta_data.push_back(StartPoint[1]*um);
        meta_data.push_back(StartPoint[2]*um);
        meta_data.push_back(EndPoint[0]*um);
        meta_data.push_back(EndPoint[1]*um);
        meta_data.push_back(EndPoint[2]*um);
        meta_data.push_back(energy*1e6*eV);

        double sigma=lateral_char;
        meta_data.push_back(sigma*um);

        StartPoint = EndPoint;
      }
    }

    H5Fclose(file_handle);
  }

  //std::cout<<meta_data.size() << std::endl;

  Parallel::broadcast(meta_data);
  for(unsigned int n=0; n<meta_data.size()/8; ++n)
  {
    track_t track;
    track.start.x()    = meta_data[8*n+0];
    track.start.y()    = meta_data[8*n+1];
    track.start.z()    = meta_data[8*n+2];
    track.end.x()      = meta_data[8*n+3];
    track.end.y()      = meta_data[8*n+4];
    track.end.z()      = meta_data[8*n+5];
    track.energy       = meta_data[8*n+6];
    track.lateral_char = meta_data[8*n+7];

    // skip zero length track
    if( track.energy <= 0.0 || (track.end - track.start).size() == 0.0 ) continue;

    _tracks.push_back(track);
  }

  Parallel::verify(_tracks.size());

}

#endif



void Particle_Source_Track::update_source()
{
  START_LOG("update_source()", "Particle_Source_Track");

  MESSAGE<< "Process carrier generation from particle track";
  RECORD();

  const double pi = 3.1415926536;
  genius_assert(_system.mesh().mesh_dimension() == 3);

  std::vector<double> region_energy(_system.n_regions(), 0.0);
  double total_energy=0.0;

  AutoPtr<NearestNodeLocator> nn_locator( new NearestNodeLocator(_system.mesh()) );

  for(unsigned int t=0; t<_tracks.size(); ++t)
  {
    if( t%(1+_tracks.size()/20) ==0 )
    {
      MESSAGE<< ".";
      RECORD();
    }

    const track_t & track = _tracks[t];
    genius_assert(track.energy > 0.0 && (track.end - track.start).size() > 0.0);

    /*
    // prevent bad tracks
    int bad_track = 0;
    if( track.energy == 0.0 || (track.end - track.start).size() == 0.0) bad_track=1;
    Parallel::broadcast(bad_track);
    if(bad_track) continue;
    */

    double track_energy = 0.0;
    std::map<const FVM_Node *, double> track_energy_density;

    const Point track_dir = (track.end - track.start).unit(); // track direction
    const double dEdx = track.energy/(track.end - track.start).size(); // linear energy density
    const double lateral_char = track.lateral_char;
    // find the nodes that near the track
    for(unsigned int r=0; r<_system.n_regions(); r++)
    {
      const SimulationRegion * region = _system.region(r);

      // fast return
      const std::pair<Point, Real> bsphere = region->boundingsphere();
      const Point cent = 0.5*(track.start+track.end);
      const double diag = 0.5*(track.start-track.end).size();
      if( (bsphere.first - cent).size() > bsphere.second + diag + 5*lateral_char ) continue;


      std::vector<const Node *> nn = nn_locator->nearest_nodes(track.start, track.end, 5*lateral_char, r);
      for(unsigned int n=0; n<nn.size(); ++n)
      {
        Point loc = *nn[n];
        const FVM_Node * fvm_node = region->region_fvm_node(nn[n]); // may be NULL, if not on local
        if(!fvm_node || !fvm_node->on_processor()) continue;

        Point loc_pp = track.start + (loc-track.start)*track_dir*track_dir;
        Real r = (loc-loc_pp).size();
        double e_r = exp(-r*r/(lateral_char*lateral_char));
        double e_z = Erf((loc_pp-track.start)*track_dir/lateral_char) - Erf((loc_pp-track.end)*track_dir/lateral_char);
        double energy_density = dEdx/(2*pi*lateral_char*lateral_char)*e_r*e_z;
        track_energy += energy_density*fvm_node->volume();
        track_energy_density[fvm_node] += energy_density;
      }
    }

    Parallel::sum(track_energy);

    if(track_energy > 0.0)
    {
      double alpha = track.energy/track_energy; //used for keep energy conservation track.energy;
      std::map<const FVM_Node *, double>::const_iterator it = track_energy_density.begin();
      for(; it != track_energy_density.end(); ++it)
      {
        const FVM_Node * fvm_node = it->first;
        if(!fvm_node) continue;

        double energy_density = it->second;

        const SimulationRegion * region = _system.region(fvm_node->subdomain_id());
        // if( region->type() != SemiconductorRegion) continue;
        double _quan_eff = quan_eff(region);

        if(fvm_node && fvm_node->on_local())
        {
          _fvm_node_particle_deposit[fvm_node] += alpha*energy_density/_quan_eff/(_t_char/2.0*sqrt(pi)*(1+Erf((_t_max-_t0)/_t_char)));
          //node_data->PatE() += alpha*energy_density;
          total_energy += alpha*energy_density*fvm_node->volume();
          region_energy[ fvm_node->subdomain_id() ] += alpha*energy_density*fvm_node->volume();
        }
      }
    }
    else
    {
      // find the nearest node of the track
      const FVM_Node * neraset_node=0;
      double distance=std::numeric_limits<double>::infinity();
      for(unsigned int r=0; r<_system.n_regions(); r++)
      {
        const SimulationRegion * region = _system.region(r);
        double dist;
        const Node * n = nn_locator->nearest_node(0.5*(track.start+track.end), r, dist);
        if(n == NULL) continue;

        const FVM_Node * fvm_node = region->region_fvm_node(n); // may be NULL, if not on local
        if(!fvm_node || !fvm_node->on_processor()) continue;

        if( dist < distance)
        {
          distance = dist;
          neraset_node = fvm_node;
        }
      }

      double min_distance=distance;
      Parallel::min(min_distance);

      if( min_distance == distance && neraset_node)
      {
        const SimulationRegion * region = _system.region(neraset_node->subdomain_id());
        double _quan_eff = quan_eff(region);

        double energy_density = track.energy/neraset_node->volume();
        _fvm_node_particle_deposit[neraset_node] += energy_density/_quan_eff/(_t_char/2.0*sqrt(pi)*(1+Erf((_t_max-_t0)/_t_char)));
        //node_data->PatE() += energy_density;
        total_energy += track.energy;
        region_energy[neraset_node->subdomain_id()] += track.energy;
      }
    }

  }

  MESSAGE<< "ok" <<std::endl;
  RECORD();

  STOP_LOG("update_source()", "Particle_Source_Track");

}



#if 0
void Particle_Source_Track::update_source()
{
  const double pi = 3.1415926536;
  genius_assert(_system.mesh().mesh_dimension() == 3);

  std::map<std::string, double> region_energy;
  double total_energy=0.0;
  int lost_energy=0;

  AutoPtr<NearestNodeLocator> nn_locator( new NearestNodeLocator(_system.mesh()) );

  // for each track
  for(unsigned int t=0; t<_tracks.size(); ++t)
  {
    const track_t & track = _tracks[t];
    // prevent bad tracks
    if( track.energy == 0.0 ) continue;
    if((track.end - track.start).size() == 0.0) continue;

    const Point track_dir = (track.end - track.start).unit(); // track direction
    const double ed = track.energy/(track.end - track.start).size(); // linear energy density
    const double lateral_char = track.lateral_char;
    // find the nodes that near the track
    std::multimap<const Node *, std::pair<FVM_Node *, double> > impacted_fvm_nodes;
    for(unsigned int r=0; r<_system.n_regions(); r++)
    {
      const SimulationRegion * region = _system.region(r);

      unsigned int try_radii = 0;
      Real radii = 5*lateral_char;
      std::vector<const Node *> nn = nn_locator->nearest_nodes(track.start, track.end, radii, r);
      while(nn.size() < 10 && try_radii++ < 3)
      {
        radii *= 3.0;
        nn = nn_locator->nearest_nodes(track.start, track.end, radii, r);
      }

      for(unsigned int n=0; n<nn.size(); ++n)
      {
        Point loc = *nn[n];
        FVM_Node * fvm_node = region->region_fvm_node(nn[n]); // may be NULL, if not on local
        if(!fvm_node || !fvm_node->on_processor()) continue;

        Point loc_pp = track.start + (loc-track.start)*track_dir*track_dir;
        Real r = (loc-loc_pp).size();
        double e_r = exp(-r*r/(lateral_char*lateral_char));
        double e_z = Erf((loc_pp-track.start)*track_dir/lateral_char) - Erf((loc_pp-track.end)*track_dir/lateral_char);
        double energy = ed/(2*pi*lateral_char*lateral_char)*e_r*e_z;
        impacted_fvm_nodes.insert(std::make_pair(nn[n], std::make_pair(fvm_node, energy)));
      }
    }


    // statistic total energy deposit
    double alpha = 1.0; //used for keep energy conservation
    {
      double energy_statistics = 0.0;
      std::multimap<const Node *, std::pair<FVM_Node *, double> >::iterator it=impacted_fvm_nodes.begin();
      for( ; it!=impacted_fvm_nodes.end(); ++it)
      {
        const FVM_Node * fvm_node = it->second.first;
        if(fvm_node && fvm_node->on_processor())
          energy_statistics += it->second.second*fvm_node->volume();
      }
      Parallel::sum(energy_statistics);
      if(energy_statistics == 0.0) lost_energy ++;
      alpha = energy_statistics > 0.0 ? track.energy/energy_statistics : 1.0; //used for keep energy conservation
    }

    std::multimap<const Node *, std::pair<FVM_Node *, double> >::iterator it=impacted_fvm_nodes.begin();
    for( ; it!=impacted_fvm_nodes.end(); ++it)
    {
      FVM_Node * fvm_node = it->second.first;
      if(!fvm_node) continue;

      SimulationRegion * region = _system.region(fvm_node->subdomain_id());
      //if( region->type() != SemiconductorRegion) continue;

      if(fvm_node && fvm_node->on_local())
      {
        FVM_NodeData * node_data = fvm_node->node_data();
        node_data->PatG() += alpha*it->second.second/_quan_eff/(_t_char/2.0*sqrt(pi)*(1+Erf((_t_max-_t0)/_t_char)));
        node_data->PatE() += alpha*it->second.second;//*fvm_node->volume();
        total_energy += alpha*it->second.second*fvm_node->volume();
        region_energy[ _system.region(fvm_node->subdomain_id())->name() ] += alpha*it->second.second*fvm_node->volume();
      }
    }
  }

  std::map<std::string, double>::const_iterator it = region_energy.begin();
  for(;it != region_energy.end(); ++it)
    std::cout<<it->first << " " << it->second<<std::endl;
  std::cout<<total_energy<< " " << lost_energy << std::endl;
}
#endif

