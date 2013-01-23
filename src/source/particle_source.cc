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
#include "interpolation_2d_csa.h"
//#include "interpolation_3d_qshep.h"
#include "interpolation_3d_nbtet.h"
#include "nearest_node_locator.h"
#include "parallel.h"
#include "mathfunc.h"
#include "log.h"

//#define DEBUG

using PhysicalUnit::s;
using PhysicalUnit::um;
using PhysicalUnit::eV;
using PhysicalUnit::cm;
using PhysicalUnit::g;




double Particle_Source::carrier_generation(double t) const
{
  if( t>= _t0 && (t-_t_max)*(t-_t_max)/(_t_char*_t_char)<30)
    return exp(-(t-_t_max)*(t-_t_max)/(_t_char*_t_char));
  return 0;
}

#if 0
double Particle_Source::limit_dt(double time, double dt) const
{
  const double total = 0.5*sqrt(3.1415926536)*_t_char*(1 - Erf((_t0-_t_max)/_t_char));
  while( 0.5*(carrier_generation(time) + carrier_generation(time+dt))*dt  > 0.05*total ) dt*=0.9;
  return dt;
}
#endif

double Particle_Source::limit_dt(double time, double dt) const
{
  const double total = 0.5*sqrt(3.1415926536)*_t_char*(1 - Erf((_t0-_t_max)/_t_char));
  if(carrier_generation(time) < 0.001*total && carrier_generation(time+dt) < 0.001*total ) return dt;

  do{
    const double interval = 0.5*sqrt(3.1415926536)*_t_char*(Erf((time+dt-_t_max)/_t_char) - Erf((time-_t_max)/_t_char));
    double error = std::abs(0.5*(carrier_generation(time) + carrier_generation(time+dt))*dt - interval);
    if( error < 0.001*total || error < 0.001*interval || dt < 0.05*_t_char) break;
    dt *= 0.9;
  }
  while(1);

  return dt;
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
  _quan_eff = c.get_real("quan.eff", 3.6)*eV;

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


void Particle_Source_DataFile::update_system()
{
  // interpolate to mesh node
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

      double E = interpolator->get_interpolated_value(*(fvm_node->root_node()), 0);
      node_data->PatG() += 2*E/_quan_eff/_t_char/sqrt(3.1415926536)/(1+Erf((_t_max-_t0)/_t_char));

    }
  }
}



//-------------------------------------------------------------------------------------------------------------------------

Particle_Source_Analytic::Particle_Source_Analytic(SimulationSystem &system, const Parser::Card &c):Particle_Source(system)
{
  MESSAGE<<"Setting Radiation Source from analytic expression..."; RECORD();

  genius_assert(c.key() == "PARTICLE");

  _start.x() = c.get_real("x", 0.0)*um;
  _start.y() = c.get_real("y", 0.0)*um;
  _start.z() = c.get_real("z", 0.0)*um;

  // ray direction
  double phi   = c.get_real("k.phi",   0.0)/180.0*3.14159265358979323846;
  double theta = c.get_real("k.theta", 0.0)/180.0*3.14159265358979323846;
  _dir = Point(sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta));

  _t0     = c.get_real("t0", 0.0)*s;
  _t_max  = c.get_real("tmax", 0.0)*s;
  _t_char = c.get_real("t.char", 2e-12)*s;
  _quan_eff = c.get_real("quan.eff", 3.6)*eV;

  _lateral_char = c.get_real("lateral.char", 0.1)*um;
  _length = c.get_real("length", 50.0)*um;

  _LET = c.get_real("let", 0.0)*1e6*eV*cm*cm/(g/1e3);

  MESSAGE<<"ok\n"<<std::endl; RECORD();
}





void Particle_Source_Analytic::update_system()
{
  const double pi = 3.1415926536;
  bool is_2d = _system.mesh().mesh_dimension() == 2;

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);

    if( region->type() != SemiconductorRegion ) continue;

    SemiconductorSimulationRegion * semi_region = dynamic_cast<SemiconductorSimulationRegion *>(region);
    double E = _LET*semi_region->material()->basic->Density(semi_region->T_external());

    double lateral_norm;
    if(is_2d)
      lateral_norm = sqrt(pi)*_lateral_char*1.0*um;
    else
      lateral_norm = pi*_lateral_char*_lateral_char;

    double G0 = E/_quan_eff/(lateral_norm)/(_t_char/2.0*sqrt(pi)*(1+Erf((_t_max-_t0)/_t_char)));

    SimulationRegion::processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = (*it);

      Point  p = (*fvm_node->root_node());
      double d = (p-_start).dot(_dir);
      double r = (p-_start - (p-_start).dot(_dir)*_dir).size();

      FVM_NodeData * node_data = fvm_node->node_data();
      if(d < _length)
        node_data->PatG() += G0*exp(-r*r/(_lateral_char*_lateral_char));
      else
        node_data->PatG() += 0;
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
  _quan_eff = c.get_real("quan.eff", 3.6)*eV;

  std::string track_file = c.get_string("profile.file", "");

  MESSAGE<<"Setting Radiation Source from particle event file " << track_file << "..."; RECORD();
  _read_particle_profile_track(track_file, c.get_real("lateral.char", 0.1));
  MESSAGE<<"ok\n"<<std::endl; RECORD();
}



void Particle_Source_Track::_read_particle_profile_track(const std::string & file, Real lateral_char)
{
  std::vector<double> meta_data;
  if(Genius::processor_id()==0)
  {
    std::ifstream in(file.c_str());
    if(!in.good())
    {
      MESSAGE<<"ERROR PARTICLE: file "<<file<<" can't be opened."<<std::endl; RECORD();
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







#if 1
void Particle_Source_Track::update_system()
{
  START_LOG("update_system()", "Particle_Source_Track");

  MESSAGE<< "  process particle generation";
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
    std::map<FVM_Node *, double> track_energy_density;

    const Point track_dir = (track.end - track.start).unit(); // track direction
    const double ed = track.energy/(track.end - track.start).size(); // linear energy density
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
        FVM_Node * fvm_node = region->region_fvm_node(nn[n]); // may be NULL, if not on local
        if(!fvm_node || !fvm_node->on_processor()) continue;

        Point loc_pp = track.start + (loc-track.start)*track_dir*track_dir;
        Real r = (loc-loc_pp).size();
        double e_r = exp(-r*r/(lateral_char*lateral_char));
        double e_z = Erf((loc_pp-track.start)*track_dir/lateral_char) - Erf((loc_pp-track.end)*track_dir/lateral_char);
        double energy_density = ed/(2*pi*lateral_char*lateral_char)*e_r*e_z;
        track_energy += energy_density*fvm_node->volume();
        track_energy_density[fvm_node] += energy_density;
      }
    }

    Parallel::sum(track_energy);

    if(track_energy > 0.0)
    {
      double alpha = track.energy/track_energy; //used for keep energy conservation track.energy;
      std::map<FVM_Node *, double>::const_iterator it = track_energy_density.begin();
      for(; it != track_energy_density.end(); ++it)
      {
        FVM_Node * fvm_node = it->first;
        if(!fvm_node) continue;

        double energy_density = it->second;

        SimulationRegion * region = _system.region(fvm_node->subdomain_id());
        //if( region->type() != SemiconductorRegion) continue;

        if(fvm_node && fvm_node->on_local())
        {
          FVM_NodeData * node_data = fvm_node->node_data();
          node_data->PatG() += alpha*energy_density/_quan_eff/(_t_char/2.0*sqrt(pi)*(1+Erf((_t_max-_t0)/_t_char)));
          node_data->PatE() += alpha*energy_density;
          total_energy += alpha*energy_density*fvm_node->volume();
          region_energy[ fvm_node->subdomain_id() ] += alpha*energy_density*fvm_node->volume();
        }
      }
    }
    else
    {
      // find the nearest node of the track
      FVM_Node * neraset_node=0;
      double distance=std::numeric_limits<double>::infinity();
      for(unsigned int r=0; r<_system.n_regions(); r++)
      {
        const SimulationRegion * region = _system.region(r);
        double dist;
        const Node * n = nn_locator->nearest_node(0.5*(track.start+track.end), r, dist);
        if(n == NULL) continue;

        FVM_Node * fvm_node = region->region_fvm_node(n); // may be NULL, if not on local
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
        double energy_density = track.energy/neraset_node->volume();
        FVM_NodeData * node_data = neraset_node->node_data();
        node_data->PatG() += energy_density/_quan_eff/(_t_char/2.0*sqrt(pi)*(1+Erf((_t_max-_t0)/_t_char)));
        node_data->PatE() += energy_density;
        total_energy += track.energy;
        region_energy[neraset_node->subdomain_id()] += track.energy;
      }
    }

  }

  MESSAGE<< "ok" <<std::endl;
  RECORD();

  STOP_LOG("update_system()", "Particle_Source_Track");

}

#endif

#if 0
void Particle_Source_Track::update_system()
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

