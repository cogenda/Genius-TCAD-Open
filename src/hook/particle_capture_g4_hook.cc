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

#include "solver_base.h"
#include "mesh_base.h"
#include "simulation_system.h"
#include "point_locator_base.h"
#include "particle_capture_g4_hook.h"
#include "nearest_node_locator.h"
#include "dose_rate.h"
#include "parallel.h"
#include "mathfunc.h"

#ifdef HAVE_AMS

using PhysicalUnit::um;
using PhysicalUnit::mm;
using PhysicalUnit::cm;
using PhysicalUnit::eV;
using PhysicalUnit::J;
using PhysicalUnit::kg;
using PhysicalUnit::s;


ParticleCaptureG4Hook::ParticleCaptureG4Hook(SolverBase & solver, const std::string & name, void *param)
  : Hook ( solver, name ), _particle("e-"), _weight(1.0), _particle_flux(1e7/(cm*cm)/s), _source_area(1*cm*cm), _particle_num_per_run(1e3),
   _resolution(1*mm), _lateral_char(10*um), _beam_num(10000),
   _const_flux(true), _dose_rate("full"), host("localhost")
{
  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for(std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
      parm_it != parm_list.end(); parm_it++)
  {
    if(parm_it->name() == "particle" && parm_it->type() == Parser::STRING)
      _particle = parm_it->get_string();
    if(parm_it->name() == "host" && parm_it->type() == Parser::STRING)
      host = parm_it->get_string();
    if(parm_it->name() == "particle.run" && parm_it->type() == Parser::REAL)
      _particle_num_per_run = parm_it->get_real();
    if(parm_it->name() == "particle.flux" && parm_it->type() == Parser::REAL)
      _particle_flux = parm_it->get_real()/(cm*cm)/s;
    if(parm_it->name() == "source.area" && parm_it->type() == Parser::REAL)
      _source_area = parm_it->get_real()*(cm*cm);
    if(parm_it->name() == "resolution" && parm_it->type() == Parser::REAL)
      _resolution = parm_it->get_real()*mm;
    if(parm_it->name() == "lateral.char" && parm_it->type() == Parser::REAL)
      _lateral_char = parm_it->get_real()*um;
    if(parm_it->name() == "beam" && parm_it->type() == Parser::INTEGER)
      _beam_num = parm_it->get_int();
    if(parm_it->name() == "const" && parm_it->type() == Parser::BOOL)
      _const_flux = parm_it->get_bool();
    if(parm_it->name() == "doserate" && parm_it->type() == Parser::STRING)
      _dose_rate = parm_it->get_string();
  }

  _first_g4_calculation_flag = false;

  dose_rate = new DoseRate(solver.get_system());
  dose_rate->set_min_distance(_resolution);
}


ParticleCaptureG4Hook::~ParticleCaptureG4Hook()
{
  delete dose_rate;
}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void ParticleCaptureG4Hook::on_init()
{
  if ( Genius::is_first_processor() )
  {
    char **comm_list;

    err = AMS_Connect(host.c_str(), -1, &comm_list);
    AMS_Check_error(err, &msg);


    /* Attach to a Communicator */
    err = AMS_Comm_attach("control", &control);
    AMS_Check_error(err, &msg);


    err = AMS_Comm_attach("ParticleEvent", &data);
    AMS_Check_error(err, &msg);
  }
}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void ParticleCaptureG4Hook::pre_solve()
{
  bool sim = !_first_g4_calculation_flag || (_first_g4_calculation_flag && !_const_flux);
  if(!sim ) return;

  double dt = SolverSpecify::dt;
  if(!SolverSpecify::TimeDependent)
    dt = 1*PhysicalUnit::s;

  SimulationSystem & system = _solver.get_system();
  const double T = system.T_external();


  // clear old dose rate information
  dose_rate->clear_energy_deposite();

  // reset particle gen Field_G/DoseRate
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);
    if(region->type() == InsulatorRegion)
    {
      SimulationRegion::local_node_iterator it = region->on_local_nodes_begin();
      SimulationRegion::local_node_iterator it_end = region->on_local_nodes_end();
      for(; it!=it_end; ++it)
      {
        FVM_Node * fvm_node = *it;
        fvm_node->node_data()->Field_G() = 0.0;
        fvm_node->node_data()->DoseRate() = 0.0;
      }
    }
  }

  int total_particle_num = std::max(1, static_cast<int>(_particle_num_per_run));
  _weight = _particle_flux*_source_area/total_particle_num; // weight in the unit of #/s
  do
  {
    int beam = total_particle_num > _beam_num ? _beam_num : total_particle_num;

    if ( Genius::is_first_processor() )
    {
      std::stringstream ss;
      ss << "/run/beamOn " << beam;
      send_command(ss.str());
      load_events();
    }

    Parallel::broadcast(EventNumber);
    Parallel::broadcast(ParticleEnergy);
    Parallel::broadcast(ParticleBeginPosition);
    Parallel::broadcast(ParticleEndPosition);
    Parallel::broadcast(EventSteps);
    Parallel::broadcast(ParticleStep);
    Parallel::broadcast(StepEnergy);

    process_particle_gen();
    process_particle_dose_octree();

    total_particle_num -= beam;

  }while(total_particle_num>0);

  // set dose rate to each insulator node
  dose_rate->sync_energy_deposite();


  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);
    //if(region->type() == InsulatorRegion)
    {
      SimulationRegion::local_node_iterator it = region->on_local_nodes_begin();
      SimulationRegion::local_node_iterator it_end = region->on_local_nodes_end();
      for(; it!=it_end; ++it)
      {
        FVM_Node * fvm_node = *it;
        const std::vector< std::pair<const Elem *, unsigned int> > & elem_has_this_node = fvm_node->elem_has_this_node();
        for(unsigned int n=0; n<elem_has_this_node.size(); ++n)
        {
          const Elem * elem = elem_has_this_node[n].first;
          fvm_node->node_data()->DoseRate() += dose_rate->energy_deposite_density(elem->centroid())/elem_has_this_node.size();
        }
      }
    }
  }




  if(!_first_g4_calculation_flag)
    _first_g4_calculation_flag = true;

  //dose_rate->export_vtk("cc.vtk");

}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void ParticleCaptureG4Hook::post_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void ParticleCaptureG4Hook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void ParticleCaptureG4Hook::on_close()
{
  if ( Genius::is_first_processor() )
  {
    send_command("exit");


    err = AMS_Disconnect();
    AMS_Check_error(err, &msg);
  }
}


//----------------------------------------------------------------------

void ParticleCaptureG4Hook::send_command(const std::string &cmd)
{
  AMS_Memory memory;
  AMS_Memory_type mtype; /* Memory type */
  AMS_Data_type dtype; /* Data type */
  AMS_Shared_type stype; /* Shared type */
  AMS_Reduction_type rtype; /* Reduction type */
  unsigned int step; /* Memory step number */

  char **mem_list, **fld_list;

  void *addr; /* Pointer to the data */
  int len; /* Data length */

  /* Attach the memory structure */
  err = AMS_Memory_attach(control, "control_buf", &memory, &step); assert(!err);
  AMS_Check_error(err, &msg);

  err = AMS_Memory_get_field_info(memory, "command", &addr, &len, &dtype, &mtype, &stype, &rtype); assert(!err);
  AMS_Check_error(err, &msg);

  char *cmd_buf = (char *)addr;
  memset(cmd_buf, 0, sizeof(char)*1024);
  strcpy( cmd_buf, cmd.c_str());

  AMS_Memory_update_send_begin(memory);
  AMS_Memory_update_send_end(memory);

  /* Receive an update */
  int changed;
  do
  {
    err = AMS_Memory_update_recv_end(memory, &changed, &step);
    AMS_Check_error(err, &msg);
  }while(!changed);

  err = AMS_Memory_detach(memory);
  AMS_Check_error(err, &msg);
}


void ParticleCaptureG4Hook::load_events()
{
  AMS_Memory memory;
  AMS_Memory_type mtype; /* Memory type */
  AMS_Data_type dtype; /* Data type */
  AMS_Shared_type stype; /* Shared type */
  AMS_Reduction_type rtype; /* Reduction type */
  unsigned int step; /* Memory step number */

  char **mem_list, **fld_list;

  void *addr; /* Pointer to the data */
  int len; /* Data length */

  /* Attach the memory structure */
  err = AMS_Memory_attach(control, "EventList", &memory, &step); assert(!err);
  AMS_Check_error(err, &msg);


  // get event number
  err = AMS_Memory_get_field_info(memory, "EventNumber", &addr, &len, &dtype, &mtype, &stype, &rtype); assert(!err);
  AMS_Check_error(err, &msg);
  assert(dtype == AMS_INT);
  EventNumber = *((int *)(addr));


  // get particle energy
  err = AMS_Memory_get_field_info(memory, "ParticleEnergy", &addr, &len, &dtype, &mtype, &stype, &rtype); assert(!err);
  AMS_Check_error(err, &msg);
  ParticleEnergy.clear();
  for(int i=0; i<len; i++)
  {
    ParticleEnergy.push_back(*((double *)(addr)+i));
  }

  // get particle start position
  err = AMS_Memory_get_field_info(memory, "ParticleBeginPosition", &addr, &len, &dtype, &mtype, &stype, &rtype); assert(!err);
  AMS_Check_error(err, &msg);
  ParticleBeginPosition.clear();
  for(int i=0; i<len; i++)
  {
    ParticleBeginPosition.push_back(*((double *)(addr)+i));
  }

  // get particle stop position
  err = AMS_Memory_get_field_info(memory, "ParticleEndPosition", &addr, &len, &dtype, &mtype, &stype, &rtype); assert(!err);
  AMS_Check_error(err, &msg);
  ParticleEndPosition.clear();
  for(int i=0; i<len; i++)
  {
    ParticleEndPosition.push_back(*((double *)(addr)+i));
  }

  // get step number for each event
  err = AMS_Memory_get_field_info(memory, "EventStep", &addr, &len, &dtype, &mtype, &stype, &rtype); assert(!err);
  AMS_Check_error(err, &msg);
  EventSteps.clear();
  for(int i=0; i<len; i++)
  {
    EventSteps.push_back(*((int *)(addr)+i));
  }

  // step end position
  err = AMS_Memory_get_field_info(memory, "StepEndPosition", &addr, &len, &dtype, &mtype, &stype, &rtype); assert(!err);
  AMS_Check_error(err, &msg);
  ParticleStep.clear();
  for(int i=0; i<len; i++)
  {
    ParticleStep.push_back(*((double *)(addr)+i));
  }

  // step energy
  err = AMS_Memory_get_field_info(memory, "StepEnergy", &addr, &len, &dtype, &mtype, &stype, &rtype); assert(!err);
  AMS_Check_error(err, &msg);
  StepEnergy.clear();
  for(int i=0; i<len; i++)
  {
    StepEnergy.push_back(*((double *)(addr)+i));
  }

  err = AMS_Memory_detach(memory);
  AMS_Check_error(err, &msg);
}



void ParticleCaptureG4Hook::process_particle_gen()
{

  double dt = SolverSpecify::dt;
  if(!SolverSpecify::TimeDependent)
    dt = 1*PhysicalUnit::s;

  double charge_at_this_step = _weight*dt;

  std::map<SimulationRegion *,  BoundaryCondition *> floating_metal_region;
  SimulationSystem & system = _solver.get_system();
  for(unsigned int n=0; n<system.n_regions(); n++)
  {
    SimulationRegion * region = system.region(n);
    BoundaryCondition * charge_integral_bc = region->floating_metal();
    if(charge_integral_bc)
      floating_metal_region.insert( std::make_pair(region, charge_integral_bc) );
  }

  std::vector<Point> particles;
  for(unsigned int i=0; i<ParticleEndPosition.size(); )
  {
    double x = ParticleEndPosition[i++]*um;
    double y = ParticleEndPosition[i++]*um;
    double z = ParticleEndPosition[i++]*um;
    particles.push_back(Point(x, y, z));
  }

  const PointLocatorBase & point_locator = system.mesh().point_locator();
  for(unsigned int i=0; i<particles.size(); i++)
  {
    const Elem * elem = point_locator(particles[i]);
    if(!elem || !elem->on_processor()) continue;

    SimulationRegion * region = system.region(elem->subdomain_id());

    //if(region->type() == InsulatorRegion)
    {
      double volume = 1e-10;
      for(unsigned int nd=0; nd<elem->n_nodes(); nd++)
      {
        volume += elem->partial_volume_truncated(nd);
      }
      for(unsigned int nd=0; nd<elem->n_nodes(); nd++)
      {
        FVM_Node * fvm_node = elem->get_fvm_node(nd);
        assert(fvm_node->node_data());

        double vol_ratio = elem->partial_volume_truncated(nd)/volume;
        fvm_node->node_data()->Field_G() += charge_at_this_step*vol_ratio/fvm_node->volume()/dt;
      }
    }


    if( floating_metal_region.find(region) != floating_metal_region.end() )
    {
      BoundaryCondition * charge_integral_bc = floating_metal_region.find(region)->second;
      charge_integral_bc->scalar("qf") += -charge_at_this_step;
    }
  }
}



void ParticleCaptureG4Hook::build_tracks()
{
  tracks.clear();

  if(_dose_rate == "simple")
  {
    // simple: for each particle, only consider begin to end
    for(int i=0; i<EventNumber; i++)
    {
      track_t track;
      track.start.x()    = ParticleBeginPosition[3*i+0]*um;
      track.start.y()    = ParticleBeginPosition[3*i+1]*um;
      track.start.z()    = ParticleBeginPosition[3*i+2]*um;
      track.end.x()      = ParticleEndPosition[3*i+0]*um;
      track.end.y()      = ParticleEndPosition[3*i+1]*um;
      track.end.z()      = ParticleEndPosition[3*i+2]*um;
      track.energy       = _weight*ParticleEnergy[i]*1e6*eV;  //weight
      track.lateral_char = _lateral_char;

      if(track.energy > 1e-6*eV && (track.start - track.end).size() > 1e-6*mm)
        tracks.push_back(track);
    }
  }
  if(_dose_rate == "full")
  {
    int step_counter = 0;
    for(int i=0; i<EventNumber; i++)
    {
      Point start;
      start.x()    = ParticleBeginPosition[3*i+0]*um;
      start.y()    = ParticleBeginPosition[3*i+1]*um;
      start.z()    = ParticleBeginPosition[3*i+2]*um;

      for(int s=0; s<EventSteps[i]; s++)
      {
        Point end;
        end.x()    = ParticleStep[3*step_counter+0]*um;
        end.y()    = ParticleStep[3*step_counter+1]*um;
        end.z()    = ParticleStep[3*step_counter+2]*um;

        track_t track;
        track.start  = start;
        track.end    = end;
        track.energy = _weight*StepEnergy[step_counter]*1e6*eV; //weight
        track.lateral_char = _lateral_char;

        if(track.energy > 1e-6*eV && (track.start - track.end).size() > 1e-6*mm)
          tracks.push_back(track);

        start = end;
        step_counter++;
      }
    }
  }

}



void ParticleCaptureG4Hook::process_particle_dose_nodal()
{
  SimulationSystem & system = _solver.get_system();
  const double T = system.T_external();

  build_tracks();
  if(tracks.empty()) return;

  double charge_at_this_step = _weight;
  double dt = SolverSpecify::dt;
  if(!SolverSpecify::TimeDependent)
    dt = 1*PhysicalUnit::s;

  const double pi = 3.1415926536;
  genius_assert(system.mesh().mesh_dimension() == 3);

  AutoPtr<NearestNodeLocator> nn_locator( new NearestNodeLocator(system.mesh()) );

  for(unsigned int t=0; t<tracks.size(); ++t)
  {
    const track_t & track = tracks[t];
    genius_assert(track.energy > 0.0 && (track.end - track.start).size() > 0.0);

    double track_energy = 0.0;
    std::map<FVM_Node *, double> track_energy_density;

    const Point track_dir = (track.end - track.start).unit(); // track direction
    const double ed = track.energy/(track.end - track.start).size(); // linear energy density
    const double lateral_char = track.lateral_char;
    // find the nodes that near the track
    for(unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);

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
        const SimulationRegion * region = system.region(fvm_node->subdomain_id());

        if(fvm_node && fvm_node->on_local())
        {
          FVM_NodeData * node_data = fvm_node->node_data();
          node_data->DoseRate() += alpha*energy_density/region->get_density(T)/dt;
        }
      }
    }
    else
    {
      // find the nearest node of the track
      FVM_Node * neraset_node=0;
      double distance=std::numeric_limits<double>::infinity();
      for(unsigned int r=0; r<system.n_regions(); r++)
      {
        const SimulationRegion * region = system.region(r);

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
        double energy_rate = track.energy/dt;
        const SimulationRegion * region = system.region(neraset_node->subdomain_id());

        FVM_NodeData * node_data = neraset_node->node_data();
        node_data->DoseRate() += energy_rate/region->get_density(T);
      }
    }
  }


  // sync doserate for all the ghost node
  for(unsigned int i=0; i<system.n_regions(); i++)
  {
    SimulationRegion * region = system.region(i);

    if(region->type() == InsulatorRegion)
      region->sync_point_variable<Real>("dose_rate");
  }
}



void ParticleCaptureG4Hook::process_particle_dose_octree()
{
  SimulationSystem & system = _solver.get_system();
  const double T = system.T_external();

  build_tracks();
  if(tracks.empty()) return;

  unsigned int bin_size = tracks.size()/Genius::n_processors();
  unsigned int begin = bin_size*Genius::processor_id();
  unsigned int end   = std::min(begin+bin_size, tracks.size());
  for(unsigned int t=begin; t<end; ++t)
  {
    const track_t & track = tracks[t];
    genius_assert(track.energy > 0.0 && (track.end - track.start).size() > 0.0);

    dose_rate->energy_deposite(track.start, track.end, track.energy);
  }

}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new ParticleCaptureG4Hook ( solver, name, fun_data );
  }

}

#endif


#endif


