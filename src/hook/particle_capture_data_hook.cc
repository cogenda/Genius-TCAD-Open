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
#include "particle_capture_data_hook.h"
#include "dose_rate.h"
#include "parallel.h"


using PhysicalUnit::um;
using PhysicalUnit::mm;
using PhysicalUnit::cm;
using PhysicalUnit::eV;
using PhysicalUnit::s;

ParticleCaptureDataHook::ParticleCaptureDataHook(SolverBase & solver, const std::string & name, void *param)
    : Hook ( solver, name ), _particle("e-"), _fraction(1.0), _weight(1.0), _resolution(1*mm)
{
  double _particle_num_per_run = 1;
  double _particle_flux = 1e7/(cm*cm)/s;
  double _source_area = 1*cm*cm;
  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for(std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
      parm_it != parm_list.end(); parm_it++)
  {
    if(parm_it->name() == "track.data" && parm_it->type() == Parser::STRING)
      _track_data_file = parm_it->get_string();

    if(parm_it->name() == "particle.data" && parm_it->type() == Parser::STRING)
      _particle_data_file = parm_it->get_string();

    if(parm_it->name() == "resolution" && parm_it->type() == Parser::REAL)
      _resolution = parm_it->get_real()*mm;

    if(parm_it->name() == "particle.run" && parm_it->type() == Parser::REAL)
      _particle_num_per_run = parm_it->get_real();

    if(parm_it->name() == "fraction" && parm_it->type() == Parser::REAL)
      _fraction = parm_it->get_real();

    if(parm_it->name() == "particle.flux" && parm_it->type() == Parser::REAL)
      _particle_flux = parm_it->get_real()/(cm*cm)/s;

    if(parm_it->name() == "source.area" && parm_it->type() == Parser::REAL)
      _source_area = parm_it->get_real()*(cm*cm);
  }

  _weight = _fraction*_particle_flux*_source_area/_particle_num_per_run; // weight in the unit of #/s
  _first_g4_calculation_flag = false;
}


ParticleCaptureDataHook::~ParticleCaptureDataHook()
{
}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void ParticleCaptureDataHook::on_init()
{}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void ParticleCaptureDataHook::pre_solve()
{
  if(_first_g4_calculation_flag) return;

  build_particles();
  process_particle_gen();
  elem_deposite.clear();


  dose_rate = new DoseRate(_solver.get_system());
  dose_rate->set_min_distance(_resolution);
  build_tracks();
  process_particle_dose_octree();
  delete dose_rate;
  dose_rate = 0;

  if(!_first_g4_calculation_flag)
    _first_g4_calculation_flag = true;

}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void ParticleCaptureDataHook::post_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void ParticleCaptureDataHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void ParticleCaptureDataHook::on_close()
{}



void ParticleCaptureDataHook::process_particle_gen()
{
  double dt = SolverSpecify::dt;
  if(!SolverSpecify::TimeDependent)
    dt = 1*PhysicalUnit::s;

  double charge_at_this_step = _weight*dt;

  SimulationSystem & system = _solver.get_system();
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);
    if(region->type() != InsulatorRegion) continue;

    SimulationRegion::element_iterator eit = region->elements_begin();
    for(; eit != region->elements_end(); ++eit)
    {
      const Elem * elem = *eit;
      double deposite = elem_deposite[elem->id()];

      double volume = 1e-10;
      for(unsigned int nd=0; nd<elem->n_nodes(); nd++)
        volume += elem->partial_volume_truncated(nd);

      for(unsigned int nd=0; nd<elem->n_nodes(); nd++)
      {
        FVM_Node * fvm_node = elem->get_fvm_node(nd);
        if(!fvm_node->on_local()) continue;
        assert(fvm_node->node_data());

        double vol_ratio = elem->partial_volume_truncated(nd)/volume;
        fvm_node->node_data()->Field_G() += deposite*charge_at_this_step*vol_ratio/fvm_node->volume()/dt;
      }
    }
  }

}



void ParticleCaptureDataHook::build_particles()
{
  if(Genius::processor_id()==0)
  {
    std::ifstream in(_particle_data_file.c_str());
    if(!in.good())
    {
      MESSAGE<<"ERROR PARTICLE: file "<<_particle_data_file<<" can't be opened."<<std::endl; RECORD();
      genius_error();
    }

    const SimulationSystem & system = _solver.get_system();
    const PointLocatorBase & point_locator = system.mesh().point_locator();

    elem_deposite.resize(system.mesh().n_elem(), 0.0);
    while(!in.eof())
    {
      std::string context;
      std::getline(in, context);
      if(context.empty()) continue;
      if(context[0] == '#') continue;

      std::stringstream ss(context);
      double x, y, z;
      ss >> x >> y >> z;

      Point p(x*um, y*um, z*um);

      const Elem * elem = point_locator(p);
      if(elem) elem_deposite[elem->id()] += 1.0;
    };
    in.close();
  }

  Parallel::broadcast(elem_deposite);
}




void ParticleCaptureDataHook::build_tracks()
{
  if(Genius::processor_id()==0)
  {
    std::ifstream in(_track_data_file.c_str());
    if(!in.good())
    {
      MESSAGE<<"ERROR PARTICLE: file "<<_track_data_file<<" can't be opened."<<std::endl; RECORD();
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

      double p1[3], p2[3], energy;
      ss >>  p1[0] >> p1[1] >> p1[2] >> p2[0] >> p2[1] >> p2[2] >> energy;

      // skip track with very little energy
      // a very low energy such as 1.38073e-315 may break some floating point compare result
      if(energy < 1e-12) continue;
      // skip very short track
      Point begin(p1[0]*um,p1[1]*um,p1[2]*um);
      Point end(p2[0]*um,p2[1]*um,p2[2]*um);
      if( (begin-end).size() < 1e-6*um ) continue;

      dose_rate->energy_deposite(begin, end, _weight*energy*1e6*eV);
    }

    in.close();
  }

}



void ParticleCaptureDataHook::process_particle_dose_octree()
{
  SimulationSystem & system = _solver.get_system();
  const double T = system.T_external();

  // set dose rate to each insulator node
  dose_rate->sync_energy_deposite();
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);
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

  //dose_rate->export_vtk("bb.vtk");

}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new ParticleCaptureDataHook ( solver, name, fun_data );
  }

}

#endif


