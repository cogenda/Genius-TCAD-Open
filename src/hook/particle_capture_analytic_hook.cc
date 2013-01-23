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
#include "particle_capture_analytic_hook.h"
#include "parallel.h"



using PhysicalUnit::um;
using PhysicalUnit::mm;
using PhysicalUnit::cm;
using PhysicalUnit::J;
using PhysicalUnit::kg;
using PhysicalUnit::s;


ParticleCaptureAnalyticHook::ParticleCaptureAnalyticHook(SolverBase & solver, const std::string & name, void *param)
  : Hook ( solver, name ), _particle("e-"), _type("uniform"), _generation(0.0), _dose_rate(0.0), _char_length(1*mm)
{
  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for(std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
      parm_it != parm_list.end(); parm_it++)
  {
    if(parm_it->name() == "particle" && parm_it->type() == Parser::STRING)
      _particle = parm_it->get_string();
    if(parm_it->name() == "type" && parm_it->type() == Parser::STRING)
      _type = parm_it->get_string();
    if(parm_it->name() == "generation" && parm_it->type() == Parser::REAL)
      _generation = parm_it->get_real()*pow(cm, -3)/s;
    if(parm_it->name() == "doserate" && parm_it->type() == Parser::REAL)
      _dose_rate = parm_it->get_real()*J/kg/s;

    if(parm_it->name() == "x" && parm_it->type() == Parser::REAL)
      _center(0)=parm_it->get_real()*um;
    if(parm_it->name() == "y" && parm_it->type() == Parser::REAL)
      _center(1)=parm_it->get_real()*um;
    if(parm_it->name() == "z" && parm_it->type() == Parser::REAL)
      _center(2)=parm_it->get_real()*um;

    if(parm_it->name() == "char.length" && parm_it->type() == Parser::REAL)
      _char_length=parm_it->get_real()*um;
  }

}


ParticleCaptureAnalyticHook::~ParticleCaptureAnalyticHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void ParticleCaptureAnalyticHook::on_init()
{

}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void ParticleCaptureAnalyticHook::pre_solve()
{
  double dt = SolverSpecify::dt;
  if(!SolverSpecify::TimeDependent)
    dt = 1*PhysicalUnit::s;

  std::map<SimulationRegion *,  BoundaryCondition *> floating_metal_region;
  SimulationSystem & system = _solver.get_system();
  for(unsigned int n=0; n<system.n_regions(); n++)
  {
    SimulationRegion * region = system.region(n);
    BoundaryCondition * charge_integral_bc = region->floating_metal();
    if(charge_integral_bc)
      floating_metal_region.insert( std::make_pair(region, charge_integral_bc) );
  }

  for(unsigned int i=0; i<system.n_regions(); i++)
  {
    SimulationRegion * region = system.region(i);

    if(region->type() == InsulatorRegion && region->material() != "Air")
    {
      SimulationRegion::local_node_iterator it = region->on_local_nodes_begin();
      SimulationRegion::local_node_iterator it_end = region->on_local_nodes_end();
      for(; it!=it_end; ++it)
      {
        FVM_Node * fvm_node = *it;
        fvm_node->node_data()->Field_G()  = gen(fvm_node->root_node());
        fvm_node->node_data()->DoseRate() = dose(fvm_node->root_node());
      }
    }


    if( floating_metal_region.find(region) != floating_metal_region.end() )
    {
      BoundaryCondition * charge_integral_bc = floating_metal_region.find(region)->second;
      SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        FVM_Node * fvm_node = *it;
        charge_integral_bc->scalar("qf") += -gen(fvm_node->root_node())*dt;
      }
    }
  }

}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void ParticleCaptureAnalyticHook::post_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void ParticleCaptureAnalyticHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void ParticleCaptureAnalyticHook::on_close()
{

}



double ParticleCaptureAnalyticHook::gen(const Point *p)
{
  if( _type == "uniform" )
    return _generation;
  if( _type == "gaussian")
  {
    double d = (*p-_center).size();
    return _generation*exp(-d*d/_char_length/_char_length);
  }
  return 0.0;
}


double ParticleCaptureAnalyticHook::dose(const Point *p)
{
  if( _type == "uniform" )
    return _dose_rate;
  if( _type == "gaussian")
  {
    double d = (*p-_center).size();
    return _dose_rate*exp(-d*d/_char_length/_char_length);
  }
  return 0.0;
}

#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new ParticleCaptureAnalyticHook ( solver, name, fun_data );
  }

}

#endif


