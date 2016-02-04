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
#include "particle_monitor_hook.h"
#include "parallel.h"


ParticleMonitorHook::ParticleMonitorHook(SolverBase & solver, const std::string & name, void *)
  :Hook(solver, name)
{}


ParticleMonitorHook::~ParticleMonitorHook()
{}



/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void ParticleMonitorHook::on_init()
{
  _carrier_number = 0.0;
}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void ParticleMonitorHook::pre_solve()
{

}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void ParticleMonitorHook::post_solve()
{
  const SimulationSystem & system = this->get_solver().get_system();
  double z_width = system.z_width();

  for(unsigned int n=0; n<system.n_regions(); n++)
  {
    const SimulationRegion * region = system.region(n);

    if( region->type()== SemiconductorRegion)
    {
      SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        const FVM_Node * fvm_node = (*it);
        const FVM_NodeData * fvm_node_data = fvm_node->node_data();

        _carrier_number += fvm_node_data->Field_G()*SolverSpecify::dt*fvm_node->volume()*z_width;
      }
    }
  }

  double carrier_number = _carrier_number;
  Parallel::sum(carrier_number);

  if ( !Genius::processor_id() )
    std::cout<<"Charge Generated: " << carrier_number << " e-h pairs, " << PhysicalUnit::e*carrier_number/(1e-15*PhysicalUnit::C) << " fC"<< std::endl;
}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void ParticleMonitorHook::post_iteration()
{

}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void ParticleMonitorHook::on_close()
{

}



#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook (SolverBase & solver, const std::string & name, void * fun_data)
  {
    return new ParticleMonitorHook(solver, name, fun_data );
  }
}

#endif

