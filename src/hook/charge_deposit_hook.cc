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

#include <string>
#include <cstdlib>
#include <ctime>
#include <iomanip>

#include "solver_base.h"
#include "spice_ckt.h"
#include "charge_deposit_hook.h"
#include "parallel.h"

using PhysicalUnit::e;


/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
ChargeDepositHook::ChargeDepositHook(SolverBase & solver, const std::string & name, void * param)
    : Hook(solver, name)
{
  _p_solver = & solver;
  _p_fvm_node   = 0;
  _min_loc = invalid_uint;
  _charge = 0;
  _apply  = false;
  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for(std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
      parm_it != parm_list.end(); parm_it++)
  {
    if(parm_it->name() == "x" && parm_it->type() == Parser::REAL)
      _pp(0)=parm_it->get_real() * PhysicalUnit::um;
    if(parm_it->name() == "y" && parm_it->type() == Parser::REAL)
      _pp(1)=parm_it->get_real() * PhysicalUnit::um;
    if(parm_it->name() == "z" && parm_it->type() == Parser::REAL)
      _pp(2)=parm_it->get_real() * PhysicalUnit::um;
    if(parm_it->name() == "charge" && parm_it->type() == Parser::REAL)
      _charge=parm_it->get_real() * PhysicalUnit::e;
    if(parm_it->name() == "region" && parm_it->type() == Parser::STRING)
      _region = parm_it->get_string();
  }

  SimulationSystem & system = _p_solver->get_system();
  if(!_region.empty() && !system.has_region(_region)) _region.clear();
}


/*----------------------------------------------------------------------
 * destructor, close file
 */
ChargeDepositHook::~ChargeDepositHook()
{

}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void ChargeDepositHook::on_init()
{

  double min_dis = 1e100;
  for( unsigned int r=0; r<_p_solver->get_system().n_regions(); r++)
  {
    SimulationRegion * region = _p_solver->get_system().region(r);
    if(!_region.empty() && region->name() != _region) continue;

    SimulationRegion::processor_node_iterator node_it = region->on_processor_nodes_begin();
    SimulationRegion::processor_node_iterator node_it_end = region->on_processor_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      FVM_Node * fvm_node = *node_it;
      const Node * node = fvm_node->root_node();

      double dis = ((*node)-_pp).size();
      if(dis<min_dis)
      {
        min_dis = dis;
        _p_fvm_node = fvm_node;
      }
    }
  }

  // after this call, the _min_loc contains processor_id with minimal min_dis
  Parallel::min_loc(min_dis, _min_loc);

  if (Genius::processor_id() == _min_loc)
  {
    SimulationRegion * region = _p_solver->get_system().region(_p_fvm_node->subdomain_id());
    assert(region->type() == InsulatorRegion);
  }
}

/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void ChargeDepositHook::pre_solve()
{
  if(_apply) return;

  if (Genius::processor_id() == _min_loc)
  {
    const Node * node = _p_fvm_node->root_node();
    FVM_NodeData * fvm_node_data = _p_fvm_node->node_data();
    if(_charge >0)
      fvm_node_data->p() += _charge/_p_fvm_node->volume();
    else
      fvm_node_data->n() += -_charge/_p_fvm_node->volume();
  }

  _apply = true;
}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void ChargeDepositHook::post_solve()
{

}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void ChargeDepositHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void ChargeDepositHook::on_close()
{}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook (SolverBase & solver, const std::string & name, void * fun_data)
  {
    return new ChargeDepositHook(solver, name, fun_data );
  }

}

#endif

