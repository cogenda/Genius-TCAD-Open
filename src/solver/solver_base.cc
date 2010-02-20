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

SolverBase::SolverBase(SimulationSystem & system)
  :_system(system), _dom_solution_root(NULL), _dom_curr_solution(NULL)
{
}

SolverBase::~SolverBase()
{
}

int SolverBase::create_solver()
{
  // call hook function on_init
  hook_list()->on_init();
  return 0;
}

int SolverBase::destroy_solver()
{
  // call (user defined) hook functions on_close
  // and then hook_list will delete all the hooks!
  hook_list()->on_close();

  return 0;
}



int SolverBase::pre_solve_process(bool /*load_solution*/)
{
  // call (user defined) hook function hook_pre_solve_process
  hook_list()->pre_solve();

  return 0;
}

int SolverBase::post_solve_process()
{
  // call (user defined) hook function hook_post_solve_process
  hook_list()->post_solve();

  return 0;
}


void SolverBase::set_solution_dom_root(mxml_node_t* root)
{
  _dom_solution_root = root;
}

mxml_node_t* SolverBase::new_dom_solution_elem() const
{
  mxml_node_t *eSolution = mxmlNewElement(_dom_solution_root, "solution");
  _dom_curr_solution = eSolution;

  mxmlNewElement(eSolution, "output");

  return eSolution;
}

mxml_node_t* SolverBase::current_dom_solution_elem() const
{
  if (!_dom_curr_solution)
    return new_dom_solution_elem();
  return _dom_curr_solution;
}

