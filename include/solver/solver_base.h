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

//  $Id: solver_base.h,v 1.15 2008/07/09 12:25:19 gdiso Exp $

#ifndef __solver_base_h__
#define __solver_base_h__

#include "genius_common.h"
#include "genius_env.h"
#include "simulation_system.h"
#include "boundary_condition_collector.h"
#include "solver_specify.h"
#include "hook_list.h"
#include "log.h"
#include "perf_log.h"

#include "mxml.h"

#ifndef CYGWIN
  #include "dlhook.h"
#endif


/**
 * the base class for all the solver.
 */
class SolverBase
{
public:

  /**
   * constructor, take one parameter as SimulationSystem
   */
  SolverBase(SimulationSystem & system);

  /**
   * virtual destructor
   */
  virtual ~SolverBase();

  /**
   * virtual functions, prepare solution and aux variables used by this solver
   */
  virtual int set_variables() { return 0; }

  /**
   * virtual function, create the solver
   */
  virtual int create_solver();

  /**
   * do pre-process before each solve action
   */
  virtual int pre_solve_process(bool load_solution=true);

  /**
   * do post-process after each solve action
   */
  virtual int post_solve_process();

  /**
   * virtual function, do the solve process
   */
  virtual int solve()=0;

  /**
   * virtual function, write solver intermediate data into system
   * It can be used to monitor the field data evolution during solve action
   */
  virtual void flush_system() {}

  /**
   * virtual function, destroy the solver
   */
  virtual int destroy_solver();

  /**
   * @return reference to system
   */
  SimulationSystem & get_system()
  { return _system; }

  /**
   * @return const reference to system
   */
  const SimulationSystem & get_system() const
  { return _system; }

  /**
   * @return the solver type
   */
  virtual SolverSpecify::SolverType solver_type()const=0;

  /**
   * @return the label of the solver
   */
  const std::string& label() const { return _label; }

  /**
   * set the label of the solver
   */
  void set_label(const std::string &label) { _label = label; }

  /**
   * add hook function to solver
   */
  void add_hook( Hook * hk )
  {
    _hooks.add_hook(hk);
  }

  /**
   * get the pointer to hook list
   */
  HookList * hook_list()
  { return & _hooks; }

  /**
   * set the root node of solution dom
   */
  void set_solution_dom_root(mxml_node_t* root);

  /**
   * return the current solution dom element
   */
  mxml_node_t* current_dom_solution_elem() const;

protected:
  /**
   * writable reference to SimulationSystem, solver
   * will read mesh/data from it and write result into it
   */
  SimulationSystem & _system;

  /**
   * the hook functions
   */
  HookList           _hooks;

  /**
   * create a solution dom element, and add to the dom document
   */
  mxml_node_t* new_dom_solution_elem() const;

private:
  /**
   * XML dom of the solution records
   */
  mxml_node_t* _dom_solution_root;

  mutable mxml_node_t* _dom_curr_solution;

  std::string _label;

};


#endif // #define __solver_base_h__
