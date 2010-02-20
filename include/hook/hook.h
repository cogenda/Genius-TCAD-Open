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

//  $Id: hook.h,v 1.3 2008/07/09 05:58:16 gdiso Exp $


#ifndef __hook_h__
#define __hook_h__

#include <string>


// predefine
class SolverBase;

/**
 * Hooks are executed previous and after each solution
 */
class Hook
{
public:

  /**
   *   Initializes the hook. this is done at the end of each solver's create() call
   *   @param SimulationSystem (input) the system
   *   @param name (input) the name of this hook
   */
  Hook(SolverBase & solver, const std::string & name)
      :_solver(solver),  _name(name)
  {}

  /**
   * virtual destroctor, do nothing here
   */
  virtual ~Hook() {}

  /**
   *   This is executed before the initialization of the solver
   */
  virtual void on_init() {}

  /**
   *   This is executed previously to each solution step.
   */
  virtual void pre_solve() {}

  /**
   *  This is executed after each solution step.
   */
  virtual void post_solve() {}

  /**
   *  This is executed after each (nonlinear) iteration
   *  i.e. for collecting convergence information or implementing various damping strategy
   */
  virtual void post_iteration() {}

  /**
   * This is executed after the finalization of the solver
   */
  virtual void on_close() {}

  /**
   * @return the name of the hook
   */
  const std::string & name() const
  { return _name; }

  /**
   * @return const reference of solver
   * NOTE although we allow hook change the solver internal data
   * however it is highly not recommend
   */
  const SolverBase  & get_solver() const
  { return _solver; }

protected:

  /**
   * reference to base class of solver, user can use dynamic_cast to convert
   * it to derived solver.
   */
  SolverBase    &    _solver;

  /**
   * the name of this hook
   */
  std::string        _name;

};


#endif //define __hook_h__
