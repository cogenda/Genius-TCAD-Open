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


#ifndef __spice_monitor_hook_h__
#define __spice_monitor_hook_h__


#include "hook.h"


/**
 * monitor spice variables during each nonlinear iteration
 * for check the convergence histroy
 */
class SpiceMonitorHook : public Hook
{

public:

  SpiceMonitorHook(SolverBase & solver, const std::string & name, void *);

  virtual ~SpiceMonitorHook();

  /**
   *   This is executed before the initialization of the solver
   */
  virtual void on_init();

  /**
   *   This is executed previously to each solution step.
   */
  virtual void pre_solve();

  /**
   *  This is executed after each solution step.
   */
  virtual void post_solve();

  /**
   *  This is executed before each (nonlinear) iteration
   */
  virtual void pre_iteration();

  /**
   *  This is executed after each (nonlinear) iteration
   */
  virtual void post_check(void * f, void * x, void * dx, void * w, bool & change_y, bool &change_w);

  /**
   * This is executed after the finalization of the solver
   */
  virtual void on_close();

private:

  /**
   * monitor mix solver
   */
  bool            _mix_solver;

  /**
   * iteration count
   */
  unsigned int iteration_count;

  /**
     * solution count
   */
  unsigned int solution_count;

};

#endif
