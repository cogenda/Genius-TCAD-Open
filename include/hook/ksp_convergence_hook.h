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


#ifndef __ksp_convergence_hook_h__
#define __ksp_convergence_hook_h__


#include "hook.h"
#include <ctime>
#include <vector>
#include <string>

/**
 * report function residual
 */
class KSPConvergenceHook : public Hook
{

public:
  KSPConvergenceHook(SolverBase & solver, const std::string & name, void *);

  virtual ~KSPConvergenceHook();

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
   *  This is executed after each (nonlinear) iteration
   */
  virtual void post_iteration();

  /**
   * This is executed after the finalization of the solver
   */
  virtual void on_close();

private:

  /**
   * monitor poisson solver
   */
  bool            _poisson_solver;

  /**
   * monitor ddm solver
   */
  bool            _ddm_solver;


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
