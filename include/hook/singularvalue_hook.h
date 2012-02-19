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


#ifndef __singularvalue_hook_h__
#define __singularvalue_hook_h__

#include <string>

#include "config.h"
#include "hook.h"



/**
 * calculate singular value for jacobian matrix on each nonlinear iteration
 * for check the condition number of device
 */
class SingularValueHook : public Hook
{

public:

  SingularValueHook(SolverBase & solver, const std::string & name, void *);

  virtual ~SingularValueHook();

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

private:


};

#endif
