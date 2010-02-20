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

//  $Id: shell_hook.h,v 1.2 2008/07/09 05:58:16 gdiso Exp $


#ifndef __shell_hook_h__
#define __shell_hook_h__


#include "hook.h"


/**
 * Executes a shell command.
 */
class ShellHook : public Hook
{

public:
  ShellHook(SolverBase & solver, const std::string & name, void * command)
  : Hook(solver, name), _command((const char *)command)
  {}

  virtual ~ShellHook()
  {}

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
   *  i.e. for collecting convergence information or implementing various damping strategy
   */
  virtual void post_iteration();

  /**
   * This is executed after the finalization of the solver
   */
  virtual void on_close();

private:
  /**
   * The command to be executed
   */
 std::string _command;

};

#endif
