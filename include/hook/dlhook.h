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

//  $Id: dlhook.h,v 1.4 2008/07/09 05:58:15 gdiso Exp $

#ifndef __dll_hook_h__
#define __dll_hook_h__

#include "config.h"

#ifdef DLLHOOK

#include "hook.h"

class DllHook : public Hook
{
public:

  DllHook(SolverBase & solver, const std::string & name, void *fun_data);

  ~DllHook();

  typedef Hook * GET_HOOK(SolverBase & solver, const std::string & name, void * fun_data);

private:

  void             * dll_handle;

  GET_HOOK         * get_hook;

  Hook             * hook;


public:
  /**
   *   This is executed before the initialization of the solver
   */
  virtual void on_init()
  { if(hook) hook->on_init();        }

  /**
   *   This is executed previously to each solution step.
   */
  virtual void pre_solve()
  { if(hook) hook->pre_solve();      }

  /**
   *  This is executed after each solution step.
   */
  virtual void post_solve()
  { if(hook) hook->post_solve();     }


  /**
   *  This is executed before each (nonlinear) iteration
   *  i.e. for analysis the condition number of jacobian matrix
   */
  virtual void pre_iteration()
  { if(hook) hook->pre_iteration(); }

  /**
   *  This is executed after each (nonlinear) iteration
   *  i.e. for collecting convergence information or implementing various damping strategy
   */
  virtual void post_iteration()
  { if(hook) hook->post_iteration(); }

  /**
   *  This is executed after each (nonlinear) iteration
   *  i.e. for collecting convergence information or implementing various damping strategy
   */
  virtual void post_iteration(void * f, void * x, void * y, void * w, bool & change_y, bool &change_w)
  { if(hook) hook->post_iteration(f, x, y, w, change_y, change_w); }

  /**
   * This is executed after the finalization of the solver
   */
  virtual void on_close()
  { if(hook) hook->on_close();       }

};

#endif

#endif
