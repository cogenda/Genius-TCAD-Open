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


#ifndef __probe_hook_h__
#define __probe_hook_h__


#include "hook.h"
#include <ctime>

/**
 * write electrode IV into file which can be plotted by probe
 */
class ProbeHook : public Hook
{

public:
  ProbeHook(SolverBase & solver, const std::string & name, void * file);

  virtual ~ProbeHook();

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


 SolverBase*     _p_solver;

 /**
  * probe point
  */
 Point _pp;

 const FVM_Node* _p_fvm_node;

 unsigned int    _min_loc;

 /**
  * the output file name
  */
 std::string     _probe_file;

 /**
  * file stream
  */
 std::ofstream   _out;


};

#endif
