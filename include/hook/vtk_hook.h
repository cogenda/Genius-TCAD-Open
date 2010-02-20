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


#ifndef __vtk_hook_h__
#define __vtk_hook_h__


#include "hook.h"
#include <ctime>
#include <vector>
#include <string>

/**
 * write vtk file
 */
class VTKHook : public Hook
{

public:
  VTKHook(SolverBase & solver, const std::string & name, void *);

  virtual ~VTKHook();

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
   * the output file name
   */
  std::string     _vtk_prefix;

  /**
  * count
  */
  unsigned int count;

  /**
   * last value
   */
  double _t_last;
  double _v_last;
  double _i_last;

  /**
   * step
   */
  double _t_step;
  double _v_step;
  double _i_step;

  /**
   * the time(voltage, current) sequence for paraview animation
   */
  std::vector< std::pair<double, std::string> > time_sequence;

  /**
   * if we are in mixA mode
   */
  bool            _mixA;
};

#endif
