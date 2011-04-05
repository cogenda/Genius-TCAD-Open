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


#ifndef __threshold_hook_h__
#define __threshold_hook_h__


#include "hook.h"
#include <ctime>
#include <vector>
#include <string>

/**
 * export specified variable at every step, stop the simulation when the value exceed given threshold
 */
class ThresholdHook : public Hook
{

public:
  ThresholdHook(SolverBase & solver, const std::string & name, void *);

  virtual ~ThresholdHook();

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
   * monitor which region
   */
  std::string _region;

  /**
   * monitor which box
   */
  Point _lower_bound;

  /**
   * monitor which box
   */
  Point _upper_bound;

  bool _is_bound_box_valid();

  bool _in_bound_box(const Point &p);

  /**
   * the field scalar variable to be monitor and their threshold
   */
  std::map<SolutionVariable, double>     _scalar_variable_threshold_map;

  /**
   * the field vector variable to be monitor and their threshold
   */
  std::map<SolutionVariable, double>     _vector_variable_threshold_map;

  void _check_T_threshold();

  void _check_E_threshold();

  bool _violate_threshold;

  bool _stop_when_violate_threshold;


  /**
   * the output file name
   */
  std::string     _threshold_prefix;

};

#endif
