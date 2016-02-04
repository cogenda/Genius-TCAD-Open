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

#ifndef __data_hook_h__
#define __data_hook_h__


#include "hook.h"


/**
 *
 */
class DataHook : public Hook
{

public:
  DataHook(SolverBase & solver, const std::string & name, void * datafile);

  virtual ~DataHook();

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
   * flag when we do 2D simulation
   */
  unsigned int _dim;

  /**
   * count
   */
  unsigned int    _count;

  /**
   * the output file name
   */
  std::string     _output_prefix;

  /**
   * file stream
   */
  std::ofstream   _out;

  /**
   * the variable to record
   */
  std::vector<std::string>  _variable_name;

  /**
   * which region to record
   */
  std::vector<std::string> _region;
};

#endif
