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

//  $Id: rawfile_hook.h,v 1.5 2008/07/09 05:58:16 gdiso Exp $


#ifndef __rawfile_hook_h__
#define __rawfile_hook_h__


#include "hook.h"
#include <time.h>

/**
 * write electrode IV into spice raw file (Ascii format).
 * then user can view the IV curve by some other program.
 * ( can we do real time display here? )
 */
class RawFileHook : public Hook
{

public:
  RawFileHook(SolverBase & solver, const std::string & name, void * rawfile);

  virtual ~RawFileHook();

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
  * current time
  */
 time_t          _time;

 /**
  * the input file name
  */
 std::string     _input_file;

 /**
  * the raw file name
  */
 std::string     _raw_file;

 /**
  * file stream
  */
 std::ofstream   _out;

 /**
  * if we are in mixA mode
  */
 bool            _mixA;

 /**
  * the variable name buffer
  */
 std::vector<std::pair<std::string, std::string> >  _variables;

 /**
  * the variable value buffer
  */
 std::vector< std::vector<double> > _values;

 /**
  * the total number of values
  */
 unsigned int _n_values;

};

#endif
