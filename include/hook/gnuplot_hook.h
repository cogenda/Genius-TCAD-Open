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

//  $Id: gnuplot_hook.h,v 1.5 2008/07/09 12:25:19 gdiso Exp $


#ifndef __gnuplot_hook_h__
#define __gnuplot_hook_h__


#include "hook.h"


/**
 * write electrode IV into file which can be plotted by gnuplot
 */
class GnuplotHook : public Hook
{

public:
  GnuplotHook(SolverBase & solver, const std::string & name, void * file);

  virtual ~GnuplotHook();

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


 /**
  * the input file name
  */
 std::string     _input_file;

 /**
  * the output file name
  */
 std::string     _gnuplot_file;

 /**
  * file stream
  */
 std::ofstream   _out;

 /**
  * write the head of file
  */
 void  _write_gnuplot_head();

 /**
  * if we are in ddm solver
  */
 bool            _ddm;

 /**
  * if we are in mixA mode
  */
 bool            _mixA;
};

#endif
