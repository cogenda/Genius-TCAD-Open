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




#ifndef __cv_hook_h__
#define __cv_hook_h__


#include "hook.h"
#include <time.h>

/**
 * write electrode IV into spice raw file (Ascii format).
 * then user can view the IV curve by some other program.
 * ( can we do real time display here? )
 */ 
class CVHook : public Hook 
{

public:
  CVHook(SolverBase & solver, const std::string & name, void * file);
  
  virtual ~CVHook();
  
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
  * the variable name buffer
  */
 std::string  _sweep_electrode;
 std::vector<std::string>  _gate_electrodes;
 
 /**
  * the variable value buffer
  */
 std::vector<double> _vsweep;
 std::vector< std::vector<double> > _gate_charge;
 
 /**
  * the total number of values 
  */
 unsigned int _n_values;
 
};

#endif 
