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


#ifndef __particle_capture_1d_hook_h__
#define __particle_capture_1d_hook_h__



#include "hook.h"
#include "plane.h"
#include "monot_cubic_interpolator.h"


/**
 * set electron deposit and dose rate by DICTAT 1D 
 */
class ParticleCapture1DHook : public Hook
{

public:

  ParticleCapture1DHook(SolverBase & solver, const std::string & name, void *);

  virtual ~ParticleCapture1DHook();

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
   * G4 particle type, only support electron at present
   */
  std::string _particle;


private:
  
  Plane _base;
 
  MonotCubicInterpolator _gen;
  
  MonotCubicInterpolator _dose;
  
  void build(const std::string &);

};


#endif

