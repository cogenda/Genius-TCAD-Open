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


#ifndef __injection_hook_h__
#define __injection_hook_h__


#include "hook.h"


class InjectionHook : public Hook
{

public:

  InjectionHook(SolverBase & solver, const std::string & name, void *);

  virtual ~InjectionHook();

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

  BoundaryCondition * _injection_bc;

  struct InjectionSurface
  {
    FVM_Node * fvm_node;
    unsigned int surface_patches;
    std::vector<Point> surface_norm;
    std::vector<Real>  surface_area;
  };

  std::vector<InjectionSurface> _injection_surfaces;


  void *dll_file;

  typedef double (*external_funxtion)(double x, double y, double z, double nx, double ny, double nz, double t);

  /**
   * calculate the current density of electron injected to semiconductor surface
   * given a surface point(x,y,z) with unit um, norm(nx, ny, nz) and time with unit s
   * return the electron current density in the unit of A/cm^2
   * @param x   the x location of surface point, unit um
   * @param y   the y location of surface point, unit um
   * @param z   the z location of surface point, unit um
   * @param nx   the norm vector
   * @param ny   the norm vector
   * @param nz   the norm vector
   * @param t   the current time, unit second
   * @return the electron current density. A per square cm
   */
  double (*electron_inject_density)(double x, double y, double z, double nx, double ny, double nz, double t);

  /**
   * calculate the current density of hole injected to semiconductor surface
   * given a surface point(x,y,z) with unit um, norm(nx, ny, nz) and time with unit s
   * return the hole current density in the unit of A/cm^2
   * @param x   the x location of surface point, unit um
   * @param y   the y location of surface point, unit um
   * @param z   the z location of surface point, unit um
   * @param nx   the norm vector
   * @param ny   the norm vector
   * @param nz   the norm vector
   * @param t   the current time, unit second
   * @return the electron current density. A per square cm
   */
  double (*hole_inject_density)(double x, double y, double z, double nx, double ny, double nz, double t);

};

#endif
