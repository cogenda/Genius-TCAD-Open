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

//  $Id: ddm_solver.h,v 1.4 2008/07/09 05:58:16 gdiso Exp $

#ifndef __ddm_solver_h__
#define __ddm_solver_h__

#include "fvm_nonlinear_solver.h"

/**
 * the common method for device drift-diffusion method solver
 */
class DDMSolverBase : public FVM_NonlinearSolver
{
public:
  /**
   * constructor
   */
  DDMSolverBase(SimulationSystem & system);

  /**
   * destructor
   */
  virtual ~DDMSolverBase()
  {}

  /**
   * do equilibrium state simulation
   */
  virtual int solve_equ();

  /**
   * compute steady-state
   */
  virtual int solve_steadystate();

  /**
   * do dcsweep, can be both voltage scan or current scan
   */
  virtual int solve_dcsweep();

  /**
   * do transient simulation
   */
  virtual int solve_transient();

  /**
   * IV curve automatically trace
   */
  virtual int solve_iv_trace();

  /**
   * virtual function, create the solver
   */
  virtual int create_solver();

  /**
   * do post-process after each solve action
   */
  virtual int post_solve_process();

  /**
   * virtual function, destroy the solver
   */
  virtual int destroy_solver();

  /**
   * load previous state into solution vector, empty here
   */
  virtual int diverged_recovery()=0;

  /**
   * compute the norm of local truncate error (LTE)
   * each derived DDM solver should override it.
   */
  virtual PetscScalar LTE_norm()=0;

  /**
   * force carrier density to be positive during projection
   */
  virtual void projection_positive_density_check(Vec , Vec )  {}

  /**
   * compute the abs and relative error norm of the solution
   * each derived DDM solver should override it.
   */
  virtual void error_norm()=0;

  /**
   * snes monitor, do nothing
   */
  virtual void petsc_snes_monitor(PetscInt , PetscReal ) {}

  /**
   * sens convergence criteria for dd solvers
   */
  virtual void petsc_snes_convergence_test(PetscInt its, PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason);

  /**
   * ksp convergence criteria for dd solvers
   */
  virtual void petsc_ksp_convergence_test(PetscInt its, PetscReal rnorm, KSPConvergedReason* reason);

protected:

  /**
   * the global privious solution vector at n step
   */
  Vec            x_n;

  /**
   * the global privious solution vector at n-1 step
   */
  Vec            x_n1;

  /**
   * the global privious solution vector at n-2 step
   */
  Vec            x_n2;

  /**
   * predict solution vector
   */
  Vec            xp;

  /**
   * local truncation error vector;
   */
  Vec            LTE;



  // extra KSP solver for Trace mode

  /**
   * vec for dI/dx, I is the current of trace electrode
   */
  Vec          pdI_pdx;

  /**
   * vec for df(x)/dV, V is the app. voltage of trace electrode
   */
  Vec          pdF_pdV;

  /**
   * vec for dx/dV. J dx/dV = df(x)/dV, J=df(x)/dx is the Jacobian matrix
   */
  Vec          pdx_pdV;

  /**
   * dI/dV = dI/dx * dx/dV
   */
  PetscScalar  dI_dV;

  /**
   * ksp solver for Trace mode
   */
  KSP          kspc;

  /**
   * PC for Trace mode
   */
  PC           pcc;

  /**
   * create ksp solver for trace mode
   */
  void solve_iv_trace_begin();

  /**
   * destroy ksp solver for trace mode
   */
  void solve_iv_trace_end();

  /**
   * virtual function for set electrode dI/dV, each ddm solver should re-implement this function
   */
  virtual void set_trace_electrode(BoundaryCondition *)
  { genius_error(); }

  /**
   * x norm of potential
   */
  PetscScalar potential_norm;

  /**
   * x norm of electron density
   */
  PetscScalar electron_norm;

  /**
   * x norm of hole density
   */
  PetscScalar hole_norm;

  /**
   * x norm of lattice temperature
   */
  PetscScalar temperature_norm;

  /**
   * x norm of electron temperature
   */
  PetscScalar elec_temperature_norm;

  /**
   * x norm of hole temperature
   */
  PetscScalar hole_temperature_norm;

  /**
   * f norm of poisson's equation
   */
  PetscScalar poisson_norm;

  /**
   * f norm of electron continuity equation
   */
  PetscScalar elec_continuity_norm;

  /**
   * f norm of hole continuity equation
   */
  PetscScalar hole_continuity_norm;

  /**
   * f norm of lattice temperature equation
   */
  PetscScalar heat_equation_norm;

  /**
   * f norm of electron energy balance equation
   */
  PetscScalar elec_energy_equation_norm;

  /**
   * f norm of hole energy balance equation
   */
  PetscScalar hole_energy_equation_norm;

  /**
   * f norm of electron IV equation
   */
  PetscScalar electrode_norm;

  /**
   * nonlinear function norm
   */
  PetscScalar function_norm;

  /**
   * nonlinear iteration
   */
  PetscInt nonlinear_iteration;

};

#endif //#define __ddm_solver_h__
