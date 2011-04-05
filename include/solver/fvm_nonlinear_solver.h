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


#ifndef __fvm_nonlinear_solver_h__
#define __fvm_nonlinear_solver_h__

#include <vector>

#include "config.h"
#include "enum_petsc_type.h"
#include "fvm_pde_solver.h"
//#include "petscis.h"
//#include "petscvec.h"
//#include "petscmat.h"
//#include "petscksp.h"
#include "petscsnes.h"



/**
 * The nonlinear solver contex.
 */
class FVM_NonlinearSolver : public FVM_PDESolver
{
public:

  /**
   * constructor
   */
  FVM_NonlinearSolver(SimulationSystem & system);


  /**
   * destructor
   */
  virtual ~FVM_NonlinearSolver();

  /**
   * set all the nonlinear solver contex:
   * solution vector, function vector, jacobian matrix
   * as well as parallel scatter
   */
  void setup_nonlinear_data();

  /**
   * dump jacobian matrix in petsc format to external file
   * for more detailed analysis of the properties of jacobian matrix
   */
  void dump_jacobian_matrix_petsc(const std::string &file) const;

  /**
   * dump jacobian matrix in Harwell-Boeing format to external file
   * for more detailed analysis of the properties of jacobian matrix
   * when rhs is true, also dump rhs vector
   * however, this must run in serial
   */
  //void dump_jacobian_matrix_hb(const std::string &file, bool rhs=false) const;

  /**
   * dump function to external file
   * for more detailed analysis of the properties
   */
  void dump_function_vector_petsc(const std::string &file) const;

  /**
   * @return the condition number of jacobian matrix
   */
  double condition_number_of_jacobian_matrix();

  /**
   * calculate eigen value of jacobian matrix
   */
  void eigen_value_of_jacobian_matrix();

  /**
   * @return the jacobian_matrix
   */
  Mat & jacobian_matrix()  { return J; }

  /**
   * @return the rhs vector
   */
  Vec & rhs_vector()  { return f; }

  /**
   * @return the solution vector
   */
  Vec & solution_vector()  { return x; }

  /**
   * do snes solve!
   */
  virtual void sens_solve();

  /**
   * clear all the nonlinear solver contex
   */
  void clear_nonlinear_data();

  /**
   * Returns the type of solver to use.
   */
  SolverSpecify::LinearSolverType linear_solver_type () const { return _linear_solver_type; }

  /**
   * Sets the type of solver to use.
   */
  void set_solver_type (const SolverSpecify::LinearSolverType st)
  { _linear_solver_type = st; }

  /**
   * Returns the type of preconditioner to use.
   */
  SolverSpecify::PreconditionerType preconditioner_type () const
  { return _preconditioner_type; }

  /**
   * Sets the type of preconditioner to use.
   */
  void set_preconditioner_type (const SolverSpecify::PreconditionerType pct)
  { _preconditioner_type = pct; }


  /**
   * pure virtual function for evaluating the residual of function f at x
   */
  virtual void build_petsc_sens_residual(Vec x, Vec r)=0;

  /**
   * pure virtual function for evaluating the Jacobian J of function f at x
   */
  virtual void build_petsc_sens_jacobian(Vec x, Mat *jac, Mat *pc)=0;

  /**
   * virtual function for snes monitor. derived class can override it as needed.
   */
  virtual void petsc_snes_monitor(PetscInt its, PetscReal fnorm);

  /**
   * virtual function for snes convergence test. derived class can override it as needed.
   */
  virtual void petsc_snes_convergence_test(PetscInt its, PetscReal xnorm, PetscReal gnorm, PetscReal fnorm, SNESConvergedReason *reason);

  /**
   * virtual function for ksp convergence test. derived class can override it as needed.
   */
  virtual void petsc_ksp_convergence_test(PetscInt its, PetscReal rnorm, KSPConvergedReason* reason);

  /**
   * virtual function for line search pre check. derived class can override it as needed.
   */
  virtual void sens_line_search_pre_check(Vec x, Vec y, PetscBool *changed_y);

  /**
   * virtual function for line search post check. derived class can override it as needed.
   */
  virtual void sens_line_search_post_check(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w);


protected:


  /**
   * incidate that the jacobian_matrix is never assembled
   */
  bool jacobian_matrix_first_assemble;

  /**
   * which type of nonlinear solver to use.
   */
  void set_petsc_nonelinear_solver_type(SolverSpecify::NonLinearSolverType t);

  /**
   * which type of linear solver to use.
   */
  void set_petsc_linear_solver_type(SolverSpecify::LinearSolverType t);

  /**
   * which type of proconditioner to use.
   */
  void set_petsc_preconditioner_type(SolverSpecify::PreconditionerType t);

  /**
   * the global solution vector
   */
  Vec            x;

  /**
   * the global nonlinear function vector
   */
  Vec            f;

  /**
   * the jacobian matrix
   */
  Mat            J;

  /**
   * the left scaling vector of J
   */
  Vec            L;

  /**
   * the local solution vector
   */
  Vec            lx;

  /**
   * the local function vector
   */
  Vec            lf;

  /**
   * the global index set of solution vector this processor needs
   */
  IS             gis;

  /**
   * the local index set of solution vector this processor needs
   */
  IS             lis;

  /**
   * data structure used to map global solution vector to local vector
   */
  VecScatter     scatter;

  /**
   * petsc nonlinear solver contex
   */
  SNES           snes;

  /**
   * the KSP context for a SNES solver
   */
  KSP            ksp;

  /**
   * the PC context for KSP solver
   */
  PC             pc;

  /**
   * Enum stating which type of nonlinear solver to use.
   */
  SolverSpecify::NonLinearSolverType _nonlinear_solver_type;


  /**
   * Enum stating which type of iterative solver to use.
   */
  SolverSpecify::LinearSolverType   _linear_solver_type;

  /**
   * a stack for previous linear solvers
   */
  std::vector<SolverSpecify::LinearSolverType>  _linear_solver_stack;

  /**
   * Enum statitng with type of preconditioner to use.
   */
  SolverSpecify::PreconditionerType _preconditioner_type;

};


#endif // #define __fvm_nonlinear_solver_h__
