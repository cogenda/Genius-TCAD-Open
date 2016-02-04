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
  virtual void dump_matrix_petsc(const Mat mat, const std::string &file) const;

  /**
   * dump jacobian matrix in Harwell-Boeing format to external file
   * for more detailed analysis of the properties of jacobian matrix
   * when rhs is true, also dump rhs vector
   * however, this must run in serial
   */
  //void dump_jacobian_matrix_hb(const std::string &file, bool rhs=false) const;

  /**
   * dump jacobian matrix in asc format to external file
   * for more detailed analysis of the properties of jacobian matrix
   */
  virtual void dump_matrix_asc(const Mat mat, const std::string &file) const;

   /**
   * dump jacobian matrix in triplet format to external file
   * for more detailed analysis of the properties of jacobian matrix
    */
  virtual void dump_matrix_triplet(const Mat mat, const std::string &file) const;

  /**
   * dump vector to external file
   * for more detailed analysis of the properties
   */
  virtual void dump_vector_petsc(const Vec vec, const std::string &file) const;

  /**
   * load vector from external file
   * for more detailed analysis of the properties
   */
  virtual void load_vector_petsc(Vec vec, const std::string &file) const;

  /**
   * @return the condition number of jacobian matrix
   */
  double condition_number_of_jacobian_matrix();

  /**
   * calculate the n largest and smallest eigen value of jacobian matrix
   * optinally get the ith smallest vec and jth largest vec
   */
  void eigen_value_of_jacobian_matrix(int n=1, int i=0, Vec = PETSC_NULL, int j=0, Vec = PETSC_NULL);

  /**
   * @return the jacobian_matrix
   */
  Mat jacobian_matrix() const  { return J; }

  /**
   * @return the rhs vector
   */
  Vec rhs_vector() const { return f; }

  /**
   * @return the solution vector
   */
  Vec solution_vector() const { return x; }

  /**
   * @return the ksp residual history
   */
  std::vector<double> ksp_residual_history() const;

  /**
   * create a new vector with (compatable parallel) pattern
   */
  void create_vector(Vec &);

  /**
   * destroy a vector
   */
  void destroy_vector(Vec &);

  /**
   * do snes solve!
   */
  virtual void snes_solve();

  /**
   * clear all the nonlinear solver contex
   */
  void clear_nonlinear_data();


  /**
   * Sets the type of nonlinear solver to use.
   */
  void set_nonlinear_solver_type (const SolverSpecify::NonLinearSolverType nst)
  {  _nonlinear_solver_type = nst; }

  /**
   * Returns the type of nonlinear solver to use.
   */
  SolverSpecify::NonLinearSolverType nonlinear_solver_type () const { return _nonlinear_solver_type; }

  /**
   * Returns the type of linear solver to use.
   */
  SolverSpecify::LinearSolverType linear_solver_type () const { return _linear_solver_type; }

  /**
   * Sets the type of solver to use.
   */
  void set_linear_solver_type (const SolverSpecify::LinearSolverType st)
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
   * PETSC SNES can have an individual prefix
   */
  virtual std::string snes_prefix() const = 0;

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

  /**
   * virtual function, write solver intermediate data into system
   * It can be used to monitor the field data evolution during solve action
   */
  virtual void flush_system(Vec ) {}

protected:


  /**ksp_residual_history
   * incidate that the jacobian_matrix is never assembled
   */
  bool jacobian_matrix_first_assemble;

  /**
   * which type of nonlinear solver to use.
   */
  void set_petsc_nonelinear_solver_type();

  /**
   * which type of linear solver to use.
   */
  void set_petsc_linear_solver_type();

  /**
   * which type of proconditioner to use.
   */
  void set_petsc_preconditioner_type();

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
   * the local scaling vector
   */
  Vec            ll;
  
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


#if PETSC_VERSION_GE(3,3,0)
  /**
   * petsc nonlinear solver linesearch contex
   */
  SNESLineSearch snesls;
#else
  /**
   * petsc < 3.3 do not have snesls object
   */
  #define  snesls snes

#endif

  /**
   * the KSP context for a SNES solver
   */
  KSP            ksp;

  /**
   * the PC context for KSP solver
   */
  PC             pc;

  /**
   * array for ksp residual history
   */
  std::vector<PetscReal> _ksp_residual_history;

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


  /**
   * all the petsc options, will be delete when this class is destroied.
   */
  std::map<std::string, std::string> petsc_options;

  /**
   * apply and buffer petsc options
   */
  int set_petsc_option(const std::string &key, const std::string &value, bool has_prefix=true);

};


#endif // #define __fvm_nonlinear_solver_h__
