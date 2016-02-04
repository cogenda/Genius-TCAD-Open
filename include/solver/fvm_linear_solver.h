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



#ifndef __fvm_linear_solver_h__
#define __fvm_linear_solver_h__

#include "enum_petsc_type.h"
#include "fvm_pde_solver.h"
#include "petscksp.h"


/**
 * The linear solver contex.
 */
class FVM_LinearSolver : public FVM_PDESolver
{
public:

  /**
   * construct the linear solver contex:
   * solution vector, right hand side (RHS) vector, matrix
   * as well as parallel scatter
   */
  FVM_LinearSolver(SimulationSystem & system);

  /**
   * free all the contex
   */
  virtual ~FVM_LinearSolver();

  /**
   * setup linear context
   */
  void setup_linear_data();

  /**
   * clear linear data
   */
  void clear_linear_data();

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
   * @return the condition number of matrix
   */
  double condition_number_of_matrix();

  /**
   * calculate the n largest and smallest eigen value of matrix
   * optinally get the ith smallest vec and jth largest vec
   */
  void eigen_value_of_matrix(int n=1, int i=0, Vec = PETSC_NULL, int j=0, Vec = PETSC_NULL);

  /**
   * @return the jacobian_matrix
   */
  Mat jacobian_matrix() const  { return A; }


  /**
   * @return the rhs vector
   */
  Vec rhs_vector() const { return b; }


  /**
   * @return the residual vector
   */
  Vec residual_vector() const { return r; }


  /**
   * @return the working vector
   */
  Vec working_vector() const  { return w; }


  /**
   * @return the solution vector
   */
  Vec solution_vector() const { return x; }


  /**
   * Returns the type of solver to use.
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
   * PETSC KSP can have an individual prefix
   */
  virtual std::string ksp_prefix() const = 0;

  /**
   * @return the number of linear iteration
   */
  PetscInt get_linear_iteration() const;

  /**
   * virtual function for ksp convergence test. derived class can override it as needed.
   */
  virtual void petsc_ksp_convergence_test(PetscInt its, PetscReal rnorm, KSPConvergedReason* reason);


protected:

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
   * the global right handside vector
   */
  Vec            b;

  /**
   * the global residual vector, defined as Ax-b
   */
  Vec            r;

  /**
   * the global working vector
   */
  Vec            w;

  /**
   * the matrix
   */
  Mat            A;

  /**
   * the left scaling vector of A
   */
  Vec            L;

  /**
   * the local solution vector
   */
  Vec            lx;

  /**
   * the local function vector
   */
  Vec            lb;

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
   * the KSP context for a SNES solver
   */
  KSP            ksp;

  /**
   * the PC context for KSP solver
   */
  PC             pc;

  /**
   * incidate that the matrix is never assembled
   */
  bool matrix_first_assemble;

  /**
   * Enum stating which type of iterative solver to use.
   */
  SolverSpecify::LinearSolverType   _linear_solver_type;

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


#endif // #define __fvm_linear_solver_h__

