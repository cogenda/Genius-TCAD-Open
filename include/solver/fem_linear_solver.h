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



#ifndef __fem_linear_solver_h__
#define __fem_linear_solver_h__

#include <vector>

#include "enum_petsc_type.h"
#include "fem_pde_solver.h"
#include "petscksp.h"


/**
 * The linear solver contex.
 */
class FEM_LinearSolver : public FEM_PDESolver
{
public:

  /**
   * construct the linear solver contex:
   * solution vector, right hand side (RHS) vector, matrix
   * as well as parallel scatter
   */
  FEM_LinearSolver(SimulationSystem & system);

  /**
   * free all the contex
   */
  virtual ~FEM_LinearSolver();

  /**
   * setup linear context
   */
  void setup_linear_data();

  /**
   * clear linear data
   */
  void clear_linear_data();

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
   * PETSC KSP can have an individual prefix
   */
  virtual std::string ksp_prefix() const = 0;


  /**
   * virtual function for building the RHS vector
   */
  virtual void build_rhs(Vec ) {}

  /**
   * virtual function for building the matrix A and precondition matrix PC
   */
  virtual void build_matrix(Mat , Mat ) {}



protected:

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
   * the global right handside vector
   */
  Vec            b;

  /**
   * the matrix
   */
  Mat            A;

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


#endif // #define __fem_linear_solver_h__

