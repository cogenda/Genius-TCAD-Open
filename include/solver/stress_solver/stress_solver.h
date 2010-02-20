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

//  $Id: linear_solver.h,v 1.3 2008/07/09 05:58:16 gdiso Exp $

#ifndef __stress_solver_h__
#define __stress_solver_h__

#include <vector>

#include "enum_petsc_type.h"
#include "fvm_linear_solver.h"
//#include "petscis.h"
//#include "petscvec.h"
//#include "petscmat.h"
#include "petscksp.h"


/**
 * The linear solver contex.
 */
class StressSolver : public FVM_LinearSolver
{
public:

  /**
   * construct the linear solver contex:
   * solution vector, right hand side (RHS) vector, matrix
   * as well as parallel scatter
   */
  StressSolver(SimulationSystem & system, Parser::InputParser & decks)
  : FVM_LinearSolver(system), _decks(decks)
  {system.record_active_solver(this->solver_type());}

  /**
   * free all the contex
   */
  virtual ~StressSolver(){}

  /**
   * @return the solver type
   */
  virtual SolverSpecify::SolverType solver_type() const
  {return SolverSpecify::STRESS;}

  /**
   * virtual function, create the solver
   */
  virtual int create_solver();

  /**
   * virtual function, do the solve process
   */
  virtual int solve();

  /**
   * virtual function, destroy the solver, release internal data
   */
  virtual int destroy_solver() ;

  /**
   * virtual function for building the RHS vector
   */
  virtual void build_rhs(Vec b);

  /**
   * virtual function for building the matrix A and precondition matrix PC
   */
  virtual void build_matrix(Mat A, Mat pc);


  /**
   * @return node's dof for each region. only 1 dof here
   */
  virtual unsigned int node_dofs(const SimulationRegion * region) const
  { assert(region!=NULL); return 1; }

  /**
   * @return the dofs of each boundary condition
   */
  virtual unsigned int bc_dofs(const BoundaryCondition * bc) const
  { assert(bc!=NULL); return 0; }

  /**
   * @return the boundary condition dofs contributed by each boundary node
   */
  virtual unsigned int bc_node_dofs(const BoundaryCondition * bc) const
  { assert(bc!=NULL); return 0; }

  virtual bool all_neighbor_elements_involved(const SimulationRegion *) const
  { return false; }

 private:

   /**
  * since reading deck involves stack operate, we can not use const here
  */
 Parser::InputParser & _decks;

};


#endif // #define __stress_solver_h__

