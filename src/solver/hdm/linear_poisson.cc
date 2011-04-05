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

#include "hdm/linear_poisson.h"
#include "parallel.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;



int LinearPoissonSolver::create_solver()
{
  // adjust default snes/ksp/pc/mat/vec settings

  FVM_Node::set_solver_index ( 1 );
  BoundaryCondition::set_solver_index(1);

#ifdef PETSC_HAVE_MUMPS
  // use MUMPS as fast linear solver
  SolverSpecify::LS = SolverSpecify::MUMPS;
#else
  SolverSpecify::LS = SolverSpecify::LU;
#endif

  // must set linear matrix/vector here!
  setup_linear_data();

  // ask mumps to do LU factorization
  PCFactorSetMatSolverPackage ( pc, "mumps" );

  // user can do further adjusment from command line
  KSPSetFromOptions ( ksp );

  // build the matrix here
  build_matrix ( A, A );

  return 0;
}



int LinearPoissonSolver::solve()
{
  build_rhs ( b );
  KSPSolve ( ksp,b,x );
  update_solution();
  return 0;
}



int LinearPoissonSolver::destroy_solver()
{
  // clear linear contex
  clear_linear_data();

  return 0;
}


void LinearPoissonSolver::build_rhs ( Vec b )
{
  VecZeroEntries ( b );

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;
  FVM_Node::set_solver_index ( 1 );
  for ( unsigned int n=0; n<_system.n_regions(); ++n )
  {
    SimulationRegion * region = _system.region ( n );
    region->LinearPoissin_RHS ( b, add_value_flag );
  }

  BoundaryCondition::set_solver_index(1);
  for ( unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n )
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc ( n );
    bc->LinearPoissin_RHS ( b, add_value_flag );
  }

  VecAssemblyBegin ( b );
  VecAssemblyEnd ( b );
}


void LinearPoissonSolver::build_matrix ( Mat A, Mat )
{
  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;
  FVM_Node::set_solver_index ( 1 );
  for ( unsigned int n=0; n<_system.n_regions(); ++n )
  {
    SimulationRegion * region = _system.region ( n );
    region->LinearPoissin_Matrix ( A, add_value_flag );
  }

  BoundaryCondition::set_solver_index(1);
  if ( !matrix_first_assemble )
    for ( unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n )
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc ( n );
      bc->LinearPoissin_Reserve ( A, add_value_flag );
    }
  matrix_first_assemble = true;

  for ( unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n )
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc ( n );
    bc->LinearPoissin_Matrix ( A, add_value_flag );
  }

  MatAssemblyBegin ( A, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd ( A, MAT_FINAL_ASSEMBLY );

}





void LinearPoissonSolver::update_solution()
{
  // scatte global solution vector x to local vector lx
  VecScatterBegin ( scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD );
  VecScatterEnd ( scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD );

  PetscScalar *lxx;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray ( lx, &lxx );

  // update phi with damping
  for ( unsigned int n=0; n<_system.n_regions(); n++ )
  {
    SimulationRegion * region = _system.region ( n );
    region->LinearPoissin_Update_Solution(lxx);
  }

  // restore array back to Vec
  VecRestoreArray ( lx, &lxx );
}

