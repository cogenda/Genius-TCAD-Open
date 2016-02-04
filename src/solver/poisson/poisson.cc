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

//  $Id: poisson.cc,v 1.36 2008/07/09 07:53:36 gdiso Exp $



#include "poisson/poisson.h"
#include "solver_specify.h"
#include "electrical_source.h"
#include "parallel.h"
#include "petsc_utils.h"

#define DEBUG


using PhysicalUnit::kb;
using PhysicalUnit::e;


/*------------------------------------------------------------------
 * create the poisson solver contex
 */
int PoissonSolver::create_solver()
{

  MESSAGE<< "\nPoisson Solver init..." << std::endl;
  RECORD();

  set_nonlinear_solver_type ( SolverSpecify::NS );
  set_linear_solver_type    ( SolverSpecify::LS );
  set_preconditioner_type   ( SolverSpecify::PC );

  // must set nonlinear matrix/vector here!
  setup_nonlinear_data();
  //distable lag pc/jacobian
  SNESSetLagPreconditioner(snes, 1);
  SNESSetLagJacobian(snes, 1);

  //abstol = 1e-12*n_global_dofs        - absolute convergence tolerance
  //rtol   = 1e-14                      - relative convergence tolerance
  //stol   = 1e-9                       - convergence tolerance in terms of the norm of the change in the solution between steps
  SNESSetTolerances(snes, 1e-12*n_global_dofs, 1e-14, 1e-9, SolverSpecify::MaxIteration, 1000);

  // rtol   = 1e-12*n_global_dofs  - the relative convergence tolerance (relative decrease in the residual norm)
  // abstol = 1e-20*n_global_dofs  - the absolute convergence tolerance (absolute size of the residual norm)
  KSPSetTolerances(ksp, 1e-12*n_global_dofs, 1e-20*n_global_dofs, PETSC_DEFAULT, std::max(50, static_cast<int>(n_global_dofs/10)) );

  // user can do further adjusment from command line
  SNESSetFromOptions (snes);

  // init (user defined) hook functions here


  return FVM_FlexNonlinearSolver::create_solver();

}



/*------------------------------------------------------------------
 * set initial value to solution vector and scaling vector
 */
int PoissonSolver::pre_solve_process(bool load_solution)
{

  // search for all the regions, call corresponding function
  for(unsigned int n=0; n<_system.n_regions(); ++n)
  {
    SimulationRegion * region = _system.region(n);
    region->Poissin_Fill_Value(x, L);
  }

  // for all the bcs
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->Poissin_Fill_Value(x, L);
  }

  VecAssemblyBegin(x);
  VecAssemblyBegin(L);

  VecAssemblyEnd(x);
  VecAssemblyEnd(L);

  return FVM_FlexNonlinearSolver::pre_solve_process(load_solution);
}




/*------------------------------------------------------------------
 *
 */
int PoissonSolver::solve()
{

  START_LOG("solve()", "PoissonSolver");

  // set each electrode with external stimulate. transient 0 value is used here
  _system.get_electrical_source()->update(0);

  // call pre_solve_process
  pre_solve_process();

  // here call Petsc to solve the nonlinear Poisson's equation
  snes_solve();

  // get the converged reason
  SNESConvergedReason reason;
  SNESGetConvergedReason(snes, &reason);

  // print convergence/divergence reason
  {
    MESSAGE
    <<"----------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[reason]<<"\n\n\n";
    RECORD();
  }

  // call post_solve_process
  post_solve_process();

  //dump_jacobian_matrix_petsc("poisson.mat");

  STOP_LOG("solve()", "PoissonSolver");

  return 0;
}


/*------------------------------------------------------------------
 * restore the solution to each region
 */
int PoissonSolver::post_solve_process()
{

  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  VecGetArray(lx, &lxx);

  // update all the regions
  for(unsigned int n=0; n<_system.n_regions(); ++n)
  {
    SimulationRegion * region = _system.region(n);
    region->Poissin_Update_Solution(lxx);
  }


  // extra work, calculate the electric field
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);

    SimulationRegion::processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = *it;
      FVM_NodeData * node_data = fvm_node->node_data();
      node_data->E() = -fvm_node->gradient(POTENTIAL, true);
    }
  }


  // update bcs, set electrode potential equal to its vapp
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); ++b)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->Poissin_Update_Solution(lxx);
  }

  VecRestoreArray(lx, &lxx);

  // call (user defined) hook function hook_post_solve_process

  return FVM_FlexNonlinearSolver::post_solve_process();
}



/*------------------------------------------------------------------
 * destroy solver
 */
int PoissonSolver::destroy_solver()
{

  // clear nonlinear contex
  clear_nonlinear_data();

  return FVM_FlexNonlinearSolver::destroy_solver();
}


void PoissonSolver::potential_damping(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{

  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;

  VecGetArray(x, &xx);  // previous iterate value
  VecGetArray(y, &yy);  // new search direction and length
  VecGetArray(w, &ww);  // current candidate iterate

  const PetscScalar T = this->get_system().T_external();
  PetscScalar dV_max = 0.0; // the max changes of psi

  // we should find dV_max;
  // first, we find in local
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = _system.region(n);
    if( region->type() != SemiconductorRegion ) continue;

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = *it;
      // we konw the fvm_node->local_offset() is psi in semiconductor region
      unsigned int local_offset = fvm_node->local_offset();
      dV_max = std::max(dV_max, std::abs(yy[local_offset]));
    }
  }

  // for parallel situation, we should find the dv_max in global.
  Parallel::max( dV_max );

  if( dV_max > 1e-6*kb*T/e )
  {
    // compute logarithmic potential damping factor f;
    PetscScalar Vut = kb*T/e * SolverSpecify::potential_update;
    PetscScalar f = log(1+dV_max/Vut)/(dV_max/Vut);

    // do newton damping here
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      const SimulationRegion * region = _system.region(n);
      // only consider semiconductor region
      //if( region->type() != SemiconductorRegion ) continue;

      SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        const FVM_Node * fvm_node = *it;
        unsigned int local_offset = fvm_node->local_offset();
        ww[local_offset] = xx[local_offset] - f*yy[local_offset];
      }
    }
  }

  VecRestoreArray(x, &xx);
  VecRestoreArray(y, &yy);
  VecRestoreArray(w, &ww);

  *changed_y = PETSC_FALSE;
  *changed_w = PETSC_TRUE;

  return;
}


/////////////////////////////////////////////////////////////////////////
// provide function and jacobian evaluation function for poisson solver//
/////////////////////////////////////////////////////////////////////////





/*------------------------------------------------------------------
 * evaluate the residual of function f at x
 */
void PoissonSolver::build_petsc_sens_residual(Vec x, Vec r)
{

  START_LOG("PoissonSolver_Residual()", "PoissonSolver");

  // scatte global solution vector x to local vector lx
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray(lx, &lxx);

  // clear old data
  VecZeroEntries (r);

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate Poisson's equation in all the regions
  for(unsigned int n=0; n<_system.n_regions(); ++n)
  {
    SimulationRegion * region = _system.region(n);
    region->Poissin_Function(lxx, r, add_value_flag);
  }

  // process hanging node here
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->Poissin_Function_Hanging_Node(lxx, r, add_value_flag);
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // preprocess each bc
  VecAssemblyBegin(r);
  VecAssemblyEnd(r);
  std::vector<PetscInt> src_row,  dst_row,  clear_row;
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->Poissin_Function_Preprocess(lxx, r, src_row, dst_row, clear_row);
  }
  //add source rows to destination rows, and clear rows
  PetscUtils::VecAddClearRow(r, src_row, dst_row, clear_row);
  add_value_flag = NOT_SET_VALUES;

  // evaluate Poisson's equation for all the boundaries
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); ++b)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->Poissin_Function(lxx, r, add_value_flag);
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);

  // final assembly the function Vec
  VecAssemblyBegin(r);
  VecAssemblyEnd(r);

  // scale the function vec
  VecPointwiseMult(r, r, L);

  STOP_LOG("PoissonSolver_Residual()", "PoissonSolver");

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}



/*------------------------------------------------------------------
 * evaluate the Jacobian J of function f at x
 */
void PoissonSolver::build_petsc_sens_jacobian(Vec x, Mat *, Mat *)
{

  START_LOG("Poissin_Jacobian()", "PoissonSolver");

  // scatte global solution vector x to local vector lx
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray(lx, &lxx);

  Jac->zero();

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

  // evaluate Jacobian matrix of Poisson's equation in all the regions
  for(unsigned int n=0; n<_system.n_regions(); ++n)
  {
    SimulationRegion * region = _system.region(n);
    region->Poissin_Jacobian(lxx, Jac, add_value_flag);
  }

  // process hanging node here
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->Poissin_Jacobian_Hanging_Node(lxx, Jac, add_value_flag);
  }


#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif


  // assembly matrix
  Jac->close(false);


  // evaluate Jacobian matrix of governing equations of DDML1 for all the boundaries
  std::vector<PetscInt> src_row,  dst_row,  clear_row;
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->Poissin_Jacobian_Preprocess(lxx, Jac, src_row, dst_row, clear_row);
  }

  //add source rows to destination rows
  Jac->add_row_to_row(src_row, dst_row);
  // clear row
  Jac->clear_row(clear_row);

  add_value_flag = NOT_SET_VALUES;
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->Poissin_Jacobian(lxx, Jac, add_value_flag);
  }

  Jac->close(true);

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);

  //scaling the matrix
  MatDiagonalScale(J, L, PETSC_NULL);

  STOP_LOG("Poissin_Jacobian()", "PoissonSolver");

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}


