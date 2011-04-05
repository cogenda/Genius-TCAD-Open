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

#include "electrical_source.h"
#include "hdm/hdm.h"
#include "hdm/linear_poisson.h"


HDMSolver::HDMSolver(SimulationSystem & system)
    :FVM_ExplicitSolver(system)
{
  // we need a linear poisson solver
  poisson_solver = new LinearPoissonSolver(system);
}



HDMSolver::~HDMSolver()
{
  delete poisson_solver;
  FVM_Node::set_solver_index(0);
  BoundaryCondition::set_solver_index(0);
}


int HDMSolver::create_solver()
{
  MESSAGE<< '\n' << "Hydrodynamic Solver init..." << std::endl;
  RECORD();

  poisson_solver->create_solver();

  // create explicit hdm context
  FVM_Node::set_solver_index(0);
  BoundaryCondition::set_solver_index(0);
  setup_explicit_data();


  return FVM_ExplicitSolver::create_solver();
}


int HDMSolver::solve()
{

  switch( SolverSpecify::Type )
  {
  case SolverSpecify::EQUILIBRIUM :
    solve_equ();
    break;

  case SolverSpecify::STEADYSTATE:
    solve_steadystate();
    break;

  case SolverSpecify::TRANSIENT:
    solve_transient();
    break;

  default:
    MESSAGE<< '\n' << "HDMSolver: Unsupported solve type.";
    RECORD();
    break;
  }
  return 0;
}


int HDMSolver::destroy_solver()
{
  poisson_solver->destroy_solver();

  clear_explicit_data();

  return FVM_ExplicitSolver::destroy_solver();
}


int HDMSolver::pre_solve_process ( bool load_solution )
{
  if(load_solution)
  {
    // fill solution vector with initial value and vol vector with 1.0/volume
    FVM_Node::set_solver_index(0);
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->HDM_Fill_Value(x, vol);
    }

    // consider ghost cell volume here. that is 1.0/(2.0*volume)
    BoundaryCondition::set_solver_index(0);
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      if(bc->bc_type() == NeumannBoundary || bc->bc_type() == IF_Insulator_Semiconductor )
        bc->HDM_Ghostcell_Volume(vol);
    }
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    VecAssemblyBegin(vol);
    VecAssemblyEnd(vol);
  }

  return FVM_ExplicitSolver::pre_solve_process();
}


int HDMSolver::post_solve_process()
{
  FVM_Node::set_solver_index(0);

  return FVM_ExplicitSolver::post_solve_process();
}


void HDMSolver::solve_equ()
{
  // call pre_solve_process
  this->pre_solve_process(true);
    // set electrode with transient time 0 value of stimulate source(s)
  _system.get_sources()->update ( 0 );

  // do local time advancing until convergence
  int steps=0;
  SolverSpecify::clock = 0.0;
  bool converg = false;
  do
  {
    local_time_advance(steps, 1.0/3.0);
    update_solution();
    // call post_solve_process
    this->post_solve_process();
    converg = convergence_criteria(steps);
  }
  while( !converg && ++steps < SolverSpecify::MaxIteration );
}


void HDMSolver::solve_steadystate()
{}


void HDMSolver::solve_transient()
{}



void HDMSolver::local_time_advance(int step, Real alpha)
{
  // save current solution to x_prev
  VecCopy(x, x_prev);

  PetscScalar *lxx, *ltt;

  // scatte global solution vector x to local vector lx
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray(lx, &lxx);

  VecZeroEntries(f);
  VecZeroEntries(t);

  // build flux, compute local time step
  FVM_Node::set_solver_index(0);
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->HDM_Flux(lxx, f, t);
  }

  // process solid wall boundary by ghost cell
  InsertMode add_value_flag = NOT_SET_VALUES;
  BoundaryCondition::set_solver_index(0);
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    if(bc->bc_type() == NeumannBoundary || bc->bc_type() == IF_Insulator_Semiconductor )
      bc->HDM_Boundary(lxx, f, add_value_flag);
  }
  VecRestoreArray(lx, &lxx);


  VecAssemblyBegin(f);
  VecAssemblyBegin(t);
  VecAssemblyEnd(f);
  VecAssemblyEnd(t);

  //global time step

  PetscScalar dt;
  PetscInt p;
  VecMin(t, &p, &dt);
  VecSet(t, dt);

  SolverSpecify::clock += dt;

  VecScale(t, alpha);

  // update flux term
  VecPointwiseMult(f, vol, f); // f <- f/vol
  VecPointwiseMult(f, t, f); // f <- t*f
  VecAXPY(x, -1, f); // x^(n+1) = x^(n) - t*f^(n)


#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // update source term
  VecScatterBegin(scatter, t, lt, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, t, lt, INSERT_VALUES, SCATTER_FORWARD);

  // solve linear poisson's equation
  //if(step%100==0)
  {
    sync_rho();
    poisson_solver->solve();
  }

  VecGetArray(lx, &lxx);
  VecGetArray(lt, &ltt);
  FVM_Node::set_solver_index(0);
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->HDM_Source(lxx, ltt, x);
  }

  VecAssemblyBegin(x);
  VecAssemblyEnd(x);



  // add boundary constrain
  add_value_flag = NOT_SET_VALUES;  // flag for indicate ADD_VALUES operator.
  BoundaryCondition::set_solver_index(0);
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    if(bc->bc_type() != NeumannBoundary && bc->bc_type() != IF_Insulator_Semiconductor)
      bc->HDM_Boundary(lxx, x, add_value_flag);
  }

  VecAssemblyBegin(x);
  VecAssemblyEnd(x);


  VecRestoreArray(lx, &lxx);
  VecRestoreArray(lt, &ltt);

  // hook
  hook_list()->post_iteration();

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}


bool HDMSolver::convergence_criteria(int step)
{
  PetscScalar rtol=0.0, atol=0.0;

  VecWAXPY( w, -1, x_prev, x); //  w = x - x_prev
  VecPointwiseDivide( w, w, t);// here w = (x - x_prev)/t
  VecNorm(w, NORM_2, &atol);

  VecAXPY( x_prev, -1,  x); //   x_prev <- x - x_prev
  VecSet(w, 1e-12);
  VecAXPY(w, 1.0, x); // w <- x+1e-12
  VecPointwiseDivide( w, x_prev, w);// here w <- (x - x_prev)/x
  VecNorm(w, NORM_2, &rtol);

  if(Genius::processor_id() == 0)
    std::cout<<"step="<<step<<" " <<atol<<" "<<rtol<<std::endl;
  return atol < 1e-4 && rtol < 1e-4;
}


void HDMSolver::update_solution()
{
  // scatte global solution vector x to local vector lx
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray(lx, &lxx);

  // for all the regions
  FVM_Node::set_solver_index(0);
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->HDM_Update_Solution(lxx);
  }

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);
}


void HDMSolver::sync_rho()
{
  // scatte global solution vector x to local vector lx
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray(lx, &lxx);

  // for all the regions
  FVM_Node::set_solver_index(0);
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    if( region->type() != SemiconductorRegion )  continue;

    SimulationRegion::local_node_iterator it = region->on_local_nodes_begin();
    SimulationRegion::local_node_iterator it_end = region->on_local_nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = (*it);
      FVM_NodeData * node_data = fvm_node->node_data(); genius_assert(node_data);

      PetscScalar n    = lxx[fvm_node->local_offset()];
      PetscScalar p    = lxx[fvm_node->local_offset()+4];
      node_data->rho() = node_data->Net_doping() - n + p;
      //node_data->rho() = node_data->Net_doping() - 0.5*(node_data->n() + n) + 0.5*(node_data->p() + p);
      node_data->n()   = n;
      node_data->p()   = p;
    }
  }

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);
}


