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

#include "fvm_flex_nonlinear_solver.h"
#include "lsmonitor_hook.h"


/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
LinearSolverMonitorHook::LinearSolverMonitorHook ( SolverBase & solver, const std::string & name, void * param )
    : Hook ( solver, name )
{
  this->_poisson_solver = false;
  this->_ddm_solver = false;
  this->solution_count=0;
  this->iteration_count=0;
}


/*----------------------------------------------------------------------
 * destructor, close file
 */
LinearSolverMonitorHook::~LinearSolverMonitorHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void LinearSolverMonitorHook::on_init()
{
  this->_poisson_solver = false;
  this->_ddm_solver = false;

  if(_solver.solver_type() == SolverSpecify::POISSON)
  {
    _poisson_solver = true;
  }

  if( _solver.solver_type() == SolverSpecify::DDML1 ||
      _solver.solver_type() == SolverSpecify::DDML2 ||
      _solver.solver_type() == SolverSpecify::EBML3 ||
      _solver.solver_type() == SolverSpecify::DDML1MIXA ||
      _solver.solver_type() == SolverSpecify::DDML2MIXA ||
      _solver.solver_type() == SolverSpecify::EBML3MIXA ||
      _solver.solver_type() == SolverSpecify::HALLDDML1
    )
  {
    _ddm_solver = true;
  }


}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void LinearSolverMonitorHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void LinearSolverMonitorHook::post_solve()
{
  this->solution_count++;
  this->iteration_count=0;
}



/*----------------------------------------------------------------------
 *  This is executed before each (nonlinear) iteration
 */
void LinearSolverMonitorHook::pre_iteration()
{}


/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void LinearSolverMonitorHook::post_iteration()
{
  FVM_FlexNonlinearSolver & nonlinear_solver = dynamic_cast<FVM_FlexNonlinearSolver &>(_solver);
  const SimulationSystem &system = nonlinear_solver.get_system();


  // save jacobian matrix
  {
    std::string matrix_prefix = SolverSpecify::out_prefix+".jacobian";
    std::ostringstream matrix_filename;
    matrix_filename << matrix_prefix << '.' << this->solution_count<< '.' << this->iteration_count << ".mat";
    nonlinear_solver.dump_matrix_petsc(nonlinear_solver.jacobian_matrix(), matrix_filename.str());
  }

  // save rhs vector
  {
    std::string vector_prefix = SolverSpecify::out_prefix+".function";
    std::ostringstream vector_filename;
    vector_filename << vector_prefix << '.' << this->solution_count<< '.' << this->iteration_count << ".vec";

    nonlinear_solver.dump_vector_petsc(nonlinear_solver.rhs_vector(), vector_filename.str());
  }

  // save solution vector
  {
    std::string vector_prefix = SolverSpecify::out_prefix+".solution";
    std::ostringstream vector_filename;
    vector_filename << vector_prefix << '.' << this->solution_count<< '.' << this->iteration_count << ".vec";

    nonlinear_solver.dump_vector_petsc(nonlinear_solver.solution_vector(), vector_filename.str());
  }

  this->iteration_count++;

}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void LinearSolverMonitorHook::on_close()
{}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new LinearSolverMonitorHook ( solver, name, fun_data );
  }
}

#endif


