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
#include "singularvalue_hook.h"


/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
SingularValueHook::SingularValueHook ( SolverBase & solver, const std::string & name, void * param )
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
SingularValueHook::~SingularValueHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void SingularValueHook::on_init()
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
void SingularValueHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void SingularValueHook::post_solve()
{
  this->solution_count++;
  this->iteration_count=0;
}



/*----------------------------------------------------------------------
 *  This is executed before each (nonlinear) iteration
 */
void SingularValueHook::pre_iteration()
{
  FVM_FlexNonlinearSolver & nonlinear_solver = dynamic_cast<FVM_FlexNonlinearSolver &>(_solver);
  const SimulationSystem &system = nonlinear_solver.get_system();

  // calculate the eigen value as well as the smallest eigen vector
  nonlinear_solver.condition_number_of_jacobian_matrix();
}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void SingularValueHook::on_close()
{}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new SingularValueHook ( solver, name, fun_data );
  }

}

#endif


