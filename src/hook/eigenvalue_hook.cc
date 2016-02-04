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

#include "genius_common.h"
#include "fvm_flex_nonlinear_solver.h"
#include "eigenvalue_hook.h"


/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
EigenValueHook::EigenValueHook ( SolverBase & solver, const std::string & name, void * param )
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
EigenValueHook::~EigenValueHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void EigenValueHook::on_init()
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
void EigenValueHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void EigenValueHook::post_solve()
{
  this->solution_count++;
  this->iteration_count=0;
}



/*----------------------------------------------------------------------
 *  This is executed before each (nonlinear) iteration
 */
void EigenValueHook::pre_iteration()
{}


/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void EigenValueHook::post_iteration()
{
  FVM_FlexNonlinearSolver & nonlinear_solver = dynamic_cast<FVM_FlexNonlinearSolver &>(_solver);
  const SimulationSystem &system = nonlinear_solver.get_system();

  // create a vector for eigenvector
  Vec Vrs;
  nonlinear_solver.create_vector(Vrs);

  // calculate the eigen value as well as the smallest eigen vector
  nonlinear_solver.eigen_value_of_jacobian_matrix(5, 0, Vrs);

#if 0
  // export the value of eigen vector to vtk file
  // NOTE: will change the state of system
  {
    std::string vtk_prefix = SolverSpecify::out_prefix+".eigenvector";
    std::ostringstream vtk_filename;
    vtk_filename << vtk_prefix << '.' << this->solution_count<< '.' << this->iteration_count << ".vtu";
    // save the value of eigen vector to mesh nodes
    nonlinear_solver.flush_system(Vrs);
    system.export_vtk ( vtk_filename.str(), false );
  }
#endif

#if 0
  // save eigen vector
  {
    std::string vector_prefix = SolverSpecify::out_prefix+".eigenvector";
    std::ostringstream vector_filename;
    vector_filename << vector_prefix << '.' << this->solution_count<< '.' << this->iteration_count << ".vec";
    nonlinear_solver.dump_vector_petsc(Vrs, vector_filename.str());
  }
#endif

  this->iteration_count++;

  nonlinear_solver.destroy_vector(Vrs);
}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void EigenValueHook::on_close()
{}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new EigenValueHook ( solver, name, fun_data );
  }

}

#endif


