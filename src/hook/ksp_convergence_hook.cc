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


#include "solver_base.h"
#include "fvm_flex_nonlinear_solver.h"
#include "ksp_convergence_hook.h"
#include "parallel.h"


KSPConvergenceHook::KSPConvergenceHook(SolverBase & solver, const std::string & name, void *)
  :Hook(solver, name)
{
  this->_poisson_solver = false;
  this->_ddm_solver = false;
  this->solution_count=0;
  this->iteration_count=0;
}


KSPConvergenceHook::~KSPConvergenceHook()
{}



/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void KSPConvergenceHook::on_init()
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
      _solver.solver_type() == SolverSpecify::DDML1MIX  ||
      _solver.solver_type() == SolverSpecify::DDML2MIX  ||
      _solver.solver_type() == SolverSpecify::EBML3MIX  ||
      _solver.solver_type() == SolverSpecify::HALLDDML1
    )
  {
    _ddm_solver = true;
  }
}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void KSPConvergenceHook::pre_solve()
{

}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void KSPConvergenceHook::post_solve()
{
  this->solution_count++;
  this->iteration_count=0;
}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void KSPConvergenceHook::post_iteration()
{
  if ( !Genius::processor_id() )
  {
    const FVM_FlexNonlinearSolver & nonlinear_solver = dynamic_cast<const FVM_FlexNonlinearSolver &>(_solver);
    std::vector<double>  residual_history = nonlinear_solver.ksp_residual_history();

    std::string ksp_prefix = SolverSpecify::out_prefix+".residual_history";
    std::ostringstream ksp_filename;
    ksp_filename << ksp_prefix << '.' << this->solution_count<< '.' << this->iteration_count << ".dat";

    std::ofstream   fout;
    fout.open(ksp_filename.str().c_str());

    fout << "# KSP residual history" << std::endl;
    for(unsigned int n=0; n<residual_history.size(); ++n)
      fout << std::setw(10) << n << std::setw(15) << residual_history[n] << std::endl;

    fout.close();
  }


  this->iteration_count++;
}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void KSPConvergenceHook::on_close()
{

}



#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook (SolverBase & solver, const std::string & name, void * fun_data)
  {
    return new KSPConvergenceHook(solver, name, fun_data );
  }
}

#endif

