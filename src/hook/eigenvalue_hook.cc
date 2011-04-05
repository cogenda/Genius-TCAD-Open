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
#include <time.h>
#include <string>
#include <cstdlib>
#include <iomanip>


#include "fvm_nonlinear_solver.h"
#include "eigenvalue_hook.h"

#ifdef HAVE_SLEPC
#include "slepceps.h"
#include "slepcsvd.h"
#endif

using PhysicalUnit::um;

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
{

}


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
{
#ifdef HAVE_SLEPC

  EPS            eps_s;
  EPS            eps_l;

  FVM_NonlinearSolver & nonlinear_solver = dynamic_cast<FVM_NonlinearSolver &>(_solver);
  Mat & J = nonlinear_solver.jacobian_matrix();

  // create eigen value problem solver
  EPSCreate(PETSC_COMM_WORLD, &eps_s);
  EPSCreate(PETSC_COMM_WORLD, &eps_l);
  // Set operator
  EPSSetOperators(eps_s, J, PETSC_NULL);
  EPSSetOperators(eps_l, J, PETSC_NULL);

  // calculate smallest and largest eigen value
  EPSSetWhichEigenpairs(eps_s, EPS_SMALLEST_MAGNITUDE);
  EPSSetWhichEigenpairs(eps_l, EPS_LARGEST_MAGNITUDE);

  // shift-and-invert spectral transformation to enhance convergence of eigenvalues near zero
  ST st_s;
  EPSGetST(eps_s, &st_s);
  STSetType(st_s, STSINVERT);
  //EPSSetTrueResidual(eps_s, PETSC_TRUE);

  //EPSSetTolerances(eps_s, 1e-6, 100);
  //EPSSetTolerances(eps_l, 1e-6, 100);

  // Set solver parameters at runtime
  EPSSetFromOptions(eps_s);
  EPSSetFromOptions(eps_l);


  PetscScalar kr_s, ki_s;
  PetscScalar kr_l, ki_l;
  PetscReal error_s;
  PetscReal error_l;
  PetscInt nconv_s;
  PetscInt nconv_l;

  // get the smallest eigen value
  EPSSolve( eps_s );
  EPSGetConverged( eps_s, &nconv_s );
  if( nconv_s > 0 )
  {
    EPSGetEigenvalue( eps_s, 0, &kr_s, &ki_s );
    EPSComputeRelativeError( eps_s, 0, &error_s );
  }

  // get the largest eigen value
  EPSSolve( eps_l );
  EPSGetConverged( eps_l, &nconv_l );
  if( nconv_l > 0 )
  {
    EPSGetEigenvalue( eps_l, 0, &kr_l, &ki_l );
    EPSComputeRelativeError( eps_l, 0, &error_l );
  }

  if(Genius::processor_id() == 0)
  {
    if( nconv_s > 0 )
      std::cout<< "Smallest eigen value: " << std::scientific << std::setprecision(6)<<std::setw(10) << kr_s << " with error " << error_s << std::endl;
    if( nconv_l > 0 )
      std::cout<< "Largest  eigen value: " << std::scientific << std::setprecision(6)<<std::setw(10) << kr_l << " with error " << error_l << std::endl;
    if( nconv_s > 0 && nconv_l > 0 )
      std::cout<< "Approx reciprocal condition number: " << std::scientific << std::setprecision(6)<<std::setw(10) << kr_s/kr_l << std::endl;
  }


  EPSDestroy(eps_s);
  EPSDestroy(eps_l);

#endif
}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void EigenValueHook::on_close()
{}


#ifndef CYGWIN

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new EigenValueHook ( solver, name, fun_data );
  }

}

#endif


