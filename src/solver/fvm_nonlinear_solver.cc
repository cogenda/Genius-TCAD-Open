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


#include <numeric>

#include "fvm_nonlinear_solver.h"
#include "parallel.h"

#ifdef HAVE_SLEPC
#include "slepceps.h"
#include "slepcsvd.h"
#endif

//--------------------------------------------------------------------
// Functions with C linkage to pass to PETSc.  PETSc will call these
// methods as needed.
//
// Since they must have C linkage they have no knowledge of a namespace.
// Give them an obscure name to avoid namespace pollution.
extern "C"
{
  //-------------------------------------------------------------------
  // this function is called by PETSc at the end of each nonlinear step
  static PetscErrorCode  __genius_petsc_snes_monitor (SNES, PetscInt its, PetscReal fnorm, void *ctx)
  {
    int ierr=0;

    // convert void* to FVM_NonlinearSolver*
    FVM_NonlinearSolver * nonlinear_solver = (FVM_NonlinearSolver *)ctx;

    nonlinear_solver->petsc_snes_monitor(its, fnorm);

    return ierr;

  }

  //---------------------------------------------------------------
  // this function is called by PETSc to test for convergence of the nonlinear iterative solution.
  static PetscErrorCode __genius_petsc_snes_convergence_test
  (SNES, PetscInt its, PetscReal xnorm, PetscReal gnorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
  {
    int ierr=0;

    // convert void* to FVM_NonlinearSolver*
    FVM_NonlinearSolver * nonlinear_solver = (FVM_NonlinearSolver *)ctx;

    nonlinear_solver->petsc_snes_convergence_test(its, xnorm, gnorm, fnorm, reason);

    return ierr;
  }

  //---------------------------------------------------------------
  // this function is called by PETSc to test for convergence of the linear iterative solution.
  static PetscErrorCode __genius_petsc_ksp_convergence_test(KSP, PetscInt its, PetscReal rnorm, KSPConvergedReason* reason, void *ctx)
  {
    int ierr=0;
    // convert void* to FVM_NonlinearSolver*
    FVM_NonlinearSolver * nonlinear_solver = (FVM_NonlinearSolver *)ctx;

    nonlinear_solver->petsc_ksp_convergence_test(its, rnorm, reason);

    return ierr;
  }

  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the residual at X
  static PetscErrorCode  __genius_petsc_snes_residual (SNES, Vec x, Vec f, void *ctx)
  {
    int ierr=0;

    // convert void* to FVM_NonlinearSolver*
    FVM_NonlinearSolver * nonlinear_solver = (FVM_NonlinearSolver *)ctx;

    nonlinear_solver->build_petsc_sens_residual(x, f);

    return ierr;
  }

  //---------------------------------------------------------------
  // this function is called by PETSc to evaluate the Jacobian at X
  static PetscErrorCode  __genius_petsc_snes_jacobian (SNES, Vec x, Mat *jac, Mat *pc, MatStructure *msflag, void *ctx)
  {
    int ierr=0;

    // convert void* to FVM_NonlinearSolver*
    FVM_NonlinearSolver * nonlinear_solver = (FVM_NonlinearSolver *)ctx;

    nonlinear_solver->build_petsc_sens_jacobian(x, jac, pc);

    *msflag = SAME_NONZERO_PATTERN;

    //*msflag = DIFFERENT_NONZERO_PATTERN;

    return ierr;
  }

  //---------------------------------------------------------------
  // this function is called by PETSc to do pre check after each line search
  static PetscErrorCode __genius_petsc_snes_pre_check(SNES , Vec x, Vec y, void *ctx, PetscBool *changed_y)
  {
    PetscErrorCode ierr=0;

    // convert void* to FVM_NonlinearSolver*
    FVM_NonlinearSolver * nonlinear_solver = (FVM_NonlinearSolver *)ctx;

    nonlinear_solver->sens_line_search_pre_check(x, y, changed_y);

    return ierr;
  }

  //---------------------------------------------------------------
  // this function is called by PETSc to do post check after each line search
  static PetscErrorCode __genius_petsc_snes_post_check(SNES , Vec x, Vec y, Vec w, void *ctx, PetscBool *changed_y, PetscBool *changed_w)
  {
    PetscErrorCode ierr=0;

    // convert void* to FVM_NonlinearSolver*
    FVM_NonlinearSolver * nonlinear_solver = (FVM_NonlinearSolver *)ctx;

    nonlinear_solver->sens_line_search_post_check(x, y, w, changed_y, changed_w);

    return ierr;
  }

} // end extern "C"
//---------------------------------------------------------------------



/*------------------------------------------------------------------
 * constructor, setup context
 */
FVM_NonlinearSolver::FVM_NonlinearSolver(SimulationSystem & system): FVM_PDESolver(system)
{
  PetscErrorCode ierr;

  // create petsc nonlinear solver context
  ierr = SNESCreate(PETSC_COMM_WORLD, &snes); genius_assert(!ierr);

}


/*------------------------------------------------------------------
 * setup nonlinear matrix/vector
 */
void FVM_NonlinearSolver::setup_nonlinear_data()
{
  // map mesh to PETSC solver
  build_dof_map();

  // set petsc routine

  PetscErrorCode ierr;

  // create the global solution vector
  ierr = VecCreateMPI(PETSC_COMM_WORLD, n_local_dofs, n_global_dofs, &x); genius_assert(!ierr);
  // use VecDuplicate to create vector with same pattern
  ierr = VecDuplicate(x, &f);genius_assert(!ierr);
  ierr = VecDuplicate(x, &L);genius_assert(!ierr);

  // set all the components of scale vector L to 1.0
  ierr = VecSet(L, 1.0); genius_assert(!ierr);

  // create local vector, which has extra room for ghost dofs! the MPI_COMM here is PETSC_COMM_SELF
  ierr = VecCreateSeq(PETSC_COMM_SELF,  local_index_array.size() , &lx); genius_assert(!ierr);
  // use VecDuplicate to create vector with same pattern
  ierr = VecDuplicate(lx, &lf);genius_assert(!ierr);

  // create the index for vector statter
#ifdef PETSC_VERSION_DEV
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, global_index_array.size(), &global_index_array[0] , PETSC_COPY_VALUES, &gis); genius_assert(!ierr);
  ierr = ISCreateGeneral(PETSC_COMM_SELF,  local_index_array.size(),  &local_index_array[0] ,  PETSC_COPY_VALUES, &lis); genius_assert(!ierr);
#else
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, global_index_array.size(), &global_index_array[0] , &gis); genius_assert(!ierr);
  ierr = ISCreateGeneral(PETSC_COMM_SELF,  local_index_array.size(),  &local_index_array[0] ,  &lis); genius_assert(!ierr);
#endif

  // it seems we can free global_index_array and local_index_array to save the memory

  // create the vector statter
  ierr = VecScatterCreate(x, gis, lx, lis, &scatter); genius_assert(!ierr);



  // create the jacobian matrix
  ierr = MatCreate(PETSC_COMM_WORLD,&J); genius_assert(!ierr);
  ierr = MatSetSizes(J, n_local_dofs, n_local_dofs, n_global_dofs, n_global_dofs); genius_assert(!ierr);


  if (Genius::n_processors()>1)
  {
    ierr = MatSetType(J,MATMPIAIJ); genius_assert(!ierr);
    // alloc memory for parallel matrix here
    ierr = MatMPIAIJSetPreallocation(J, 0, &n_nz[0], 0, &n_oz[0]); genius_assert(!ierr);
  }
  else
  {
    ierr = MatSetType(J,MATSEQAIJ); genius_assert(!ierr);
    // alloc memory for sequence matrix here
    ierr = MatSeqAIJSetPreallocation(J, 0, &n_nz[0]); genius_assert(!ierr);
  }



  // indicates when PetscUtils::MatZeroRows() is called the zeroed entries are kept in the nonzero structure
#if PETSC_VERSION_GE(3,1,0)
  ierr = MatSetOption(J, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); genius_assert(!ierr);
#endif

#if PETSC_VERSION_EQ(3,0,0)
  ierr = MatSetOption(J, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE); genius_assert(!ierr);
#endif


  //ierr = MatSetOption(J, MAT_YES_NEW_NONZERO_LOCATIONS); genius_assert(!ierr);

  // indicates when MatSetValue with ADD_VALUES mode, the 0 entries will be ignored
  // since I uas ADD_VALUES 0 for keeping nonzero pattern, do not use this parameter!
  //ierr = MatSetOption(J, MAT_IGNORE_ZERO_ENTRIES); genius_assert(!ierr);

  // the jacobian matrix is not assembled yet.
  jacobian_matrix_first_assemble = false;

  ierr = MatSetFromOptions(J); genius_assert(!ierr);



  // set petsc nonlinear solver type here
  set_petsc_nonelinear_solver_type( SolverSpecify::NS );

  // set the nonlinear function
  ierr = SNESSetFunction (snes, f, __genius_petsc_snes_residual, this);genius_assert(!ierr);

  // set the nonlinear Jacobian
  ierr = SNESSetJacobian (snes, J, J, __genius_petsc_snes_jacobian, this);genius_assert(!ierr);

  // set nonlinear solver monitor
  ierr = SNESMonitorSet (snes, __genius_petsc_snes_monitor, this, PETSC_NULL); genius_assert(!ierr);

  // set nonlinear solver convergence test routine
  ierr = SNESSetConvergenceTest(snes, __genius_petsc_snes_convergence_test, this, PETSC_NULL); genius_assert(!ierr);


  // get petsc linear solver context
  // Have the Krylov subspace method use our good initial guess rather than 0
  ierr = SNESGetKSP (snes, &ksp); genius_assert(!ierr);

  // set user defined ksy convergence criterion
  ierr = KSPSetConvergenceTest (ksp, __genius_petsc_ksp_convergence_test, this, PETSC_NULL); genius_assert(!ierr);


  // Have the Krylov subspace method use our good initial guess rather than 0
  //ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE); genius_assert(!ierr);
  //ierr = KSPSetInitialGuessKnoll(ksp, PETSC_TRUE); genius_assert(!ierr);

  // get petsc preconditional context
  ierr = KSPGetPC(ksp, &pc); genius_assert(!ierr);

  // Set user-specified linear solver and preconditioner types
  set_petsc_linear_solver_type    ( SolverSpecify::LS );
  set_petsc_preconditioner_type   ( SolverSpecify::PC );



}




/*------------------------------------------------------------------
 * dump jacobian matrix to external file for analysis
 */
void FVM_NonlinearSolver::dump_jacobian_matrix_petsc(const std::string &file) const
{
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, file.c_str(), FILE_MODE_WRITE, &viewer);
  //PetscViewerSetFormat(viewer,PetscViewerFormat format)
  MatView(J, viewer);
  PetscViewerDestroy(viewer);
}


/*------------------------------------------------------------------
 * dump function to external file
 */
void FVM_NonlinearSolver::dump_function_vector_petsc(const std::string &file) const
{
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, file.c_str(), FILE_MODE_WRITE, &viewer);
  VecView(f, viewer);
  PetscViewerDestroy(viewer);
}




/*------------------------------------------------------------------
 * destroy nonlinear data
 */
void FVM_NonlinearSolver::clear_nonlinear_data()
{
  PetscErrorCode ierr;
  // free everything
  ierr = VecDestroy(x);              genius_assert(!ierr);
  ierr = VecDestroy(f);              genius_assert(!ierr);
  ierr = VecDestroy(L);              genius_assert(!ierr);
  ierr = VecDestroy(lx);             genius_assert(!ierr);
  ierr = VecDestroy(lf);             genius_assert(!ierr);
  ierr = ISDestroy(gis);             genius_assert(!ierr);
  ierr = ISDestroy(lis);             genius_assert(!ierr);
  ierr = VecScatterDestroy(scatter); genius_assert(!ierr);
  ierr = MatDestroy(J);              genius_assert(!ierr);
}

/*------------------------------------------------------------------
 * destructor: destroy context
 */
FVM_NonlinearSolver::~FVM_NonlinearSolver()
{
  PetscErrorCode ierr;
  ierr = SNESDestroy(snes);          genius_assert(!ierr);
}


/*------------------------------------------------------------------
 * default snes monitor
 */
void FVM_NonlinearSolver::petsc_snes_monitor(PetscInt its, PetscReal fnorm)
{
  static PetscInt pre_lits;

  if( its==0 ) pre_lits=0;

  // get the iterations of linear solver

  PetscInt lits;
  SNESGetLinearSolveIterations(snes, &lits);

  MESSAGE<< " its " << its << '\t'
  << std::scientific
  << " |residual|_2 = " << fnorm
  << "  linear iter = "  << lits - pre_lits
  << std::endl;
  RECORD();

  pre_lits = lits;

}


/*------------------------------------------------------------------
 * default snes convergence test
 */
void FVM_NonlinearSolver::petsc_snes_convergence_test(PetscInt its, PetscReal xnorm, PetscReal gnorm, PetscReal fnorm, SNESConvergedReason *reason)
{
  SNESDefaultConverged(snes, its, xnorm, gnorm, fnorm, reason, this);
}


/*------------------------------------------------------------------
 * default ksp convergence test
 */
void FVM_NonlinearSolver::petsc_ksp_convergence_test(PetscInt its, PetscReal rnorm, KSPConvergedReason* reason)
{
  KSPDefaultConverged(ksp, its, rnorm, reason, this);
}


/*------------------------------------------------------------------
 * default line search pre check, do nothing
 */
void FVM_NonlinearSolver::sens_line_search_pre_check(Vec , Vec , PetscBool *)
{
  hook_list()->pre_iteration();
  return;
}

/*------------------------------------------------------------------
 * default line search post check, call each post_iteration for the hooks
 */
void FVM_NonlinearSolver::sens_line_search_post_check(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{
  bool _changed_y = false, _changed_w = false;
  hook_list()->post_iteration((void*)f, (void*)x, (void*)y, (void*)w, _changed_y, _changed_w);
  *changed_y = _changed_y ? PETSC_TRUE : *changed_y;
  *changed_w = _changed_w ? PETSC_TRUE : *changed_w;
  return;
}



//--------------------------------------------------------------------------------------------------------------------

void FVM_NonlinearSolver::set_petsc_nonelinear_solver_type(SolverSpecify::NonLinearSolverType nonlinear_solver_type)
{
  int ierr = 0;

  _nonlinear_solver_type = nonlinear_solver_type;

  switch (nonlinear_solver_type)
  {
      case SolverSpecify::Newton:
      ierr = SNESSetType(snes,SNESLS); genius_assert(!ierr);
      ierr = SNESLineSearchSet(snes,SNESLineSearchNo,PETSC_NULL); genius_assert(!ierr);
      // set the LineSearch pre/post check routine
      ierr = SNESLineSearchSetPreCheck(snes, __genius_petsc_snes_pre_check, this); genius_assert(!ierr);
      ierr = SNESLineSearchSetPostCheck(snes, __genius_petsc_snes_post_check, this); genius_assert(!ierr);
      return;

      case SolverSpecify::LineSearch:
      ierr = SNESSetType(snes,SNESLS); genius_assert(!ierr);
      ierr = SNESLineSearchSet(snes,SNESLineSearchCubic,PETSC_NULL); genius_assert(!ierr);
      // set the LineSearch pre/post check routine
      ierr = SNESLineSearchSetPreCheck(snes, __genius_petsc_snes_pre_check, this); genius_assert(!ierr);
      ierr = SNESLineSearchSetPostCheck(snes, __genius_petsc_snes_post_check, this); genius_assert(!ierr);
      return;

      case SolverSpecify::TrustRegion:
      ierr = SNESSetType(snes, SNESTR); genius_assert(!ierr);
      ierr = SNESSetTrustRegionTolerance(snes, 1e-30); genius_assert(!ierr);
      return;

      default :
      std::cerr << "ERROR:  Unsupported PETSC Nonlinear Solver: "
      << nonlinear_solver_type     << std::endl
      << "Continuing with PETSC defaults" << std::endl;
  }

}



void FVM_NonlinearSolver::set_petsc_linear_solver_type(SolverSpecify::LinearSolverType linear_solver_type)
{
  int ierr = 0;

  _linear_solver_type = linear_solver_type;

  switch (linear_solver_type)
  {

      case SolverSpecify::CG:
      MESSAGE<< "Using CG linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPCG);         genius_assert(!ierr); return;

      case SolverSpecify::CR:
      MESSAGE<< "Using CR linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPCR);         genius_assert(!ierr); return;

      case SolverSpecify::CGS:
      MESSAGE<< "Using CGS linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPCGS);        genius_assert(!ierr); return;

      case SolverSpecify::BICG:
      MESSAGE<< "Using BICG linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPBICG);       genius_assert(!ierr); return;

      case SolverSpecify::TCQMR:
      MESSAGE<< "Using TCQMR linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPTCQMR);      genius_assert(!ierr); return;

      case SolverSpecify::TFQMR:
      MESSAGE<< "Using TFQMR linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPTFQMR);      genius_assert(!ierr); return;

      case SolverSpecify::LSQR:
      MESSAGE<< "Using LSQR linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPLSQR);       genius_assert(!ierr); return;

      case SolverSpecify::BICGSTAB:
      MESSAGE<< "Using BCGS linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPBCGS);       genius_assert(!ierr); return;

      case SolverSpecify::BCGSL:
      MESSAGE<< "Using BCGSL linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPBCGSL);      genius_assert(!ierr); return;

      case SolverSpecify::MINRES:
      MESSAGE<< "Using MINRES linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPMINRES);     genius_assert(!ierr); return;

      case SolverSpecify::GMRES:
      MESSAGE<< "Using GMRES linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPGMRES);      genius_assert(!ierr);
      // for GMRES method, we need to enlarge restart step (default is 30)
      ierr = KSPGMRESSetRestart(ksp, 100);            genius_assert(!ierr);
      return;

      case SolverSpecify::FGMRES:
      MESSAGE<< "Using FGMRES linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPFGMRES);      genius_assert(!ierr);
      // for GMRES method, we need to enlarge restart step (default is 30)
      ierr = KSPGMRESSetRestart(ksp, 100);            genius_assert(!ierr);
      return;

      case SolverSpecify::RICHARDSON:
      MESSAGE<< "Using RICHARDSON linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPRICHARDSON); genius_assert(!ierr); return;

      case SolverSpecify::CHEBYSHEV:
      MESSAGE<< "Using CHEBYSHEV linear solver..."<<std::endl;  RECORD();
      ierr = KSPSetType (ksp, (char*) KSPCHEBYCHEV);  genius_assert(!ierr); return;

      case SolverSpecify::LU:
      case SolverSpecify::UMFPACK:
      case SolverSpecify::SuperLU:
      case SolverSpecify::MUMPS:
      case SolverSpecify::PASTIX:
      case SolverSpecify::SuperLU_DIST:
      if (Genius::n_processors()>1)
      {
        switch(linear_solver_type)
        {
            case   SolverSpecify::LU :
            case   SolverSpecify::MUMPS :
            // if LU method is required, we should check if SuperLU_DIST/MUMPS is installed
            // the default parallel LU solver is set to MUMPS
#ifdef PETSC_HAVE_MUMPS
            MESSAGE<< "Using MUMPS linear solver..."<<std::endl;
            RECORD();
            ierr = KSPSetType (ksp, (char*) KSPPREONLY); genius_assert(!ierr);
            ierr = PCSetType  (pc, (char*) PCLU); genius_assert(!ierr);
            ierr = PCFactorSetMatSolverPackage (pc, "mumps"); genius_assert(!ierr);
            ierr = PetscOptionsSetValue("-mat_mumps_icntl_14", "80");  genius_assert(!ierr);
            //ierr = PetscOptionsSetValue("-mat_mumps_icntl_23","4000");
#else
            MESSAGE<< "Warning:  no MUMPS solver configured, use BCGS instead!" << std::endl;
            RECORD();
            ierr = KSPSetType (ksp, (char*) KSPBCGS);       genius_assert(!ierr); return;
#endif
            break;

            case SolverSpecify::PASTIX:
#ifdef PETSC_HAVE_PASTIX
            MESSAGE<< "Using PaStiX linear solver..."<<std::endl;
            RECORD();
            ierr = KSPSetType (ksp, (char*) KSPPREONLY); genius_assert(!ierr);
            ierr = PCSetType  (pc, (char*) PCLU); genius_assert(!ierr);
            ierr = PCFactorSetMatSolverPackage (pc, "pastix"); genius_assert(!ierr);
#else
            MESSAGE<< "Warning:  no PaStiX solver configured, use BCGS instead!" << std::endl;
            RECORD();
            ierr = KSPSetType (ksp, (char*) KSPBCGS);       genius_assert(!ierr); return;
#endif
            break;

            case SolverSpecify::SuperLU_DIST:
#ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
            MESSAGE<< "Using SuperLU_DIST linear solver..."<<std::endl;
            RECORD();
            ierr = KSPSetType (ksp, (char*) KSPPREONLY); genius_assert(!ierr);
            ierr = PCSetType  (pc, (char*) PCLU); genius_assert(!ierr);
            ierr = PCFactorSetMatSolverPackage (pc, "superlu_dist"); genius_assert(!ierr);
#else
            MESSAGE<< "Warning:  no SuperLU_DIST solver configured, use BCGS instead!" << std::endl;
            RECORD();
            ierr = KSPSetType (ksp, (char*) KSPBCGS);       genius_assert(!ierr); return;
#endif
            break;
            default:
            ierr = KSPSetType (ksp, (char*) KSPBCGS);       genius_assert(!ierr); return;
        }
      }
      else
      {
        ierr = KSPSetType (ksp, (char*) KSPPREONLY); genius_assert(!ierr);
        ierr = PCSetType  (pc, (char*) PCLU); genius_assert(!ierr);
        switch (linear_solver_type)
        {
            case SolverSpecify::LU :
            MESSAGE<< "Using LAPACK-LU linear solver..."<<std::endl;  RECORD();
            break;
            case SolverSpecify::UMFPACK :
#ifdef PETSC_HAVE_UMFPACK
            MESSAGE<< "Using UMFPACK linear solver..."<<std::endl;
            RECORD();
            ierr = PCFactorSetMatSolverPackage (pc, "umfpack"); genius_assert(!ierr);
#else
            MESSAGE << "Warning:  no UMFPACK solver configured, use default LU solver instead!" << std::endl;
            RECORD();
#endif
            break;
            case SolverSpecify::SuperLU :
#ifdef PETSC_HAVE_SUPERLU
            MESSAGE<< "Using SuperLU linear solver..."<<std::endl;
            RECORD();
            ierr = PCFactorSetMatSolverPackage (pc, "superlu"); genius_assert(!ierr);
#else
            MESSAGE << "Warning:  no SuperLU solver configured, use default LU solver instead!" << std::endl;
            RECORD();
#endif
            break;
            case SolverSpecify::MUMPS :
#ifdef PETSC_HAVE_MUMPS
            MESSAGE<< "Using MUMPS linear solver..."<<std::endl;
            RECORD();
            ierr = PCFactorSetMatSolverPackage (pc, "mumps"); genius_assert(!ierr);
#else
            MESSAGE << "Warning:  no MUMPS solver configured, use default LU solver instead!" << std::endl;
            RECORD();
#endif
            break;
            case SolverSpecify::PASTIX:
#ifdef PETSC_HAVE_PASTIX
            MESSAGE<< "Using PaStiX linear solver..."<<std::endl;
            RECORD();
            ierr = PCFactorSetMatSolverPackage (pc, "pastix"); genius_assert(!ierr);
#else
            MESSAGE<< "Warning:  no PaStiX solver configured, use default LU solver instead!" << std::endl;
            RECORD();
#endif
            break;
            case   SolverSpecify::SuperLU_DIST:
#ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
            MESSAGE<< "Using SuperLU_DIST linear solver..."<<std::endl;
            RECORD();
            ierr = PCFactorSetMatSolverPackage (pc, "superlu_dist"); genius_assert(!ierr);
#else
            MESSAGE << "Warning:  no SuperLU_DIST solver configured, use default LU solver instead!" << std::endl;
            RECORD();
#endif
            break;

            default:
            // should never reach here
            genius_error();
        }
      }


      ierr = PCFactorSetReuseFill(pc, PETSC_TRUE);genius_assert(!ierr);
      ierr = PCFactorSetReuseOrdering(pc, PETSC_TRUE); genius_assert(!ierr);
      // prevent zero pivot in LU factorization
#if PETSC_VERSION_LE(3,0,0)
      ierr = PCFactorSetPivoting(pc, 1.0); genius_assert(!ierr);
      //ierr = PCFactorReorderForNonzeroDiagonal(pc, 1e-20); genius_assert(!ierr);<-- Caught signal number 11 SEGV error will occure when diag value < 1e-20
      ierr = PCFactorSetShiftNonzero(pc, 1e-20); genius_assert(!ierr);
#endif

#if PETSC_VERSION_GE(3,1,0)
      ierr = PCFactorSetColumnPivot(pc, 1.0); genius_assert(!ierr);
      //ierr = PCFactorReorderForNonzeroDiagonal(pc, 1e-20); genius_assert(!ierr);<-- Caught signal number 11 SEGV error will occure when diag value < 1e-20
      ierr = PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO);genius_assert(!ierr);
#endif
      return;

      default:
      std::cerr << "ERROR:  Unsupported PETSC Solver: "
      << linear_solver_type        << std::endl
      << "Continuing with PETSC defaults" << std::endl;
  }
}



void FVM_NonlinearSolver::set_petsc_preconditioner_type(SolverSpecify::PreconditionerType preconditioner_type)
{
  int ierr = 0;

  // skip direct methods
  if (_linear_solver_type == SolverSpecify::LU ||
      _linear_solver_type == SolverSpecify::UMFPACK ||
      _linear_solver_type == SolverSpecify::SuperLU ||
      _linear_solver_type == SolverSpecify::MUMPS   ||
      _linear_solver_type == SolverSpecify::PASTIX  ||
      _linear_solver_type == SolverSpecify::SuperLU_DIST
     )
  {
    ierr = KSPSetType (ksp, (char*) KSPPREONLY);
    genius_assert(!ierr);
    return;
  }

  _preconditioner_type = preconditioner_type;

  switch (preconditioner_type)
  {
      case SolverSpecify::IDENTITY_PRECOND:
      ierr = PCSetType (pc, (char*) PCNONE);      genius_assert(!ierr); return;

      case SolverSpecify::CHOLESKY_PRECOND:
      ierr = PCSetType (pc, (char*) PCCHOLESKY);  genius_assert(!ierr); return;

      case SolverSpecify::ICC_PRECOND:
      ierr = PCSetType (pc, (char*) PCICC);       genius_assert(!ierr); return;

      case SolverSpecify::ILU_PRECOND:
      if (Genius::n_processors()>1)
      {
#ifdef PETSC_HAVE_LIBHYPRE
        MESSAGE<< "Using Hypre/Euclid ILU preconditioner..."<<std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCHYPRE);     genius_assert(!ierr);
        ierr = PCHYPRESetType (pc, "euclid");       genius_assert(!ierr);
        return;
#else
        MESSAGE << "Warning:  no parallel ILU preconditioner configured, use ASM instead!" << std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCASM);       genius_assert(!ierr);
        return;
#endif

      }
      else
      {
#ifdef PETSC_HAVE_LIBHYPRE
        MESSAGE<< "Using Hypre/Euclid ILU preconditioner..."<<std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCHYPRE);     genius_assert(!ierr);
        ierr = PCHYPRESetType (pc, "euclid");       genius_assert(!ierr);
        return;
#else
        ierr = PCSetType (pc, (char*) PCILU);            genius_assert(!ierr);
        ierr = PCFactorSetReuseFill(pc, PETSC_TRUE);     genius_assert(!ierr);
        ierr = PCFactorSetReuseOrdering(pc, PETSC_TRUE); genius_assert(!ierr);
#if PETSC_VERSION_LE(3,0,0)
        ierr = PCFactorSetPivoting(pc, 1.0); genius_assert(!ierr);
        ierr = PCFactorSetShiftNonzero(pc, 1e-12);  genius_assert(!ierr);
#endif
#if PETSC_VERSION_GE(3,1,0)
        ierr = PCFactorSetColumnPivot(pc, 1.0); genius_assert(!ierr);
        ierr = PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO);genius_assert(!ierr);
#endif
        //ierr = PCFactorSetMatOrderingType(pc, MATORDERING_ND); genius_assert(!ierr);
        //ierr = PCFactorSetMatOrderingType(pc, MATORDERING_RCM);
        //ierr = PetscOptionsSetValue("-pc_factor_nonzeros_along_diagonal", 0); genius_assert(!ierr);
        ierr = PCFactorSetAllowDiagonalFill(pc);genius_assert(!ierr);
        //ierr = PetscOptionsSetValue("-pc_factor_diagonal_fill", 0);
        return;
#endif

      }

      case SolverSpecify::ILUT_PRECOND:
      if (Genius::n_processors()==1)
      {
#ifdef PETSC_HAVE_SUPERLU
        MESSAGE<< "Using SuperLU ILUT preconditioner..."<<std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCILU);     genius_assert(!ierr);
        ierr = PCFactorSetMatSolverPackage (pc, "superlu"); genius_assert(!ierr);
        //ierr = PetscOptionsSetValue("-mat_superlu_ilu_droptol","1e-7"); genius_assert(!ierr);
        //ierr = PetscOptionsSetValue("-mat_superlu_ilu_filltol","1e-2"); genius_assert(!ierr);
        //ierr = PetscOptionsSetValue("-mat_superlu_ilu_fillfactor","200"); genius_assert(!ierr);
        return;
#else
        MESSAGE << "Warning:  no ILUT preconditioner configured, use ASM instead!" << std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCASM);       genius_assert(!ierr);
        return;
#endif

      }
      else
      {
        MESSAGE << "Warning:  no parallel ILUT preconditioner configured, use ASM instead!" << std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCASM);       genius_assert(!ierr);
        return;
      }

      // some times, we still need a LU solver as strong preconditioner
      case SolverSpecify::LU_PRECOND :
      {
        if (Genius::n_processors()==1)
        {
          ierr = PCSetType (pc, (char*) PCLU);       genius_assert(!ierr);
#ifdef PETSC_HAVE_MUMPS
          MESSAGE<< "Using MUMPS as LU preconditioner..."<<std::endl;    RECORD();
          ierr = PCFactorSetMatSolverPackage (pc, "mumps"); genius_assert(!ierr);
#endif
          ierr = PCFactorSetReuseFill(pc, PETSC_TRUE);genius_assert(!ierr);
          ierr = PCFactorSetReuseOrdering(pc, PETSC_TRUE); genius_assert(!ierr);
#if PETSC_VERSION_LE(3,0,0)
          ierr = PCFactorSetPivoting(pc, 1.0); genius_assert(!ierr);
          ierr = PCFactorSetShiftNonzero(pc, 1e-12);  genius_assert(!ierr);
#endif
#if PETSC_VERSION_GE(3,1,0)
          ierr = PCFactorSetColumnPivot(pc, 1.0); genius_assert(!ierr);
          ierr = PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);genius_assert(!ierr);
#endif
          return;
        }
        else
        {
#ifdef PETSC_HAVE_MUMPS
          ierr = PCSetType (pc, (char*) PCLU);       genius_assert(!ierr);
          MESSAGE<< "Using MUMPS as parallel LU preconditioner..."<<std::endl;    RECORD();
          ierr = PCFactorSetMatSolverPackage (pc, "mumps"); genius_assert(!ierr);
          ierr = PCFactorSetReuseFill(pc, PETSC_TRUE);genius_assert(!ierr);
          ierr = PCFactorSetReuseOrdering(pc, PETSC_TRUE); genius_assert(!ierr);
#if PETSC_VERSION_LE(3,0,0)
          ierr = PCFactorSetPivoting(pc, 1.0); genius_assert(!ierr);
          ierr = PCFactorSetShiftNonzero(pc, 1e-12);  genius_assert(!ierr);
#endif
  #if PETSC_VERSION_GE(3,1,0)
          ierr = PCFactorSetColumnPivot(pc, 1.0); genius_assert(!ierr);
          ierr = PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO);genius_assert(!ierr);
#endif
#else
          MESSAGE << "Warning:  no parallel LU preconditioner configured, use ASM instead!" << std::endl;
          RECORD();
          ierr = PCSetType (pc, (char*) PCASM);       genius_assert(!ierr);
#endif
          return;
        }
      }
      case SolverSpecify::ASM_PRECOND:
      {
        if (Genius::n_processors() > 1)
        {
          ierr = PCSetType (pc, (char*) PCASM);       genius_assert(!ierr);
          //#ifdef PETSC_HAVE_MUMPS
          //ierr = SNESSetLagPreconditioner(snes, 2); genius_assert(!ierr);
          //ierr = PetscOptionsSetValue("-sub_ksp_type","preonly"); genius_assert(!ierr);
          //ierr = PetscOptionsSetValue("-sub_pc_type","lu"); genius_assert(!ierr);
          //ierr = PetscOptionsSetValue("-sub_pc_factor_mat_solver_package","mumps"); genius_assert(!ierr);
          //#else
          ierr = PetscOptionsSetValue("-sub_pc_type","ilu"); genius_assert(!ierr);
          //#endif
          ierr = PetscOptionsSetValue("-sub_pc_factor_shift_nonzero","1e-12"); genius_assert(!ierr);
          ierr = PCFactorSetReuseFill(pc, PETSC_TRUE);genius_assert(!ierr);
          ierr = PCFactorSetReuseOrdering(pc, PETSC_TRUE); genius_assert(!ierr);
        }
        else
        {
          ierr = PCSetType (pc, (char*) PCILU);       genius_assert(!ierr);
          ierr = PCFactorSetReuseFill(pc, PETSC_TRUE);genius_assert(!ierr);
          ierr = PCFactorSetReuseOrdering(pc, PETSC_TRUE); genius_assert(!ierr);
#if PETSC_VERSION_LE(3,0,0)
          ierr = PCFactorSetPivoting(pc, 1.0); genius_assert(!ierr);
          ierr = PCFactorSetShiftNonzero(pc, 1e-12);  genius_assert(!ierr);
#endif
#if PETSC_VERSION_GE(3,1,0)
          ierr = PCFactorSetColumnPivot(pc, 1.0); genius_assert(!ierr);
          ierr = PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO);genius_assert(!ierr);
#endif
          //ierr = PetscOptionsSetValue("-pc_factor_nonzeros_along_diagonal", 0); genius_assert(!ierr);
          ierr = PetscOptionsSetValue("-pc_factor_diagonal_fill", 0); genius_assert(!ierr);
        }
        return;
      }

      case SolverSpecify::PARMS_PRECOND:
      {
#ifdef PETSC_HAVE_PARMS
        MESSAGE<< "Using pARMS preconditioner..."<<std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCPARMS);     genius_assert(!ierr);
        return;
#else
        MESSAGE << "Warning:  no pARMS preconditioner configured, use ASM instead!" << std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCASM);       genius_assert(!ierr);
        return;
#endif

      }

      case SolverSpecify::BOOMERAMG_PRECOND:
      {
#ifdef PETSC_HAVE_LIBHYPRE
        MESSAGE<< "Using Hypre/BoomerAMG preconditioner..."<<std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCHYPRE);     genius_assert(!ierr);
        ierr = PCHYPRESetType (pc, "boomeramg");    genius_assert(!ierr);
        return;
#else
        MESSAGE << "Warning:  no Hypre/BoomerAMG preconditioner configured, use ASM instead!" << std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCASM);       genius_assert(!ierr);
        return;
#endif

      }

      case SolverSpecify::JACOBI_PRECOND:
      ierr = PCSetType (pc, (char*) PCJACOBI);    genius_assert(!ierr); return;

      case SolverSpecify::BLOCK_JACOBI_PRECOND:
      ierr = PCSetType (pc, (char*) PCBJACOBI);   genius_assert(!ierr); return;

      case SolverSpecify::SOR_PRECOND:
      ierr = PCSetType (pc, (char*) PCSOR);       genius_assert(!ierr); return;

      case SolverSpecify::EISENSTAT_PRECOND:
      ierr = PCSetType (pc, (char*) PCEISENSTAT); genius_assert(!ierr); return;


      case SolverSpecify::USER_PRECOND:
      ierr = PCSetType (pc, (char*) PCMAT);       genius_assert(!ierr); return;


      case SolverSpecify::SHELL_PRECOND:
      ierr = PCSetType (pc, (char*) PCSHELL);     genius_assert(!ierr); return;

      default:
      std::cerr
      << "ERROR:  Unsupported PETSC Preconditioner: "
      << preconditioner_type       << std::endl
      << "Continuing with PETSC defaults" << std::endl;
  }
}

#ifndef PETSC_VERSION_DEV
#define SNES_DIVERGED_LINE_SEARCH SNES_DIVERGED_LS_FAILURE
#endif
void FVM_NonlinearSolver::sens_solve()
{
  // do snes solve
  SNESSolve ( snes, PETSC_NULL, x );

  // get the converged reason
  SNESConvergedReason reason;
  SNESGetConvergedReason ( snes,&reason );

  // if Line search failed, disable Line search
  if ( reason == SNES_DIVERGED_LINE_SEARCH || reason == SNES_DIVERGED_LOCAL_MIN )
  {
    MESSAGE <<"------> nonlinear solver " << SNESConvergedReasons[reason] <<". Disable Line Search.\n\n\n";
    RECORD();
    SNESLineSearchSet ( snes,SNESLineSearchNo,PETSC_NULL );
    SNESSolve ( snes, PETSC_NULL, x );
  }

}



double FVM_NonlinearSolver::condition_number_of_jacobian_matrix()
{
#ifdef HAVE_SLEPC

  // SVD solver for largest singular value
  SVD            svd_l;

  // SVD solver for smallest singular value
  SVD            svd_s;
  EPS            eps_s;
  ST             st_s;
  KSP            ksp_s;
  PC             pc_s;

  PetscErrorCode ierr;

  // Create singular value solver context
  ierr = SVDCreate(PETSC_COMM_WORLD, &svd_s);
  ierr = SVDCreate(PETSC_COMM_WORLD, &svd_l);
  // Set operator
  ierr = SVDSetOperator(svd_s, J);
  ierr = SVDSetOperator(svd_l, J);

  // thick-restart Lanczos as default solver for large singular value
  ierr = SVDSetWhichSingularTriplets(svd_l, SVD_LARGEST);
  ierr = SVDSetType(svd_l, SVDTRLANCZOS);

  // small singular value use eigen value solver on Cyclic Matrix
  ierr = SVDSetWhichSingularTriplets(svd_s, SVD_SMALLEST);
  ierr = SVDSetType(svd_s, SVDCYCLIC);
  ierr = SVDCyclicSetExplicitMatrix(svd_s, PETSC_TRUE);
  // shift-and-invert spectral transformation to enhance convergence of eigenvalues near zero
  ierr = SVDCyclicGetEPS(svd_s, &eps_s);
  ierr = EPSSetType(eps_s, EPSKRYLOVSCHUR);
  ierr = EPSGetST(eps_s, &st_s);
  ierr = STSetType(st_s, STSINVERT);

  ierr = STGetKSP(st_s, &ksp_s);
  ierr = KSPGetPC(ksp_s, &pc_s);
  // since we have to deal with bad conditioned problem, we choose direct solver whenever possible
  if (Genius::n_processors() == 1)
  {
    // direct solver as preconditioner
    ierr = KSPSetType (ksp_s, (char*) KSPGMRES); assert(!ierr);
    // superlu which use static pivot seems very stable
#ifdef PETSC_HAVE_SUPERLU
    ierr = PCSetType  (pc_s, (char*) PCLU); assert(!ierr);
    ierr = PCFactorSetMatSolverPackage (pc_s, "superlu"); assert(!ierr);
#else
    ierr = PCSetType  (pc_s, (char*) PCLU); assert(!ierr);
    ierr = PCFactorSetShiftType(pc_s,MAT_SHIFT_NONZERO); assert(!ierr);
#endif
  }
  else
  {
    // direct solver as preconditioner
    ierr = KSPSetType (ksp_s, (char*) KSPGMRES); assert(!ierr);
    // superlu_dist which also use static pivot (right?). anyway, mumps use partial pivot does not work
#ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
    ierr = PCSetType  (pc_s, (char*) PCLU); assert(!ierr);
    ierr = PCFactorSetMatSolverPackage (pc_s, "superlu_dist"); assert(!ierr);
#else
    // no direct solver? use ASM. seems never work
    ierr = PCSetType (pc_s, (char*) PCASM);       assert(!ierr);
#endif
  }

  // Set solver parameters at runtime
  ierr = SVDSetFromOptions(svd_s);  assert(!ierr);
  ierr = SVDSetFromOptions(svd_l);  assert(!ierr);

  ierr = SVDSetUp(svd_s); assert(!ierr);
  ierr = SVDSetUp(svd_l); assert(!ierr);

  PetscReal sigma_large=1, sigma_small=1;
  PetscInt nconv;
  PetscReal error;

  // find the largest singular value

  SVDSolve(svd_l);
  SVDGetConverged(svd_l, &nconv);  assert(nconv>0);
  SVDGetSingularTriplet(svd_l, 0, &sigma_large, PETSC_NULL, PETSC_NULL);
  SVDComputeRelativeError(svd_l, 0, &error);

  //std::cout<<sigma_large << " " << error << std::endl;
  // find the smallest singular value

  SVDSolve(svd_s);
  SVDGetConverged(svd_s, &nconv);  assert(nconv>0);
  SVDGetSingularTriplet(svd_s, 0, &sigma_small, PETSC_NULL, PETSC_NULL);
  SVDComputeRelativeError(svd_s, 0, &error);

  //std::cout<<sigma_small << " " << error << std::endl;

  ierr = SVDDestroy(svd_s);            assert(!ierr);
  ierr = SVDDestroy(svd_l);            assert(!ierr);

  return sigma_large/sigma_small;
#else
  return 1.0;
#endif
}



void FVM_NonlinearSolver::eigen_value_of_jacobian_matrix()
{
#ifdef HAVE_SLEPC

  EPS            eps_s;
  EPS            eps_l;

  // create eigen value problem solver
  EPSCreate(PETSC_COMM_WORLD, &eps_s);
  EPSCreate(PETSC_COMM_WORLD, &eps_l);
  // Set operator
  EPSSetOperators(eps_s, J, PETSC_NULL);
  EPSSetOperators(eps_l, J, PETSC_NULL);

  EPSSetWhichEigenpairs(eps_s, EPS_SMALLEST_MAGNITUDE);
  EPSSetWhichEigenpairs(eps_l, EPS_LARGEST_MAGNITUDE);

  //EPSSetType(eps_l, EPSARPACK);

  // shift-and-invert spectral transformation to enhance convergence of eigenvalues near zero
  ST st_s;
  EPSGetST(eps_s, &st_s);
  STSetType(st_s, STSINVERT);

  // Set solver parameters at runtime
  EPSSetFromOptions(eps_s);
  EPSSetFromOptions(eps_l);

  PetscScalar kr, ki; /* eigenvalue, k */
  PetscInt nconv;
  PetscReal error;

  EPSSolve( eps_s );
  EPSGetConverged( eps_s, &nconv );
  genius_assert(nconv>0);
  EPSGetEigenvalue( eps_s, 0, &kr, &ki );
  EPSComputeRelativeError( eps_s, 0, &error );
  std::cout<<std::complex<PetscScalar>(kr, ki) << " " << error << std::endl;

  EPSSolve( eps_l );
  EPSGetConverged( eps_l, &nconv );
  genius_assert(nconv>0);
  EPSGetEigenvalue( eps_l, 0, &kr, &ki );
  EPSComputeRelativeError( eps_l, 0, &error );
  std::cout<<std::complex<PetscScalar>(kr, ki) << " " << error << std::endl;

  EPSDestroy(eps_s);
  EPSDestroy(eps_l);

#endif
}


