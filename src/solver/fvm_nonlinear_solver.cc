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

    //*msflag = SAME_NONZERO_PATTERN;

    *msflag = DIFFERENT_NONZERO_PATTERN;

    return ierr;
  }

  //---------------------------------------------------------------
  // this function is called by PETSc to do pre check after each line search
  static PetscErrorCode __genius_petsc_snes_pre_check(SNES , Vec x, Vec y, void *ctx, PetscTruth *changed_y)
  {
    PetscErrorCode ierr=0;

    // convert void* to FVM_NonlinearSolver*
    FVM_NonlinearSolver * nonlinear_solver = (FVM_NonlinearSolver *)ctx;

    nonlinear_solver->sens_line_search_pre_check(x, y, changed_y);

    return ierr;
  }

  //---------------------------------------------------------------
  // this function is called by PETSc to do post check after each line search
  static PetscErrorCode __genius_petsc_snes_post_check(SNES , Vec x, Vec y, Vec w, void *ctx, PetscTruth *changed_y, PetscTruth *changed_w)
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
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, global_index_array.size() , &global_index_array[0] , &gis); genius_assert(!ierr);
  ierr = ISCreateGeneral(PETSC_COMM_SELF,  local_index_array.size() ,  &local_index_array[0] ,  &lis); genius_assert(!ierr);

  // it seems we can free global_index_array and local_index_array to save the memory

  // create the vector statter
  ierr = VecScatterCreate(x, gis, lx, lis, &scatter); genius_assert(!ierr);



  // create the jacobian matrix
  ierr = MatCreate(PETSC_COMM_WORLD,&J); genius_assert(!ierr);
  ierr = MatSetSizes(J, n_local_dofs, n_local_dofs, n_global_dofs, n_global_dofs); genius_assert(!ierr);

#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
   // we are using petsc-3 or petsc-devel
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
#endif


#if PETSC_VERSION_LE(2,3,3)
  if (Genius::n_processors()>1)
  {
    // if LU method is required, we should check if SuperLU_DIST/MUMPS is installed
    // the deault parallel LU solver is set to MUMPS
    switch(SolverSpecify::LS)
    {
    case   SolverSpecify::LU :
    case   SolverSpecify::MUMPS :
      {
#ifdef PETSC_HAVE_MUMPS
        MESSAGE<< "Using MUMPS linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(J,MATAIJMUMPS); genius_assert(!ierr);
        ierr = PetscOptionsSetValue("-mat_mumps_icntl_14","80");
        //ierr = PetscOptionsSetValue("-mat_mumps_icntl_23","4000");
#else
        MESSAGE<< "Warning:  no MUMPS solver configured, use BCGS instead!" << std::endl;
        RECORD();
        ierr = MatSetType(J,MATMPIAIJ); genius_assert(!ierr);
        SolverSpecify::LS = SolverSpecify::BICGSTAB;
#endif
        break;
      }
    case   SolverSpecify::SuperLU_DIST:
      {
#ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
        MESSAGE<< "Using SuperLU_DIST linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(J,MATSUPERLU_DIST); genius_assert(!ierr);
        ierr = PetscOptionsSetValue("-mat_superlu_dist_rowperm","NATURAL");
#else
        MESSAGE << "Warning:  no SuperLU_DIST solver configured, use BCGS instead!" << std::endl;
        RECORD();
        ierr = MatSetType(J,MATMPIAIJ); genius_assert(!ierr);
        SolverSpecify::LS = SolverSpecify::BICGSTAB;
#endif
        break;
      }

      // no SuperLU_DIST/MUMPS? we return to KSP method
    default:
      ierr = MatSetType(J,MATMPIAIJ); genius_assert(!ierr);
    }
    // alloc memory for parallel matrix here
    ierr = MatMPIAIJSetPreallocation(J, 0, &n_nz[0], 0, &n_oz[0]); genius_assert(!ierr);
  }
  else
  {
    switch (SolverSpecify::LS)
    {
      // check if we can use UMFPACK
    case SolverSpecify::UMFPACK  :
      {
#ifdef PETSC_HAVE_UMFPACK
        MESSAGE<< "Using UMFPACK linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(J,MATUMFPACK); genius_assert(!ierr);
#else
        MESSAGE << "Warning:  no UMFPACK solver configured, use default LU solver instead!" << std::endl;
        RECORD();
        ierr = MatSetType(J,MATSEQAIJ); genius_assert(!ierr);
#endif
        break;
      }

      // check if we can use SuperLU
    case SolverSpecify::SuperLU :
      {
#ifdef PETSC_HAVE_SUPERLU
        MESSAGE<< "Using SuperLU linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(J,MATSUPERLU);genius_assert(!ierr);
#else
        MESSAGE << "Warning:  no SuperLU solver configured, use default LU solver instead!" << std::endl;
        RECORD();
        ierr = MatSetType(J,MATSEQAIJ); genius_assert(!ierr);
#endif
        break;
      }
      //use MUMPS for serial problem
    case   SolverSpecify::MUMPS :
      {
#ifdef PETSC_HAVE_MUMPS
        MESSAGE<< "Using MUMPS linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(J,MATAIJMUMPS); genius_assert(!ierr);
#else
        MESSAGE << "Warning:  no MUMPS solver configured, use default LU solver instead!" << std::endl;
        RECORD();
        ierr = MatSetType(J,MATSEQAIJ); genius_assert(!ierr);
#endif
        break;
      }
      //use SuperLU_DIST for serial problem
    case   SolverSpecify::SuperLU_DIST:
      {
#ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
        MESSAGE<< "Using SuperLU_DIST linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(J,MATSUPERLU_DIST); genius_assert(!ierr);
        ierr = PetscOptionsSetValue("-mat_superlu_dist_rowperm","NATURAL"); genius_assert(!ierr);
#else
        MESSAGE << "Warning:  no SuperLU_DIST solver configured, use default LU solver instead!" << std::endl;
        RECORD();
        ierr = MatSetType(J,MATSEQAIJ); genius_assert(!ierr);
#endif
        break;
      }
      // we have to use the default Mat type MATSEQAIJ
    default:
      ierr = MatSetType(J,MATSEQAIJ); genius_assert(!ierr);

    }
    // alloc memory for sequence matrix here
    ierr = MatSeqAIJSetPreallocation(J, 0, &n_nz[0]); genius_assert(!ierr);

  }
#endif

  // indicates when MatZeroRows() is called the zeroed entries are kept in the nonzero structure
#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
  ierr = MatSetOption(J, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE); genius_assert(!ierr);
#endif

#if PETSC_VERSION_LE(2,3,3)
  ierr = MatSetOption(J, MAT_KEEP_ZEROED_ROWS); genius_assert(!ierr);
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
#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
  ierr = SNESSetConvergenceTest(snes, __genius_petsc_snes_convergence_test, this, PETSC_NULL); genius_assert(!ierr);
#endif

#if PETSC_VERSION_LE(2,3,3)
  ierr = SNESSetConvergenceTest(snes, __genius_petsc_snes_convergence_test, this); genius_assert(!ierr);
#endif

  // get petsc linear solver context
  // Have the Krylov subspace method use our good initial guess rather than 0
  ierr = SNESGetKSP (snes, &ksp); genius_assert(!ierr);

  // set user defined ksy convergence criterion
#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
  ierr = KSPSetConvergenceTest (ksp, __genius_petsc_ksp_convergence_test, this, PETSC_NULL); genius_assert(!ierr);
#endif

#if PETSC_VERSION_LE(2,3,3)
  ierr = KSPSetConvergenceTest (ksp, __genius_petsc_ksp_convergence_test, this); genius_assert(!ierr);
#endif

  // Have the Krylov subspace method use our good initial guess rather than 0
  //ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  // get petsc preconditional context
  ierr = KSPGetPC(ksp, &pc); genius_assert(!ierr);

  // Set user-specified linear solver and preconditioner types
  set_petsc_linear_solver_type    ( SolverSpecify::LS );
  set_petsc_preconditioner_type   ( SolverSpecify::PC );

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
void FVM_NonlinearSolver::sens_line_search_pre_check(Vec , Vec , PetscTruth *)
{
  return;
}

/*------------------------------------------------------------------
 * default line search post check, do nothing
 */
void FVM_NonlinearSolver::sens_line_search_post_check(Vec , Vec , Vec , PetscTruth *, PetscTruth *)
{
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
    ierr = KSPSetType (ksp, (char*) KSPCG);         genius_assert(!ierr); return;

  case SolverSpecify::CR:
    ierr = KSPSetType (ksp, (char*) KSPCR);         genius_assert(!ierr); return;

  case SolverSpecify::CGS:
    ierr = KSPSetType (ksp, (char*) KSPCGS);        genius_assert(!ierr); return;

  case SolverSpecify::BICG:
    ierr = KSPSetType (ksp, (char*) KSPBICG);       genius_assert(!ierr); return;

  case SolverSpecify::TCQMR:
    ierr = KSPSetType (ksp, (char*) KSPTCQMR);      genius_assert(!ierr); return;

  case SolverSpecify::TFQMR:
    ierr = KSPSetType (ksp, (char*) KSPTFQMR);      genius_assert(!ierr); return;

  case SolverSpecify::LSQR:
    ierr = KSPSetType (ksp, (char*) KSPLSQR);       genius_assert(!ierr); return;

  case SolverSpecify::BICGSTAB:
    ierr = KSPSetType (ksp, (char*) KSPBCGS);       genius_assert(!ierr); return;

  case SolverSpecify::MINRES:
    ierr = KSPSetType (ksp, (char*) KSPMINRES);     genius_assert(!ierr); return;

  case SolverSpecify::GMRES:
    ierr = KSPSetType (ksp, (char*) KSPGMRES);      genius_assert(!ierr);
    // for GMRES method, we need to enlarge restart step (default is 30)
    ierr = KSPGMRESSetRestart(ksp, 150);            genius_assert(!ierr);
    return;

  case SolverSpecify::RICHARDSON:
    ierr = KSPSetType (ksp, (char*) KSPRICHARDSON); genius_assert(!ierr); return;

  case SolverSpecify::CHEBYSHEV:
    ierr = KSPSetType (ksp, (char*) KSPCHEBYCHEV);  genius_assert(!ierr); return;

  case SolverSpecify::LU:
  case SolverSpecify::UMFPACK:
  case SolverSpecify::SuperLU:
  case SolverSpecify::MUMPS:
  case SolverSpecify::PASTIX:
  case SolverSpecify::SuperLU_DIST:
#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
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
          ierr = PCFactorSetMatSolverPackage (pc, MAT_SOLVER_MUMPS); genius_assert(!ierr);
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
      ierr = PCFactorSetMatSolverPackage (pc, MAT_SOLVER_PASTIX); genius_assert(!ierr);
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
          ierr = PCFactorSetMatSolverPackage (pc, MAT_SOLVER_SUPERLU_DIST); genius_assert(!ierr);
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
        	break;
        case SolverSpecify::UMFPACK :
#ifdef PETSC_HAVE_UMFPACK
          MESSAGE<< "Using UMFPACK linear solver..."<<std::endl;
          RECORD();
          ierr = PCFactorSetMatSolverPackage (pc, MAT_SOLVER_UMFPACK); genius_assert(!ierr);
#else
          MESSAGE << "Warning:  no UMFPACK solver configured, use default LU solver instead!" << std::endl;
          RECORD();
#endif
          break;
        case SolverSpecify::SuperLU :
#ifdef PETSC_HAVE_SUPERLU
          MESSAGE<< "Using SuperLU linear solver..."<<std::endl;
          RECORD();
          ierr = PCFactorSetMatSolverPackage (pc, MAT_SOLVER_SUPERLU); genius_assert(!ierr);
#else
          MESSAGE << "Warning:  no SuperLU solver configured, use default LU solver instead!" << std::endl;
          RECORD();
#endif
          break;
        case SolverSpecify::MUMPS :
#ifdef PETSC_HAVE_MUMPS
          MESSAGE<< "Using MUMPS linear solver..."<<std::endl;
          RECORD();
          ierr = PCFactorSetMatSolverPackage (pc, MAT_SOLVER_MUMPS); genius_assert(!ierr);
#else
          MESSAGE << "Warning:  no MUMPS solver configured, use default LU solver instead!" << std::endl;
          RECORD();
#endif
          break;
        case SolverSpecify::PASTIX:
#ifdef PETSC_HAVE_PASTIX
          MESSAGE<< "Using PaStiX linear solver..."<<std::endl;
          RECORD();
      ierr = KSPSetType (ksp, (char*) KSPPREONLY); genius_assert(!ierr);
      ierr = PCSetType  (pc, (char*) PCLU); genius_assert(!ierr);
      ierr = PCFactorSetMatSolverPackage (pc, MAT_SOLVER_PASTIX); genius_assert(!ierr);
#else
          MESSAGE<< "Warning:  no PaStiX solver configured, use BCGS instead!" << std::endl;
          RECORD();
      ierr = KSPSetType (ksp, (char*) KSPBCGS);       genius_assert(!ierr); return;
#endif
          break;
        case   SolverSpecify::SuperLU_DIST:
#ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
          MESSAGE<< "Using SuperLU_DIST linear solver..."<<std::endl;
          RECORD();
          ierr = PCFactorSetMatSolverPackage (pc, MAT_SOLVER_SUPERLU_DIST); genius_assert(!ierr);
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
#endif

#if PETSC_VERSION_LE(2,3,3)
    ierr = KSPSetType (ksp, (char*) KSPPREONLY); genius_assert(!ierr);
    ierr = PCSetType  (pc, (char*) PCLU); genius_assert(!ierr);
#endif
    // prevent zero pivot in LU factorization
    ierr = PCFactorSetReuseOrdering(pc, PETSC_TRUE); genius_assert(!ierr);
    ierr = PCFactorSetPivoting(pc, 1.0); genius_assert(!ierr);
    //ierr = PCFactorReorderForNonzeroDiagonal(pc, 1e-20); genius_assert(!ierr); //<-- Caught signal number 11 SEGV error will occure when diag value < 1e-20
    ierr = PCFactorSetShiftNonzero(pc, 1e-20); genius_assert(!ierr);
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
      _linear_solver_type == SolverSpecify::SuperLU_DIST
     )
    return;

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
      MESSAGE<< "Using Hypre/Euclid preconditioner..."<<std::endl;
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
      ierr = PCSetType (pc, (char*) PCILU);       genius_assert(!ierr);
      ierr = PCFactorSetReuseFill(pc, PETSC_TRUE);genius_assert(!ierr);
      ierr = PCFactorSetPivoting(pc, 1.0); genius_assert(!ierr);
      ierr = PCFactorReorderForNonzeroDiagonal(pc, 1e-20); genius_assert(!ierr);
      ierr = PCFactorSetShiftNonzero(pc, 1e-12);  genius_assert(!ierr);
      //ierr = PetscOptionsSetValue("-pc_factor_nonzeros_along_diagonal", 0); genius_assert(!ierr);
      ierr = PetscOptionsSetValue("-pc_factor_diagonal_fill", 0); genius_assert(!ierr);
      return;
    }
    // some times, we still need a LU solver as strong preconditioner
  case SolverSpecify::LU_PRECOND :
    {
      if (Genius::n_processors()==1)
      {
        ierr = PCSetType (pc, (char*) PCLU);       genius_assert(!ierr);
        ierr = PCFactorSetReuseFill(pc, PETSC_TRUE);genius_assert(!ierr);
        ierr = PCFactorSetShiftNonzero(pc, 1e-12);  genius_assert(!ierr);
        return;
      }
      else
      {
        MESSAGE << "Warning:  no parallel LU preconditioner configured, use ASM instead!" << std::endl;
        RECORD();
        ierr = PCSetType (pc, (char*) PCASM);       genius_assert(!ierr);
        return;
      }
    }
  case SolverSpecify::ASM_PRECOND:
    if (Genius::n_processors() > 1)
    {
      ierr = PCSetType (pc, (char*) PCASM);       genius_assert(!ierr);
      ierr = PetscOptionsSetValue("-sub_pc_type","ilu"); genius_assert(!ierr);
      ierr = PetscOptionsSetValue("-sub_pc_factor_reuse_fill","1"); genius_assert(!ierr);
      ierr = PetscOptionsSetValue("-sub_pc_factor_shift_nonzero","1e-12"); genius_assert(!ierr);
    }
    else
    {
      ierr = PCSetType (pc, (char*) PCILU);       genius_assert(!ierr);
      ierr = PCFactorSetReuseFill(pc, PETSC_TRUE);genius_assert(!ierr);
      ierr = PCFactorSetPivoting(pc, 1.0); genius_assert(!ierr);
      ierr = PCFactorReorderForNonzeroDiagonal(pc, 1e-20); genius_assert(!ierr);
      ierr = PCFactorSetShiftNonzero(pc, 1e-12);  genius_assert(!ierr);
      //ierr = PetscOptionsSetValue("-pc_factor_nonzeros_along_diagonal", 0); genius_assert(!ierr);
      ierr = PetscOptionsSetValue("-pc_factor_diagonal_fill", 0); genius_assert(!ierr);
    }
    return;

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


void FVM_NonlinearSolver::sens_solve()
{
  // do snes solve
  SNESSolve ( snes, PETSC_NULL, x );

  // get the converged reason
  SNESConvergedReason reason;
  SNESGetConvergedReason ( snes,&reason );

  // if Line search failed, disable Line search
  if ( reason == SNES_DIVERGED_LS_FAILURE || reason == SNES_DIVERGED_LOCAL_MIN )
  {
    MESSAGE <<"------> nonlinear solver " << SNESConvergedReasons[reason] <<". Disable Line Search.\n\n\n";
    RECORD();
    SNESLineSearchSet ( snes,SNESLineSearchNo,PETSC_NULL );
    SNESSolve ( snes, PETSC_NULL, x );
  }

#if 0
  if ( reason == SNES_DIVERGED_LINEAR_SOLVE )
  {
  MESSAGE <<"------> nonlinear solver " << SNESConvergedReasons[reason] <<". Set to use direct solver.\n\n\n";
  RECORD();
#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
  KSPSetType (ksp, (char*) KSPPREONLY);
  PCSetType  (pc, (char*) PCLU);
  PCFactorSetMatSolverPackage (pc, MAT_SOLVER_MUMPS);
  PetscOptionsSetValue("-mat_mumps_icntl_14", "80");
#endif
#if PETSC_VERSION_LE(2,3,3)
  MatSetType(J,MATAIJMUMPS);
  PetscOptionsSetValue("-mat_mumps_icntl_14","80");
#endif
  SNESSolve ( snes, PETSC_NULL, x );
}
#endif


}

