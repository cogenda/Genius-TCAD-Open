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

#include "fvm_linear_solver.h"
#include "parallel.h"


/*------------------------------------------------------------------
 * constructor, setup context
 */
FVM_LinearSolver::FVM_LinearSolver(SimulationSystem & system): FVM_PDESolver(system)
{}


/*------------------------------------------------------------------
 * setup nonlinear matrix/vector
 */
void FVM_LinearSolver::setup_linear_data()
{
  // map mesh to PETSC solver
  build_dof_map();

  // set petsc routine

  PetscErrorCode ierr;

  // create the global solution vector
  ierr = VecCreateMPI(PETSC_COMM_WORLD, n_local_dofs, n_global_dofs, &x); genius_assert(!ierr);
  // use VecDuplicate to create vector with same pattern
  ierr = VecDuplicate(x, &b);genius_assert(!ierr);
  ierr = VecDuplicate(x, &L);genius_assert(!ierr);

  // set all the components of scale vector L to 1.0
  ierr = VecSet(L, 1.0); genius_assert(!ierr);

  // create local vector, which has extra room for ghost dofs! the MPI_COMM here is PETSC_COMM_SELF
  ierr = VecCreateSeq(PETSC_COMM_SELF,  local_index_array.size() , &lx); genius_assert(!ierr);
  // use VecDuplicate to create vector with same pattern
  ierr = VecDuplicate(lx, &lb);genius_assert(!ierr);

  // create the index for vector statter
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, global_index_array.size() , &global_index_array[0] , &gis); genius_assert(!ierr);
  ierr = ISCreateGeneral(PETSC_COMM_SELF,  local_index_array.size() ,  &local_index_array[0] ,  &lis); genius_assert(!ierr);

  // it seems we can free global_index_array and local_index_array to save the memory

  // create the vector statter
  ierr = VecScatterCreate(x, gis, lx, lis, &scatter); genius_assert(!ierr);


  // create the matrix
  ierr = MatCreate(PETSC_COMM_WORLD, &A); genius_assert(!ierr);
  ierr = MatSetSizes(A, n_local_dofs, n_local_dofs, n_global_dofs, n_global_dofs); genius_assert(!ierr);

#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
  // we are using petsc-devel
  if (Genius::n_processors()>1)
  {
    ierr = MatSetType(A,MATMPIAIJ); genius_assert(!ierr);
    // alloc memory for parallel matrix here
    ierr = MatMPIAIJSetPreallocation(A, 0, &n_nz[0], 0, &n_oz[0]); genius_assert(!ierr);
  }
  else
  {
    ierr = MatSetType(A,MATSEQAIJ); genius_assert(!ierr);
    // alloc memory for sequence matrix here
    ierr = MatSeqAIJSetPreallocation(A, 0, &n_nz[0]); genius_assert(!ierr);
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
      {
        MESSAGE<< "Warning:  Default LU solver can not be used in parallel, use BCGS instead!" << std::endl;
        RECORD();
        ierr = KSPSetType (ksp, (char*) KSPBCGS);       genius_assert(!ierr); return;
        break;
      }
    case   SolverSpecify::MUMPS :
      {
#ifdef PETSC_HAVE_MUMPS
        MESSAGE<< "Using MUMPS linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(A,MATAIJMUMPS); genius_assert(!ierr);
#else
        MESSAGE<< "Warning:  no MUMPS solver configured, use BCGS instead!" << std::endl;
        RECORD();
        ierr = MatSetType(A,MATMPIAIJ); genius_assert(!ierr);
        SolverSpecify::LS = SolverSpecify::BICGSTAB;
#endif
        break;
      }
    case   SolverSpecify::SuperLU_DIST:
      {
#ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
        MESSAGE<< "Using SuperLU_DIST linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(A,MATSUPERLU_DIST); genius_assert(!ierr);
        ierr = PetscOptionsSetValue("-mat_superlu_dist_rowperm","NATURAL");
#else
        MESSAGE << "Warning:  no SuperLU_DIST solver configured, use BCGS instead!" << std::endl;
        RECORD();
        ierr = MatSetType(A,MATMPIAIJ); genius_assert(!ierr);
        SolverSpecify::LS = SolverSpecify::BICGSTAB;
#endif
        break;
      }

      // no SuperLU_DIST/MUMPS? we return to KSP method
    default:
      ierr = MatSetType(A,MATMPIAIJ); genius_assert(!ierr);
    }
    // alloc memory for parallel matrix here
    ierr = MatMPIAIJSetPreallocation(A, 0, &n_nz[0], 0, &n_oz[0]); genius_assert(!ierr);
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
        ierr = MatSetType(A,MATUMFPACK); genius_assert(!ierr);
#else
        MESSAGE << "Warning:  no UMFPACK solver configured, use default LU solver instead!" << std::endl;
        RECORD();
        ierr = MatSetType(A,MATSEQAIJ); genius_assert(!ierr);
#endif
        break;
      }

      // check if we can use SuperLU
    case SolverSpecify::SuperLU :
      {
#ifdef PETSC_HAVE_SUPERLU
        MESSAGE<< "Using SuperLU linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(A,MATSUPERLU);genius_assert(!ierr);
#else
        MESSAGE << "Warning:  no SuperLU solver configured, use default LU solver instead!" << std::endl;
        RECORD();
        ierr = MatSetType(A,MATSEQAIJ); genius_assert(!ierr);
#endif
        break;
      }
      //use MUMPS for serial problem
    case   SolverSpecify::MUMPS :
      {
#ifdef PETSC_HAVE_MUMPS
        MESSAGE<< "Using MUMPS linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(A,MATAIJMUMPS); genius_assert(!ierr);
#else
        MESSAGE << "Warning:  no MUMPS solver configured, use default LU solver instead!" << std::endl;
        RECORD();
        ierr = MatSetType(A,MATSEQAIJ); genius_assert(!ierr);
#endif
        break;
      }
      //use SuperLU_DIST for serial problem
    case   SolverSpecify::SuperLU_DIST:
      {
#ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
        MESSAGE<< "Using SuperLU_DIST linear solver..."<<std::endl;
        RECORD();
        ierr = MatSetType(A,MATSUPERLU_DIST); genius_assert(!ierr);
        ierr = PetscOptionsSetValue("-mat_superlu_dist_rowperm","NATURAL"); genius_assert(!ierr);
#else
        MESSAGE << "Warning:  no SuperLU_DIST solver configured, use default LU solver instead!" << std::endl;
        RECORD();
        ierr = MatSetType(A,MATSEQAIJ); genius_assert(!ierr);
#endif
        break;
      }
      // we have to use the default Mat type MATSEQAIJ
    case SolverSpecify::LU :
    default:
      ierr = MatSetType(A,MATSEQAIJ); genius_assert(!ierr);

    }
    // alloc memory for sequence matrix here
    ierr = MatSeqAIJSetPreallocation(A, 0, &n_nz[0]); genius_assert(!ierr);

  }
#endif

  // indicates when MatZeroRows() is called the zeroed entries are kept in the nonzero structure
#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
  ierr = MatSetOption(A, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE); genius_assert(!ierr);
#endif

#if PETSC_VERSION_LE(2,3,3)
  ierr = MatSetOption(A, MAT_KEEP_ZEROED_ROWS); genius_assert(!ierr);
#endif
  //ierr = MatSetOption(A, MAT_YES_NEW_NONZERO_LOCATIONS); genius_assert(!ierr);

  // indicates when MatSetValue with ADD_VALUES mode, the 0 entries will be ignored
  // since I uas ADD_VALUES 0 for keeping nonzero pattern, do not use this parameter!
  //ierr = MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES); genius_assert(!ierr);

  // the matrix is not assembled yet.
  matrix_first_assemble = false;

  ierr = MatSetFromOptions(A); genius_assert(!ierr);

  // create petsc linear solver context
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); genius_assert(!ierr);

  // set corresponding matrix
  ierr = KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN); genius_assert(!ierr);

  // get petsc preconditional context
  ierr = KSPGetPC(ksp, &pc); genius_assert(!ierr);

  // Set user-specified linear solver and preconditioner types
  set_petsc_linear_solver_type    ( SolverSpecify::LS );
  set_petsc_preconditioner_type   ( SolverSpecify::PC );

}


/*------------------------------------------------------------------
 * destroy nonlinear data
 */
void FVM_LinearSolver::clear_linear_data()
{
  PetscErrorCode ierr;
  // free everything
  ierr = VecDestroy(x);              genius_assert(!ierr);
  ierr = VecDestroy(b);              genius_assert(!ierr);
  ierr = VecDestroy(L);              genius_assert(!ierr);
  ierr = VecDestroy(lx);             genius_assert(!ierr);
  ierr = VecDestroy(lb);             genius_assert(!ierr);
  ierr = ISDestroy(gis);             genius_assert(!ierr);
  ierr = ISDestroy(lis);             genius_assert(!ierr);
  ierr = VecScatterDestroy(scatter); genius_assert(!ierr);
  ierr = MatDestroy(A);              genius_assert(!ierr);
}


/*------------------------------------------------------------------
 * destructor: destroy context
 */
FVM_LinearSolver::~FVM_LinearSolver()
{
  PetscErrorCode ierr;
  ierr = KSPDestroy(ksp);          genius_assert(!ierr);
}






void FVM_LinearSolver::set_petsc_linear_solver_type(SolverSpecify::LinearSolverType linear_solver_type)
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
        MESSAGE<< "Warning:  Default LU solver can not be used in parallel, use BCGS instead!" << std::endl;
        RECORD();
        ierr = KSPSetType (ksp, (char*) KSPBCGS);       genius_assert(!ierr); return;
        break;

      case   SolverSpecify::MUMPS :
        // if LU method is required, we should check if SuperLU_DIST/MUMPS is installed
        // the default parallel LU solver is set to MUMPS
#ifdef PETSC_HAVE_MUMPS
        MESSAGE<< "Using MUMPS linear solver..."<<std::endl;
        RECORD();
        ierr = KSPSetType (ksp, (char*) KSPPREONLY); genius_assert(!ierr);
        ierr = PCSetType  (pc, (char*) PCLU); genius_assert(!ierr);
        ierr = PCFactorSetMatSolverPackage (pc, MAT_SOLVER_MUMPS); genius_assert(!ierr);
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
    //ierr = PCFactorReorderForNonzeroDiagonal(pc, 1e-20); genius_assert(!ierr); <-- Caught signal number 11 SEGV error will occure when diag value < 1e-20
    ierr = PCFactorSetShiftNonzero(pc, 1e-20); genius_assert(!ierr);
    return;

  default:
    std::cerr << "ERROR:  Unsupported PETSC Solver: "
    << linear_solver_type        << std::endl
    << "Continuing with PETSC defaults" << std::endl;
  }

}



void FVM_LinearSolver::set_petsc_preconditioner_type(SolverSpecify::PreconditionerType preconditioner_type)
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
        ierr = PCFactorSetPivoting(pc, 1.0); genius_assert(!ierr);
        ierr = PCFactorReorderForNonzeroDiagonal(pc, 1e-20); genius_assert(!ierr);
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
    std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
    << preconditioner_type       << std::endl
    << "Continuing with PETSC defaults" << std::endl;
  }
}


