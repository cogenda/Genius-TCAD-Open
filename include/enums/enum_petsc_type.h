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

//  $Id: enum_petsc_type.h,v 1.5 2008/07/09 05:58:15 gdiso Exp $



#ifndef __enum_petsc_type_h__
#define __enum_petsc_type_h__



// ------------------------------------------------------------
// enum SolverType definition
namespace SolverSpecify
{

  /**
   * Defines an \p enum for iterative solver types
   */
  enum NonLinearSolverType {
    Newton=0,
    LineSearch,
    TrustRegion,
    INVALID_NONLINEAR_SOLVER};



  /**
   * Defines an \p enum for linear solver types
   */
  enum LinearSolverType {CG=0,         // iterative solvers
                         CGN,
                         CGS,
                         CR,
                         QMR,
                         TCQMR,
                         TFQMR,
                         BICG,
                         BICGSTAB,
                         BCGSL,
                         MINRES,
                         GMRES,
                         DGMRES,
                         FGMRES,
                         LSQR,
                         JACOBI,
                         SOR_FORWARD,
                         SOR_BACKWARD,
                         SSOR,
                         RICHARDSON,
                         CHEBYSHEV,
                         LU,            // direct solvers
                         UMFPACK,
                         SuperLU,
                         PASTIX,
                         MUMPS,
                         SuperLU_DIST,
                         GSS,
                         INVALID_LINEAR_SOLVER};

 /**
  * Defines an \p enum for category of linear solvers
  */
  enum LinearSolverCategory {ITERATIVE=0,
                             DIRECT,
                             HYBRID};

  /**
   * Defines an \p enum for preconditioner types
   */
  enum PreconditionerType {IDENTITY_PRECOND =0,
                           JACOBI_PRECOND,
                           BLOCK_JACOBI_PRECOND,
                           SOR_PRECOND,
                           SSOR_PRECOND,
                           EISENSTAT_PRECOND,
                           BOOMERAMG_PRECOND,
                           ASM_PRECOND,
                           ASMILU0_PRECOND,
                           ASMILU1_PRECOND,
                           ASMILU2_PRECOND,
                           ASMILU3_PRECOND,
                           ASMLU_PRECOND,
                           CHOLESKY_PRECOND,
                           ICC_PRECOND,
                           ILU_PRECOND,
                           ILUT_PRECOND,
                           LU_PRECOND,
                           PARMS_PRECOND,
                           USER_PRECOND,
                           SHELL_PRECOND,
                           INVALID_PRECONDITIONER};


}




#endif // #define __enum_petsc_type_h__




