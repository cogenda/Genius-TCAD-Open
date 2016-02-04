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

#include <iomanip>
#include "petsc_matrix.h"
#include "ddm_ac/ddm_ac.h"
#include "parallel.h"
#include "mathfunc.h"  // for PI


using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::cm;
using PhysicalUnit::K;



/**
 * free all the contex
 */
DDMACSolver::~DDMACSolver()
{
  FVM_Node::set_solver_index(0);
  BoundaryCondition::set_solver_index(0);
}


/*------------------------------------------------------------------
 * create linear solver contex and adjust some parameters
 */
int DDMACSolver::create_solver()
{
  int ierr=0;

  MESSAGE<< '\n' << "AC Small Signal Solver init..." << std::endl;
  RECORD();


  if( SolverSpecify::Electrode_ACScan.empty() )
  {
    MESSAGE<<"ERROR: You must specify one electrode for AC scan."<<std::endl; RECORD();
    genius_error();
  }

  // set ac variables for each region
  set_variables();

  set_solver_index(3);

  set_linear_solver_type    ( SolverSpecify::LS );
  set_preconditioner_type   ( SolverSpecify::PC );

  // must setup linear contex here!
  setup_linear_data();

  // rtol    - the relative convergence tolerance (relative decrease in the residual norm)
  // abstol  - the absolute convergence tolerance (absolute size of the residual norm)
  KSPSetTolerances ( ksp, 1e-15, SolverSpecify::ksp_atol, PETSC_DEFAULT, std::max ( 50, static_cast<int> ( n_global_dofs/10 ) ) );

  // user can do further adjusment from command line
  KSPSetFromOptions ( ksp );

  // set extra 2 vecs we needed here
  ierr = VecDuplicate ( x, &s );  genius_assert ( !ierr );
  ierr = VecDuplicate ( lx, &ls );  genius_assert ( !ierr );

  // extra matrix for store Jacobian
  ierr = MatCreate ( PETSC_COMM_WORLD, &J_ );  genius_assert ( !ierr );
  ierr = MatSetSizes ( J_, n_local_dofs, n_local_dofs, n_global_dofs, n_global_dofs );genius_assert ( !ierr );
  if ( Genius::n_processors() >1 )
  {
    ierr = MatSetType ( J_, MATMPIAIJ );  genius_assert ( !ierr );
    ierr = MatMPIAIJSetPreallocation ( J_, 0, &n_nz[0], 0, &n_oz[0] );  genius_assert ( !ierr );
  }
  else
  {
    ierr = MatSetType ( J_, MATSEQAIJ );    genius_assert ( !ierr );
    // alloc memory for sequence matrix here
    ierr = MatSeqAIJSetPreallocation ( J_, 0, &n_nz[0] );    genius_assert ( !ierr );
  }
  // we have to set this flag since preallocation is not exact
  ierr = MatSetOption(J_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); genius_assert(!ierr);


  // extra matrix for store A
  ierr = MatCreate ( PETSC_COMM_WORLD, &A_ );  genius_assert ( !ierr );
  ierr = MatSetSizes ( A_, n_local_dofs, n_local_dofs, n_global_dofs, n_global_dofs );genius_assert ( !ierr );
  if ( Genius::n_processors() >1 )
  {
    ierr = MatSetType ( A_, MATMPIAIJ );  genius_assert ( !ierr );
    ierr = MatMPIAIJSetPreallocation ( A_, 0, &n_nz[0], 0, &n_oz[0] );  genius_assert ( !ierr );
  }
  else
  {
    ierr = MatSetType ( A_, MATSEQAIJ );    genius_assert ( !ierr );
    // alloc memory for sequence matrix here
    ierr = MatSeqAIJSetPreallocation ( A_, 0, &n_nz[0] );    genius_assert ( !ierr );
  }
  // we have to set this flag since preallocation is not exact
  ierr = MatSetOption(A_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); genius_assert(!ierr);


  // extra matrix for transformation matrix, each row has only 2 entry
  ierr = MatCreate ( PETSC_COMM_WORLD, &T_ );  genius_assert ( !ierr );
  ierr = MatSetSizes ( T_, n_local_dofs, n_local_dofs, n_global_dofs, n_global_dofs );  genius_assert ( !ierr );
  if ( Genius::n_processors() >1 )
  {
    ierr = MatSetType ( T_, MATMPIAIJ );    genius_assert ( !ierr );
    ierr = MatMPIAIJSetPreallocation ( T_, 2, PETSC_NULL, 0, PETSC_NULL );    genius_assert ( !ierr );
  }
  else
  {
    ierr = MatSetType ( T_, MATSEQAIJ );    genius_assert ( !ierr );
    // alloc memory for sequence matrix here
    ierr = MatSeqAIJSetPreallocation ( T_, 2, PETSC_NULL );    genius_assert ( !ierr );
  }


  // extra vector for store T*b
  VecDuplicate ( b, &b_ );

  MESSAGE<< "AC Small Signal Solver init finished." << std::endl;
  RECORD();

  return FVM_LinearSolver::create_solver();

}



/*------------------------------------------------------------------
 * prepare solution and aux variables used by this solver
 */
int DDMACSolver::set_variables()
{
  for ( unsigned int n=0; n<_system.n_regions(); n++ )
  {
    SimulationRegion * region = _system.region ( n );

    switch ( region->type() )
    {
      case SemiconductorRegion :
      {
        region->add_variable("electron.ac", POINT_CENTER);
        region->add_variable("hole.ac", POINT_CENTER);
        region->add_variable("potential.ac", POINT_CENTER);
        region->add_variable("temperature.ac", POINT_CENTER);
        region->add_variable("elec_temperature.ac", POINT_CENTER);
        region->add_variable("hole_temperature.ac", POINT_CENTER);
        break;
      }
      case InsulatorRegion :
      case ElectrodeRegion :
      case MetalRegion :
      {
        region->add_variable("potential.ac", POINT_CENTER);
        region->add_variable("temperature.ac", POINT_CENTER);
        break;
      }
      default: break;
    }
  }

  return 0;
}




/*------------------------------------------------------------------
 * call this function before each solution process
 */
int DDMACSolver::pre_solve_process ( bool /*load_solution*/ )
{
  /*
   * fill previous system (DC) solution into Vec s.
   * Then Vec s is used for pre-build Jacobian matrix J.
   *
   * load_solution is ignored, we always need to fill solution
   */

  // for all the regions
  for ( unsigned int n=0; n<_system.n_regions(); n++ )
  {
    SimulationRegion * region = _system.region ( n );
    region->DDMAC_Fill_Value ( s, L );
  }

  VecAssemblyBegin ( s );
  VecAssemblyEnd ( s );

  VecAssemblyBegin ( L );
  VecAssemblyEnd ( L );

  /*
   * fill Matrix J with function EBM3_Jacobian
   */

  // scatte global solution vector s to local vector ls
  VecScatterBegin ( scatter, s, ls, INSERT_VALUES, SCATTER_FORWARD );
  VecScatterEnd ( scatter, s, ls, INSERT_VALUES, SCATTER_FORWARD );

  PetscScalar *lss;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray ( ls, &lss );

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

  // evaluate Jacobian matrix of governing equations of EBM in all the regions
  PetscMatrix<PetscScalar> *Jac = new PetscMatrix<PetscScalar>(n_global_dofs, n_global_dofs, n_local_dofs, n_local_dofs);
  for ( unsigned int n=0; n<_system.n_regions(); n++ )
  {
    SimulationRegion * region = _system.region(n);
    region->EBM3_Jacobian(lss, Jac, add_value_flag);
  }
  Jac->close(true);
  MatCopy(Jac->mat(), J_, DIFFERENT_NONZERO_PATTERN);
  delete Jac;

  // restore array back to Vec
  VecRestoreArray ( ls, &lss );



  /*
   * assign VAC to corresponding electrode
   */

  // clear VAC for all the electrode
  for ( unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n )
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc ( n );
    if ( bc->is_electrode() )
      bc->ext_circuit()->Vac() = 0.0;
  }

  for ( unsigned int n=0; n<SolverSpecify::Electrode_ACScan.size(); ++n )
  {
    std::vector<BoundaryCondition *> bcs = _system.get_bcs()->get_bcs_by_electrode_label( SolverSpecify::Electrode_ACScan[n] );
    for(size_t b=0; b<bcs.size(); b++)
    {
      BoundaryCondition * bc = bcs[b];
      if(bc && bc->is_electrode())
        bc->ext_circuit()->Vac() = SolverSpecify::VAC;
    }
  }

  return FVM_LinearSolver::pre_solve_process ( true );

}





/*------------------------------------------------------------------
 * solve implimention!
 */
int DDMACSolver::solve()
{
  START_LOG ( "solve()", "DDMACSolver" );

  set_solver_index(3);

  this->pre_solve_process();

  for ( SolverSpecify::Freq = SolverSpecify::FStart; SolverSpecify::Freq <= SolverSpecify::FStop;  )
  {

    double omega = 2*PI*SolverSpecify::Freq;

    MESSAGE
    <<"AC Scan: f("<<SolverSpecify::Electrode_ACScan[0]<<") = "
    << std::scientific
    <<SolverSpecify::Freq*PhysicalUnit::s/1e6<<" MHz "<<"\n";
    RECORD();

    build_ddm_ac ( omega );

    KSPSolve ( ksp, b, x );

    KSPConvergedReason reason;
    KSPGetConvergedReason ( ksp, &reason );

    PetscInt   its;
    KSPGetIterationNumber ( ksp, &its );

    PetscReal  rnorm;
    KSPGetResidualNorm ( ksp, &rnorm );

    MESSAGE<<"------> residual norm = "<<rnorm<<" its = "<<its<<" with "<<KSPConvergedReasons[reason]<<"\n\n";
    RECORD();

    this->post_solve_process();

    if( SolverSpecify::Freq  < SolverSpecify::FStop && SolverSpecify::Freq*SolverSpecify::FMultiple > SolverSpecify::FStop)
      SolverSpecify::Freq  = SolverSpecify::FStop;
    else
      SolverSpecify::Freq*=SolverSpecify::FMultiple;
  }


  STOP_LOG ( "solve()", "DDMACSolver" );

  return 0;
}




/*------------------------------------------------------------------
 * call this function after each solution process
 */
int DDMACSolver::post_solve_process()
{

  PetscScalar omega = 2*PI*SolverSpecify::Freq;

  VecScatterBegin ( scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD );
  VecScatterEnd ( scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD );

  PetscScalar *lxx;
  VecGetArray ( lx, &lxx );

  //update solution for all regions
  for ( unsigned int n=0; n<_system.n_regions(); n++ )
  {
    SimulationRegion * region = _system.region ( n );
    region->DDMAC_Update_Solution ( lxx );
  }

  //update solution for all bcs
  for ( unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n )
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc ( n );
    bc->DDMAC_Update_Solution ( lxx, J_, omega );
  }

  VecRestoreArray ( lx, &lxx );

  return FVM_LinearSolver::post_solve_process();
}








/*------------------------------------------------------------------
 * restore the solution to each region
 */
int DDMACSolver::destroy_solver()
{
  int ierr = 0;

  // clear linear matrix/vector
  clear_linear_data();

  // destroy the Vec and Mat we defined in class DDMACSolver
  ierr = VecDestroy ( PetscDestroyObject(s) );
  genius_assert ( !ierr );
  ierr = VecDestroy ( PetscDestroyObject(ls) );
  genius_assert ( !ierr );
  ierr = MatDestroy ( PetscDestroyObject(J_) );
  genius_assert ( !ierr );
  ierr = MatDestroy ( PetscDestroyObject(A_) );
  genius_assert ( !ierr );
  ierr = MatDestroy ( PetscDestroyObject(T_) );
  genius_assert ( !ierr );
  ierr = VecDestroy ( PetscDestroyObject(b_) );
  genius_assert ( !ierr );

  if ( !_first_create ) MatDestroy ( PetscDestroyObject(C_) );

  return FVM_LinearSolver::destroy_solver();
}






///////////////////////////////////////////////////////////////
//        Provide matrix evaluation for DDM AC solver        //
///////////////////////////////////////////////////////////////







/*------------------------------------------------------------------
 * build the matrix and right hand side vector b with certain freq omega
 */
void DDMACSolver::build_ddm_ac ( PetscScalar omega )
{

  START_LOG ( "build_ddm_ac()", "DDMACSolver" );

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

  MatZeroEntries ( A_ );
  VecZeroEntries ( b_ );

  // evaluate Jacobian matrix of governing equations of EBM for all the regions
  for ( unsigned int n=0; n<_system.n_regions(); n++ )
  {
    SimulationRegion * region = _system.region ( n );
    region->DDMAC_Fill_Matrix_Vector ( A_, b_, J_, omega, add_value_flag );
  }

  // evaluate Jacobian matrix of governing equations of EBM for all the boundaries
  for ( unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n )
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc ( n );
    bc->DDMAC_Fill_Matrix_Vector ( A_, b_, J_, omega, add_value_flag );
  }

  // assembly the matrix A
  MatAssemblyBegin ( A_, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd ( A_, MAT_FINAL_ASSEMBLY );

  // assembly the vec b
  VecAssemblyBegin ( b_ );
  VecAssemblyEnd ( b_ );


  // process transformation matrix
  {
    MatZeroEntries ( T_ );
    add_value_flag = NOT_SET_VALUES;
    for ( unsigned int n=0; n<_system.n_regions(); n++ )
    {
      SimulationRegion * region = _system.region ( n );
      region->DDMAC_Fill_Transformation_Matrix ( T_, J_, omega, add_value_flag );
    }

    if(Genius::processor_id() == Genius::n_processors() -1)
    {
      for ( unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n )
      {
        BoundaryCondition * bc = _system.get_bcs()->get_bc ( n );
        if ( !bc->is_electrode() ) continue;
        MatSetValue ( T_, bc->global_offset(), bc->global_offset(), 1.0, ADD_VALUES );
        MatSetValue ( T_, bc->global_offset() +1, bc->global_offset() +1, 1.0, ADD_VALUES );
      }
    }

    // assembly the transformation matrix
    MatAssemblyBegin ( T_, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd ( T_, MAT_FINAL_ASSEMBLY );


    // do transport
    if ( _first_create )
    {
      MatMatMult ( T_, A_, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C_ );
      _first_create = false;
    }
    else
      MatMatMult ( T_, A_, MAT_REUSE_MATRIX , PETSC_DEFAULT, &C_ );
    MatCopy(C_, A, DIFFERENT_NONZERO_PATTERN);

    MatMult ( T_, b_, b );
  }


  //MatView(A, PETSC_VIEWER_DRAW_WORLD);
  //getchar();

  STOP_LOG ( "build_ddm_ac()", "DDMACSolver" );

}


