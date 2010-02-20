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

#include "ddm_ac/ddm_ac.h"
#include "parallel.h"
#include "mathfunc.h"  // for PI


using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::cm;
using PhysicalUnit::K;




/*------------------------------------------------------------------
 * create linear solver contex and adjust some parameters
 */
int DDMACSolver::create_solver()
{
  int ierr=0;

  MESSAGE<< '\n' << "AC Small Signal Solver init..." << std::endl;
  RECORD();

  // must setup linear contex here!
  setup_linear_data();

  // rtol    - the relative convergence tolerance (relative decrease in the residual norm)
  // abstol  - the absolute convergence tolerance (absolute size of the residual norm)
  KSPSetTolerances(ksp, SolverSpecify::ksp_rtol, SolverSpecify::ksp_atol, PETSC_DEFAULT, std::max(50, static_cast<int>(n_global_dofs/10)));

  // user can do further adjusment from command line
  KSPSetFromOptions (ksp);

  // set extra 2 vecs we needed here
  ierr = VecDuplicate(x, &s);   genius_assert(!ierr);
  ierr = VecDuplicate(lx, &ls); genius_assert(!ierr);

  // extra matrix for store Jacobian
  ierr = MatCreate(PETSC_COMM_WORLD, &J); genius_assert(!ierr);
  ierr = MatSetSizes(J, n_local_dofs, n_local_dofs, n_global_dofs, n_global_dofs); genius_assert(!ierr);
  ierr = MatSetType(J, MATMPIAIJ); genius_assert(!ierr);
  ierr = MatMPIAIJSetPreallocation(J, 0, &n_nz[0], 0, &n_oz[0]); genius_assert(!ierr);

  // call hook function on_init
  hook_list()->on_init();

  return 0;

}




/*------------------------------------------------------------------
 * call this function before each solution process
 */
int DDMACSolver::pre_solve_process(bool /*load_solution*/)
{
  /*
   * fill previous system (DC) solution into Vec s.
   * Then Vec s is used for pre-build Jacobian matrix J.
   *
   * load_solution is ignored, we always need to fill solution
   */

  // for all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDMAC_Fill_Value(s, L);
  }

  VecAssemblyBegin(s);
  VecAssemblyEnd(s);

  VecAssemblyBegin(L);
  VecAssemblyEnd(L);

  /*
   * fill Matrix J with function EBM3_Jacobian
   */

  // scatte global solution vector s to local vector ls
  VecScatterBegin(scatter, s, ls, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, s, ls, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lss;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray(ls, &lss);

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

  // evaluate Jacobian matrix of governing equations of EBM in all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->EBM3_Jacobian(lss, &J, add_value_flag);
  }

  // restore array back to Vec
  VecRestoreArray(ls, &lss);

  // assembly the matrix J
  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (J, MAT_FINAL_ASSEMBLY);


  /*
   * assign VAC to corresponding electrode
   */

  // clear VAC for all the electrode
  for(unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(n);
    if(bc->is_electrode())
      bc->ext_circuit()->Vac() = 0.0;
  }

  for(unsigned int n=0; n<SolverSpecify::Electrode_ACScan.size(); ++n)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(SolverSpecify::Electrode_ACScan[n]);
    genius_assert(bc!=NULL);
    bc->ext_circuit()->Vac() = SolverSpecify::VAC;
  }

  return FVM_LinearSolver::pre_solve_process(true);

}





/*------------------------------------------------------------------
 * solve implimention!
 */
int DDMACSolver::solve()
{
  START_LOG("ACSolver()", "KSP_Solver");

  this->pre_solve_process();


  for(SolverSpecify::Freq = SolverSpecify::FStart; SolverSpecify::Freq <= SolverSpecify::FStop; SolverSpecify::Freq*=SolverSpecify::FMultiple)
  {

    double omega = 2*PI*SolverSpecify::Freq;

    MESSAGE
      <<"AC Scan: f("<<SolverSpecify::Electrode_ACScan[0]<<") = "
      << std::setiosflags(std::ios::fixed)
      <<SolverSpecify::Freq*PhysicalUnit::s/1e6<<" MHz "<<"\n";
    RECORD();

    build_ddm_ac(A, b, omega);

    KSPSolve(ksp, b, x);

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);

    PetscInt   its;
    KSPGetIterationNumber(ksp, &its);

    PetscReal  rnorm;
    KSPGetResidualNorm(ksp, &rnorm);

    MESSAGE<<"------> residual norm = "<<rnorm<<" its = "<<its<<" with "<<KSPConvergedReasons[reason]<<"\n\n";
    RECORD();

    this->post_solve_process();

  }


  STOP_LOG("ACSolver()", "KSP_Solver");

  return 0;
}




/*------------------------------------------------------------------
 * call this function after each solution process
 */
int DDMACSolver::post_solve_process()
{

  double omega = 2*PI*SolverSpecify::Freq;

  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  VecGetArray(lx, &lxx);

  //update solution for all regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDMAC_Update_Solution(lxx);
  }

  //update solution for all bcs
  for(unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(n);
    bc->DDMAC_Update_Solution(lxx, J, omega);
  }

  VecRestoreArray(lx, &lxx);

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
  ierr = VecDestroy(s);              genius_assert(!ierr);
  ierr = VecDestroy(ls);             genius_assert(!ierr);
  ierr = MatDestroy(J);              genius_assert(!ierr);

  return FVM_LinearSolver::destroy_solver();
}






///////////////////////////////////////////////////////////////
//        Provide matrix evaluation for DDM AC solver        //
///////////////////////////////////////////////////////////////







/*------------------------------------------------------------------
 * build the matrix and right hand side vector b with certain freq omega
 */
void DDMACSolver::build_ddm_ac(Mat A, Vec b, double omega)
{
  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

  MatZeroEntries(A);
  VecZeroEntries(b);

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDMAC_Fill_Matrix_Vector(A, b, J, omega, add_value_flag);
  }

  // evaluate Jacobian matrix of governing equations of EBM for all the boundaries
  for(unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(n);
    bc->DDMAC_Fill_Matrix_Vector(A, b, J, omega, add_value_flag);
  }

  // assembly the matrix
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);

  // assembly the vec
  VecAssemblyBegin(b);
  VecAssemblyEnd  (b);

  //MatView(A, PETSC_VIEWER_DRAW_WORLD);
  //getchar();

}


