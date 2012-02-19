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

#include "fvm_explicit_solver.h"
#include "parallel.h"


/*------------------------------------------------------------------
 * constructor, setup context
 */
FVM_ExplicitSolver::FVM_ExplicitSolver(SimulationSystem & system): FVM_PDESolver(system)
{}

/*------------------------------------------------------------------
 * destructor: destroy context
 */
FVM_ExplicitSolver::~FVM_ExplicitSolver()
{}


/*------------------------------------------------------------------
 * setup nonlinear matrix/vector
 */
void FVM_ExplicitSolver::setup_explicit_data()
{
  // map mesh to PETSC solver
  build_dof_map();

  // set petsc routine

  PetscErrorCode ierr;

  // create the global solution vector
  ierr = VecCreateMPI(PETSC_COMM_WORLD, n_local_dofs, n_global_dofs, &x); genius_assert(!ierr);
  // use VecDuplicate to create vector with same pattern
  ierr = VecDuplicate(x, &x_prev);genius_assert(!ierr);
  ierr = VecDuplicate(x, &vol);genius_assert(!ierr);
  ierr = VecDuplicate(x, &f);genius_assert(!ierr);
  ierr = VecDuplicate(x, &t);genius_assert(!ierr);
  ierr = VecDuplicate(x, &w);genius_assert(!ierr);


  // create local vector, which has extra room for ghost dofs! the MPI_COMM here is PETSC_COMM_SELF
  ierr = VecCreateSeq(PETSC_COMM_SELF,  local_index_array.size() , &lx); genius_assert(!ierr);
  ierr = VecDuplicate(lx, &lt);genius_assert(!ierr);

  // create the index for vector statter
#if PETSC_VERSION_GE(3,2,0)
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, global_index_array.size(), &global_index_array[0] , PETSC_COPY_VALUES, &gis); genius_assert(!ierr);
  ierr = ISCreateGeneral(PETSC_COMM_SELF,  local_index_array.size(),  &local_index_array[0] ,  PETSC_COPY_VALUES, &lis); genius_assert(!ierr);
#else
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, global_index_array.size(), &global_index_array[0] , &gis); genius_assert(!ierr);
  ierr = ISCreateGeneral(PETSC_COMM_SELF,  local_index_array.size(),  &local_index_array[0] ,  &lis); genius_assert(!ierr);
#endif

  // it seems we can free global_index_array and local_index_array to save the memory

  // create the vector statter
  ierr = VecScatterCreate(x, gis, lx, lis, &scatter); genius_assert(!ierr);

}


/*------------------------------------------------------------------
 * destroy explicit data
 */
void FVM_ExplicitSolver::clear_explicit_data()
{
  PetscErrorCode ierr;
  // free everything
  ierr = VecDestroy(PetscDestroyObject(x));              genius_assert(!ierr);
  ierr = VecDestroy(PetscDestroyObject(x_prev));         genius_assert(!ierr);
  ierr = VecDestroy(PetscDestroyObject(lx));             genius_assert(!ierr);
  ierr = VecDestroy(PetscDestroyObject(lt));             genius_assert(!ierr);
  ierr = VecDestroy(PetscDestroyObject(vol));            genius_assert(!ierr);
  ierr = VecDestroy(PetscDestroyObject(f));              genius_assert(!ierr);
  ierr = VecDestroy(PetscDestroyObject(t));              genius_assert(!ierr);
  ierr = VecDestroy(PetscDestroyObject(w));              genius_assert(!ierr);
  ierr = ISDestroy(PetscDestroyObject(gis));             genius_assert(!ierr);
  ierr = ISDestroy(PetscDestroyObject(lis));             genius_assert(!ierr);
  ierr = VecScatterDestroy(PetscDestroyObject(scatter)); genius_assert(!ierr);
}

