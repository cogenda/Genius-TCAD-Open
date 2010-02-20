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


#ifndef __fem_pde_solver_h__
#define __fem_pde_solver_h__

#include "solver_base.h"



/**
 * The pde discretization usually involves not only local node but also its neighbors
 * we should know more information to solve this type of system, especially in parallel
 */
class FEM_PDESolver : public  SolverBase
{

public:
  /**
   * constructor
   */
  FEM_PDESolver(SimulationSystem & system)
  : SolverBase(system),
  n_global_dofs(0),
  n_local_dofs(0),
  global_offset(0)
  {}

  /**
   * destructor
   */
  virtual ~FEM_PDESolver()
  {
    local_index_array.clear();
    global_index_array.clear();
  }

  /**
   * set data structure for parallel vector and matrix layout
   * this function implements the mapping from mesh structure
   * to PETSC parallel linear/nonlinear solver.
   */
  void set_parallel_dof_map();

  /**
   * set data structure for serial vector and matrix layout
   * this function implements the mapping from mesh structure
   * to PETSC serial linear/nonlinear solver. only for limited
   * serial version of Genius.
   */
  void set_serial_dof_map();

  /**
   * wrap function for call set_parallel_dof_map() or set_serial_dof_map()
   */
  void build_dof_map();

  /**
   * @return the (exact) nodal dofs of each simulation region
   */
  virtual unsigned int node_dofs() const
  { return 1; }

  /**
   * Fills the vector \p di with the global degree of freedom indices
   * for the element.
   */
  void build_dof_indices (const Elem* const elem, std::vector<PetscInt>& di) const;

protected:

  /**
   * the summary of node's dof, in global
   */
  unsigned int n_global_dofs;

  /**
   * total number of freedom on this processor,
   */
  unsigned int n_local_dofs;

  /**
   * the offset of local block of dofs at global dofs
   */
  unsigned int global_offset;

  /**
   * this array contains the local dof index as well as ghost dofs!
   * the ghost dof are located after the position of n_local_dofs
   */
  std::vector<PetscInt> local_index_array;

  /**
   * this array contains the global dof index as well as ghost dofs!
   * the ghost dof are located after the position of n_local_dofs
   */
  std::vector<PetscInt> global_index_array;

  /**
   * The number of on-processor nonzeros in my portion of the
   * global matrix.
   */
  std::vector<PetscInt> n_nz;

  /**
   * The number of off-processor nonzeros in my portion of the
   * global matrix.
   */
  std::vector<PetscInt> n_oz;

};

#endif //define __fem_pde_solver_h__
