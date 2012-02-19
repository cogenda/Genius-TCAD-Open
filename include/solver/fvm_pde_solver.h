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


#ifndef __fvm_pde_solver_h__
#define __fvm_pde_solver_h__

#include "solver_base.h"



/**
 * The pde discretization usually involves not only local node but also its neighbors
 * we should know more information to solve this type of system, especially in parallel
 */
class FVM_PDESolver : public  SolverBase
{

public:
  /**
   * constructor
   */
  FVM_PDESolver(SimulationSystem & system)
  : SolverBase(system),
  n_global_node_dofs(0),
  n_global_bc_dofs(0),
  n_global_dofs(0),
  n_local_dofs(0),
  global_offset(0)
  {}

  /**
   * destructor
   */
  virtual ~FVM_PDESolver()
  {
    local_index_array.clear();
    global_index_array.clear();
  }


  /**
   * set the index used by global_offset
   */
  void set_solver_index(unsigned int index)
  {
    FVM_Node::set_solver_index(index);
    BoundaryCondition::set_solver_index(index);
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
   * virtual function indicates if PDE involves all neighbor elements.
   * when it is true, the matrix bandwidth will include all the nodes belongs to neighbor elements, i.e. DDM solver
   * when it is false, only neighbor nodes (link local node by edge) are appeared in matrix bandwidth, i.e. poisson solver.
   */
  virtual bool all_neighbor_elements_involved(const SimulationRegion * ) const
  { return false; }

  /**
   * @return the (exact) nodal dofs of each simulation region
   */
  virtual unsigned int node_dofs(const SimulationRegion * region) const
  { genius_assert(region!=NULL); return 1; }

  /**
   * @return the (exact) dofs of each boundary condition
   */
  virtual unsigned int bc_dofs(const BoundaryCondition * bc) const
  { genius_assert(bc!=NULL); return 0; }

  /**
   * @return the matrix bandwidth of each boundary condition which owns extra dofs
   */
  virtual unsigned int bc_bandwidth(const BoundaryCondition * bc) const
  { genius_assert(bc!=NULL); return 0; }

  /**
   * @return the (approximate) dofs contributed by each boundary node
   * it is used to statistic bandwidth of boundary condition equation
   */
  virtual unsigned int bc_node_dofs(const BoundaryCondition * bc) const
  { genius_assert(bc!=NULL); return 0; }

  /**
   * @return extra dofs if some solver has
   */
  virtual unsigned int extra_dofs() const
  { return 0; }

  /**
   * set the matrix nonzero pattern for extra dofs
   */
  virtual void set_extra_matrix_nonzero_pattern()
  { return; }

protected:

  /**
   * the summary of node's dof, in global
   */
  unsigned int n_global_node_dofs;

  /**
   * dofs of boundary condition, in gloabl
   */
  unsigned int n_global_bc_dofs;

  /**
   * total number of freedom in global
   * equals to n_global_node_dofs + n_global_bc_dofs + extra_dofs
   */
  unsigned int n_global_dofs;

  /**
   * total number of freedom on this processor,
   * bc_dofs and extra_dofs are processed with the last processor.
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

#endif //define __fvm_pde_solver_h__
