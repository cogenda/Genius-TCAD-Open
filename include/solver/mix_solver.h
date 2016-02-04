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


#ifndef __mix_solver_h__
#define __mix_solver_h__

#include "ddm_solver.h"

/**
 * common functiuons for advanced mixed mode simulation
 * and redefine some methods in DDMSolverBase
 */
class MixSolverBase : public DDMSolverBase
{
public:
  /**
   * constructor
   */
  MixSolverBase(SimulationSystem & system): DDMSolverBase(system), _circuit(system.get_circuit())
  {}

  /**
   * destructor
   */
  virtual ~MixSolverBase()
  {}

  /**
   * virtual function, create the solver
   */
  virtual int create_solver();

  /**
   * virtual function, destroy the solver
   */
  virtual int destroy_solver();


  /**
   * do pre-process before each solve action
   */
  virtual int pre_solve_process(bool load_solution=true);


  /**
   * @return the bandwidth contributed by each boundary node
   */
  virtual unsigned int n_bc_node_dofs(const BoundaryCondition * ) const
  { return 0; }


  /**
   * @return the extra dofs of each boundary condition.
   * each electrode has an extra dof
   */
  virtual unsigned int bc_dofs(const BoundaryCondition * bc) const
  {
    if(bc->is_electrode()) return 1;
    return 0;
  }

  /**
   * set the global dof in the electrode, which point to spice node
   */
  void link_electrode_to_spice_node();

  /**
   * set the matrix nonzero pattern for spice circuit
   */
  virtual void set_extra_matrix_nonzero_pattern();


  /**
   * call ckt_load to build spice rhs/matrix, then
   * add spice rhs to petsc vector
   */
  void build_spice_function(PetscScalar *lxx, Vec f, InsertMode &add_value_flag);

  /**
   * add spice matrix to petsc matrix
   */
  void build_spice_jacobian(PetscScalar *lxx, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * compute the common mode voltage of all the electrode,
   * should subtract this value before simulation
   * and add it back after simulation
   */
  virtual PetscReal common_model_voltage() const;

  /**
   * output information about spice node
   */
  void print_spice_node() const;

  /**
   * compute steady-state
   */
  virtual int solve_dcop(bool tran_op=false);

  /**
   * do dcsweep, can be both voltage scan or current scan
   */
  virtual int solve_dcsweep();

  /**
   * do transient simulation
   */
  virtual int solve_transient();

  /**
   * IV curve automatically trace
   */
  //virtual int solve_iv_trace();

  /**
   * function for convergence test of pseudo time step method
   */
  virtual bool pseudo_time_step_convergence_test() { return true; }


  /**
   * function for line search post check. update SPICE solution here
   */
  virtual void sens_line_search_post_check(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w);

  /**
   * sens convergence criteria for mixed type solvers
   */
  virtual void petsc_snes_convergence_test(PetscInt its, PetscReal , PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason);

  /**
   * ksp convergence criteria for dd solvers
   */
  virtual void petsc_ksp_convergence_test(PetscInt its, PetscReal rnorm, KSPConvergedReason* reason);

protected:

  /**
   * hold the pointer to spice circuit
   */
  SPICE_CKT * _circuit;

  /**
   * bc to spice node map
   */
  std::vector< std::pair<BoundaryCondition *, unsigned int> > _bc_to_ckt_node_link;

  /**
   * f norm of spice equation
   */
  PetscScalar spice_norm;

};

#endif //#define __mix_solver_h__
