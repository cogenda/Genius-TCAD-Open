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


#ifndef __mixA_solver_h__
#define __mixA_solver_h__

#include "ddm_solver.h"

/**
 * common functiuons for advanced mixed mode simulation
 * and redefine some methods in DDMSolverBase
 */
class MixASolverBase : public DDMSolverBase
{
public:
  /**
   * constructor
   */
  MixASolverBase(SimulationSystem & system): DDMSolverBase(system), _circuit(system.get_circuit())
  {}

  /**
   * destructor
   */
  virtual ~MixASolverBase()
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
   */
  virtual unsigned int bc_dofs(const BoundaryCondition * bc) const
  {
    if(bc->is_electrode() && !bc->is_spice_electrode()) return 1;
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
   * @return the extra dofs of spice circuit
   */
  virtual unsigned int extra_dofs() const;

  /**
   * set initial guess of SNES solver
   */
  void spice_fill_value(Vec x, Vec L);

  /**
   * call ckt_load to build spice rhs/matrix, then
   * add spice rhs to petsc vector
   */
  void build_spice_function(PetscScalar *lxx, Vec f, InsertMode &add_value_flag);

  /**
   * add spice matrix to petsc matrix
   */
  void build_spice_jacobian(PetscScalar *lxx, Mat *jac, InsertMode &add_value_flag);

  /**
   * dump spice matrix in petsc format to external file
   * for more detailed analysis of the properties of jacobian matrix
   */
  void dump_spice_matrix_petsc(const std::string &file) const;

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
   * f norm of spice equation
   */
  PetscScalar spice_norm;

};

#endif //#define __mixA_solver_h__
