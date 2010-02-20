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


#ifndef __mix1_h__
#define __mix1_h__



#include "mix_common/mix_common.h"
#include "solver_specify.h"

#ifndef CYGWIN

/**
 * Mixed type simulation with Level 1 Drift-Diffusion Model equations
 */
class Mix1Solver : public MixSolverBase
{
public:
  Mix1Solver(SimulationSystem & system): MixSolverBase(system)
  {system.record_active_solver(this->solver_type());}


  ~Mix1Solver()
  {}

  /**
   * @return the solver type
   */
  virtual SolverSpecify::SolverType solver_type() const
  {return SolverSpecify::DDML1MIX;}

  /**
   * virtual function, create the solver
   */
  virtual int create_solver();

  /**
   * virtual function, do the solve process
   */
  virtual int solve();

  /**
   * do pre-process before each solve action
   */
  virtual int pre_solve_process(bool load_solution=true);

  /**
   * do post-process after each solve action
   */
  virtual int post_solve_process();

  /**
   * get pdI/pdw, pdI/pdV and pdF/pdV for each electrode
   */
  virtual int get_electrode_load();

  /**
   * load previous state into solution vector
   */
  virtual int diverged_recovery();

  /**
   * Save steady-state solution as the previous solution data in spice
   */
  virtual int init_spice_data();

  /**
   * load solution data previously accepted by spice
   */
  virtual int load_spice_data();

  /**
   * since spice accept the solution, save solution data
   */
  virtual int save_spice_data();

  /**
   * @return node's dof for each region.
   */
  virtual unsigned int node_dofs(const SimulationRegion * region) const
  {
    switch(region->type())
    {
      case SemiconductorRegion : return 3; //semiconductor node has 3 dof
      case InsulatorRegion     : return 1; //insulator node has 1 dof
      case ConductorRegion     : return 1; //conductor node has 1 dof
      default : return 0;
    }
  }

  /**
   * @return the dofs of each boundary condition.
   */
  virtual unsigned int bc_dofs(const BoundaryCondition * bc) const
  {
    switch (bc->bc_type())
    {
      case ChargedContact    :
      // the above bcs has one extra equation
      return 1;
      // others, no extra bc equation
      default: return 0;
    }
  }


  /**
   * @return the matrix bandwidth of each boundary condition which owns extra dofs
   */
  virtual unsigned int bc_bandwidth(const BoundaryCondition * bc) const
  {
    switch (bc->bc_type())
    {
      case ChargedContact    :
      // the above bc equation has a bandwidth of 1
      return 1;
      // others, bandwidth is zero
      default: return 0;
    }
  }

  /**
   * @return the boundary condition dofs contributed by each boundary node
   */
  virtual unsigned int bc_node_dofs(const BoundaryCondition * bc) const
  {
    switch (bc->bc_type())
    {
      case ChargedContact    : return 1;// for electrostatic Gauss's law
      default: return 0;
    }
  }


  /**
   * indicates if PDE involves all neighbor elements.
   * when it is true, the matrix bandwidth will include all the nodes belongs to neighbor elements, i.e. DDM solver
   * when it is false, only neighbor nodes (link local node by edge) are appeared in matrix bandwidth, i.e. poisson solver.
   */
  virtual bool all_neighbor_elements_involved(const SimulationRegion * region) const
  {
    switch(region->type())
    {
      case SemiconductorRegion : return true;
      case InsulatorRegion     : return false;
      case ConductorRegion     : return false;
      default : return false;
    }
  }

  /**
   * wrap function for evaluating the residual of function f at x
   */
  virtual void build_petsc_sens_residual(Vec x, Vec r);

  /**
   * wrap function for evaluating the Jacobian J of function f at x
   */
  virtual void build_petsc_sens_jacobian(Vec x, Mat *jac, Mat *pc);

  /**
   * function for line search post check. do Newton damping here
   */
  virtual void sens_line_search_post_check(Vec x, Vec y, Vec w, PetscTruth *changed_y, PetscTruth *changed_w)
  {
    switch( SolverSpecify::Damping )
    {
      case SolverSpecify::DampingPotential : potential_damping(x, y, w, changed_y, changed_w); return;
      case SolverSpecify::DampingBankRose  : bank_rose_damping(x, y, w, changed_y, changed_w); return;
      case SolverSpecify::DampingNo        : positive_density_damping(x, y, w, changed_y, changed_w); return;
      default: positive_density_damping(x, y, w, changed_y, changed_w);
    }
    return;
  }

  /**
   * check carrier density after projection
   */
  virtual void projection_positive_density_check(Vec x, Vec xo);

  /**
   * compute the abs and relative error norm of the solution
   */
  virtual void error_norm();

private:

  /**
   * Potential Newton damping scheme
   */
  void potential_damping(Vec x, Vec y, Vec w, PetscTruth *changed_y, PetscTruth *changed_w);

  /**
   * Bank-Rose Newton damping scheme
   */
  void bank_rose_damping(Vec x, Vec y, Vec w, PetscTruth *changed_y, PetscTruth *changed_w);

  /**
   * Positive carrier density Newton damping scheme
   */
  void positive_density_damping(Vec x, Vec y, Vec w, PetscTruth *changed_y, PetscTruth *changed_w);
};

#endif //#ifndef CYGWIN

#endif // #define __mix1_h__
