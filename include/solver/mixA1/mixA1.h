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


#ifndef __mixA1_h__
#define __mixA1_h__


#include "mixA_solver.h"
#include "solver_specify.h"

class SPICE_CKT;

/**
 * Advanced Mixed type simulation with Level 1 Drift-Diffusion Model equations
 */
class MixA1Solver : public MixASolverBase
{
public:
  MixA1Solver(SimulationSystem & system): MixASolverBase(system)
  {
    system.record_active_solver(this->solver_type());
  }


  ~MixA1Solver()
  {}

  /**
   * @return the solver type
   */
  virtual SolverSpecify::SolverType solver_type() const
  {return SolverSpecify::DDML1MIXA;}

  /**
   * virtual function, create the solver
   */
  virtual int create_solver();

  /**
   * do pre-process before each solve action
   */
  virtual int pre_solve_process(bool load_solution=true);

  /**
   * virtual function, do the solve process
   */
  virtual int solve();

  /**
   * do post-process after each solve action
   */
  virtual int post_solve_process();

  /**
   * load previous state into solution vector
   */
  virtual int diverged_recovery();

  /**
   * @return the extra dofs of spice circuit
   */
  virtual unsigned int extra_dofs() const;

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
   * compute the norm of local truncate error (LTE)
   */
  virtual PetscScalar LTE_norm();

  /**
   * check carrier density after projection
   */
  virtual void projection_positive_density_check(Vec x, Vec xo);

  /**
   * compute the abs and relative error norm of the solution
   */
  virtual void error_norm();

  /**
   * @return the bandwidth contributed by each boundary node
   */
  virtual unsigned int n_bc_node_dofs(const BoundaryCondition * bc) const
  {
    switch (bc->bc_type())
    {
      case OhmicContact      : return 3; // ohmic electrode current
      case SchottkyContact   : return 1; // displacement current
      case SimpleGateContact : return 1; // displacement current
      case GateContact       : return 1; // displacement current
      case ChargedContact    : return 1; // for electrostatic Gauss's law
      default: return 0;
    }
  }

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



#endif // #define __mixA1_h__
