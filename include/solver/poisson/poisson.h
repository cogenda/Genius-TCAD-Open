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

//  $Id: poisson.h,v 1.15 2008/07/09 05:58:16 gdiso Exp $

#ifndef __poisson_h__
#define __poisson_h__

#include "fvm_flex_nonlinear_solver.h"
#include "enum_solver_specify.h"

/**
 * Solve the nonlinear Poisson's equation of
 * \nabla \eps \nabla \Phi = - \rho
 * where \rho is the net charge dengsity in all the region.
 * For semiconductor region, electron and hole density n and p
 * are the function of \Phi ( actually, related with fermi potential  )
 */
class PoissonSolver : public FVM_FlexNonlinearSolver
{
public:
  /**
   * the constructor of PoissonSolver, take system as parameter
   */
  PoissonSolver(SimulationSystem & system): FVM_FlexNonlinearSolver(system)
  {system.record_active_solver(this->solver_type());}

  /**
   * destructor, do nothing
   */
  ~PoissonSolver()
  { }

  /**
   * @return the solver type
   */
  virtual SolverSpecify::SolverType solver_type() const
  {return SolverSpecify::POISSON;}

  /**
   * virtual function, create the solver
   */
  virtual int create_solver();

  /**
   * virtual function, do the solve process
   */
  virtual int solve();

  /**
   * virtual function, destroy the solver, release internal data
   */
  virtual int destroy_solver() ;

  /**
   * do pre-process before each solve action
   */
  virtual int pre_solve_process(bool load_solution=true);

  /**
   * do post-process after each solve action
   */
  virtual int post_solve_process();

  /**
   * @return node's dof for each region. only 1 dof here
   */
  virtual unsigned int node_dofs(const SimulationRegion * region) const
  {
    assert(region!=NULL);
    switch(region->type())
    {
      case SemiconductorRegion : return 1;
      case InsulatorRegion     : return 1;
      case ElectrodeRegion     : return 1;
      case MetalRegion         : return 1;
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
      case ChargeIntegral    : return 1;// the ChargeIntegral bc has one extra equation
      default: return 0;   // others, no extra bc equation
    }
  }


  /**
   * PETSC SNES can have an individual prefix
   */
  virtual std::string snes_prefix() const { return "poisson_"; }

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
  virtual void sens_line_search_post_check(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
  {
    this->potential_damping(x, y, w, changed_y, changed_w);
    FVM_FlexNonlinearSolver::sens_line_search_post_check(x, y, w, changed_y, changed_w);
  }

private:

  void potential_damping(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w);

};


#endif // #define __poisson_h__
