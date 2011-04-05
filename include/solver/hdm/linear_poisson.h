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

#ifndef __linear_poisson_h__
#define __linear_poisson_h__

#include "fvm_linear_solver.h"
#include "enum_solver_specify.h"

/**
 * Solve the hydrodynamic equations with explicit scheme.
 * suitable for very large problem
 */
class LinearPoissonSolver : public FVM_LinearSolver
{
public:
  /**
   * the constructor of HDMSolver, take system as parameter
   */
  LinearPoissonSolver(SimulationSystem & system): FVM_LinearSolver(system)
  {system.record_active_solver(this->solver_type());}

  /**
   * destructor, do nothing
   */
  ~LinearPoissonSolver()
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
    case ChargedContact    : return 1;   // the ChargedContact bc has one extra equation
    default: return 0;   // others, no extra bc equation
    }
  }

  /**
   * @return the matrix bandwidth of each boundary condition which owns extra dofs
   */
  virtual unsigned int bc_bandwidth(const BoundaryCondition * bc) const
  {
    switch (bc->bc_type())
    {
      case ChargedContact   : return 1;  // the ChargedContact bc equation has a bandwidth of 1
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
      case ChargedContact    : return 1; // node on ChargedContact has one dof to bc equation
      default: return 0;
    }
  }


  /**
   * indicates if PDE involves all neighbor elements.
   * when it is true, the matrix bandwidth will include all the nodes belongs to neighbor elements, i.e. DDM solver
   * when it is false, only neighbor nodes (link local node by edge) are appeared in matrix bandwidth, i.e. poisson solver.
   */
  virtual bool all_neighbor_elements_involved(const SimulationRegion * ) const
  { return false; }



  /**
   * virtual function for building the RHS vector
   */
  virtual void build_rhs(Vec );

  /**
   * virtual function for building the matrix A and precondition matrix PC
   */
  virtual void build_matrix(Mat , Mat );


private:
  /**
   * write solution to system
   */
  void update_solution();
};


#endif // #define __linear_poisson_h__
