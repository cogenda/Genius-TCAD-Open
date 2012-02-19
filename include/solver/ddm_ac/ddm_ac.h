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


#ifndef __ddm_ac_solver_h__
#define __ddm_ac_solver_h__

#include "enum_petsc_type.h"
#include "fvm_linear_solver.h"
#include "petscksp.h"
#include "solver_specify.h"

/**
 * The AC small signal solver contex.
 */
class DDMACSolver : public FVM_LinearSolver
{
public:

  /**
   * construct the linear solver contex:
   * solution vector, right hand side (RHS) vector, matrix
   * as well as parallel scatter
   */
  DDMACSolver(SimulationSystem & system)
  : FVM_LinearSolver(system),_first_create(true)
  {
    system.record_active_solver(this->solver_type());
  }

  /**
   * free all the contex
   */
  virtual ~DDMACSolver();

  /**
   * @return the solver type
   */
  virtual SolverSpecify::SolverType solver_type() const
  {return SolverSpecify::DDMAC;}

  /**
   * virtual function, create the solver
   */
  virtual int create_solver();

  /**
   * virtual functions, prepare solution and aux variables used by this solver
   */
  virtual int set_variables();

  /**
   * call this function before each solve process
   */
  virtual int pre_solve_process(bool load_solution=true);

  /**
   * virtual function, do the solve process
   */
  virtual int solve();

  /**
   * call this function after each solve process
   */
  virtual int post_solve_process();

  /**
   * virtual function, destroy the solver, release internal data
   */
  virtual int destroy_solver() ;

  /**
   * @return node's dof for each region.
   * since they are all complex numbers, the required dofs are doubled for real solver
   */
  virtual unsigned int node_dofs(const SimulationRegion * region) const
  {
    switch(region->type())
    {
    case SemiconductorRegion :
      switch(region->get_advanced_model()->EB_Level)
      {
      case ModelSpecify::NONE : return 6; // 3 basic dofs, psi, n and p.
      case ModelSpecify::Tl   :
      case ModelSpecify::Tn   :
      case ModelSpecify::Tp   : return 8; // 4 dofs for one of the extra Tl, Tn and Tp equation
      case ModelSpecify::TnTp :
      case ModelSpecify::TnTl :
      case ModelSpecify::TpTl : return 10; // 5 dofs for Tn/Tp, Tn/Tl or Tp/Tl pair
      case ModelSpecify::ALL  : return 12; // 6 dofs for all of the extra Tl, Tn and Tp equations
      }
      //dofs for insulator/conductor node
    case InsulatorRegion     :
    case ElectrodeRegion     :
    case MetalRegion         :
      switch(region->get_advanced_model()->EB_Level)
      {
      case ModelSpecify::NONE :
      case ModelSpecify::Tn   :
      case ModelSpecify::Tp   :
      case ModelSpecify::TnTp : return 2; // only potential equation
      case ModelSpecify::Tl   :
      case ModelSpecify::TnTl :
      case ModelSpecify::TpTl :
      case ModelSpecify::ALL  : return 4;// potential equation plus Tl equation
      }
    default : return 0;
    }
  }

  /**
   * @return the dofs of each boundary condition.
   * since they are all complex numbers, the required dofs are doubled for real solver
   */
  virtual unsigned int bc_dofs(const BoundaryCondition * bc) const
  {
    switch (bc->bc_type())
    {
    case OhmicContact      :
    case SchottkyContact   :
    case SimpleGateContact :
    case GateContact       :
    case SolderPad         :
    case ChargedContact    :
      // the above bcs has one extra equation
      return 2;
      // others, no extra bc equation
    default: return 0;
    }
  }

  /**
   * @return the (approximate) boundary condition dofs contributed by each boundary node
   */
  virtual unsigned int bc_node_dofs(const BoundaryCondition * bc) const
  {
    // return maximum possiable dofs of all the EBM level
    switch (bc->bc_type())
    {
    case OhmicContact      : return 12; // ohmic electrode current
    case SchottkyContact   : return 2;  // displacement current
    case SimpleGateContact : return 2;  // displacement current
    case GateContact       : return 2;  // displacement current
    case SolderPad         : return 2;  // conductance current
    case ChargedContact    : return 2;
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
    case ElectrodeRegion     : return false;
    case MetalRegion         : return false;
    default : return false;
    }
  }


private:

  /**
   * the global system solution vector
   */
  Vec            s;

  /**
   * the local system solution vector
   */
  Vec            ls;

  /**
   * the Jacobian matrix
   */
  Mat            J_;

  /**
   * the unpreconditioned AC matrix
   */
  Mat            A_;

  /**
   * the unpreconditioned rhs vector
   */
  Vec            b_;

  /**
   * the transformation matrix, used for precondition
   * Z.-Y. Wang, K.-C. Wu, and R. W. Dutton, “An approach to construct pre-conditioning matrices
   * for block iteration of linear equations,” IEEE Trans. CAD, Vol. 11, No. 11, pp. 1334-1343,
   * 1992.
   */
  Mat            T_;

  /**
   * temp matrix
   */
  Mat            C_;

  /**
   * flag to show if A_ is created
   */
  bool           _first_create;

  /**
   * building the Matrix A, RHS vector b under certain freq omega
   */
  void build_ddm_ac(double omega);
};


#endif // #define __ddm_ac_solver_h__

