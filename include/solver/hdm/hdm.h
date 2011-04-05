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

#ifndef __hdm_solver_h__
#define __hdm_solver_h__

#include "fvm_explicit_solver.h"
#include "enum_solver_specify.h"

/**
 * Solve the hydrodynamic equations with explicit scheme.
 * suitable for very large problem
 */
class HDMSolver : public FVM_ExplicitSolver
{
  public:
    HDMSolver ( SimulationSystem & system );

    virtual ~HDMSolver();

    /**
     * @return the solver type
     */
    virtual SolverSpecify::SolverType solver_type() const
    {return SolverSpecify::HDM;}

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
    virtual int pre_solve_process ( bool load_solution=true );

    /**
     * do post-process after each solve action
     */
    virtual int post_solve_process();

    /**
     * @return node's dof for each region. only Semiconductor Region has 8 dof here
     */
    virtual unsigned int node_dofs ( const SimulationRegion * region ) const
    {
      assert ( region!=NULL );
      switch ( region->type() )
      {
        case SemiconductorRegion : return 8;
        default : return 0;
      }
    }

    /**
     * indicates if PDE involves all neighbor elements.
     * when it is true, the matrix bandwidth will include all the nodes belongs to neighbor elements, i.e. DDM solver
     * when it is false, only neighbor nodes (link local node by edge) are appeared in matrix bandwidth, i.e. poisson solver.
     */
    virtual bool all_neighbor_elements_involved ( const SimulationRegion * ) const
    { return false; }


  private:

    /**
     * internal poisson solver
     */
    SolverBase * poisson_solver;

    /**
     * solve
     */
    void solve_equ();

    /**
     * solve steady-state, this is nearly the same as
     */
    void solve_steadystate();

    /**
     * transient simulation, dual time step method is used here
     */
    void solve_transient();

    /**
     * write solution to system
     */
    void update_solution();

  private:

    void local_time_advance(int , Real alpha);

    bool convergence_criteria(int );

    /**
     * write charge density to system variable rho
     * as data exchange with poisson solver
     */
    void sync_rho();

};


#endif // #define __hdm_solver_h__
