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



#ifndef __emfem2d_h__
#define __emfem2d_h__

#include <vector>

#include "enum_petsc_type.h"
#include "fem_linear_solver.h"
#include "petscksp.h"


/**
 * derived from fem linear solver contex.
 */
class EMFEM2DSolver : public FEM_LinearSolver
{
public:

  /**
   * construct the linear solver contex:
   * solution vector, right hand side (RHS) vector, matrix
   * as well as parallel scatter
   */
  EMFEM2DSolver(SimulationSystem & system, const Parser::Card & c)
  : FEM_LinearSolver(system), _card(c)
  {system.record_active_solver(this->solver_type());}

  /**
   * free all the contex
   */
  virtual ~EMFEM2DSolver() {}

  /**
   * @return the solver type
   */
  virtual SolverSpecify::SolverType solver_type() const
  {return SolverSpecify::EM_FEM_2D;}

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
   * @return nodal dof number.
   * because we use complex equtions, set 2 for each node
   */
  virtual unsigned int node_dofs() const
  { return 2; }


private:

  /**
   * parse the command card
   */
  void setup_solver_parameters();

  /**
   * read wave spectrum from file
   */
  void parse_spectrum_file(const std::string &);

  /**
   * hold the reference of command card
   */
  const Parser::Card & _card;

  /**
   * the parameters of optical source
   */
  struct OpticalSource
  {
     double  wave_length;// wave length
     double  power;      // wave power
     double  TE_weight;  // weight of TE mode in incidence light
     double  TM_weight;  // weight of TM mode in incidence light
     double  phi_TE;     // the initial phase angle for TE
     double  phi_TM;     // the initial phase angle for TM
     double  eta;        // quantum efficiency at this wave length
     bool    eta_auto;   // if it is true, determine quantum efficiency by photon energy and material bandgap
  };

  /**
   * all the optical sources are listed here
   */
  std::vector<OpticalSource> _optical_sources;

  /**
   * the angle of incident wave
   */
  double _incidence_angle;

  /**
   * solve TE scatter problem
   * @param lamda  wave length
   * @param power  wave power
   * @param phase0 wave phase
   */
  void solve_TM_scatter_problem(double lamda, double power, double phase0);

  /**
   * solve TM scatter problem
   * @param lamda  wave length
   * @param power  wave power
   * @param phase0 wave phase
   */
  void solve_TE_scatter_problem(double lamda, double power, double phase0);

  /**
   * build the matrix A, precondition matrix PC and RHS for TE scatter problem
   */
  void build_TE_matrix_rhs(double lamda, double power, double phase0);

  /**
   * build the matrix A, precondition matrix PC and RHS for TM scatter problem
   */
  void build_TM_matrix_rhs(double lamda, double power, double phase0);

  /**
   * save nodal solution to fvm node data structure
   * when append is true, the new solution will be added to previous solution
   */
  void save_TE_solution(double lamda, double power, double phase0, double eta, bool eta_auto, bool append=true);

  /**
   * save nodal solution to fvm node data structure
   * when append is true, the new solution will be added to previous solution
   */
  void save_TM_solution(double lamda, double power, double phase0, double eta, bool eta_auto, bool append=true);

  /**
   * the edge chain of Absorbing Boundary
   */
  std::vector< std::pair<const Elem *, unsigned int> > absorb_edge_chain;

  /**
   * the node chain of Absorbing Boundary
   */
  std::vector< std::pair<const Node *, const Node *> > absorb_node_chain;

  /**
   * set the absorb boundary
   */
  void build_absorb_chain();

  /**
   * aux function to build looped vector index
   */
  unsigned int loop_count(int size, int current, int step)
  {
     int index = current + step;
     if(index<0) index = size+index;
     if(index>=size) index = index-size;
     return static_cast<unsigned int>(index);
  }

  /**
   * @return the curvature at two nodes of absorbing edge of ith element in Absorbing Boundary
   */
  std::vector<double> curvature_at_edge(unsigned int i, const Elem * edge);

  /**
   * @return curvature of the circle (we know Absorbing Boundary is circle)
   */
  std::vector<double> curvature_of_circle(const Point &p1, const Point &p2, const Point &p3);

  /**
   * the Absorbing Boundary type
   */
  enum ABC_TYPE {FirstOrder, SecondOrder, PML};

  ABC_TYPE _abc_type;

  /**
   * when we use Mu Absorbing Boundary, we have to build the curvature of the boundary
   * if we know the Absorbing Boundary is circle, things are much easier.
   */
  enum ABC_SHAPE{Circle,  Ellipse, UNKNOWN_SHAPE};

  ABC_SHAPE _abc_shape;

};


#endif // #define __emfem2d_h__

