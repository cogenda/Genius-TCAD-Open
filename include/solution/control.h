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

//  $Id: control.h,v 1.9 2008/07/09 05:58:16 gdiso Exp $

#ifndef __device_solver_control_h__
#define __device_solver_control_h__


#include "auto_ptr.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "simulation_system.h"
#include "solver_base.h"
#include "mxml.h"


// purify namespace
// undefine "REAL", which was defined by the tetgen and triangle mesh generator
#ifdef REAL
#undef REAL
#endif


/**
 * control the genius simulation solvers
 */
class SolverControl
{
public:

  /**
   * constructor
   */
  SolverControl();

  /**
   * destructor
   */
  ~SolverControl();

  /**
   * get the solution records as a dom
   */
  mxml_node_t* get_dom_solution() const;

  /**
   * get count of solution records
   */
  int get_dom_solution_count() const;

  /**
   * we loop here, parse each input card and do corresponding operator
   */
  int mainloop();

  /**
   * set the input deck
   */
  void setDecks(Parser::InputParser *input);

  /**
   * set result output file (xml format)
   */
  void setSolutionFile(const std::string &fname_result) { _fname_solution=fname_result; }

  /**
   * get result output file name
   */
  std::string getSolutionFile() { return _fname_solution; }

  Parser::InputParser& decks() { return *_decks;}
  const Parser::InputParser& decks() const { return *_decks; }

  Mesh& mesh() { return *_mesh; }
  const Mesh& mesh() const { return *_mesh; }

  SimulationSystem& system() { return *_system; }
  const SimulationSystem& system() const { return *_system; }

  //----------------------------------------------
  // device simulation
  //----------------------------------------------

  /**
   * re-create the simulation system
   * the current deck must contain GLOBAL command, external electrical and
   * magnetic source, and circuit net-list commands
   */
  int reset_simulation_system();

  /**
   * do process simulation. set doping profile and mole fraction to semiconductor region
   */
  int do_process();

  /**
   * call mesh generator to build the mesh
   */
  int  do_mesh ();

  /**
   * process "METHOD" card.
   * @note only restore solver parameters.
   */
  int  set_method ( const Parser::Card & c );

   /**
   * process "MODEL" card.
   * @note only restore physical model parameters.
   */
  int  set_model ( const Parser::Card & c );

  /**
   * process "HOOK" card
   */
  int do_hook ( const Parser::Card & c);

  /**
   * process and do "SOLVE" card
   */
  int  do_solve   ( const Parser::Card & c );

  /**
   * process and do "EXPORT" card
   */
  int  do_export  ( const Parser::Card & c );

  /**
   * process and do "IMPORT" card
   */
  int  do_import  ( const Parser::Card & c );

  /**
   * create 2d fem solver to do electromagnetic calculation
   * mainly for optical device simulation
   */
  int do_em_fem2d_solve ( const Parser::Card & c );

  /**
   * create ray trace solver to do electromagnetic calculation
   * mainly for optical device simulation
   */
  int do_ray_trace( const Parser::Card & c );

  /**
   * process and do "NODESET" card
   */
  int  set_initial_node_voltage  ( const Parser::Card & c );

  /**
   * process and do "ATTACH" card
   */
  int  set_electrode_source  ( const Parser::Card & c );

  /**
   * process and do "PMI" card
   */
  int  set_physical_model  ( const Parser::Card & c );

  /**
   * process and do "REFINE.CONFORM" card
   */
  int  do_refine_conform ( const Parser::Card & c );

  /**
   * process and do "REFINE.HIERARCHICAL" card
   */
  int do_refine_hierarchical ( const Parser::Card & c );

  /**
   * process and do "REFINE.UNIFORM" card
   */
  int  do_refine_uniform  ( const Parser::Card & c );

  /**
   * extend 2d mesh to 3d mesh
   */
  int  extend_to_3d ( const Parser::Card & c );

  /**
   * plot 2d mesh?
   */
  int  plot_mesh ( const Parser::Card & c );

private:

  Parser::InputParser *_decks;

  /**
   * the main mesh structure
   */
  AutoPtr<Mesh> _mesh;

  /**
   * the main simulation system
   */
  AutoPtr<SimulationSystem> _system;

  /**
   * define the meshgen pointer,
   * we may use it if user want to generate mesh
   * At the same time, meshgen is useful when AMR is required
   */
  AutoPtr<MeshGeneratorBase> meshgen;

  /**
   * import mesh from CAD system
   */
  int  import_mesh ();

  /**
   * also, a solver which provides doping profile is required in most case.
   * it can be a serious process simulator,
   * a interpolation method reading doping information from text file
   * or just simple analytic functions
   * when AMR is required, we can easily build new doping profile.
   * however, if user only import an existing model into Genius,
   * no enough information to build meshgen and DopingSolver structure.
   * when AMR requires, only do an general purpose refinement
   * and interpolation doping from exist mesh.
   */
  AutoPtr<SolverBase> DopingSolver;

  /**
   * mole fraction for compound material. set it as a solver
   */
  AutoPtr<SolverBase> MoleSolver;

  /**
   * solution records as a XML document
   */
  mxml_node_t *_dom_solution;

  std::string _fname_solution;
};

class SolverControlHook : public Hook
{
public:
  SolverControlHook(SolverBase & solver, const std::string & name, SolverControl& control, const std::string& fname);

  virtual ~SolverControlHook();

  /**
   *   This is executed before the initialization of the solver
   */
  virtual void on_init();

  /**
   *   This is executed previously to each solution step.
   */
  virtual void pre_solve();

  /**
   *  This is executed after each solution step.
   */
  virtual void post_solve();

  /**
   *  This is executed after each (nonlinear) iteration
   */
  virtual void post_iteration();

  /**
   *  This is executed after each (nonlinear) iteration
   */
  virtual void post_iteration(void * , void * , void * , bool &, bool &);

  /**
   * This is executed after the finalization of the solver
   */
  virtual void on_close();

private:
  SolverControl& _control;

  std::string _fname;

};

#endif
