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


#ifndef __mob_monitor_hook_h__
#define __mob_monitor_hook_h__


#include "hook.h"
#include <ctime>
#include <vector>
#include <map>
#include <string>

class MeshBase;

#ifdef HAVE_VTK
class vtkUnstructuredGrid;
#endif

/**
 * write vtk file for each nonlinear iteration
 * for check the convergence histroy
 */
class MobMonitorHook : public Hook
{

public:

  MobMonitorHook(SolverBase & solver, const std::string & name, void *);

  virtual ~MobMonitorHook();

  /**
   *   This is executed before the initialization of the solver
   */
  virtual void on_init();

  /**
   *  This is executed after each solution step.
   */
  virtual void post_solve();

  /**
   * This is executed after the finalization of the solver
   */
  virtual void on_close();

private:

  /**
   * monitor ddm solver
   */
  bool _ddm_solver;

  /**
   * solution count
   */
  unsigned int solution_count;

private:

  const MeshBase & mesh;

  double _t_start;
  double _t_end;

  bool   _surface;

#ifdef HAVE_VTK

  /**
   * mesh nodes
   */
  std::vector<double> points;

  /**
   * mesh edges
   */
  std::map< std::pair<unsigned int, unsigned int>, unsigned int > edges;

  /**
   * write the nodes from the mesh into a vtkUnstructuredGrid
   */
  void nodes_to_vtk();

  /**
   * write the edges from the mesh into a vtkUnstructuredGrid
   */
  void edges_to_vtk();

  /**
   * write mesh relative solution data to vtk file
   */
  void solution_to_vtk();

  /**
   * write node id to vtk file
   */
  void node_id_to_vtk();

  /**
   * aux function to write a cell based scaler solution to vtkUnstructuredGrid
   */
  void write_cell_scaler_solution( const std::vector<float> &sol, const std::string & sol_name);

  /**
   * pointer to the VTK grid
   */
  vtkUnstructuredGrid* grid;

#endif
};

#endif
