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


#ifndef __gummel_monitor_hook_h__
#define __gummel_monitor_hook_h__


#include "config.h"
#ifdef COGENDA_COMMERCIAL_PRODUCT

#include "hook.h"
#include <ctime>
#include <vector>
#include <string>

class MeshBase;

#ifdef HAVE_VTK
class vtkUnstructuredGrid;
#endif

/**
 * write vtk file for each nonlinear iteration
 * for check the convergence histroy
 */
class GummelMonitorHook : public Hook
{

public:

  GummelMonitorHook(SolverBase & solver, const std::string & name, void *);

  virtual ~GummelMonitorHook();

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
   * This is executed after the finalization of the solver
   */
  virtual void on_close();

private:

  /**
   * monitor poisson solver
   */
  bool            _gummel_solver;

  /**
   * solution count
   */
  unsigned int     solution_count;

private:

  /**
   * monitor total carrier in the system
   */
  std::pair<double, double> carrier_count();

private:

  const MeshBase & mesh;

  // boundary info
  std::vector<unsigned int>       el;
  std::vector<unsigned short int> sl;
  std::vector<short int>          il;

  /**
   * export vtk after this begin time
   */
  double          _t_start;

  /**
   * export vtk before this end time
   */
  double          _t_end;

  /**
   * export vtk every this time step
   */
  double          _t_step;

  /**
   * last time we export vtk
   */
  double          _t_last;

  /**
   * flag to indicate vtk export
   */
  bool            _vtk_out;


#ifdef HAVE_VTK
  /**
   * write the nodes from the mesh into a vtkUnstructuredGrid
   */
  void nodes_to_vtk();

  /**
 * write the cells from the mesh into a vtkUnstructuredGrid
   */
  void cells_to_vtk();

  /**
   * write out mesh data to the VTK file, this might come in handy to display
   * material data and partition information
   */
  void meshinfo_to_vtk();

  /**
   * write mesh relative solution data to vtk file
   */
  void current_solution_to_vtk(void * _x, void * _dx);

  /**
   * write mesh relative solution data to vtk file
   */
  void poisson_correction_solution_to_vtk(void * _x, void * _dx);


  /**
   * write extra region name/material information and boundary information into xml vtk file
   */
  std::string export_extra_info();

  /**
   * aux function to write a scaler solution to vtkUnstructuredGrid
   */
  void write_node_scaler_solution(const std::vector<unsigned int> & order, std::vector<float> &sol, const std::string & sol_name);

  /**
   * aux function to write a cell based scaler solution to vtkUnstructuredGrid
   */
  void write_cell_scaler_solution(const std::vector<unsigned int> & order, std::vector<float> &sol, const std::string & sol_name);

  /**
   * pointer to the VTK grid
   */
  vtkUnstructuredGrid* grid;

  class XMLUnstructuredGridWriter;

#endif // #ifdef HAVE_VTK
};

#endif // #ifdef COGENDA_COMMERCIAL_PRODUCT

#endif // #ifndef __gummel_monitor_hook_h__
