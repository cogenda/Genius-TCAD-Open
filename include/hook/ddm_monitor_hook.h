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


#ifndef __ddm_monitor_hook_h__
#define __ddm_monitor_hook_h__


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
class DDMMonitorHook : public Hook
{

public:

  DDMMonitorHook(SolverBase & solver, const std::string & name, void *);

  virtual ~DDMMonitorHook();

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
   *  This is executed before each (nonlinear) iteration
   */
  virtual void pre_iteration();

  /**
   *  This is executed after each (nonlinear) iteration
   */
  virtual void post_check(void * f, void * x, void * dx, void * w, bool & change_y, bool &change_w);

  /**
   * This is executed after the finalization of the solver
   */
  virtual void on_close();

private:

  /**
   * monitor poisson solver
   */
  bool            _poisson_solver;

  /**
   * monitor ddm solver
   */
  bool            _ddm_solver;

  /**
   * iteration count
   */
  unsigned int iteration_count;

  /**
   * solution count
   */
  unsigned int solution_count;

private:

  const MeshBase & mesh;

  // boundary info
  std::vector<unsigned int>       el;
  std::vector<unsigned short int> sl;
  std::vector<short int>          il;

  double _t_start;
  double _t_end;

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
  void solution_to_vtk_ddm(void * _f, void * _x, void * _dx, void * _w);

  /**
   * write mesh relative solution data to vtk file
   */
  void solution_to_vtk_poisson(void * _f, void * _x, void * _dx, void * _w);

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
#endif
};

#endif
