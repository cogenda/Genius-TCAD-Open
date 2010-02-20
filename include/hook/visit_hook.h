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



#ifndef __visit_hook_h__
#define __visit_hook_h__



#include "hook.h"

// visit include
#include <VisItControlInterface_V1.h>
#include <VisItDataInterface_V1.h>


enum VisitState
{
  UNKNOWN=0,
  CONNECTED,
  DISCONNECTED
};


/**
 * Use LLNL Visit to do real time display
 * NOTE since we can only connect to one visit process, we set many member static here.
 */
class VisitHook : public Hook
{

public:

  /**
   * constructor.
   */
  VisitHook(SolverBase & solver, const std::string & name, void * dummy);

  /**
   * destructor.
   */
  virtual ~VisitHook();


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
   *  i.e. for collecting convergence information or implementing various damping strategy
   */
  virtual void post_iteration();

  /**
   * This is executed after the finalization of the solver
   */
  virtual void on_close();

public:

  void connect_to_visit();

  void comm_to_visit();

  void disconnect_to_visit();

public:

  static VisIt_SimulationMetaData * get_meta_data(void);

  /**
   *
   */
  static VisIt_CurveData * get_electrode_iv(const std::string &name);

  static VisIt_MeshData *  get_mesh(int domain, const std::string &name);

  static VisIt_ScalarData * get_scalar(int domain, const std::string &name);

private:

  /**
   * the visit sim state
   */
  static VisitState visit_state;

  /*
   * static pointer to solver
   */
  static SolverBase  *  _p_solver;

  /**
  * the variable name buffer
  */
 static std::vector<std::pair<std::string, std::string> >  _variables;

 /**
  * the variable value buffer
  */
 static std::vector< std::vector<double> > _values;

 /**
  * the total number of values
  */
 static unsigned int _n_values;

 static std::map<unsigned int, int> global_id_to_index;

 /**
  * get the cell type that can be recognized by visit
  */
 static int visit_cell_type(const Elem *);


};



#endif //#define __visit_hook_h__
