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


#ifndef __mix_solver_h__
#define __mix_solver_h__




#include "fvm_nonlinear_solver.h"

#ifndef CYGWIN

//data interface with ngspice
#include "ndevexch.h"

/* network function */
#include <errno.h>
#include <netinet/in.h>  /* IPv4 socket address structres. */
#include <netdb.h>       /* Access to DNS lookups. */
#include <arpa/inet.h>   /* inet_ntop function. */
#include <sys/socket.h>  /* Socket functions. */
#include <setjmp.h>
#include <signal.h>

#ifdef CYGWIN
#define MSG_FLAG 0
#else
#define MSG_FLAG MSG_WAITALL
#endif

typedef struct
{
  PetscScalar pdI_pdV;
  Vec         pdI_pdw;
  Vec         pdF_pdV;
  Vec         pdw_pdV;
} sPINcond;



/**
 * the base class for device/circuit mixed type simulation
 */
class MixSolverBase : public FVM_NonlinearSolver
{
public:
  /**
   * constructor
   */
  MixSolverBase(SimulationSystem & system): FVM_NonlinearSolver(system)
  { port = SolverSpecify::ServerPort;  }

  /**
   * destructor
   */
  ~MixSolverBase()
  {}

  /**
   * virtual function, create the solver
   */
  virtual int create_solver();

  /**
   * do post-process after each solve action
   */
  virtual int post_solve_process();

  /**
   * virtual function, destroy the solver
   */
  virtual int destroy_solver();

  /**
   * setup network connection to ngspice with TCP port \p port
   */
  virtual int connect_to_ngspice();

  /**
   * do simulation under the control of ngspice
   */
  virtual int run_under_ngspice();

  /**
   * free the connection to ngspice
   */
  virtual int disconnet_to_ngspice();

  /**
   * get pdI/pdw, pdI/pdV and pdF/pdV for each electrode
   */
  virtual int get_electrode_load()=0;

  /**
   * load previous state into solution vector, empty here
   */
  virtual int diverged_recovery()=0;

  /**
   * Save steady-state solution as the previous solution data in spice
   */
  virtual int init_spice_data()=0;

  /**
   * load solution data previously accepted by spice
   */
  virtual int load_spice_data()=0;

  /**
   * since spice accept the solution, save solution data
   */
  virtual int save_spice_data()=0;

  /**
   * get the pin index by its label
   */
  int get_pin_index(const std::string & label)
  { return pin_name_to_pin_index_map[label]; }

  /**
   * check carrier density after projection
   */
  virtual void projection_positive_density_check(Vec x, Vec xo)=0;

  /**
   * compute the abs and relative error norm of the solution
   * each derived DDM solver should override it.
   */
  virtual void error_norm()=0;

  /**
   * snes monitor, do nothing
   */
  virtual void petsc_snes_monitor(PetscInt , PetscReal ) {}

  /**
   * sens convergence criteria for dd solvers
   */
  virtual void petsc_snes_convergence_test(PetscInt its, PetscReal xnorm, PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason);

  /**
   * ksp convergence criteria for dd solvers
   */
  virtual void petsc_ksp_convergence_test(PetscInt its, PetscReal rnorm, KSPConvergedReason* reason);

protected:

  /**
   * the global solution vector accepted by spice at n-th step
   */
  Vec            xs_n;

  /**
   * the global solution vector accepted by spice at n-1 th step
   */
  Vec            xs_n1;

  /**
   * time step from n-1 to n th step
   */
  PetscScalar    spice_dt_last;

  /**
   * My listening socket.
   */
  int listener;

  /**
   * my socket
   */
  int client;

  /**
   * my TCP port number
   */
  int port;

  /**
   * set TCP/IP socket and connect to NGSPICE.
   */
  int NDEVServer();

  /**
   * process ngspice DEV_LOAD call
   */
  int DEV_LOAD();

  /**
   * process ngspice DEV_ACCEPT call
   */
  int DEV_ACCEPT();

  /**
   * process ngspice DEV_CONV_TEST call
   */
  int DEV_CONV_TEST();

  /**
   * do transient solve
   */
  int tran_solve();

  /**
   * do dc solve
   */
  int dc_solve();

  /**
   * circuit information passed from ngspice
   */
  sCKTinfo      CKTInfo;

  /**
   * the numerical device information passed from ngspice
   */
  sDeviceinfo   Deviceinfo;

  /**
   * the pin voltage and conduction matrix enties. will be send back to ngspice.
   * max allowed 7 pins at present.
   */
  sPINinfo      PINinfos[7];

  /**
   * structure for calculate pin conductance matrix entries.
   */
  sPINcond      PINconds[7];

  /**
   * map pin name to the index
   */
  std::map<std::string, int> pin_name_to_pin_index_map;

  /**
   * special linear solver for sloving electrode conductance.
   * since this equ should be solved many times with the same matrix and different rhs
   * we should use LU method here -- just one time LU decompose needed.
   */
  KSP  linear_lu;

  /**
   * use direct method
   */
  PC   pc_direct;

  /**
   * matrix for sloving electrode conductance
   */
  Mat  G;

  /**
   * incidate that the g_matrix is never assembled
   */
  bool g_matrix_first_assemble;

  /**
   * x norm of potential
   */
  PetscScalar potential_norm;

  /**
   * x norm of electron density
   */
  PetscScalar electron_norm;

  /**
   * x norm of hole density
   */
  PetscScalar hole_norm;

  /**
   * x norm of lattice temperature
   */
  PetscScalar temperature_norm;

  /**
   * x norm of electron temperature
   */
  PetscScalar elec_temperature_norm;

  /**
   * x norm of hole temperature
   */
  PetscScalar hole_temperature_norm;

  /**
   * f norm of poisson's equation
   */
  PetscScalar poisson_norm;

  /**
   * f norm of electron continuity equation
   */
  PetscScalar elec_continuity_norm;

  /**
   * f norm of hole continuity equation
   */
  PetscScalar hole_continuity_norm;

  /**
   * f norm of lattice temperature equation
   */
  PetscScalar heat_equation_norm;

  /**
   * f norm of electron energy balance equation
   */
  PetscScalar elec_energy_equation_norm;

  /**
   * f norm of hole energy balance equation
   */
  PetscScalar hole_energy_equation_norm;

  /**
   * nonlinear function norm
   */
  PetscScalar function_norm;

};

#endif //#ifndef CYGWIN

#endif //#define __ddm_solver_h__
