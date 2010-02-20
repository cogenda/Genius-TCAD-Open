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

//  $Id: solver_specify.cc,v 1.5 2008/07/09 12:25:19 gdiso Exp $


#include <deque>
#include <string>
#include <vector>
#include <string>

#include "solver_specify.h"


/**
 * this namespace stores information for controlling the solver.
 * genius fill these information by user's input deck and pass it to each solver.
 */
namespace SolverSpecify
{
  /**
   * enum value for specify which solver to be used, can be POISSON, DDML1 ..
   */
  SolverType     Solver;


  /**
   * enum value for specify which solution this solver will do. i.e. transient , steadystate...
   */
  SolutionType     Type;

  /**
   * label to identify this solve step
   */
  std::string      label;

  /**
   * the prefix string for output file
   */
  std::string      out_prefix;

  /**
   * hooks to be installed
   */
  std::deque<std::string> Hooks;

  /**
   * parameters for the hooks
   */
  std::map<std::string, std::vector<Parser::Parameter> > Hook_Parameters;

  /**
   * nonlinear solver scheme: basic, line search, trust region...
   */
  NonLinearSolverType     NS;

  /**
   * linear solver scheme: LU, BCGS, GMRES ...
   */
  LinearSolverType        LS;

  /**
   * preconditioner scheme: ASM ILU ...
   */
  PreconditionerType      PC;

  /**
   * Newton damping
   */
  DampingScheme     Damping;


  //--------------------------------------------
  // linear solver convergence criteria
  //--------------------------------------------

  /**
   * relative error tolerance
   */
  double   ksp_rtol;

  /**
   * abs error tolerance
   */
  double   ksp_atol;

  /**
   * abs error tolerance is set as max(ksp_atol_fnorm*fnorm, ksp_atol)
   * where fnorm is the nonlinear function norm
   */
  double   ksp_atol_fnorm;

  //--------------------------------------------
  // nonlinear solver convergence criteria
  //--------------------------------------------

  /**
   * Max nonlinear iteration number
   */
  unsigned int      MaxIteration;

  /**
   * When relative error of solution variable less
   * than this value, solution is considered converged.
   */
  double      relative_toler;

  /**
   * When relative error is used as converged criteria,
   * the equation norm should satisfy the absolute
   * converged criteria with a relaxation of
   * this value.
   */
  double      toler_relax;

  /**
   * The absolute converged criteria for the Poisson equation.
   */
  double      poisson_abs_toler;

  /**
   * The absolute converged criteria for the electron continuity equation.
   */
  double      elec_continuity_abs_toler;

  /**
   * The absolute converged criteria for the hole continuity equation.
   */
  double      hole_continuity_abs_toler;

  /**
   * The absolute converged criteria for the lattice heat equation equation.
   */
  double      heat_equation_abs_toler;

  /**
   * The absolute converged criteria for the electron energy balance equation.
   */
  double      elec_energy_abs_toler;

  /**
   * The absolute converged criteria for the hole energy balance equation.
   */
  double      hole_energy_abs_toler;

  /**
   * The absolute converged criteria for the electron quantum potential equation.
   */
  double      elec_quantum_abs_toler;

  /**
   * The absolute converged criteria for the hole quantum potential equation.
   */
  double      hole_quantum_abs_toler;

  /**
   * The absolute converged criteria for the electrode bias equation.
   */
  double      electrode_abs_toler;

  //--------------------------------------------
  // TS (transient solver)
  //--------------------------------------------

  /**
   * TS indicator
   */
  bool      TimeDependent;

  /**
   * transient scheme
   */
  TSType    TS_type;

  /**
   * start time of transient simulation
   */
  double    TStart;

  /**
   * user defined time step of transient simulation
   * it is just a reference value
   */
  double    TStep;

  /**
   * the maximum time step. TStep will not exceed this value
   */
  double    TStepMax;

  /**
   * stop time of transient simulation
   */
  double    TStop;

  /**
   * indicate if auto step control should be used
   */
  bool      AutoStep;

  /**
   * indicate if predict of next solution value should be used
   */
  bool      Predict;

  /**
   * relative tol of TS truncate error, used in AutoStep
   */
  double    TS_rtol;

  /**
   * absolute tol of TS truncate error, used in AutoStep
   */
  double    TS_atol;

  /**
   * indicate BDF2 can be started.
   */
  bool      BDF2_restart;

  /**
   * use initial condition, only for mixA solver
   */
  bool      UIC;

  /**
   * current time
   */
  double    clock;

  /**
   * current time step
   */
  double    dt;

  /**
   * last time step
   */
  double    dt_last;

  /**
   * previous time step
   */
  double    dt_last_last;

  /**
   *  the simulation cycles
   */
  int       T_Cycles;

  //------------------------------------------------------
  // parameters for DC and TRACE simulation
  //------------------------------------------------------

  /**
   * electrode(s) the voltage DC sweep will be performanced
   */
  std::vector<std::string>    Electrode_VScan;

  /**
   * start voltage of DC sweep
   */
  double    VStart;

  /**
   * voltage step
   */
  double    VStep;

  /**
   * max voltage step
   */
  double    VStepMax;

  /**
   * stop voltage of DC sweep
   */
  double    VStop;

  /**
   * electrode the current DC sweep will be performanced
   */
  std::vector<std::string>    Electrode_IScan;

  /**
   * start current of DC sweep
   */
  double    IStart;

  /**
   * current step
   */
  double    IStep;

  /**
   * max current step
   */
  double    IStepMax;

  /**
   * stop current
   */
  double    IStop;

  /**
   *  the simulation cycles
   */
  int       DC_Cycles;

  /**
   * ramp up the voltage/current sources in circuit, only for mixA solver
   */
  int       RampUpSteps;

  /**
   * the voltage step for ramp up
   */
  double    RampUpVStep;

  /**
   * the current step for ramp up
   */
  double    RampUpIStep;

  //------------------------------------------------------
  // parameters for MIX simulation
  //------------------------------------------------------

  /**
   * TCP port number
   */
  unsigned short int ServerPort;

  //------------------------------------------------------
  // parameters for AC simulation
  //------------------------------------------------------
  /**
   * electrode for AC small signal sweep
   */
  std::vector<std::string>    Electrode_ACScan;

  /**
   * the amplitude of small signal for AC sweep
   */
  double    VAC;

  /**
   * start frequency
   */
  double    FStart;

  /**
   * frequency multiple factor
   */
  double    FMultiple;

  /**
   * stop frequency
   */
  double    FStop;

  /**
   * current frequency
   */
  double    Freq;


  //------------------------------------------------------
  // parameters for MIX simulation
  //------------------------------------------------------

  /**
   * when OptG is true,
   * the optical carrier generation is considered in the simulation
   */
  bool      OptG;

  /**
   * when PatG is true,
   * the particle carrier generation is considered in the simulation
   */
  bool      PatG;

  /**
   * set default values
   */
  void set_default_parameter()
  {
    Solver       = DDML1;
    Type         = EQUILIBRIUM;
    NS           = LineSearch;
    LS           = BICGSTAB;
    PC           = ASM_PRECOND;
    Damping      = DampingNo;
    MaxIteration = 20;

    TimeDependent = false;
    TS_type       = BDF2;
    BDF2_restart  = true;
    UIC           = false;
    AutoStep      = true;
    Predict       = true;
    clock         = 0.0;

    RampUpSteps   = 1;
    RampUpVStep   = 0.25;
    RampUpIStep   = 0.1;

    VAC           = 0.0;
    OptG          = false;
    PatG          = false;
  }

  SolutionType type_string_to_enum(const std::string s)
  {
    if (s == "equilibrium")                   return EQUILIBRIUM;
    if (s == "steadystate")                   return STEADYSTATE;
    if (s == "op" )                           return OP;
    if (s == "dcsweep" )                      return DCSWEEP;
    if (s == "trace" )                        return TRACE;
    if (s == "acsweep")                       return ACSWEEP;
    if (s == "transient")                     return TRANSIENT;

    return INVALID_SolutionType;
  }

}

