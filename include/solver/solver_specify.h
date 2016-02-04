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

//  $Id: solver_specify.h,v 1.4 2008/07/09 12:25:19 gdiso Exp $

#ifndef __solver_specify_h__
#define __solver_specify_h__


#include <string>
#include <vector>
#include <map>
#include <deque>
#include <string>

#include "enum_petsc_type.h"
#include "enum_solver_specify.h"

#include "parser_parameter.h"


/**
 * this namespace stores information for controlling the solver.
 * genius fill these information by user's input deck and pass it to each solver.
 */
namespace SolverSpecify
{
  /**
   * enum value for specify which solver to be used, can be POISSON, DDML1 ..
   */
  extern SolverType     Solver;

  /**
   * enum value for specify which solution this solver will do. i.e. transient , steadystate...
   */
  extern SolutionType     Type;

  /**
   * label to identify this solve step
   */
  extern std::string      label;

  /**
   * the prefix string for output file
   */
  extern std::string      out_prefix;

  /**
   * the append model for output file
   */
  extern bool      out_append;

  /**
   * hooks to be installed, \<id \<hook_name, hook_parameters\> \>
   */
  extern std::map<std::string, std::pair<std::string, std::vector<Parser::Parameter> > > Hooks;

  /**
   * nonlinear solver scheme: basic, line search, trust region...
   */
  extern NonLinearSolverType     NS;

  /**
   * Determines when the preconditioner is rebuilt in the nonlinear solve
   */
  extern int     NSLagPCLU;


  /**
   * Determines when the Jacobian is rebuilt in the nonlinear solve
   */
  extern int     NSLagJacobian;

  /**
   * linear solver scheme: LU, BCGS, GMRES ...
   */
  extern LinearSolverType        LS;

  /**
   * preconditioner scheme: ASM ILU ...
   */
  extern PreconditionerType      PC;

  /**
   * Newton damping
   */
  extern DampingScheme     Damping;

  /**
   * indicate whether voronoi truncation will be used
   */
  extern VoronoiTruncationFlag VoronoiTruncation;


  //--------------------------------------------
  // half implicit method
  //--------------------------------------------

  /**
   * linear solver scheme: LU, BCGS, GMRES ... for carrier continuity equation
   */
  extern LinearSolverType        LS_CARRIER;

  /**
   * linear solver scheme: LU, BCGS, GMRES ... for full current equation
   */
  extern LinearSolverType        LS_CURRENT;

  /**
   * linear solver scheme: LU, BCGS, GMRES ... for poisson correction equation
   */
  extern LinearSolverType        LS_POISSON;

  /**
   * preconditioner scheme: ASM ILU ... for carrier continuity equation
   */
  extern PreconditionerType      PC_CARRIER;

  /**
   * preconditioner scheme: ASM ILU ... for full current equation
   */
  extern PreconditionerType      PC_CURRENT;

  /**
   * preconditioner scheme: ASM ILU ... for poisson correction equation
   */
  extern PreconditionerType      PC_POISSON;

  /**
   * linearize error threshold of half implicit method
   */
  extern double    LinearizeErrorThreshold;


  /**
   * steady state threshold (ratio to Vt) of half implicit method
   */
  extern double    SteadyStateThreshold;

  /**
   * solve carrier distribution again after poisson correction in half implicit method
   */
  extern bool      ReSolveCarrier;

  /**
   * allow artificial carrier generation,
   * which makes poisson's equation self-consistant
   */
  extern bool      ArtificialCarrier;

  /**
   * Poisson correction parameter
   */
  extern double    PoissonCorrectionParameter;

  /**
   * dump matrix and rhs vector
   */
  extern bool      DumpMatrixVector;

  //--------------------------------------------
  // nonlinear/linear solver convergence criteria
  //--------------------------------------------

 /**
  * snes relative error tolerance
  */
  extern double   snes_rtol;

  /**
   * relative error tolerance
   */
  extern double   ksp_rtol;

  /**
   * abs error tolerance
   */
  extern double   ksp_atol;

  /**
   * abs error tolerance is set as max(ksp_atol_fnorm*fnorm, ksp_atol)
   * where fnorm is the nonlinear function norm
   */
  extern double   ksp_atol_fnorm;


  /**
   * consider the system is singular
   */
  extern bool     ksp_singular;

  //--------------------------------------------
  // nonlinear solver convergence criteria
  //--------------------------------------------

  /**
   * Max nonlinear iteration number
   */
  extern unsigned int      MaxIteration;

  /**
   * potential damping factor, in kT/q
   */
  extern double   potential_update;


  /**
   * potential damping for spice
   */
  extern bool   damping_spice;

  /**
   * voltage damping factor for spice, in V
   */
  extern double   spice_voltage_update;

  /**
   * When absolute error of equation less
   * than this value, solution is considered converged.
   */
  extern double      absolute_toler;

  /**
   * When relative error of solution variable less
   * than this value, solution is considered converged.
   */
  extern double      relative_toler;

  /**
   * When relative error is used as converged criteria,
   * the equation norm should satisfy the absolute
   * converged criteria with a relaxation of
   * this value.
   */
  extern double      toler_relax;

  /**
   * The absolute converged criteria for the Poisson equation.
   */
  extern double      poisson_abs_toler;

  /**
   * The absolute converged criteria for the electron continuity equation.
   */
  extern double      elec_continuity_abs_toler;

  /**
   * The absolute converged criteria for the hole continuity equation.
   */
  extern double      hole_continuity_abs_toler;

  /**
   * The absolute converged criteria for the trap continuity equation.
   */
  extern double      trap_continuity_abs_toler;

  /**
   * The absolute converged criteria for the lattice heat equation equation.
   */
  extern double      heat_equation_abs_toler;

  /**
   * The absolute converged criteria for the electron energy balance equation.
   */
  extern double      elec_energy_abs_toler;

  /**
   * The absolute converged criteria for the hole energy balance equation.
   */
  extern double      hole_energy_abs_toler;

  /**
   * The absolute converged criteria for the electron quantum potential equation.
   */
  extern double      elec_quantum_abs_toler;

  /**
   * The absolute converged criteria for the hole quantum potential equation.
   */
  extern double      hole_quantum_abs_toler;

  /**
   * The absolute converged criteria for the electrode bias equation.
   */
  extern double      electrode_abs_toler;

  /**
   * The absolute converged criteria for the spice circuit.
   */
  extern double      spice_abs_toler;

  /**
   * The function value larger than divergence_factor*xx_abs_toler is considered as divergence
   */
  extern double      divergence_factor;

  //--------------------------------------------
  // TS (transient solver)
  //--------------------------------------------

  /**
   * TS indicator
   */
  extern bool      TimeDependent;

  /**
   * transient scheme
   */
  extern TemporalScheme    TS_type;

  /**
   * start time of transient simulation
   */
  extern double    TStart;

  /**
   * user defined time step of transient simulation
   * it is just a reference value
   */
  extern double    TStep;

  /**
   * the minimal time step. TStep will not exceed this value
   */
  extern double    TStepMin;

  /**
   * the maximum time step. TStep will not exceed this value
   */
  extern double    TStepMax;

  /**
   * stop time of transient simulation
   */
  extern double    TStop;

  /**
   * indicate if auto step control should be used
   */
  extern bool      AutoStep;

  /**
   * reject time steps that LTE not satisfied
   */
  extern bool      RejectStep;

  /**
   * indicate if predict of next solution value should be used
   */
  extern bool      Predict;

  /**
   * relative tol of TS truncate error, used in AutoStep
   */
  extern double    TS_rtol;

  /**
   * absolute tol of TS truncate error, used in AutoStep
   */
  extern double    TS_atol;

  /**
   * indicate BDF2 can be started.
   */
  extern bool      BDF2_LowerOrder;

  /**
   * use initial condition, only for mixA solver
   */
  extern bool      UIC;

  /**
   * do operator point calculation before transient simulation, only for mixA solver
   */
  extern bool      tran_op;

  /**
   * false indicate TR based on DC state
   * true indicate a previous TR simulation
   */
  extern bool      tran_histroy;

  /**
   * current time
   */
  extern double    clock;

  /**
   * current time step
   */
  extern double    dt;

  /**
   * last time step
   */
  extern double    dt_last;

  /**
   * previous time step
   */
  extern double    dt_last_last;

  /**
   *  the simulation cycles
   */
  extern int       T_Cycles;


  //------------------------------------------------------
  // parameters for DC and TRACE simulation
  //------------------------------------------------------

  /**
   * electrode(s) the voltage DC sweep will be performanced
   */
  extern std::vector<std::string>    Electrode_VScan;

  /**
   * The voltage of electrode(s) to do DC sweep
   */
  extern double         Electrode_VScan_Voltage;

  /**
   * start voltage of DC sweep
   */
  extern double    VStart;

  /**
   * voltage step
   */
  extern double    VStep;

  /**
   * max voltage step
   */
  extern double    VStepMax;

  /**
   * stop voltage of DC sweep
   */
  extern double    VStop;

  /**
   * electrode the current DC sweep will be performanced
   */
  extern std::vector<std::string>    Electrode_IScan;

  /**
   * The current of electrode(s) to do DC sweep
   */
  extern double         Electrode_IScan_Current;

  /**
   * start current of DC sweep
   */
  extern double    IStart;

  /**
   * current step
   */
  extern double    IStep;

  /**
   * max current step
   */
  extern double    IStepMax;

  /**
   * stop current
   */
  extern double    IStop;

  /**
   *  the simulation cycles
   */
  extern int       DC_Cycles;


  /**
   * use node set, only for mixA solver
   */
  extern bool      NodeSet;

  /**
   * ramp up the voltage/current sources in circuit, only for mixA solver
   */
  extern int       RampUpSteps;

  /**
   * the voltage step for ramp up
   */
  extern double    RampUpVStep;

  /**
   * the current step for ramp up
   */
  extern double    RampUpIStep;

  /**
   * the initial value of gmin
   */
  extern double    GminInit;

  /**
   * the final value of gmin
   */
  extern double    Gmin;

  /**
   * drive the system to steadystate
   */
  extern bool OpToSteady;


  //------------------------------------------------------
  // parameters for AC simulation
  //------------------------------------------------------
  /**
   * electrode for AC small signal sweep
   */
  extern std::vector<std::string>     Electrode_ACScan;

  /**
   * the amplitude of small signal for AC sweep
   */
  extern double    VAC;

  /**
   * start frequency
   */
  extern double    FStart;

  /**
   * frequency multiple factor
   */
  extern double    FMultiple;

  /**
   * stop frequency
   */
  extern double    FStop;

  /**
   * current frequency
   */
  extern double    Freq;

  //------------------------------------------------------
  // parameters for pseudo time stepping method
  //------------------------------------------------------

  /**
   * if pseudo time stepping method enabled
   */
  extern bool PseudoTimeMethod;

  /**
   * if pseudo time method for CMOS enabled
   * in this mode, sigma of metal region will be modified
   */
  extern bool PseudoTimeCMOS;

  /**
   * characteristic length of CMOS device, default value is 0.1um
   */
  extern double PseudoTimeCMOSLambda;

  /**
   * characteristic resistance of CMOS device, default value is 1K
   */
  extern double PseudoTimeCMOSRes;

  /**
   * characteristic capatance of CMOS device, default value is 1f
   */
  extern double PseudoTimeCMOSCap;

  /**
   * characteristic time of CMOS device, default value is 0.1ns
   */
  extern double PseudoTimeCMOSTime;

  /**
   * pseudo time step for carrier in semiconductor region
   */
  extern double PseudoTimeStepCarrier;

  /**
   * pseudo time step for potential in semiconductor region
   */
  extern double PseudoTimeStepPotential;

  /**
   * pseudo time step for metal region
   */
  extern double PseudoTimeStepMetal;

  /**
   * maximum pseudo time step
   */
  extern double PseudoTimeStepMax;

  /**
   * relative X convergence criteria in pseudo time step mode
   */
  extern double PseudoTimeMethodRXTol;

  /**
   * relative function convergence criteria in pseudo time step mode
   */
  extern double PseudoTimeMethodRFTol;

  /**
   * convergence relax to absolute tol in pseudo time step mode
   */
  extern double PseudoTimeTolRelax;

  /**
   * max pseudo time step
   */
  extern int PseudoTimeSteps;

  //------------------------------------------------------
  // parameters for optical / particle effect
  //------------------------------------------------------

  /**
   * when OptG is true,
   * the optical carrier generation is considered in the simulation
   */
  extern bool      OptG;

  /**
   * when PatG is true,
   * the particle carrier generation is considered in the simulation
   */
  extern bool      PatG;


  /**
   * coupled source effect
   */
  extern bool      SourceCoupled;


  /**
   * TID total dose
   */
  extern double    TID_TotalDose;


  /**
   * TID dose rate
   */
  extern double    TID_DoseRate;

  /**
   * TID dose step
   */
  extern double    TID_DoseStep;

  /**
   * TID OP step
   */
  extern double    TID_OPStep;

  /**
   * after TID solver, only fixed charge remains, carriers in oxide are set to zero
   */
  extern bool      TID_FixedCharge;

  //------------------------------------------------------
  // initial parameters
  //------------------------------------------------------

  /**
   * constructor, set default values
   */
  extern void set_default_parameter();



}



#endif //#define __solver_specify_h__
