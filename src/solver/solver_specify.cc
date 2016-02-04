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

#include <limits>
#include <deque>
#include <string>
#include <vector>
#include <string>

#include "solver_specify.h"
#include "physical_unit.h"

using PhysicalUnit::A;
using PhysicalUnit::V;
using PhysicalUnit::W;
using PhysicalUnit::C;
using PhysicalUnit::s;
using PhysicalUnit::um;
using PhysicalUnit::rad;

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
   * the append model for output file
   */
  bool      out_append;

  /**
   * hooks to be installed \<id \<hook_name, hook_parameters\> \>
   */
  std::map<std::string, std::pair<std::string, std::vector<Parser::Parameter> > > Hooks;

  /**
   * nonlinear solver scheme: basic, line search, trust region...
   */
  NonLinearSolverType     NS;

  /**
   * Determines when the preconditioner is rebuilt in the nonlinear solve
   */
  int     NSLagPCLU;

  /**
   * Determines when the Jacobian is rebuilt in the nonlinear solve
   */
  int     NSLagJacobian;

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


  /**
   * indicate whether voronoi truncation will be used
   */
  VoronoiTruncationFlag VoronoiTruncation;

  //--------------------------------------------
  // half implicit method
  //--------------------------------------------

  /**
   * linear solver scheme: LU, BCGS, GMRES ... for carrier continuity equation
   */
  LinearSolverType        LS_CARRIER;

  /**
   * linear solver scheme: LU, BCGS, GMRES ... for full current equation
   */
  LinearSolverType        LS_CURRENT;

  /**
   * linear solver scheme: LU, BCGS, GMRES ... for poisson correction equation
   */
  LinearSolverType        LS_POISSON;

  /**
   * preconditioner scheme: ASM ILU ... for carrier continuity equation
   */
  PreconditionerType      PC_CARRIER;

  /**
   * preconditioner scheme: ASM ILU ... for full current equation
   */
  PreconditionerType      PC_CURRENT;

  /**
   * preconditioner scheme: ASM ILU ... for poisson correction equation
   */
  PreconditionerType      PC_POISSON;

  /**
   * linearize error threshold of half implicit method
   */
  double    LinearizeErrorThreshold;


  /**
   * steady state threshold (ratio to Vt) of half implicit method
   */
  double    SteadyStateThreshold;

  /**
   * solve carrier distribution again after poisson correction in half implicit method
   */
  bool      ReSolveCarrier;

  /**
   * allow artificial carrier generation,
   * which makes poisson's equation self-consistant
   */
  bool      ArtificialCarrier;

  /**
   * Poisson correction parameter
   */
  double    PoissonCorrectionParameter;


  /**
   * dump matrix and rhs vector
   */
  bool      DumpMatrixVector;


  //--------------------------------------------
  // nonlinear/linear solver convergence criteria
  //--------------------------------------------

  /**
   * snes relative error tolerance
   */
  double   snes_rtol;


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

  /**
   * consider the system is singular
   */
  bool     ksp_singular;

  //--------------------------------------------
  // nonlinear solver convergence criteria
  //--------------------------------------------

  /**
   * Max nonlinear iteration number
   */
  unsigned int      MaxIteration;

  /**
   * potential damping factor, in kT/q
   */
  double   potential_update;

  /**
   * potential damping for spice
   */
  bool   damping_spice;

  /**
   * voltage damping factor for spice, in V
   */
  double   spice_voltage_update;

  /**
   * When absolute error of equation less
   * than this value, solution is considered converged.
   */
  double      absolute_toler;

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
   * The absolute converged criteria for the trap continuity equation.
   */
  double      trap_continuity_abs_toler;

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

  /**
   * The absolute converged criteria for the spice circuit.
   */
  double      spice_abs_toler;

  /**
   * The function value larger than divergence_factor*xx_abs_toler is considered as divergence
   */
  double      divergence_factor;

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
  TemporalScheme    TS_type;

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
   * the minimal time step. TStep will not exceed this value
   */
  double    TStepMin;

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
   * reject time steps that LTE not satisfied
   */
  bool      RejectStep;

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
  bool      BDF2_LowerOrder;

  /**
   * use initial condition, only for mixA solver
   */
  bool      UIC;

  /**
   * do operator point calculation before transient simulation, only for mixA solver
   */
  bool      tran_op;

  /**
   * false indicate TR based on DC state
   * true indicate a previous TR simulation
   */
  bool      tran_histroy;

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
   * The voltage of electrode(s) to do DC sweep
   */
  double         Electrode_VScan_Voltage;

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
   * The current of electrode(s) to do DC sweep
   */
  double         Electrode_IScan_Current;


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
   * use node set, only for mixA solver
   */
  bool      NodeSet;

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


  /**
   * the initial value of gmin
   */
  double    GminInit;

  /**
   * the final value of gmin
   */
  double    Gmin;

  /**
   * drive the system to steadystate
   */
  bool     OpToSteady;


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
  // parameters for pseudo time stepping method
  //------------------------------------------------------

  /**
   * if pseudo time stepping method enabled
   */
  bool PseudoTimeMethod;

  /**
   * if pseudo time method for CMOS enabled
   * in this mode, sigma of metal region will be modified
   */
  bool PseudoTimeCMOS;

  /**
   * characteristic length of CMOS device, default value is 0.1um
   */
  double PseudoTimeCMOSLambda;

  /**
   * characteristic resistance of CMOS device, default value is 1K
   */
  double PseudoTimeCMOSRes;

  /**
   * characteristic capatance of CMOS device, default value is 1f
   */
  double PseudoTimeCMOSCap;

  /**
   * characteristic time of CMOS device, default value is 0.1ns
   */
  double PseudoTimeCMOSTime;

  /**
   * pseudo time step for carrier in semiconductor region
   */
  double PseudoTimeStepCarrier;

  /**
   * pseudo time step for potential in semiconductor region
   */
  double PseudoTimeStepPotential;

  /**
   * pseudo time step for metal region
   */
  double PseudoTimeStepMetal;

  /**
   * maximum pseudo time step
   */
  double PseudoTimeStepMax;

  /**
   * relative X convergence criteria in pseudo time step mode
   */
  double PseudoTimeMethodRXTol;

  /**
   * relative function convergence criteria in pseudo time step mode
   */
  double PseudoTimeMethodRFTol;

  /**
   * convergence relax to absolute tol in pseudo time step mode
   */
  double PseudoTimeTolRelax;

  /**
   * max pseudo time step
   */
  int PseudoTimeSteps;


  //------------------------------------------------------
  // parameters for optical / particle effect
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
   * coupled source effect
   */
  bool      SourceCoupled;


  /**
   * TID total dose
   */
  double    TID_TotalDose;


  /**
   * TID dose rate
   */
  double    TID_DoseRate;

  /**
   * TID dose step
   */
  double    TID_DoseStep;

  /**
   * TID OP step
   */
  double    TID_OPStep;

  /**
   * after TID solver, only fixed charge remains, carriers in oxide are set to zero
   */
  bool      TID_FixedCharge;

  //------------------------------------------------------
  // initial parameters
  //------------------------------------------------------

  /**
   * set default values
   */
  void set_default_parameter()
  {
    Solver            = DDML1;
    Type              = EQUILIBRIUM;
    NS                = LineSearch;
#ifdef PETSC_HAVE_MUMPS
    LS                = GMRES;
    PC                = LU_PRECOND;
    NSLagPCLU         = 5;
    NSLagJacobian     = 1;
#else
    LS                = BCGSL;
    PC                = ASM_PRECOND;
    NSLagPCLU         = 1;
    NSLagJacobian     = 1;
#endif

    out_append        = false;

    Damping           = DampingPotential;
    VoronoiTruncation = VoronoiTruncationAlways;

    LS_POISSON        = GMRES;
    PC_POISSON        = ASM_PRECOND;
    LinearizeErrorThreshold   = 1.0;
    SteadyStateThreshold      = 1e-5;
    ReSolveCarrier    = false;
    ArtificialCarrier = true;
    PoissonCorrectionParameter= 0.0;
    DumpMatrixVector  = false;

    MaxIteration              = 30;
    potential_update          = 1.0;

    damping_spice             = false;
    spice_voltage_update      = 10;

    snes_rtol                 = 1e-5;

    ksp_rtol                  = 1e-8;
    ksp_atol                  = 1e-15;
    ksp_atol_fnorm            = 1e-7;
    ksp_singular              = false;

    absolute_toler            = 1e-12;
    relative_toler            = 1e-5;
    toler_relax               = 1e5;
    poisson_abs_toler         = 1e-26*C;
    elec_continuity_abs_toler = 5e-18*A;
    hole_continuity_abs_toler = 5e-18*A;
    trap_continuity_abs_toler = 5e-18*A;
    heat_equation_abs_toler   = 1e-11*W;
    elec_energy_abs_toler     = 1e-18*W;
    hole_energy_abs_toler     = 1e-18*W;
    electrode_abs_toler       = 1e-14*A;
    spice_abs_toler           = 1e-12*A;
    elec_quantum_abs_toler    = 1e-26*C;
    hole_quantum_abs_toler    = 1e-26*C;
    divergence_factor         = 1e20;

    TimeDependent             = false;
    TStepMin                  = 1e-14*s;
    TS_type                   = BDF2;
    BDF2_LowerOrder           = true;
    UIC                       = false;
    tran_op                   = true;
    tran_histroy              = false;
    AutoStep                  = true;
    RejectStep                = true;
    Predict                   = true;
    TS_rtol                   = 1e-3;
    TS_atol                   = 1e-7;
    clock                     = 0.0;
    dt                        = 1e100;


    VStepMax          = 1.0;
    IStepMax          = 1.0;

    Electrode_VScan_Voltage = 0.0;
    Electrode_IScan_Current = 0.0;

    NodeSet           = true;
    RampUpSteps       = 0;
    RampUpVStep       = std::numeric_limits<double>::infinity();
    RampUpIStep       = std::numeric_limits<double>::infinity();

    GminInit          = 1e-12;
    Gmin              = 1e-12;

    VAC               = 0.0;

    OpToSteady        = true;

    OptG              = false;
    PatG              = false;
    SourceCoupled     = false;


    TID_TotalDose     = 0;
    TID_DoseRate      = 0;
    TID_DoseStep      = 1000*rad;
    TID_OPStep        = 1000*rad;
    TID_FixedCharge   = false;

    PseudoTimeMethod            = false;
    PseudoTimeCMOS              = true;
    PseudoTimeCMOSLambda        = 0.1*um;
    PseudoTimeCMOSRes           = 1e3*V/A;
    PseudoTimeCMOSCap           = 1e-15*C/V;
    PseudoTimeCMOSTime          = 1e-10*s;
    PseudoTimeStepPotential     = 1e-6*s;
    PseudoTimeStepCarrier       = 1e-8*s;
    PseudoTimeStepMetal         = 1e-10*s;
    PseudoTimeStepMax           = 1e-6*s;
    PseudoTimeMethodRXTol       = 1e-2;
    PseudoTimeMethodRFTol       = 1e-2;
    PseudoTimeTolRelax          = 1e8;
    PseudoTimeSteps             = 50;
  }

  SolutionType solution_type_string_to_enum(const std::string s)
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

  SolverType solver_type_string_to_enum(const std::string &s)
  {
    if (s == "poisson")            return POISSON;
    if (s == "ddml1")              return DDML1;
    if (s == "ddml1m")             return DDML1MIXA;
    if (s == "ddml1ms")            return DDML1MIX;
    if (s == "hall")               return HALLDDML1;
    if (s == "ddml2")              return DDML2;
    if (s == "ddml2m")             return DDML2MIXA;
    if (s == "ddml2ms")            return DDML2MIX;
    if (s == "ebml3")              return EBML3;
    if (s == "ebml3m")             return EBML3MIXA;
    if (s == "ebml3ms")            return EBML3MIX;
    if (s == "ddmac")              return DDMAC;
    if (s == "qddm")               return DENSITY_GRADIENT;
    if (s == "halfimplicit")       return HALF_IMPLICIT;
    if (s == "ric")                return RIC;
    if (s == "dictat")             return DICTAT;

    return INVALID_SOLVER;
  }
}

