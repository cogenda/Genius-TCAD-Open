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

//  $Id: enum_solver_specify.h,v 1.7 2008/07/09 05:58:15 gdiso Exp $



#ifndef __enum_solver_specify_h__
#define __enum_solver_specify_h__



// ------------------------------------------------------------
// enum Solver Specify

namespace SolverSpecify
{

  /**
   * define solver type
   */
  enum SolverType
  {
    POISSON = 0,    // nonlinear poisson solver
    HDM,            // hydrodynamic solver
    DDML1,          // level 1 drift-diffusion solver
    DDML1MIX,       // mixed type level 1 drift-diffusion solver (old)
    DDML1MIXA,      // mixed type level 1 drift-diffusion solver
    HALLDDML1,      // level 1 drift-diffusion solver with hall effect
    DDML1R,         // level 1 drift-diffusion solver with electron PDE in resistive metal region
    DDML2,          // level 2 drift-diffusion solver
    DDML2MIX,       // mixed type level 2 drift-diffusion solver (old)
    DDML2MIXA,      // mixed type level 2 drift-diffusion solver
    EBML3,          // level 3 energy balance solver
    EBML3MIX,       // mixed type level 3 energy balance solver (old)
    EBML3MIXA,      // mixed type level 3 energy balance solver
    GUMMEL,         // gummel steady-state solver
    HALF_IMPLICIT,  // half-implicit transient solver
    CARRIER_TRANSIENT,    // gummel carrier solver
    CARRIER_TRANSIENT_IR, // gummel carrier solver with implicit recombination term
    MOCK_CURRENT,   // Mock's full current continuity solver (and its correction by Polsky)
    POLSKY_POISSON, // Polsky's poisson correction
    DDMAC,          // ac sweep solver
    DENSITY_GRADIENT,
    MC,
    RAY_TRACE,      // ray tracing optical solver
    EM_FEM_2D,
    EM_FEM_3D,
    FDTD,
    FVTD,
    DOPING_ANALYTIC,  // doping distribution solver
    MOLE_ANALYTIC,    // mole distribution solver
    OPTG_ANALYTIC,
    STRESS,
    RIC,              // solver of radiation induced conductivity model
    SOLVER_BASE,
    INVALID_SOLVER
  };


  /**
   * define solution type
   */
  enum SolutionType
  {
    EQUILIBRIUM=0,
    STEADYSTATE,
    OP,
    DCSWEEP,
    DCSWEEP_VSCAN,
    DCSWEEP_ISCAN,
    TRANSIENT,
    ACSWEEP,
    TRACE,
    INVALID_SolutionType
  };


  /**
   * define Newton damping
   */
  enum DampingScheme
  {
    DampingNo=0,
    DampingBankRose,
    DampingPotential,
    DampingSuperPotential
  };


  /**
   * enum whether to use the truncated voronoi box
   */
  enum  VoronoiTruncationFlag
  {
    VoronoiTruncationNo=0,
    VoronoiTruncationBoundary,
    VoronoiTruncationAlways
  };


  /**
   * define order for ODE solver
   */
  enum TemporalScheme
  {
    BDF1=0,
    BDF2,
    TRBDF2
  };



  /**
   * define carriers should be considered
   */
  enum Carrier
  {
    N_Type=0,
    P_Type,
    PN_Both,
    PN_None
  };


  /**
   * define how to evaluate high field mobility force
   */
  enum MobilityForce
  {
    EJ   = 0,
    ESimple,
    CarrierTemp
  };


  /**
   * define how to evaluate Impact Ionization force
   */
  enum IIForce
  {
    EdotJ   = 0,
    GradQf,
    EVector,
    ESide,
    SoftII,
    HardII,
    TempII
  };

  /**
   * convert string to enum
   */
  extern SolverType solver_type_string_to_enum(const std::string &);

  /**
   * convert string to enum
   */
  extern SolutionType solution_type_string_to_enum(const std::string s);
}




#endif // #define __enum_solver_specify_h__




