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
    POISSON = 0,
    DDML1,
    DDML1MIX,
    DDML1MIXA,
    HALLDDML1,
    QDDML1,
    DDML2,
    DDML2MIX,
    DDML2MIXA,
    EBML3,
    EBML3MIX,
    EBML3MIXA,
    DDMAC,
    MC,
    RAY_TRACE,
    EM_FEM_2D,
    EM_FEM_3D,
    FDTD,
    FVTD,
    DOPING_ANALYTIC,
    MOLE_ANALYTIC,
    OPTG_ANALYTIC,
    STRESS,
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
    DampingPotential
  };



  /**
   * define order for ODE solver
   */
  enum TSType
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



}




#endif // #define __enum_solver_specify_h__




