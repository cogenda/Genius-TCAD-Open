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

//  $Id: physical_unit.cc,v 1.4 2008/07/09 05:58:16 gdiso Exp $

#include <cmath>
#include "physical_unit.h"


namespace PhysicalUnit
{

  /**
   *  the basic length unit
   */
  double   m;

  /**
   * the basic time unit
   */
  double   s;

  /**
   * potential unit
   */
  double   V;

  /**
   * the charge unit
   */
  double   C;

  /**
   * the temperature unit
   */
  double   K;


  // set as induced physical unit

  /**
   *  define as 1e-9*m
   */
  double   nm;

  /**
   *  define as 1e-6*m
   */
  double   um;

  /**
   *  define as 1e-2*m
   */
  double   cm;

  /**
   * the mass unit, define as J/(m*m)*s*s
   */
  double   kg;

  /**
   * the mass unit, define 1e-3*kg
   */
  double   g;

  /**
   * energy unit, define as C dot V
   */
  double   J;

  /**
   * energy unit, define as e dot V, e = C/1.602176462e-19
   */
  double   eV;

  /**
   * power unit, define as J/s
   */
  double   W;

  /**
   *  define as 1e-6*s
   */
  double   us;

  /**
   * define as 1e-12*s
   */
  double   ps;

  /**
   * define as C/s
   */
  double   A;

  /**
   * define as 1e-3*A
   */
  double   mA;

  // Fundamental Physical Constants

  /**
   * Boltzmann constant
   */
  double   kb;

  /**
   * elementary charge
   */
  double   e;

  /**
   * electron mass
   */
  double   me;

  /**
   * electric constant
   */
  double   eps0;

  /**
   * magnetic constant
   */
  double   mu0;

  /**
   * Planck constant
   */
  double   h;

  /**
   * Planck constant over 2*PI
   */
  double   hbar;

  void set_unit()
  {
    cm = 1e6;
    s  = 1e12;
    V  = 1.0;
    C  = 1.0/1.602176462e-19;
    K  = 1.0/300;

    m  = 1e2*cm;
    nm = 1e-7*cm;
    um = 1e-4*cm;
    J  = C*V;
    W  = J/s;
    kg = J/(m*m)*s*s;
    g  = 1e-3*kg;
    eV = 1.602176462e-19*J;
    us = 1e-6*s;
    ps = 1e-12*s;
    A  = C/s;
    mA = 1e-3*A;

    kb   = 1.3806503e-23*J/K;
    e    = 1.602176462e-19*C;
    me   = 9.10938188e-31*kg;
    eps0 = 8.854187818e-12*C/V/m;
    mu0  = 12.56637061e-7*std::pow(s,2)/C*V/m;
    h    = 6.62606876e-34*J*s;
    hbar = 1.054571596e-34*J*s;
  }

  //user defined scale value : length only
  void set_unit(double length)
  {

    cm = length;
    s  = 1e12;
    V  = 1.0;
    C  = 1.0/1.602176462e-19;
    K  = 1.0/300;

    m  = 1e2*cm;
    nm = 1e-7*cm;
    um = 1e-4*cm;
    J  = C*V;
    W  = J/s;
    kg = J/(m*m)*s*s;
    g  = 1e-3*kg;
    eV = 1.602176462e-19*J;
    us = 1e-6*s;
    ps = 1e-12*s;
    A  = C/s;
    mA = 1e-3*A;

    kb   = 1.3806503e-23*J/K;
    e    = 1.602176462e-19*C;
    me   = 9.10938188e-31*kg;
    eps0 = 8.854187818e-12*C/V/m;
    mu0  = 12.56637061e-7*std::pow(s,2)/C*V/m;
    h    = 6.62606876e-34*J*s;
    hbar = 1.054571596e-34*J*s;
  }

  //user defined scale value : length and potential
  void set_unit(double length,double potential)
  {
    cm = length;
    s  = 1e12;
    V  = potential;
    C  = 1.0/1.602176462e-19;
    K  = 1.0/300;

    m  = 1e2*cm;
    nm = 1e-7*cm;
    um = 1e-4*cm;
    J  = C*V;
    W  = J/s;
    kg = J/(m*m)*s*s;
    g  = 1e-3*kg;
    eV = 1.602176462e-19*J;
    us = 1e-6*s;
    ps = 1e-12*s;
    A  = C/s;
    mA = 1e-3*A;

    kb   = 1.3806503e-23*J/K;
    e    = 1.602176462e-19*C;
    me   = 9.10938188e-31*kg;
    eps0 = 8.854187818e-12*C/V/m;
    mu0  = 12.56637061e-7*std::pow(s,2)/C*V/m;
    h    = 6.62606876e-34*J*s;
    hbar = 1.054571596e-34*J*s;
  }


  //user defined scale value
  void set_unit(double length, double time,double voltage,
                double charge, double temperature)
  {
    cm = length;
    s  = time;
    V  = voltage;
    C  = charge;
    K  = temperature;

    m  = 1e2*cm;
    nm = 1e-7*cm;
    um = 1e-4*cm;
    J  = C*V;
    W  = J/s;
    kg = J/(m*m)*s*s;
    g  = 1e-3*kg;
    eV = 1.602176462e-19*J;
    us = 1e-6*s;
    ps = 1e-12*s;
    A  = C/s;
    mA = 1e-3*A;

    kb   = 1.3806503e-23*J/K;
    e    = 1.602176462e-19*C;
    me   = 9.10938188e-31*kg;
    eps0 = 8.854187818e-12*C/V/m;
    mu0  = 12.56637061e-7*std::pow(s,2)/C*V/m;
    h    = 6.62606876e-34*J*s;
    hbar = 1.054571596e-34*J*s;
  }




}
