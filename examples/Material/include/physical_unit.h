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

//  $Id: physical_unit.h,v 1.7 2008/07/09 05:58:16 gdiso Exp $

#ifndef __physical_unit_h__
#define __physical_unit_h__

/**
 * this namespace contains scaled physical unit.
 * it should keeps the same in main genius code and
 * parameter library
 */
namespace PhysicalUnit
{

  // use as fundamental physical unit

  /**
   *  the basic length unit
   */
  extern double   m;

  /**
   * the basic time unit
   */
  extern double   s;

  /**
   * potential unit
   */
  extern double   V;

  /**
   * the charge unit
   */
  extern double   C;

  /**
   * the temperature unit
   */
  extern double   K;


  // set as induced physical unit

  /**
   *  define as 1e-9*m
   */
  extern double   nm;

  /**
   *  define as 1e-6*m
   */
  extern double   um;

  /**
   *  define as 1e-2*m
   */
  extern double   cm;

  /**
   * the mass unit, define as J/(m*m)*s*s
   */
  extern double   kg;

  /**
   * the mass unit, define 1e-3*kg
   */
  extern double   g;

  /**
   * energy unit, define as C dot V
   */
  extern double   J;

  /**
   * energy unit, define as e dot V, e = C/1.602176462e-19
   */
  extern double   eV;

  /**
   * power unit, define as J/s
   */
  extern double   W;

  /**
   *  define as 1e-6*s
   */
  extern double   us;

  /**
   *  define as 1e-9*s
   */
  extern double   ns;

  /**
   * define as 1e-12*s
   */
  extern double   ps;

  /**
   * define as C/s
   */
  extern double   A;

  /**
   * define as 1e-3*A
   */
  extern double   mA;

  /**
   *  define as V/A
   */
  extern double   Ohm;

  // Fundamental Physical Constants

  /**
   * Boltzmann constant
   */
  extern double   kb;

  /**
   * elementary charge
   */
  extern double   e;

  /**
   * electron mass
   */
  extern double   me;

  /**
   * electric constant
   */
  extern double   eps0;

  /**
   * magnetic constant
   */
  extern double   mu0;

  /**
   * Planck constant
   */
  extern double   h;

  /**
   * Planck constant over 2*PI
   */
  extern double   hbar;

  /**
   * set every unit by default value
   */
  extern void set_unit();

  /**
   * set length unit by parameter length, others use default value
   */
  extern void set_unit(double length);

  /**
   * set length unit by parameter length and potential unit by parameter potential,
   * others use default value
   */
  extern void set_unit(double lenght, double potential);

  /**
   * set each unit by user provide value
   */
  extern void set_unit(double length, double time,double voltage, double charge, double temperature);

}

#endif
