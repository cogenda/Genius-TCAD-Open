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


#ifndef __enum_solution_h__
#define __enum_solution_h__

#include <string>

// ------------------------------------------------------------
// enum solution variable


  /**
   * \enum SolutionVariable defines an \p enum for variable used in Genius solution
   */
  enum SolutionVariable
  {
    POTENTIAL=0, /* potential */
    ELECTRON,    /* electron concentration */
    HOLE,        /* hole concentration */
    TEMPERATURE, /* lattice temperature */
    E_TEMP,      /* electron temperature */
    H_TEMP,      /* hole temperature */
    QFN,         /* electron quasi-Fermi level */
    QFP,         /* hole quasi-Fermi level */

    E_FIELD,     /* electric field */
    DOPING,      /* net doping */
    DOPING_Na,   /* acceptor  */
    DOPING_Nd,   /* donor */
    MIN_CARRIER, /* minority carrier concentration */
    NET_CARRIER, /* net carrier concentration */
    NET_CHARGE,  /* net charge */

    MOLE_X,
    MOLE_Y,
    OPTICAL_GEN, /* charge genetated by optical ray */
    OPTICAL_HEAT,/* heat genetated by optical ray */
    PARTICLE_GEN,/* charge genetated by particle ray */
    INVALID_Variable
  };

extern SolutionVariable solution_string_to_enum(const std::string &);

#endif // #define __enum_solution_h__
