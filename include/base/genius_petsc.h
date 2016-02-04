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

#ifndef GENIUS_PETSC_H_
#define GENIUS_PETSC_H_

#include "config.h"

#include "petscsys.h"
#include "petsc_macro.h"

#if PETSC_VERSION_LE(3,1,0)
#define PetscBool PetscTruth
#endif

#if PETSC_VERSION_GE(3,2,0)
  #define PetscDestroyObject(A) (&A)
#else
  #define PetscDestroyObject(A) (A)
#endif

// when PETSC_HAVE_MPIUNI is defined by PETSC, no MPI exist!
#ifdef PETSC_HAVE_MPIUNI
  #undef HAVE_MPI
#endif

#ifdef WINDOWS
  #define MUMPS_ICNTL_14 "50"
  #define MUMPS_ICNTL_23 "0"
#endif

#ifdef LINUX
  #define MUMPS_ICNTL_14 "50"
  #define MUMPS_ICNTL_23 "0"
#endif

#endif /* GENIUS_PETSC_H_ */
