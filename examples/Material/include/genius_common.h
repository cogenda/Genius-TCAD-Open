/*
 * genius_common.h
 *
 *  Created on: Nov 12, 2009
 *      Author: hash
 */

#ifndef GENIUS_COMMON_H_
#define GENIUS_COMMON_H_

#include "config.h"

#define genius_assert(x)    {assert(x); }

#define PetscScalar double

#ifdef CYGWIN
  // something required for building windows dll
  #define DLL_EXPORT_DECLARE  __declspec  (dllexport)
#else
  #define DLL_EXPORT_DECLARE
#endif

#define DIM 3
#define  TOLERANCE 1.e-8

typedef double Real;

#endif /* GENIUS_COMMON_H_ */
