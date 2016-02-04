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

#ifndef _genius_common_h_
#define _genius_common_h_

#define TCAD_SOLVERS
#define IDC_SOLVERS

#undef COGENDA_COMMERCIAL_PRODUCT


// C/C++ includes everyone should know about
#include <iostream> // needed for std::cout, std::cerr
#include <iomanip>  // needed for std::setw
#include <cassert>
#include <cstdlib>  // std::abort();
#include <complex>

#include "config.h"



// enable AMR
#define ENABLE_AMR


// float exception test
#ifdef HAVE_FENV_H
  #include <fenv.h>
#endif


// Undefine any existing macros
#ifdef Real
#  undef Real
#endif

#ifdef Complex
#  undef Complex
#endif

#ifdef COMPLEX
#  undef COMPLEX
#endif






// Define a corresponding tolerance.  This is what should be
// considered "good enough" when doing floating point comparisons.
// For example, v == 0 is changed to std::abs(v) < TOLERANCE.

typedef double Real;
typedef std::complex<double> Complex;

// Helper functions for complex/real numbers
// to clean up #ifdef USE_COMPLEX_NUMBERS elsewhere
template<typename T> inline T genius_real(T a) { return a; }
template<typename T> inline T genius_norm(T a) { return a*a; }

template<typename T>
inline T genius_real(std::complex<T> a) { return std::real(a); }

template<typename T>
inline T genius_norm(std::complex<T> a) { return std::norm(a); }


// Check to see if TOLERANCE has been defined by another
// package, if so we might want to change the name...
#ifdef TOLERANCE
   DIE A HORRIBLE DEATH HERE...
#  undef TOLERANCE
#endif

#define  TOLERANCE 1.e-8

const int invalid_id = -1;
const unsigned int invalid_uint = static_cast<unsigned int>(-1); // very large value: 4294967295


// Define the value type for unknowns in simulations.
// This may be double or long double, dependent on PETSC
typedef Real Number;


// Define the value type for error estimates.
// Since AMR/C decisions don't have to be precise,
// we default to float for memory efficiency.
typedef float ErrorVectorReal;


// 3D spatial dimension unless otherwise specified
#ifndef DIM
#  define DIM 3
#endif


// The genius_error() macro prints a message and aborts the code
#undef genius_error
#  define genius_error()    {  std::abort(); }

// The untested macro warns that you are using untested code
#undef genius_untested
#  define genius_untested() { std::cout << "*** Warning, This code is untested, experimental, or likely to see future API changes: " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << " ***" << std::endl; }

// The deprecated macro warns that you are using deprecated code
#undef genius_deprecated
#  define genius_deprecated() { std::cout << "*** Warning, This Code is Deprecated! " << __FILE__ << ", line " << __LINE__ << ", compiled " << __DATE__ << " at " << __TIME__ << " ***" << std::endl; }



#endif // #define _genius_common_h_
