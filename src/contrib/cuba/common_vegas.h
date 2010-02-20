/*
	common.c
		include most of the modules
		this file is part of Vegas
		last modified 14 Feb 05 th
*/


#include "Random.h"
#include "ChiSquare.h"
#include "Grid_vegas.h"
#include "Integrate_vegas.h"


static inline bool BadDimension(cint ndim, cint flags)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  return ndim < SOBOL_MINDIM || (!PSEUDORNG && ndim > SOBOL_MAXDIM);
}


static inline bool BadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}

