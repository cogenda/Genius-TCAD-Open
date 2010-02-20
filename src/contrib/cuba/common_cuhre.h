/*
	common.c
		includes most of the modules
		this file is part of Cuhre
		last modified 14 Feb 05 th
*/


#include "ChiSquare.h"
#include "Rule_cuhre.h"
#include "Integrate_cuhre.h"


static inline bool BadDimension(ccount ndim)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  return ndim < 2;
}


static inline bool BadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}

