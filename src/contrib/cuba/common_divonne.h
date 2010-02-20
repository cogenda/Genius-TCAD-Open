/*
	common.c
		includes most of the modules
		this file is part of Divonne
		last modified 14 Feb 05 th
*/


static void Explore(void *voidregion, cSamples *samples, cint depth, cint flags);

static void Split(void *voidregion, int depth);

#include "Random.h"
#include "ChiSquare.h"
#include "Rule_divonne.h"
#include "Sample_divonne.h"
#include "FindMinimum_divonne.h"
#include "Explore_divonne.h"
#include "Split_divonne.h"
#include "Integrate_divonne.h"


static inline bool BadDimension(ccount ndim, cint flags, ccount key)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  if( IsSobol(key) ) return
    ndim < SOBOL_MINDIM || (!PSEUDORNG && ndim > SOBOL_MAXDIM);
  if( IsRule(key, ndim) ) return ndim < 1;
  return ndim < KOROBOV_MINDIM || ndim > KOROBOV_MAXDIM;
}


static inline bool BadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}

