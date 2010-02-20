/*
	Cuhre.c
		Adaptive integration using cubature rules
		by Thomas Hahn
		last modified 2 Mar 06 th
*/


#include "util_cuhre.h"

#define Print(s) puts(s); fflush(stdout)

static Integrand integrand_;

/*********************************************************************/

static inline void DoSample(count n, creal *x, real *f)
{
  neval_ += n;
  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_, f);
    x += ndim_;
    f += ncomp_;
  }
}

/*********************************************************************/

#include "common_cuhre.h"

Extern void EXPORT(Cuhre)(ccount ndim, ccount ncomp,
  Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  ccount key,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( BadComponent(ncomp) || BadDimension(ndim) ) *pfail = -1;
  else {
    neval_ = 0;
    integrand_ = integrand;

    *pfail = Integrate(epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, key,
      integral, error, prob);
    *pnregions = nregions_;
    *pneval = neval_;
  }
}

/*********************************************************************/

Extern void EXPORT(cuhre)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  ccount *pkey,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  EXPORT(Cuhre)(*pndim, *pncomp, integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pkey,
    pnregions, pneval, pfail,
    integral, error, prob);
}

