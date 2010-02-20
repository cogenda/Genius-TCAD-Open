/*
	Suave.c
		Subregion-adaptive Vegas Monte-Carlo integration
		by Thomas Hahn
		last modified 30 Aug 07 th
*/


#include "util_suave.h"

#define Print(s) puts(s); fflush(stdout)

static Integrand integrand_;

/*********************************************************************/

static inline void DoSample(number n, creal *w, creal *x, real *f)
{
  neval_ += n;
  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_, f, w++);
    x += ndim_;
    f += ncomp_;
  }
}

/*********************************************************************/

#include "common_suave.h"

Extern void EXPORT(Suave)(ccount ndim, ccount ncomp,
  Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nnew, creal flatness,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( BadComponent(ncomp) || BadDimension(ndim, flags) ) *pfail = -1;
  else {
    neval_ = 0;
    integrand_ = integrand;

    *pfail = Integrate(epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, nnew, flatness,
      integral, error, prob);
    *pnregions = nregions_;
    *pneval = neval_;
  }
}

/*********************************************************************/

Extern void EXPORT(suave)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  cnumber *pnnew, creal *pflatness,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  EXPORT(Suave)(*pndim, *pncomp,
    integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pnnew, *pflatness,
    pnregions, pneval, pfail,
    integral, error, prob);
}

