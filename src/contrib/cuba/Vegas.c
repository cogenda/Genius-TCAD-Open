/*
	Vegas.c
		Vegas Monte-Carlo integration
		by Thomas Hahn
		last modified 30 Aug 07 th
*/


#include "util_vegas.h"

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

#include "common_vegas.h"

Extern void EXPORT(Vegas)(ccount ndim, ccount ncomp,
  Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
  number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( BadComponent(ncomp) || BadDimension(ndim, flags) ) *pfail = -1;
  else {
    neval_ = 0;
    integrand_ = integrand;

    *pfail = Integrate(epsrel, epsabs,
      flags, mineval, maxeval, nstart, nincrease,
      integral, error, prob);

    *pneval = neval_;
  }
}

/*********************************************************************/

Extern void EXPORT(vegas)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  cnumber *pnstart, cnumber *pnincrease, 
  number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  /* make sure the filename is null-terminated */
  if( *EXPORT(vegasstate) ) {
    char *p;
    EXPORT(vegasstate)[sizeof(EXPORT(vegasstate)) - 1] = 0;
    if( (p = strchr(EXPORT(vegasstate), ' ')) ) *p = 0;
  }

  EXPORT(Vegas)(*pndim, *pncomp,
    integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pnstart, *pnincrease,
    pneval, pfail,
    integral, error, prob);
}

