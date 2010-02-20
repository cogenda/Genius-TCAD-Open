/*
	Divonne.c
		Multidimensional integration by partitioning
		originally by J.H. Friedman and M.H. Wright
		(CERNLIB subroutine D151)
		this version by Thomas Hahn
		last modified 2 Mar 06 th
*/

#include "util_divonne.h"

#define Print(s) puts(s); fflush(stdout)

static Integrand integrand_;
static PeakFinder peakfinder_;

/*********************************************************************/

static inline void DoSample(number n, ccount ldx, creal *x, real *f)
{
  neval_ += n;
  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_, f, &phase_);
    x += ldx;
    f += ncomp_;
  }
}

/*********************************************************************/

static inline count SampleExtra(cBounds *b)
{
  number n = nextra_;
  peakfinder_(&ndim_, b, &n, xextra_);
  DoSample(n, ldxgiven_, xextra_, fextra_);
  return n;
}

/*********************************************************************/

#include "common_divonne.h"

Extern void EXPORT(Divonne)(ccount ndim, ccount ncomp,
  Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cint key1, cint key2, cint key3, ccount maxpass,
  creal border, creal maxchisq, creal mindeviation,
  cnumber ngiven, ccount ldxgiven, real *xgiven,
  cnumber nextra, PeakFinder peakfinder,
  int *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( BadComponent(ncomp) ||
      BadDimension(ndim, flags, key1) ||
      BadDimension(ndim, flags, key2) ||
      ((key3 & -2) && BadDimension(ndim, flags, key3)) ) *pfail = -1;
  else {
    neval_ = neval_opt_ = neval_cut_ = 0;
    integrand_ = integrand;
    peakfinder_ = peakfinder;
    border_.lower = border;
    border_.upper = 1 - border_.lower;
    ngiven_ = ngiven;
    xgiven_ = NULL;
    ldxgiven_ = IMax(ldxgiven, ndim_);
    nextra_ = nextra;

    if( ngiven + nextra ) {
      cnumber nxgiven = ngiven*ldxgiven;
      cnumber nxextra = nextra*ldxgiven;
      cnumber nfgiven = ngiven*ncomp;
      cnumber nfextra = nextra*ncomp;

      Alloc(xgiven_, nxgiven + nxextra + nfgiven + nfextra);
      xextra_ = xgiven_ + nxgiven;
      fgiven_ = xextra_ + nxextra;
      fextra_ = fgiven_ + nfgiven;

      if( nxgiven ) {
        phase_ = 0;
        Copy(xgiven_, xgiven, nxgiven);
        DoSample(ngiven_, ldxgiven_, xgiven_, fgiven_);
      }
    }

    *pfail = Integrate(epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, key1, key2, key3, maxpass,
      maxchisq, mindeviation,
      integral, error, prob);
    *pnregions = nregions_;
    *pneval = neval_;

    if( xgiven_ ) free(xgiven_);
  }
}

/*********************************************************************/

Extern void EXPORT(divonne)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  cint *pkey1, cint *pkey2, cint *pkey3, ccount *pmaxpass,
  creal *pborder, creal *pmaxchisq, creal *pmindeviation,
  cnumber *pngiven, ccount *pldxgiven, real *xgiven,
  cnumber *pnextra, PeakFinder peakfinder,
  int *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  EXPORT(Divonne)(*pndim, *pncomp,
    integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pkey1, *pkey2, *pkey3, *pmaxpass,
    *pborder, *pmaxchisq, *pmindeviation,
    *pngiven, *pldxgiven, xgiven,
    *pnextra, peakfinder,
    pnregions, pneval, pfail,
    integral, error, prob);
}

