/*
	Sample.c
		most of what is related to sampling
		this file is part of Divonne
		last modified 4 Mar 05 th
*/


#define MARKMASK 0xfffffff
#define Marked(x) ((x) & ~MARKMASK)
#define Unmark(x) ((x) & MARKMASK)

#define MEM(samples) (samples)->x

/*********************************************************************/

static inline void SamplesIni(Samples *samples)
{
  MEM(samples) = NULL;
}

/*********************************************************************/

static inline void SamplesFree(cSamples *samples)
{
  if( MEM(samples) ) free(MEM(samples));
}

/*********************************************************************/

static void SampleSobol(cSamples *samples, cBounds *b, creal vol)
{
  creal norm = vol*samples->weight;
  real *x = samples->x, *f = samples->f, *avg = samples->avg;
  cnumber n = samples->n;
  number i;
  count dim, comp;

  for( i = 0; i < n; ++i ) {
    GetRandom(x);
    for( dim = 0; dim < ndim_; ++x, ++dim )
      *x = b[dim].lower + *x*(b[dim].upper - b[dim].lower);
  }

  DoSample(n, ndim_, samples->x, f);

  ResCopy(avg, f);
  f += ncomp_;
  for( i = 1; i < n; ++i )
    for( comp = 0; comp < ncomp_; ++comp )
      avg[comp] += *f++;

  for( comp = 0; comp < ncomp_; ++comp )
    avg[comp] *= norm;
}

/*********************************************************************/

static void SampleKorobov(cSamples *samples, cBounds *b, creal vol)
{
  creal norm = vol*samples->weight;
  real *x = samples->x, *xlast = x + ndim_;
  real *f = samples->f, *flast = f + ncomp_;
  real *avg = samples->avg;
  cnumber n = samples->n, neff = samples->neff;
  number nextra = n, i;
  real dist = 0;
  count dim, comp;

  for( i = 1; i < n; ++i ) {
    number c = i;
    for( dim = 0; dim < ndim_; ++dim ) {
      creal dx = abs(2*c - neff)*samples->weight;
      *xlast++ = b[dim].lower + dx*(b[dim].upper - b[dim].lower);
      c = c*samples->coeff % neff;
    }
  }

  for( dim = 0; dim < ndim_; ++dim ) {
    creal dx = (x[dim] = b[dim].upper) - border_.upper;
    if( dx > 0 ) dist += Sq(dx);
  }

  if( dist > 0 ) {
    dist = sqrt(dist)/EXTRAPOLATE_EPS;
    for( dim = 0; dim < ndim_; ++dim ) {
      real x2 = x[dim], dx = x2 - border_.upper;
      if( dx > 0 ) {
        x[dim] = border_.upper;
        x2 = border_.upper - dx/dist;
      }
      xlast[dim] = x2;
    }
    ++nextra;
  }

  DoSample(nextra, ndim_, x, f);

  ResCopy(avg, flast);
  flast += ncomp_;
  for( i = 2; i < n; ++i )
    for( comp = 0; comp < ncomp_; ++comp )
      avg[comp] += *flast++;

  if( nextra > n ) {
    for( comp = 0; comp < ncomp_; ++comp )
      f[comp] += dist*(f[comp] - flast[comp]);
    for( dim = 0; dim < ndim_; ++dim )
      x[dim] = b[dim].upper;
  }

  for( comp = 0; comp < ncomp_; ++comp )
    avg[comp] = (avg[comp] + avg[comp] + f[comp])*norm;
}

/*********************************************************************/

#define IsSobol(k) NegQ(k)
#define IsRule(k, d) (k == 9 || k == 7 || (k == 11 && d == 3) || (k == 13 && d == 2))

/* The following coding is used for key1, key2, key3:
             0 = for key1, key2: use default,
                 for key3: do nothing,
             1 = for key3: split region again,
             7 = degree-7 cubature rule,
             9 = degree-9 cubature rule,
            11 = degree-11 cubature rule (only in 3 dims),
            13 = degree-13 cubature rule (only in 2 dims),
     -inf..-40 = absolute # of points, Sobol numbers,
       -39..-1 = multiplicator,        Sobol numbers,
         1..39 = multiplicator,        Korobov numbers,
       40..inf = absolute # of points, Korobov numbers.  */

static count SamplesLookup(Samples *samples, cint key,
  cnumber nwant, cnumber nmax, number nmin)
{
  number n;

  if( key == 13 && ndim_ == 2 ) {
    if( rule13_.first == NULL ) Rule13Alloc(&rule13_);
    samples->rule = &rule13_;
    samples->n = n = nmin = rule13_.n;
    samples->sampler = SampleRule;
  }
  else if( key == 11 && ndim_ == 3 ) {
    if( rule11_.first == NULL ) Rule11Alloc(&rule11_);
    samples->rule = &rule11_;
    samples->n = n = nmin = rule11_.n;
    samples->sampler = SampleRule;
  }
  else if( key == 9 ) {
    if( rule9_.first == NULL ) Rule9Alloc(&rule9_);
    samples->rule = &rule9_;
    samples->n = n = nmin = rule9_.n;
    samples->sampler = SampleRule;
  }
  else if( key == 7 ) {
    if( rule7_.first == NULL ) Rule7Alloc(&rule7_);
    samples->rule = &rule7_;
    samples->n = n = nmin = rule7_.n;
    samples->sampler = SampleRule;
  }
  else {
    n = Abs1(key);
    if( n < 40 ) n *= nwant;
    samples->sampler = (key < 0) ? SampleSobol :
      (n = n/2 + 1, SampleKorobov);
    samples->n = IMin(n, nmax);
  }

  samples->neff = samples->n;

  return IDim(n - nmax) | Marked(nmax - nmin);
}

/*********************************************************************/

static void SamplesAlloc(Samples *samples)
{
#define FIRST -INT_MAX
#define MarkLast(x) (x | Marked(INT_MAX))

#include "KorobovCoeff_divonne.h"

  number nx, nf;

  if( samples->sampler == SampleKorobov ) {
    enum { max = Elements(prime) - 2 };
    cint n = IMin(2*samples->n - 1, MAXPRIME);
    int i = Hash(n), p;
    count shift = 2 + NegQ(n - 1000);

    while( i = IMin(IDim(i), max),
           n > (p = prime[i + 1]) || n <= prime[i] ) {
      cint d = (n - Unmark(p)) >> ++shift;
      i += Min1(d);
    }

    samples->coeff = coeff[i][ndim_ - KOROBOV_MINDIM];
    samples->neff = p = Unmark(p);
    samples->n = p/2 + 1;
  }

  nx = ndim_*(samples->n + 1);		/* need 1 for extrapolation */
  nf = ncomp_*(samples->n + 1);

  Alloc(samples->x, nx + nf + ncomp_ + ncomp_);
  samples->f = samples->x + nx;
  samples->avg = samples->f + nf;
  samples->err = samples->avg + ncomp_;
  ResClear(samples->err);

  samples->weight = 1./samples->neff;
}

/*********************************************************************/

static real Sample(creal *x0)
{
  real xtmp[2*NDIM], ftmp[2*NCOMP], *xlast = xtmp, f;
  real dist = 0;
  count dim;
  number nextra = 1;

  for( dim = 0; dim < ndim_; ++dim ) {
    creal x1 = *xlast++ = Min(Max(*x0++, 0.), 1.);
    real dx;
    if( (dx = x1 - border_.lower) < 0 ||
        (dx = x1 - border_.upper) > 0 ) dist += Sq(dx);
  }

  if( dist > 0 ) {
    dist = sqrt(dist)/EXTRAPOLATE_EPS;
    for( dim = 0; dim < ndim_; ++dim ) {
      real x2 = xtmp[dim], dx, b;
      if( (dx = x2 - (b = border_.lower)) < 0 ||
          (dx = x2 - (b = border_.upper)) > 0 ) {
        xtmp[dim] = b;
        x2 = b - dx/dist;
      }
      *xlast++ = x2;
    }
    nextra = 2;
  }

  DoSample(nextra, ndim_, xtmp, ftmp);

  f = ftmp[selectedcomp_];
  if( nextra > 1 ) f += dist*(f - ftmp[selectedcomp_ + ncomp_]);

  return f;
}

