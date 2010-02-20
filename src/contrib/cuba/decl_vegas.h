/*
	decl.h
		Type declarations
		this file is part of Vegas
		last modified 30 Aug 07 th
*/


#include "stddecl.h"

#define MAXGRIDS 10

#define MAXSTATESIZE 128

#define NBINS 128

typedef unsigned char bin_t;
/* Note: bin_t must be wide enough to hold the numbers 0..NBINS */

typedef const bin_t cbin_t;

typedef real Grid[NBINS];

typedef struct {
  real sum, sqsum;
  real weightsum, avgsum;
  real chisum, chisqsum, guess;
  real avg, err, chisq;
} Cumulants;

typedef const Cumulants cCumulants;

typedef void (*Integrand)(ccount *, creal *, ccount *, real *, creal *);

