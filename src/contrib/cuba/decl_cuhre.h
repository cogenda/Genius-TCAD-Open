/*
	decl.h
		Type declarations
		this file is part of Cuhre
		last modified 19 Jan 05 th
*/


#include "stddecl.h"


typedef struct {
  real avg, err;
  count bisectdim;
} Result;

typedef const Result cResult;


typedef struct {
  real avg, err, lastavg, lasterr;
  real weightsum, avgsum;
  real guess, chisum, chisqsum, chisq;
} Totals;

typedef const Totals cTotals;


typedef struct {
  real lower, upper;
} Bounds;

typedef const Bounds cBounds;


typedef struct {
  real *x, *f;
  void *first, *last;
  real errcoeff[3];
  count n;
} Rule;

typedef const Rule cRule;


#define TYPEDEFREGION \
  typedef struct region { \
    struct region *next; \
    count div; \
    Result result[NCOMP]; \
    Bounds bounds[NDIM]; \
  } Region


typedef void (*Integrand)(ccount *, creal *, ccount *, real *);

