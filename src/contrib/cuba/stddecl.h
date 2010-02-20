/*
	stddecl.h
		Type declarations common to all Cuba routines
		last modified 17 Dec 07 th
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>


//#if   __GNUC__
// support variable-size array
//#else

  #define NDIM  3
  #define NCOMP 4

//#endif

#ifndef NDIM
#define NDIM ndim_
#endif
#ifndef NCOMP
#define NCOMP ncomp_
#endif


#define VERBOSE (flags & 3)
#define LAST (flags & 4)
#define PSEUDORNG (flags & 8)
#define SHARPEDGES (flags & 16)
#define REGIONS (flags & 256)

#define INFTY DBL_MAX

#define NOTZERO 0x1p-104


#define Elements(x) (sizeof(x)/sizeof(*x))

#define Copy(d, s, n) memcpy(d, s, (n)*sizeof(*(d)))

#define VecCopy(d, s) Copy(d, s, ndim_)

#define ResCopy(d, s) Copy(d, s, ncomp_)

#define Clear(d, n) memset(d, 0, (n)*sizeof(*(d)))

#define VecClear(d) Clear(d, ndim_)

#define ResClear(d) Clear(d, ncomp_)

#define Zap(d) memset(d, 0, sizeof(d))

#define MaxErr(avg) Max(epsrel*fabs(avg), epsabs)


#ifdef __cplusplus
#define Malloc(p, n) (*(void **)&p = malloc(n))
#else
#define Malloc(p, n) (p = malloc(n))
#endif

#define MemAlloc(p, n) if( Malloc(p, n) == NULL ) { \
  fprintf(stderr, "Out of memory in " __FILE__ " line %d.\n", __LINE__); \
  exit(1); \
}

#define Alloc(p, n) MemAlloc(p, (n)*sizeof(*p))


#ifdef __cplusplus
#define Extern extern "C"
#else
#define Extern extern
typedef enum { false, true } bool;
#endif

typedef const bool cbool;

typedef const int cint;

typedef const long clong;

#define COUNT "%d"
typedef /*unsigned*/ int count;
typedef const count ccount;

#ifdef LONGLONGINT
#define PREFIX(s) ll##s
#define NUMBER "%lld"
#define NUMBER7 "%7lld"
typedef long long int number;
#else
#define PREFIX(s) s
#define NUMBER "%d"
#define NUMBER7 "%7d"
typedef int number;
#endif
typedef const number cnumber;

#define REAL "%g"
#define REALF "%f"
typedef /*long*/ double real;
	/* Switching to long double is not as trivial as it
	   might seem here.  sqrt, erf, exp, pow need to be
	   replaced by their long double versions (sqrtl, ...),
	   printf formats need to be updated similarly, and
	   ferrying long doubles to Mathematica is of course
	   quite another matter, too. */

typedef const real creal;


#ifdef UNDERSCORE
#define EXPORT(s) PREFIX(s##_)
#else
#define EXPORT(s) PREFIX(s)
#endif


static inline real Sq(creal x)
{
  return x*x;
}

static inline real Min(creal a, creal b)
{
  return (a < b) ? a : b;
}

static inline real Max(creal a, creal b)
{
  return (a > b) ? a : b;
}

static inline real Weight(creal sum, creal sqsum, cnumber n)
{
  creal w = sqrt(sqsum*n);
  return (n - 1)/Max((w + sum)*(w - sum), NOTZERO);
}


/* (a < 0) ? -1 : 0 */
#define NegQ(a) ((a) >> (sizeof(a)*8 - 1))

/* (a < 0) ? 0 : a */
#define IDim(a) ((a) & NegQ(-(a)))

/* (a < b) ? a : b */
#define IMin(a, b) ((a) - IDim((a) - (b)))

/* (a > b) ? a : b */
#define IMax(a, b) ((b) + IDim((a) - (b)))

/* (a == 0) ? 0 : -1 */
#define TrueQ(a) NegQ((a) | (-a))

/* a + (a == 0) */
#define Min1(a) ((a) + 1 + TrueQ(a))

/* abs(a) + (a == 0) */
#define Abs1(a) (((a) ^ NegQ(a)) - NegQ((a) - 1))

