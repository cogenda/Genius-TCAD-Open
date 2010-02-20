/*
	Random.c
		quasi- and pseudo-random-number generation
		last modified 5 Jan 06 th
*/


/*
	PART 1: Sobol quasi-random-number generator
	adapted from ACM TOMS algorithm 659
*/

#define SOBOL_MINDIM 1
#define SOBOL_MAXDIM 40

static struct {
  real norm;
  number v[SOBOL_MAXDIM][30], prev[SOBOL_MAXDIM];
  number seq;
} sobol_;


static inline void SobolIni(cnumber n)
{
  static number ini[9*40] = {
      3,   1,   0,   0,   0,   0,   0,   0,   0,
      7,   1,   1,   0,   0,   0,   0,   0,   0,
     11,   1,   3,   7,   0,   0,   0,   0,   0,
     13,   1,   1,   5,   0,   0,   0,   0,   0,
     19,   1,   3,   1,   1,   0,   0,   0,   0,
     25,   1,   1,   3,   7,   0,   0,   0,   0,
     37,   1,   3,   3,   9,   9,   0,   0,   0,
     59,   1,   3,   7,  13,   3,   0,   0,   0,
     47,   1,   1,   5,  11,  27,   0,   0,   0,
     61,   1,   3,   5,   1,  15,   0,   0,   0,
     55,   1,   1,   7,   3,  29,   0,   0,   0,
     41,   1,   3,   7,   7,  21,   0,   0,   0,
     67,   1,   1,   1,   9,  23,  37,   0,   0,
     97,   1,   3,   3,   5,  19,  33,   0,   0,
     91,   1,   1,   3,  13,  11,   7,   0,   0,
    109,   1,   1,   7,  13,  25,   5,   0,   0,
    103,   1,   3,   5,  11,   7,  11,   0,   0,
    115,   1,   1,   1,   3,  13,  39,   0,   0,
    131,   1,   3,   1,  15,  17,  63,  13,   0,
    193,   1,   1,   5,   5,   1,  27,  33,   0,
    137,   1,   3,   3,   3,  25,  17, 115,   0,
    145,   1,   1,   3,  15,  29,  15,  41,   0,
    143,   1,   3,   1,   7,   3,  23,  79,   0,
    241,   1,   3,   7,   9,  31,  29,  17,   0,
    157,   1,   1,   5,  13,  11,   3,  29,   0,
    185,   1,   3,   1,   9,   5,  21, 119,   0,
    167,   1,   1,   3,   1,  23,  13,  75,   0,
    229,   1,   3,   3,  11,  27,  31,  73,   0,
    171,   1,   1,   7,   7,  19,  25, 105,   0,
    213,   1,   3,   5,   5,  21,   9,   7,   0,
    191,   1,   1,   1,  15,   5,  49,  59,   0,
    253,   1,   1,   1,   1,   1,  33,  65,   0,
    203,   1,   3,   5,  15,  17,  19,  21,   0,
    211,   1,   1,   7,  11,  13,  29,   3,   0,
    239,   1,   3,   7,   5,   7,  11, 113,   0,
    247,   1,   1,   5,   3,  15,  19,  61,   0,
    285,   1,   3,   1,   1,   9,  27,  89,   7,
    369,   1,   1,   3,   7,  31,  15,  45,  23,
    299,   1,   3,   3,   9,   9,  25, 107,  39 };

  count dim, bit, nbits;
  number max, *pini = ini;

  for( nbits = 0, max = 1; max <= n; max <<= 1 ) ++nbits;
  sobol_.norm = 1./max;

  for( bit = 0; bit < nbits; ++bit )
    sobol_.v[0][bit] = (max >>= 1);

  for( dim = 1; dim < ndim_; ++dim ) {
    number *pv = sobol_.v[dim], *pvv = pv;
    number powers = *pini++, j;
    int inibits = -1, bit;
    for( j = powers; j; j >>= 1 ) ++inibits;

    memcpy(pv, pini, inibits*sizeof(*pini));
    pini += 8;

    for( bit = inibits; bit < nbits; ++bit ) {
      number newv = *pvv, j = powers;
      int b;
      for( b = 0; b < inibits; ++b ) {
        if( j & 1 ) newv ^= pvv[b] << (inibits - b);
        j >>= 1;
      }
      pvv[inibits] = newv;
      ++pvv;
    }

    for( bit = 0; bit < nbits - 1; ++bit )
      pv[bit] <<= nbits - bit - 1;
  }

  sobol_.seq = 0;
  VecClear(sobol_.prev);
}


static inline void SobolGet(real *x)
{
  number seq = sobol_.seq++;
  count zerobit = 0, dim;

  while( seq & 1 ) {
    ++zerobit;
    seq >>= 1;
  }

  for( dim = 0; dim < ndim_; ++dim ) {
    sobol_.prev[dim] ^= sobol_.v[dim][zerobit];
    x[dim] = sobol_.prev[dim]*sobol_.norm;
  }
}


static inline void SobolSkip(number n)
{
  while( n-- ) {
    number seq = sobol_.seq++;
    count zerobit = 0, dim;

    while( seq & 1 ) {
      ++zerobit;
      seq >>= 1;
    }

    for( dim = 0; dim < ndim_; ++dim )
      sobol_.prev[dim] ^= sobol_.v[dim][zerobit];
  }
}


/*
	PART 2: Mersenne Twister pseudo-random-number generator
	adapted from T. Nishimura's and M. Matsumoto's C code at
	http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
*/

/* length of state vector */
#define MERSENNE_N 624

/* period parameter */
#define MERSENNE_M 397

#define DEFAULT_SEED 5489

/* 32 or 53 random bits */
#define RANDOM_BITS 32

typedef unsigned int state_t;

static struct {
  state_t state[MERSENNE_N];
  count next;
} mersenne_;


static inline state_t Twist(state_t a, state_t b)
{
  state_t mixbits = (a & 0x80000000) | (b & 0x7fffffff);
  state_t matrixA = (-(b & 1)) & 0x9908b0df;
  return (mixbits >> 1) ^ matrixA;
}


static inline void MersenneReload()
{
  state_t *s = mersenne_.state;
  int j;

  for( j = MERSENNE_N - MERSENNE_M + 1; --j; ++s )
    *s = s[MERSENNE_M] ^ Twist(s[0], s[1]);
  for( j = MERSENNE_M; --j; ++s )
    *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], s[1]);
  *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], mersenne_.state[0]);
}


static inline void MersenneIni(state_t seed)
{
  state_t *next = mersenne_.state;
  int j;

  for( j = 1; j <= MERSENNE_N; ++j ) {
    *next++ = seed;
    seed = 0x6c078965*(seed ^ (seed >> 30)) + j;
    /* see Knuth TAOCP Vol 2, 3rd Ed, p. 106 for multiplier */
  }

  MersenneReload();
  mersenne_.next = 0;
}


static inline state_t MersenneInt(count next)
{
  state_t s = mersenne_.state[next];
  s ^= s >> 11;
  s ^= (s << 7) & 0x9d2c5680;
  s ^= (s << 15) & 0xefc60000;
  return s ^ (s >> 18);
}


static inline void MersenneGet(real *x)
{
  count next = mersenne_.next, dim;

  for( dim = 0; dim < ndim_; ++dim ) {
#if RANDOM_BITS == 53
    state_t a, b;
#endif

    if( next >= MERSENNE_N ) {
      MersenneReload();
      next = 0;
    }

#if RANDOM_BITS == 53
    a = MersenneInt(next++) >> 5;
    b = MersenneInt(next++) >> 6;
    x[dim] = (67108864.*a + b)/9007199254740992.;
#else
    x[dim] = MersenneInt(next++)/4294967295.;
#endif
  }

  mersenne_.next = next;
}


static inline void MersenneSkip(number n)
{
#if RANDOM_BITS == 53
  n = 2*n*ndim_ + mersenne_.next;
#else
  n = n*ndim_ + mersenne_.next;
#endif
  mersenne_.next = n % MERSENNE_N;
  n /= MERSENNE_N;
  while( n-- ) MersenneReload();
}


/*
	PART 3: User routines:

	- IniRandom sets up the random-number generator to produce a
	  sequence of at least n ndim_-dimensional random vectors.

	- GetRandom retrieves one random vector.

	- SkipRandom skips over n random vectors.
*/

static void IniRandom(cnumber n, cint flags)
{
  if( PSEUDORNG ) {
    sobol_.seq = -1;
    MersenneIni(DEFAULT_SEED);
  }
  else SobolIni(n);
}

static inline void GetRandom(real *x)
{
  if( sobol_.seq == -1 ) MersenneGet(x);
  else SobolGet(x);
}

static inline void SkipRandom(cnumber n)
{
  if( sobol_.seq == -1 ) MersenneSkip(n);
  else SobolSkip(n);
}

