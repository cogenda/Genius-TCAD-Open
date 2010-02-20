/*
	Split.c
		determine optimal cuts for splitting a region
		this file is part of Divonne
		last modified 7 Mar 05 th
*/


#define BNDTOL .05
#define FRACT .5
#define BIG 1e10
#define SINGTOL 1e-4

#define LHSTOL .1
#define GAMMATOL .1

/* the next four macros must be in sync with the typedef of Bounds! */
#define Lower(d) (2*d)
#define Upper(d) (2*d + 1)
#define Dim(i) ((i) >> 1)
#define SignedDelta(i) ((i & 1) ? delta[i] : -delta[i])

typedef struct {
  count i;
  real save, delta;
  real f, df, fold;
  real row, lhs, sol;
} Cut;

typedef struct {
  real diff, err, spread;
} Errors;

typedef const Errors cErrors;

/*********************************************************************/

static void SomeCut(Cut *cut, Bounds *b)
{
  count dim, maxdim;
  static count nextdim = 0;
  real xmid[NDIM], ymid, maxdev;

  for( dim = 0; dim < ndim_; ++dim )
    xmid[dim] = .5*(b[dim].upper + b[dim].lower);
  ymid = Sample(xmid);

  maxdev = 0;
  maxdim = 0;
  for( dim = 0; dim < ndim_; ++dim ) {
    real ylower, yupper, dev;
    creal x = xmid[dim];
    xmid[dim] = b[dim].lower;
    ylower = Sample(xmid);
    xmid[dim] = b[dim].upper;
    yupper = Sample(xmid);
    xmid[dim] = x;

    dev = fabs(ymid - .5*(ylower + yupper));
    if( dev >= maxdev ) {
      maxdev = dev;
      maxdim = dim;
    }
  }

  if( maxdev > 0 ) nextdim = 0;
  else maxdim = nextdim++ % ndim_;

  cut->i = Upper(maxdim);
  cut->save = b[maxdim].upper;
  b[maxdim].upper = xmid[maxdim];
}

/*********************************************************************/

static inline real Volume(creal *delta)
{
  real vol = 1;
  count dim;
  for( dim = 0; dim < ndim_; ++dim )
    vol *= delta[Lower(dim)] + delta[Upper(dim)];
  return vol;
}

/*********************************************************************/

static inline real SetupEqs(Cut *cutfirst, Cut *cutlast, real f)
{
  real sqsum = Sq(cutlast->lhs = f - cutfirst->f);
  while( cutlast > cutfirst ) {
    f = cutlast->f;
    --cutlast;
    sqsum += Sq(cutlast->lhs = f - cutlast->f);
  }
  return sqsum;
}

/*********************************************************************/

static inline void SolveEqs(Cut *cutfirst, Cut *cutlast,
  creal *delta, creal diff)
{
  real last = cutlast->lhs;
  real r = 1;
  Cut *cut;

  for( cut = cutfirst; ; ++cut ) {
    ccount dim = Dim(cut->i);
    cut->row = r*cut->df - diff/(delta[Lower(dim)] + delta[Upper(dim)]);
    if( cut == cutlast ) break;
    r = (fabs(cut->row) < BIG*fabs(cut->df)) ? cut->row/cut->df : 0;
    last -= r*cut->lhs;
  }

  if( cut->row != 0 && fabs(cut->row) < BIG*fabs(last) )
    last /= cut->row;

  for( ; ; ) {
    creal delmin = -(cut->delta = delta[cut->i]);
    creal delmax = FRACT*(delmin + cut->save);
    cut->sol = last;
    if( cut->sol > delmax ) cut->sol = .75*delmax;
    if( cut->sol < delmin ) cut->sol = .75*delmin;
    if( cut == cutfirst ) break;
    last *= cut->df;
    --cut;
    last += cut->lhs;
    if( cut->df != 0 && fabs(cut->df) < BIG*fabs(last) )
      last /= cut->df;
  }
}

/*********************************************************************/

static void Split(void *voidregion, int depth)
{
  TYPEDEFREGION;

  Region *const region = (Region *)voidregion, *reg, *next;
  creal fdiff = region->fmajor - region->fminor;
  cint sign = (fdiff < 0) ? -1 : 1;

  Cut cutfirst[2*NDIM], *cutlast = cutfirst, *cut;
  real delta[2*NDIM];
  real gamma, fgamma, lhssq;
  real *b0, *b1, tmp;
  Errors errors[NCOMP];
  count dim, comp, div, nsplit;

  selectedcomp_ = region->cutcomp;
  neval_cut_ -= neval_;

  for( dim = 0; dim < ndim_; ++dim ) {
    cBounds *b = &region->bounds[dim];
    real *x = &region->xmajor[dim];
    creal xsave = *x;
    real dist = b->upper - xsave;
    if( dist >= BNDTOL*(b->upper - b->lower) ) {
      cutlast->i = Upper(dim);
      cutlast->save = dist;
      *x += dist *= FRACT;
      cutlast->f = Sample(region->xmajor);
      *x = xsave;
      ++cutlast;
    }
    delta[Upper(dim)] = dist;
  }

  for( dim = 0; dim < ndim_; ++dim ) {
    cBounds *b = &region->bounds[dim];
    real *x = &region->xmajor[dim];
    creal xsave = *x;
    real dist = xsave - b->lower;
    if( dist >= BNDTOL*(b->upper - b->lower) ) {
      cutlast->i = Lower(dim);
      cutlast->save = dist;
      *x -= dist *= FRACT;
      cutlast->f = Sample(region->xmajor);
      *x = xsave;
      ++cutlast;
    }
    delta[Lower(dim)] = dist;
  }

  if( cutlast == cutfirst ) {
    SomeCut(cutlast, region->bounds);
    goto dissect;
  }

  for( ; ; ) {
    real mindiff = INFTY;
    Cut *mincut = cutfirst;

    --cutlast;
    for( cut = cutfirst; cut <= cutlast; ++cut ) {
      creal diff = fabs(region->fmajor - cut->f);
      if( diff <= mindiff ) {
        mindiff = diff;
        mincut = cut;
      }
    }

    gamma = Volume(delta)/region->vol;
    fgamma = gamma*region->fmajor + (1 - gamma)*region->fminor;

    if( sign*(mincut->f - fgamma) < 0 ) break;

    if( cutlast == cutfirst ) {
      SomeCut(cutlast, region->bounds);
      goto dissect;
    }

    delta[mincut->i] = mincut->save;
    memcpy(mincut, mincut + 1, (char *)cutlast - (char *)mincut);
  }

  for( cut = cutfirst; cut <= cutlast; ++cut ) {
    cut->fold = cut->f;
    cut->df = (cut->f - region->fmajor)/delta[cut->i];
  }

  lhssq = SetupEqs(cutfirst, cutlast, fgamma);

repeat:
  SolveEqs(cutfirst, cutlast, delta, gamma*fdiff);

  for( div = 1; div <= 16; div *= 4 ) {
    real gammanew, lhssqnew;

    for( cut = cutfirst; cut <= cutlast; ++cut ) {
      real *x = &region->xmajor[Dim(cut->i)];
      creal xsave = *x;
      delta[cut->i] = cut->delta + cut->sol/div;
      *x += SignedDelta(cut->i);
      cut->f = Sample(region->xmajor);
      *x = xsave;
    }

    gammanew = Volume(delta)/region->vol;
    fgamma = gammanew*region->fmajor + (1 - gammanew)*region->fminor;
    lhssqnew = SetupEqs(cutfirst, cutlast, fgamma);

    if( lhssqnew <= lhssq ) {
      real fmax;

      if( fabs(gammanew - gamma) < GAMMATOL*gamma ) break;
      gamma = gammanew;

      fmax = fabs(fgamma);
      for( cut = cutfirst; cut <= cutlast; ++cut ) {
        creal dfmin = SINGTOL*cut->df;
        creal sol = cut->sol/div;
        real df = cut->f - cut->fold;
        df = (fabs(sol) < BIG*fabs(df)) ? df/sol : 1;
        cut->df = (fabs(df) < fabs(dfmin)) ? dfmin : df;
        fmax = Max(fmax, fabs(cut->f));
        cut->fold = cut->f;
      }

      if( lhssqnew < Sq((1 + fmax)*LHSTOL) ) break;
      lhssq = lhssqnew;
      goto repeat;
    }
  }

#if 0
  if( region->depth == 0 && cutlast > cutfirst ) {
    real mindev = INFTY;
    Cut *mincut;
    for( cut = cutfirst; cut <= cutlast; ++cut ) {
      creal dev = fabs(delta[cut->i]);
      if( dev < mindev ) {
        mindev = dev;
        mincut = cut;
      }
    }
    if( mincut != cutfirst ) Copy(cutfirst, mincut, 1);
    cutlast = cutfirst;
  }
#endif

  for( cut = cutfirst; cut <= cutlast; ++cut ) {
    creal x = region->xmajor[Dim(cut->i)];
    real *b = (real *)region->bounds + cut->i;
    cut->save = *b;
    *b = x + SignedDelta(cut->i);
  }

dissect:
  neval_cut_ += neval_;

  next = region->next;
  for( comp = 0; comp < ncomp_; ++comp ) {
    Errors *e = &errors[comp];
    e->diff = region->result[comp].avg;
    e->spread = e->err = 0;
  }

  Explore(region, &samples_[0], depth -= cutlast - cutfirst + 1, 1);

  b1 = &tmp;
  for( cut = cutfirst; cut <= cutlast; ++cut ) {
    *b1 = tmp;
    b0 = (real *)region->bounds + cut->i;
    b1 = (real *)region->bounds + (cut->i ^ 1);
    tmp = *b1;
    *b1 = *b0;
    *b0 = cut->save;
    Explore(region, &samples_[0], depth++, cut != cutlast);
  }

  nsplit = 0;
  for( reg = region; reg != next; reg = reg->next ) {
    for( comp = 0; comp < ncomp_; ++comp ) {
      cResult *r = &reg->result[comp];
      Errors *e = &errors[comp];
      e->diff -= r->avg;
      e->err += r->err;
      e->spread += Sq(r->spread);
    }
    ++nsplit;
  }

  tmp = 1./nsplit;
  for( comp = 0; comp < ncomp_; ++comp ) {
    Errors *e = &errors[comp];
    e->diff = tmp*fabs(e->diff);
    e->err = (e->err == 0) ? 1 : 1 + e->diff/e->err;
    e->spread = (e->spread == 0) ? 1 : 1 + e->diff/sqrt(e->spread);
  }

  tmp = 1 - tmp;
  for( reg = region; reg != next; reg = reg->next ) {
    for( comp = 0; comp < ncomp_; ++comp ) {
      Result *r = &reg->result[comp];
      cErrors *e = &errors[comp];
      creal c = tmp*e->diff;
      if( r->err > 0 ) r->err = r->err*e->err + c;
      r->spread = r->spread*e->spread + c*samples_[0].neff;
    }
  }
}

