/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Suave
		last modified 2 Jan 08 th
*/


static int Integrate(creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nnew, creal flatness,
  real *integral, real *error, real *prob)
{
  TYPEDEFREGION;

  count dim, comp, df;
  int fail = 1;
  Result totals[NCOMP];
  Region *anchor = NULL, *region = NULL;

  if( VERBOSE > 1 ) {
    char s[256];
    sprintf(s, "Suave input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  nnew " NUMBER "\n  flatness " REAL,
      ndim_, ncomp_,
      epsrel, epsabs,
      flags, mineval, maxeval,
      nnew, flatness);
    Print(s);
  }

#ifdef MLVERSION
  if( setjmp(abort_) ) goto abort;
#endif

  IniRandom(2*maxeval, flags);

  RegionAlloc(anchor, nnew, nnew);
  anchor->next = NULL;
  anchor->div = 0;

  for( dim = 0; dim < ndim_; ++dim ) {
    Bounds *b = &anchor->bounds[dim];
    b->lower = 0;
    b->upper = 1;
    b->mid = .5;

    if( dim == 0 ) {
      count bin;
      /* define the initial distribution of bins */
      for( bin = 0; bin < NBINS; ++bin )
        b->grid[bin] = (bin + 1)/(real)NBINS;
    }
    else Copy(b->grid, anchor->bounds[0].grid, NBINS);
  }

  Sample(nnew, anchor, anchor->w,
    anchor->w + nnew, anchor->w + (ndim_ + 1)*nnew, flags);
  df = anchor->df;
  ResCopy(totals, anchor->result);

  for( nregions_ = 1; ; ++nregions_ ) {
    Var var[NDIM][2], *vLR;
    real maxratio, maxerr, minfluct, bias, mid;
    Region *regionL, *regionR, *reg, **parent, **par;
    Bounds *bounds, *boundsL, *boundsR;
    count maxcomp, bisectdim;
    number n, nL, nR, nnewL, nnewR;
    real *w, *wL, *wR, *x, *xL, *xR, *f, *fL, *fR, *wlast, *flast;

    if( VERBOSE ) {
      char s[128 + 128*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        nregions_, neval_);

      for( comp = 0; comp < ncomp_; ++comp ) {
        cResult *tot = &totals[comp];
        p += sprintf(p, "\n[" COUNT "] " 
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          comp + 1, tot->avg, tot->err, tot->chisq, df);
      }

      Print(s);
    }

    maxratio = -INFTY;
    maxcomp = 0;
    for( comp = 0; comp < ncomp_; ++comp ) {
      creal ratio = totals[comp].err/MaxErr(totals[comp].avg);
      if( ratio > maxratio ) {
        maxratio = ratio;
        maxcomp = comp;
      }
    }

    if( maxratio <= 1 && neval_ >= mineval ) {
      fail = 0;
      break;
    }

    if( neval_ >= maxeval ) break;

    maxerr = -INFTY;
    parent = &anchor;
    region = anchor;
    for( par = &anchor; (reg = *par); par = &reg->next ) {
      creal err = reg->result[maxcomp].err;
      if( err > maxerr ) {
        maxerr = err;
        parent = par;
        region = reg;
      }
    }

    Fluct(var[0], flatness,
      region->bounds, region->w, region->n, maxcomp,
      region->result[maxcomp].avg, Max(maxerr, epsabs));

    bias = (epsrel < 1e-50) ? 2 :
      Max(pow(2., -(real)region->div/ndim_)/epsrel, 2.);
    minfluct = INFTY;
    bisectdim = 0;
    for( dim = 0; dim < ndim_; ++dim ) {
      cBounds *b = &region->bounds[dim];
      creal fluct = (var[dim][0].fluct + var[dim][1].fluct)*
        (bias - b->upper + b->lower);
      if( fluct < minfluct ) {
        minfluct = fluct;
        bisectdim = dim;
      }
    }

    vLR = var[bisectdim];
    minfluct = vLR[0].fluct + vLR[1].fluct;
    nnewL = IMax(
      (minfluct == 0) ? nnew/2 : (count)(vLR[0].fluct/minfluct*nnew),
      MINSAMPLES );
    nL = vLR[0].n + nnewL;
    nnewR = IMax(nnew - nnewL, MINSAMPLES);
    nR = vLR[1].n + nnewR;

    RegionAlloc(regionL, nL, nnewL);
    RegionAlloc(regionR, nR, nnewR);

    *parent = regionL;
    regionL->next = regionR;
    regionR->next = region->next;
    regionL->div = regionR->div = region->div + 1;

    bounds = &region->bounds[bisectdim];
    mid = bounds->mid;
    n = region->n;
    w = wlast = region->w;  x = w + n;     f = flast = x + n*ndim_;
    wL = regionL->w;        xL = wL + nL;  fL = xL + nL*ndim_;
    wR = regionR->w;        xR = wR + nR;  fR = xR + nR*ndim_;

    while( n-- ) {
      cbool final = (*w < 0);
      if( x[bisectdim] < mid ) {
        if( final && wR > regionR->w ) *(wR - 1) = -fabs(*(wR - 1));
        *wL++ = *w++;
        VecCopy(xL, x);
        xL += ndim_;
        ResCopy(fL, f);
        fL += ncomp_;
      }
      else {
        if( final && wL > regionL->w ) *(wL - 1) = -fabs(*(wL - 1));
        *wR++ = *w++;
        VecCopy(xR, x);
        xR += ndim_;
        ResCopy(fR, f);
        fR += ncomp_;
      }
      x += ndim_;
      f += ncomp_;
      if( n && final ) wlast = w, flast = f;
    }

    Reweight(region->bounds, wlast, flast, f, totals, flags);
    VecCopy(regionL->bounds, region->bounds);
    VecCopy(regionR->bounds, region->bounds);

    boundsL = &regionL->bounds[bisectdim];
    boundsR = &regionR->bounds[bisectdim];
    boundsL->mid = .5*(boundsL->lower + (boundsL->upper = mid));
    boundsR->mid = .5*((boundsR->lower = mid) + boundsR->upper);

    StretchGrid(bounds->grid, boundsL->grid, boundsR->grid);

    Sample(nnewL, regionL, wL, xL, fL, flags);
    Sample(nnewR, regionR, wR, xR, fR, flags);

    df += regionL->df + regionR->df - region->df;

    for( comp = 0; comp < ncomp_; ++comp ) {
      cResult *r = &region->result[comp];
      Result *rL = &regionL->result[comp];
      Result *rR = &regionR->result[comp];
      Result *tot = &totals[comp];
      real diff, sigsq;

      tot->avg += diff = rL->avg + rR->avg - r->avg;

      diff = Sq(.25*diff);
      sigsq = rL->sigsq + rR->sigsq;
      if( sigsq > 0 ) {
        creal c = Sq(1 + sqrt(diff/sigsq));
        rL->sigsq *= c;
        rR->sigsq *= c;
      }
      rL->err = sqrt(rL->sigsq += diff);
      rR->err = sqrt(rR->sigsq += diff);

      tot->sigsq += rL->sigsq + rR->sigsq - r->sigsq;
      tot->err = sqrt(tot->sigsq);

      tot->chisq += rL->chisq + rR->chisq - r->chisq;
    }

    free(region);
    region = NULL;
  }

  for( comp = 0; comp < ncomp_; ++comp ) {
    cResult *tot = &totals[comp];
    integral[comp] = tot->avg;
    error[comp] = tot->err;
    prob[comp] = ChiSquare(tot->chisq, df);
  }

#ifdef MLVERSION
  if( REGIONS ) {
    MLPutFunction(stdlink, "List", 2);
    MLPutFunction(stdlink, "List", nregions_);
    for( region = anchor; region; region = region->next ) {
      real lower[NDIM], upper[NDIM];

      for( dim = 0; dim < ndim_; ++dim ) {
        cBounds *b = &region->bounds[dim];
        lower[dim] = b->lower;
        upper[dim] = b->upper;
      }

      MLPutFunction(stdlink, "Cuba`Suave`region", 4);
      MLPutRealList(stdlink, lower, ndim_);
      MLPutRealList(stdlink, upper, ndim_);

      MLPutFunction(stdlink, "List", ncomp_);
      for( comp = 0; comp < ncomp_; ++comp ) {
        cResult *r = &region->result[comp];
        real res[] = {r->avg, r->err, r->chisq};
        MLPutRealList(stdlink, res, Elements(res));
      }

      MLPutInteger(stdlink, region->df);
    }
  }
#endif

#ifdef MLVERSION
abort:
#endif

  if( region ) free(region);

  while( (region = anchor) ) {
    anchor = anchor->next;
    free(region);
  }

  return fail;
}

