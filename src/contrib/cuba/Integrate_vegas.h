/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Vegas
		last modified 17 Dec 07 th
*/


static int Integrate(creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
  real *integral, real *error, real *prob)
{
  real *sample;
  count dim, comp;
  int fail = 1;
  struct {
    count niter;
    number nsamples, neval;
    Cumulants cumul[NCOMP];
    Grid grid[NDIM];
  } state;
  int statemsg = VERBOSE;
  struct stat st;

  if( VERBOSE > 1 ) {
    char s[256];
    sprintf(s, "Vegas input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  nstart " NUMBER "\n  nincrease " NUMBER "\n"
      "  vegasgridno %d\n  vegasstate \"%s\"\n",
      ndim_, ncomp_,
      epsrel, epsabs,
      flags, mineval, maxeval,
      nstart, nincrease,
      EXPORT(vegasgridno), EXPORT(vegasstate));
    Print(s);
  }

#ifdef MLVERSION
  if( setjmp(abort_) ) goto abort;
#endif

  IniRandom(2*maxeval, flags);

  if( *EXPORT(vegasstate) && stat(EXPORT(vegasstate), &st) == 0 &&
      st.st_size == sizeof(state) && (st.st_mode & 0400) ) {
    cint h = open(EXPORT(vegasstate), O_RDONLY);
    read(h, &state, sizeof(state));
    close(h);
    SkipRandom(neval_ = state.neval);

    if( VERBOSE ) {
      char s[256];
      sprintf(s, "\nRestoring state from %s.", EXPORT(vegasstate));
      Print(s);
    }
  }
  else {
    state.niter = 0;
    state.nsamples = nstart;
    Zap(state.cumul);
    GetGrid(state.grid);
  }

  SamplesAlloc(sample, EXPORT(vegasnbatch));

  /* main iteration loop */

  for( ; ; ) {
    number nsamples = state.nsamples;
    creal jacobian = 1./nsamples;
    Grid margsum[NCOMP][NDIM];

    Zap(margsum);

    for( ; nsamples > 0; nsamples -= EXPORT(vegasnbatch) ) {
      cnumber nbatch = IMin(EXPORT(vegasnbatch), nsamples);
      real *w = sample;
      real *x = w + nbatch;
      real *f = x + nbatch*ndim_;
      real *lastf = f + nbatch*ncomp_;
      bin_t *bin = (bin_t *)lastf;

      while( x < f ) {
        real weight = jacobian;

        GetRandom(x);

        for( dim = 0; dim < ndim_; ++dim ) {
          creal pos = *x*NBINS;
          ccount ipos = (count)pos;
          creal prev = (ipos == 0) ? 0 : state.grid[dim][ipos - 1];
          creal diff = state.grid[dim][ipos] - prev; 
          *x++ = prev + (pos - ipos)*diff;
          *bin++ = ipos;
          weight *= diff*NBINS;
        }

        *w++ = weight;
      }

      DoSample(nbatch, sample, w, f);

      w = sample;
      bin = (bin_t *)lastf;

      while( f < lastf ) {
        creal weight = *w++;

        for( comp = 0; comp < ncomp_; ++comp ) {
          real wfun = weight*(*f++);
          if( wfun ) {
            Cumulants *c = &state.cumul[comp];
            Grid *m = margsum[comp];

            c->sum += wfun;
            c->sqsum += wfun *= wfun;
            for( dim = 0; dim < ndim_; ++dim )
              m[dim][bin[dim]] += wfun;
          }
        }

        bin += ndim_;
      }
    }

    fail = 0;

    /* compute the integral and error values */

    for( comp = 0; comp < ncomp_; ++comp ) {
      Cumulants *c = &state.cumul[comp];
      real avg, sigsq;
      real w = Weight(c->sum, c->sqsum, state.nsamples);

      sigsq = 1/(c->weightsum += w);
      avg = sigsq*(c->avgsum += w*c->sum);

      c->avg = LAST ? (sigsq = 1/w, c->sum) : avg;
      c->err = sqrt(sigsq);
      fail |= (c->err > MaxErr(c->avg));

      if( state.niter == 0 ) c->guess = c->sum;
      else {
        c->chisum += w *= c->sum - c->guess;
        c->chisqsum += w*c->sum;
      }
      c->chisq = c->chisqsum - avg*c->chisum;

      c->sum = c->sqsum = 0;
    }

    if( VERBOSE ) {
      char s[128 + 128*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        state.niter + 1, neval_);

      for( comp = 0; comp < ncomp_; ++comp ) {
        cCumulants *c = &state.cumul[comp];
        p += sprintf(p, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          comp + 1, c->avg, c->err, c->chisq, state.niter);
      }

      Print(s);
    }

    if( fail == 0 && neval_ >= mineval ) {
      if( *EXPORT(vegasstate) ) unlink(EXPORT(vegasstate));
      break;
    }

    if( neval_ >= maxeval && *EXPORT(vegasstate) == 0 ) break;

    if( ncomp_ == 1 )
      for( dim = 0; dim < ndim_; ++dim )
        RefineGrid(state.grid[dim], margsum[0][dim], flags);
    else {
      for( dim = 0; dim < ndim_; ++dim ) {
        Grid wmargsum;
        Zap(wmargsum);
        for( comp = 0; comp < ncomp_; ++comp ) {
          real w = state.cumul[comp].avg;
          if( w != 0 ) {
            creal *m = margsum[comp][dim];
            count bin;
            w = 1/Sq(w);
            for( bin = 0; bin < NBINS; ++bin )
              wmargsum[bin] += w*m[bin];
          }
        }
        RefineGrid(state.grid[dim], wmargsum, flags);
      }
    }

    ++state.niter;
    state.nsamples += nincrease;

    if( *EXPORT(vegasstate) ) {
      cint h = creat(EXPORT(vegasstate), 0666);
      if( h != -1 ) {
        state.neval = neval_;
        write(h, &state, sizeof(state));
        close(h);

        if( statemsg ) {
          char s[256];
          sprintf(s, "\nSaving state to %s.", EXPORT(vegasstate));
          Print(s);
          statemsg = false;
        }
      }
      if( neval_ >= maxeval ) break;
    }
  }

  for( comp = 0; comp < ncomp_; ++comp ) {
    cCumulants *c = &state.cumul[comp];
    integral[comp] = c->avg;
    error[comp] = c->err;
    prob[comp] = ChiSquare(c->chisq, state.niter);
  }

#ifdef MLVERSION
abort:
#endif

  free(sample);
  PutGrid(state.grid);

  return fail;
}

