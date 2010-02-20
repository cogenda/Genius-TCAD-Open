/*
	util.c
		Utility functions
		this file is part of Vegas
		last modified 2 Mar 06 th
*/


#include "decl_vegas.h"

static count ndim_, ncomp_;
static number neval_;
static Grid *gridptr_[MAXGRIDS];
static count griddim_[MAXGRIDS];
int EXPORT(vegasnbatch) = 1000;
int EXPORT(vegasgridno) = 0;
char EXPORT(vegasstate)[MAXSTATESIZE] = "";


#define SamplesAlloc(p, n) \
  MemAlloc(p, (n)*((ndim_ + ncomp_ + 1)*sizeof(real) + ndim_*sizeof(bin_t)))


#ifdef DEBUG
#include "debug.h"
#endif

