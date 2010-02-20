/*
	cuba.h
		Prototypes for the Cuba library
		this file is part of Cuba
		last modified 2 Mar 06 th
*/

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*integrand_t)(const int *, const double [],
  const int *, double []);

/* Note: Divonne actually passes a fifth argument, a const int *
   which points to the integration phase.  This is used only rarely
   and most users are confused by the warnings the compiler emits
   if the `correct' prototype is used.  Thus, if you need to access
   this argument, use an explicit cast to integrand_t when invoking 
   Divonne. */


void Vegas(const int ndim, const int ncomp, integrand_t integrand,
  const double epsrel, const double epsabs,
  const int flags, const int mineval, const int maxeval,
  const int nstart, const int nincrease,
  int *neval, int *fail,
  double integral[], double error[], double prob[]);


void Suave(const int ndim, const int ncomp, integrand_t integrand,
  const double epsrel, const double epsabs,
  const int flags, const int mineval, const int maxeval,
  const int nnew, const double flatness,
  int *nregions, int *neval, int *fail,
  double integral[], double error[], double prob[]);


void Divonne(const int ndim, const int ncomp, integrand_t integrand,
  const double epsrel, const double epsabs,
  const int flags, const int mineval, const int maxeval,
  const int key1, const int key2, const int key3, const int maxpass,
  const double border, const double maxchisq, const double mindeviation,
  const int ngiven, const int ldxgiven, double xgiven[],
  const int nextra,
  void (*peakfinder)(const int *, const double [], int *, double []),
  int *nregions, int *neval, int *fail,
  double integral[], double error[], double prob[]);


void Cuhre(const int ndim, const int ncomp, integrand_t integrand,
  const double epsrel, const double epsabs,
  const int flags, const int mineval, const int maxeval,
  const int key,
  int *nregions, int *neval, int *fail,
  double integral[], double error[], double prob[]);


extern int vegasnbatch;
extern int vegasgridno;
extern char vegasstate[128];


#ifdef __cplusplus
}
#endif

