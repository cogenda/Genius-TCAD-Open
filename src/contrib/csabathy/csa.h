/******************************************************************************
 *
 * File:           csa.h
 *
 * Created:        16/10/2002
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        A header for csa library (2D data approximation with
 *                 bivariate C1 cubic spline)
 *
 * Revisions:      None
 *
 *****************************************************************************/

#if !defined(_CSA_H)
#define _CSA_H

namespace CSA
{

  struct point
  {
    double x;
    double y;
    double z;
  };

  extern int csa_verbose;

  // predefine
  struct square;
  struct csa;
  
  struct triangle
  {
    square* parent;
    double xc, yc;
    double r;                   /* data visibility radius */

    int nallocated;
    int npoints;
    point** points;
    double** std;
  };

  struct square
  {
    csa* parent;
    int i, j;                   /* indices */

    double xmin, ymin;
    double xc, yc;

    int nallocated;
    int npoints;
    point** points;
    double** std;

    triangle* t;

    int primary;                /* flag -- whether this square contains a
                                     * primary triangle */
    int order;                  /* spline order */

    int hascoeffs[4];           /* flag -- whether there are no NaNs among
                                     * the spline coefficients */

    double coeffs[25];
  };

  struct csa
  {
    double xmin;
    double xmax;
    double ymin;
    double ymax;

    int npoints;
    point** points;
    int npointsallocated;

    int nstd;
    double** std;
    int nstdallocated;

    /*
     * squarization 
     */
    int ni;
    int nj;
    double h;
    square*** squares;          /* square* [j][i] */

    int npt;                    /* Number of Primary Triangles */
    triangle** pt;              /* Primary Triangles -- triangle* [npt] */
    int nincreased;             /* Number of sub-datasets thinned */
    int nthinned;               /* Number of sub-datasets increased */
    int norder[4];              /* Number of fittings of given interpolation
                                     * order */

    /*
     * algorithm parameters 
     */
    int npmin;                  /* minimal number of points locally involved
                                     * in spline calculation (normally = 3) */
    int npmax;                  /* maximal number of points locally involved
                                     * in spline calculation (required > 10,
                                     * recommended 20 < npmax < 60) */
    double k;                   /* relative tolerance multiple in fitting
                                     * spline coefficients: the higher this
                                     * value, the higher degree of the locally
                                     * fitted spline (recommended 80 < k < 200) */
    int nppc;                   /* average number of points per cell */
  };


  csa* csa_create();
  void csa_destroy(csa* a);
  void csa_addpoints(csa* a, int n, point points[]);
  void csa_addstd(csa* a, int n, double variance[]);
  void csa_calculatespline(csa* a);
  void csa_approximatepoint(csa* a, point* p);
  void csa_approximatepoints(csa* a, int n, point* points);

  void csa_setnpmin(csa* a, int npmin);
  void csa_setnpmax(csa* a, int npmax);
  void csa_setk(csa* a, int k);
  void csa_setnppc(csa* a, int nppc);

  void points_read(char* fname, int dim,  int* n, point** points, double** std);
  void parse_commandline(char* arg, int *invY, int *logz, int* invariant, int* square, int* nppc, int* k);

}

#endif
