/******************************************************************************
 *
 * File:           csa.c
 *
 * Created:        16/10/2002
 *
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *
 * Purpose:        2D data approximation with bivariate C1 cubic spline.
 *                 A set of library functions + standalone utility.
 *
 * Description:    See J. Haber, F. Zeilfelder, O.Davydov and H.-P. Seidel,
 *                 Smooth approximation and rendering of large scattered data
 *                 sets, in  ``Proceedings of IEEE Visualization 2001''
 *                 (Th.Ertl, K.Joy and A.Varshney, Eds.), pp.341-347, 571,
 *                 IEEE Computer Society, 2001.
 *                 http://www.uni-giessen.de/www-Numerische-Mathematik/
 *                        davydov/VIS2001.ps.gz
 *                 http://www.math.uni-mannheim.de/~lsmath4/paper/
 *                        VIS2001.pdf.gz
 *
 * Revisions:      09/04/2003 PS: Modified points_read() to read from a
 *                   file specified by name, not by handle.
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <errno.h>
#include "csa.h"

extern "C"
{
#include "svd.h"
}

namespace CSA
{

  int csa_verbose = 0;

#define NPASTART_S 10           /* Number of Points Allocated at Start for a square */
#define NPASTART_T 100          /* Number of Points Allocated at Start for a triangle */

  /* default algorithm parameters */
#define NPMIN_DEF 3
#define NPMAX_DEF 40
#define K_DEF 150
#define NPPC_DEF 5

#define EPS 1.0e-15

  static double __makeNaN()
  {
#ifdef NAN
    // this is defined in GCC.
    // What to do under win32/msvc?
    return NAN;
#else
    // FIXME: this will raise an divide by zero exception, and cause problem else where.
    double a,b;
    a=0.0;
    b=0.0;
    return a/b;
#endif
  }

  static bool __isNaN(double v)
  {
    return v!=v;
  }

  static void quit(char* format, ...)
  {
    va_list args;

    fflush(stdout);             /* just in case -- to have the exit message
                                     * last */

    fprintf(stderr, "\nerror: csa: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    exit(1);
  }

  /* Allocates n1xn2 matrix of something. Note that it will be accessed as
   * [n2][n1].
   * @param n1 Number of columns
   * @param n2 Number of rows
   * @return Matrix
   */
  static void* alloc2d(int n1, int n2, size_t unitsize)
  {
    unsigned int size;
    char* p;
    char** pp;
    int i;

    if (n1 <= 0 || n2 <= 0)
      quit("alloc2d(): invalid size (n1 = %d, n2 = %d)\n", n1, n2);

    size = n1 * n2;
    if ((p = (char*)calloc(size, unitsize)) == NULL)
      quit("alloc2d(): %s\n", strerror(errno));

    size = n2 * sizeof(void*);
    if ((pp = (char**)malloc(size)) == NULL)
      quit("alloc2d(): %s\n", strerror(errno));
    for (i = 0; i < n2; i++)
      pp[i] = &p[i * n1 * unitsize];

    return pp;
  }

  /* Destroys a matrix.
   * @param pp Matrix
   */
  static void free2d(void* pp)
  {
    void* p;

    p = ((void**) pp)[0];
    free(pp);
    free(p);
  }

  static int str2double(char* token, double* value)
  {
    char* end = NULL;

    if (token == NULL)
    {
      *value = __makeNaN();
      return 0;
    }

    *value = strtod(token, &end);

    if (end == token)
    {
      *value = __makeNaN();
      return 0;
    }

    return 1;
  }

  /*
   * `index' changes from 0 to 3, indicating position of the triangle within
   * the parent square:
   *    -----
   *   |  1  |
   *   |0   2|
   *   |  3  |
   *    -----
   */
  static triangle* triangle_create(square* s)
  {
    triangle* t = (triangle*)malloc(sizeof(triangle));
    double h = s->parent->h;

    t->parent = s;

    t->xc = s->xmin + h / 6.0;
    t->yc = s->ymin + h / 2.0;

    t->r = 0.0;
    t->points = NULL;
    t->std = NULL;
    t->nallocated = 0;
    t->npoints = 0;

    return t;
  }

  static void triangle_addpoint(triangle* t, point* p, double* std)
  {
    if (t->nallocated == t->npoints)
    {
      if (t->nallocated == 0)
      {
        t->points = (point **)malloc(NPASTART_T * sizeof(point*));
        t->nallocated = NPASTART_T;
        if (std != NULL)
          t->std = (double **)malloc(NPASTART_T * sizeof(double*));
      }
      else
      {
        t->nallocated *= 2;
        t->points = (point **)realloc(t->points, t->nallocated * sizeof(point*));
        if (std != NULL)
          t->std = (double **)realloc(t->std, t->nallocated * sizeof(double*));
      }
    }

    t->points[t->npoints] = p;
    if (t->std != NULL)
      t->std[t->npoints] = std;
    t->npoints++;
  }

  static void triangle_destroy(triangle* t)
  {
    if (t->points != NULL)
      free(t->points);
    if (t->std != NULL)
      free(t->std);
    free(t);
  }

  /* Calculates barycentric coordinates of a point.
   * Takes into account that possible triangles are rectangular, with the right
   * angle at t->vertices[0], the vertices[1] vertex being in
   * (-3*PI/4) + (PI/2) * t->index direction from vertices[0], and
   * vertices[2] being at (-5*PI/4) + (PI/2) * t->index.
   */
  static void triangle_calculatebc(square* s, int tindex, point* p, double bc[])
  {
    double dx = p->x - s->xc;
    double dy = p->y - s->yc;
    double h = s->parent->h;

    if (tindex == 0)
    {
      bc[1] = (dy - dx) / h;
      bc[2] = -(dx + dy) / h;
    }
    else if (tindex == 1)
    {
      bc[1] = (dx + dy) / h;
      bc[2] = (dy - dx) / h;
    }
    else if (tindex == 2)
    {
      bc[1] = (dx - dy) / h;
      bc[2] = (dx + dy) / h;
    }
    else
    {
      bc[1] = -(dx + dy) / h;
      bc[2] = (dx - dy) / h;
    }
    bc[0] = 1.0 - bc[1] - bc[2];
  }

  static square* square_create(csa* parent, double xmin, double ymin, int i, int j)
  {
    int ii;

    square* s = (square*)malloc(sizeof(square));
    double h = parent->h;

    s->parent = parent;
    s->i = i;
    s->j = j;

    s->xmin = xmin;
    s->ymin = ymin;
    s->xc = xmin + h / 2.0;
    s->yc = ymin + h / 2.0;

    s->points = NULL;
    s->std = NULL;
    s->nallocated = 0;
    s->npoints = 0;

    s->t = triangle_create(s);

    s->primary = 0;
    s->order = -1;

    for (ii = 0; ii < 4; ++ii)
      s->hascoeffs[ii] = 0;

    for (ii = 0; ii < 25; ++ii)
      s->coeffs[ii] = __makeNaN();

    return s;
  }

  static void square_destroy(square* s)
  {
    if (s->t != NULL)
      triangle_destroy(s->t);
    if (s->points != NULL)
      free(s->points);
    if (s->std != NULL)
      free(s->std);
    free(s);
  }

  static void square_addpoint(square* s, point* p, double* std)
  {
    if (s->nallocated == s->npoints)
    {
      if (s->nallocated == 0)
      {
        s->points = (point**)malloc(NPASTART_S * sizeof(point*));
        if (std != NULL)
          s->std = (double**)malloc(NPASTART_S * sizeof(double*));
        s->nallocated = NPASTART_S;
      }
      else
      {
        s->nallocated *= 2;
        s->points = (point**)realloc(s->points, s->nallocated * sizeof(point*));
        if (std != NULL)
          s->std = (double**)realloc(s->std, s->nallocated * sizeof(double*));
      }
    }

    s->points[s->npoints] = p;
    if (std != NULL)
      s->std[s->npoints] = std;
    s->npoints++;
  }

  csa* csa_create()
  {
    csa* a = (csa*)malloc(sizeof(csa));

    a->xmin = DBL_MAX;
    a->xmax = -DBL_MAX;
    a->ymin = DBL_MAX;
    a->ymax = -DBL_MAX;

    a->npoints = 0;
    a->points = (point**)malloc(NPASTART_T * sizeof(point*));
    a->npointsallocated = NPASTART_T;

    a->std = NULL;
    a->nstd = 0;
    a->nstdallocated = 0;

    a->ni = 0;
    a->nj = 0;
    a->h = __makeNaN();
    a->squares = NULL;

    a->npt = 0;
    a->pt = NULL;
    a->nincreased = 0;
    a->nthinned = 0;
    a->norder[0] = 0;
    a->norder[1] = 0;
    a->norder[2] = 0;
    a->norder[3] = 0;

    a->npmin = NPMIN_DEF;
    a->npmax = NPMAX_DEF;
    a->k = K_DEF;
    a->nppc = NPPC_DEF;

    svd_verbose = (csa_verbose > 1) ? 1 : 0;

    return a;
  }

  void csa_destroy(csa* a)
  {
    int i, j;

    if (a->squares != NULL)
    {
      for (j = 0; j < a->nj; ++j)
        for (i = 0; i < a->ni; ++i)
          square_destroy(a->squares[j][i]);
      free2d(a->squares);
    }
    if (a->pt != NULL)
      free(a->pt);
    if (a->points != NULL)
      free(a->points);
    if (a->std != NULL)
      free(a->std);
    free(a);
  }

  void csa_addpoints(csa* a, int n, point points[])
  {
    int na = a->npointsallocated;
    int i;

    assert(a->squares == NULL);
    /*
     * (can be called prior to squarization only)
     */

    while (na < a->npoints + n)
      na *= 2;
    if (na != a->npointsallocated)
    {
      a->points = (point**)realloc(a->points, na * sizeof(point*));
      a->npointsallocated = na;
    }

    for (i = 0; i < n; ++i)
    {
      point* p = &points[i];

      a->points[a->npoints] = p;
      a->npoints++;

      if (p->x < a->xmin)
        a->xmin = p->x;
      if (p->x > a->xmax)
        a->xmax = p->x;
      if (p->y < a->ymin)
        a->ymin = p->y;
      if (p->y > a->ymax)
        a->ymax = p->y;
    }
  }

  /* Adds std data.
   */
  void csa_addstd(csa* a, int n, double std[])
  {
    int na = a->nstdallocated;
    int i;

    if (std == NULL)
      return;

    if (a->std == NULL)
    {
      na = (n < NPASTART_S) ? NPASTART_S : n;
      a->std = (double**)malloc(na * sizeof(double*));
      a->nstdallocated = na;
    }

    while (na < a->nstd + n)
      na *= 2;
    if (na != a->nstdallocated)
    {
      a->std = (double**)realloc(a->std, na * sizeof(double*));
      a->nstdallocated = na;
    }

    for (i = 0; i < n; ++i)
    {
      assert(std[i] > 0.0);
      a->std[a->nstd] = &std[i];
      a->nstd++;
    }
  }

  /* Marks the squares containing "primary" triangles by setting "primary" flag
   * to 1.
   */
  static void csa_setprimaryflag(csa* a)
  {
    square*** squares = a->squares;
    int nj1 = a->nj - 1;
    int ni1 = a->ni - 1;
    int i, j;

    for (j = 1; j < nj1; ++j)
    {
      for (i = 1; i < ni1; ++i)
      {
        if (squares[j][i]->npoints > 0)
        {
          if ((i + j) % 2 == 0)
          {
            squares[j][i]->primary = 1;
            squares[j - 1][i - 1]->primary = 1;
            squares[j + 1][i - 1]->primary = 1;
            squares[j - 1][i + 1]->primary = 1;
            squares[j + 1][i + 1]->primary = 1;
          }
          else
          {
            squares[j - 1][i]->primary = 1;
            squares[j + 1][i]->primary = 1;
            squares[j][i - 1]->primary = 1;
            squares[j][i + 1]->primary = 1;
          }
        }
      }
    }
  }

  /* Splits the data domain in a number of squares.
   */
  static void csa_squarize(csa* a)
  {
    int nps[7] = { 0, 0, 0, 0, 0, 0 };  /* stats on number of points per
                                             * square */
    double dx = a->xmax - a->xmin;
    double dy = a->ymax - a->ymin;
    int npoints = a->npoints;
    double h;
    int i, j, ii, nadj;

    if (csa_verbose)
    {
      fprintf(stderr, "squarizing:\n");
      fflush(stderr);
    }

    if (npoints == 0)
      return;

    assert(a->squares == NULL);
    /*
     * (can be done only once)
     */

    h = sqrt(dx * dy * a->nppc / npoints);      /* square edge size */
    if (dx < h)
      h = dy * a->nppc / npoints;
    if (dy < h)
      h = dx * a->nppc / npoints;
    a->h = h;

    a->ni = (int) ceil(dx / h) + 2;
    a->nj = (int) ceil(dy / h) + 2;

    if (csa_verbose)
    {
      fprintf(stderr, "  %d x %d squares\n", a->ni, a->nj);
      fflush(stderr);
    }

    /*
     * create squares
     */
    a->squares = (square ***)alloc2d(a->ni, a->nj, sizeof(void*));
    for (j = 0; j < a->nj; ++j)
      for (i = 0; i < a->ni; ++i)
        a->squares[j][i] = square_create(a, a->xmin + h * (i - 1), a->ymin + h * (j - 1), i, j);

    /*
     * map points to squares
     */
    for (ii = 0; ii < npoints; ++ii)
    {
      point* p = a->points[ii];

      i = (int) floor((p->x - a->xmin) / h) + 1;
      j = (int) floor((p->y - a->ymin) / h) + 1;
      square_addpoint(a->squares[j][i], p, (a->std == NULL) ? NULL : a->std[ii]);
    }

    /*
     * mark relevant squares with no points
     */
    csa_setprimaryflag(a);

    /*
     * Create a list of "primary" triangles, for which spline coefficients
     * will be calculated directy (by least squares method), without using
     * C1 smoothness conditions.
     */
    a->pt = (triangle **)malloc((a->ni / 2 + 1) * a->nj * sizeof(triangle*));
    for (j = 0, ii = 0, nadj = 0; j < a->nj; ++j)
    {
      for (i = 0; i < a->ni; ++i)
      {
        square* s = a->squares[j][i];

        if (s->npoints > 0)
        {
          int nn = s->npoints / 5;

          if (nn > 6)
            nn = 6;
          nps[nn]++;
          ii++;
        }
        if (s->primary && s->npoints == 0)
          nadj++;
        if (s->primary)
        {
          a->pt[a->npt] = s->t;
          a->npt++;
        }
      }
    }
    assert(a->npt > 0);

    if (csa_verbose)
    {
      fprintf(stderr, "  %d non-empty squares\n", ii);
      fprintf(stderr, "  %d primary squares\n", a->npt);
      fprintf(stderr, "  %d primary squares with no data\n", nadj);
      fprintf(stderr, "  %.2f points per square \n", (double) a->npoints / ii);
    }

    if (csa_verbose == 2)
    {
      for (i = 0; i < 6; ++i)
        fprintf(stderr, "  %d-%d points -- %d squares\n", i * 5, i * 5 + 4, nps[i]);
      fprintf(stderr, "  %d or more points -- %d squares\n", i * 5, nps[i]);
    }

    if (csa_verbose == 2)
    {
      fprintf(stderr, " j\\i");
      for (i = 0; i < a->ni; ++i)
        fprintf(stderr, "%3d ", i);
      fprintf(stderr, "\n");
      for (j = a->nj - 1; j >= 0; --j)
      {
        fprintf(stderr, "%3d ", j);
        for (i = 0; i < a->ni; ++i)
        {
          square* s = a->squares[j][i];

          if (s->npoints > 0)
            fprintf(stderr, "%3d ", s->npoints);
          else
            fprintf(stderr, "  . ");
        }
        fprintf(stderr, "\n");
      }
    }

    /*
     * all necessary data is now copied to squares;
     * release redundant memory in the host structure
     */
    free(a->points);
    a->points = NULL;
    a->npoints = 0;
    if (a->std != NULL)
    {
      free(a->std);
      a->std = NULL;
    }

    if (csa_verbose)
      fflush(stderr);
  }

  /* Returns all squares intersecting with a square with center in t->middle
   * and edges of length 2*t->r parallel to X and Y axes.
   */
  static void getsquares(csa* a, triangle* t, int* n, square*** squares)
  {
    double h = t->parent->parent->h;
    int imin = (int) floor((t->xc - t->r - a->xmin) / h);
    int imax = (int) ceil((t->xc + t->r - a->xmin) / h);
    int jmin = (int) floor((t->yc - t->r - a->ymin) / h);
    int jmax = (int) ceil((t->yc + t->r - a->ymin) / h);
    int i, j;

    if (imin < 0)
      imin = 0;
    if (imax >= a->ni)
      imax = a->ni - 1;
    if (jmin < 0)
      jmin = 0;
    if (jmax >= a->nj)
      jmax = a->nj - 1;

    *n = 0;
    (*squares) = (square **)malloc((imax - imin + 1) * (jmax - jmin + 1) * sizeof(square*));

    for (j = jmin; j <= jmax; ++j)
    {
      for (i = imin; i <= imax; ++i)
      {
        square* s = a->squares[j][i];

        if (s->npoints > 0)
        {
          (*squares)[*n] = a->squares[j][i];
          (*n)++;
        }
      }
    }
  }

  static double distance(point* p1, point* p2)
  {
    return hypot(p1->x - p2->x, p1->y - p2->y);
  }

  static double distance_xy(double x1, double y1, double x2, double y2)
  {
    return hypot(x1 - x2, y1 - y2);
  }

  /* Thins data by creating an auxiliary regular grid and leaving only the most
   * central point within each grid cell.
   * (I follow the paper here. It is possible that taking average -- in terms of
   * both value and position -- of all points within a cell would be a bit more
   * robust. However, because of keeping only shallow copies of input points,
   * this would require quite a bit of structural changes. So, leaving it as is
   * for now.)
   */
  static void thindata(triangle* t, int npmax)
  {
    csa* a = t->parent->parent;
    int imax = (int) ceil(sqrt((double) (npmax * 3 / 2)));
    square*** squares = (square ***)alloc2d(imax, imax, sizeof(void*));
    double h = t->r * 2.0 / imax;
    double h2 = h / 2.0;
    double xmin = t->xc - t->r;
    double ymin = t->yc - t->r;
    int i, j, ii;

    for (j = 0; j < imax; ++j)
      for (i = 0; i < imax; ++i)
        squares[j][i] = square_create(a, xmin + h * i, ymin + h * j, i, j);

    for (ii = 0; ii < t->npoints; ++ii)
    {
      point* p = t->points[ii];
      int i = (int) floor((p->x - xmin) / h);
      int j = (int) floor((p->y - ymin) / h);
      square* s = squares[j][i];

      if (s->npoints == 0)
        square_addpoint(s, p, (t->std == NULL) ? NULL : t->std[ii]);
      else
      {                  /* npoints == 1 */

        point pmiddle;

        pmiddle.x = xmin + h * i + h2;
        pmiddle.y = ymin + h * j + h2;

        if (s->std == NULL)
        {
          if (distance(s->points[0], &pmiddle) > distance(p, &pmiddle))
            s->points[0] = p;
        }
        else
        {
          if ((*t->std[ii] < *s->std[0]) || (*t->std[ii] == *s->std[0] && distance(s->points[0], &pmiddle) > distance(p, &pmiddle)))
          {
            s->points[0] = p;
            s->std[0] = t->std[ii];
          }
        }
      }
    }

    t->npoints = 0;
    for (j = 0; j < imax; ++j)
    {
      for (i = 0; i < imax; ++i)
      {
        square* s = squares[j][i];

        if (squares[j][i]->npoints != 0)
          triangle_addpoint(t, s->points[0], (s->std == NULL) ? NULL : s->std[0]);
        square_destroy(s);
      }
    }

    free2d(squares);
    imax++;
  }

  /* Finds data points to be used in calculating spline coefficients for each
   * primary triangle.
   */
  static void csa_attachpointstriangle(csa* a, triangle* t)
  {
    int increased = 0;

    if (csa_verbose)
    {
      fprintf(stderr, ".");
      fflush(stderr);
    }

    t->r = a->h * 1.25;
    while (1)
    {
      int nsquares = 0;
      square** squares = NULL;
      int ii;

      getsquares(a, t, &nsquares, &squares);
      for (ii = 0; ii < nsquares; ++ii)
      {
        square* s = squares[ii];
        int iii;

        for (iii = 0; iii < s->npoints; ++iii)
        {
          point* p = s->points[iii];

          if (distance_xy(p->x, p->y, t->xc, t->yc) <= t->r)
            triangle_addpoint(t, p, (s->std == NULL) ? NULL : s->std[iii]);
        }
      }

      free(squares);

      if (t->npoints < a->npmin)
      {
        if (!increased)
        {
          increased = 1;
          a->nincreased++;
        }
        t->r *= 1.25;
        t->npoints = 0;
      }
      else if (t->npoints > a->npmax)
      {
        a->nthinned++;
        thindata(t, a->npmax);
        if (t->npoints > a->npmin)
          break;
        else
        {
          /*
           * Sometimes you have too much data, you thin it and --
           * oops -- you have too little. This is not a frequent
           * event, so let us not bother to put a new subdivision.
           */
          t->r *= 1.25;
          t->npoints = 0;
        }
      }
      else
        break;
    }
  }

  static int n2q(int n)
  {
    assert(n >= 3);

    if (n >= 10)
      return 3;
    else if (n >= 6)
      return 2;
    else                        /* n == 3 */
      return 1;
  }

  /*
   *  square->coeffs[]:
   *
   *   ---------------------
   *  | 3    10    17    24 |
   *  |    6    13    20    |
   *  | 2     9    16    23 |
   *  |    5    12    19    |
   *  | 1     8    15    22 |
   *  |    4    11    18    |
   *  | 0     7    14    21 |
   *   ---------------------
   */

  static void csa_findprimarycoeffstriangle(csa* a, triangle* t)
  {
    square* s = t->parent;
    int npoints = t->npoints;
    point** points = t->points;
    double* z = (double *)malloc(npoints * sizeof(double));
    double* std = NULL;
    int q = n2q(t->npoints);
    int ok = 1;
    double b[10];
    double b1[6];
    int ii;

    if (csa_verbose)
    {
      fprintf(stderr, ".");
      fflush(stderr);
    }

    for (ii = 0; ii < npoints; ++ii)
      z[ii] = points[ii]->z;

    if (t->std != NULL)
    {
      std = (double *)malloc(npoints * sizeof(double));
      for (ii = 0; ii < npoints; ++ii)
        std[ii] = *t->std[ii];
    }

    do
    {
      double bc[3];
      double wmin, wmax;

      if (!ok)
        q--;

      assert(q >= 0);

      if (q == 3)
      {
        double** A = (double **)alloc2d(10, npoints, sizeof(double));
        double w[10];

        for (ii = 0; ii < npoints; ++ii)
        {
          double* aii = A[ii];
          double tmp;

          triangle_calculatebc(s, 0, points[ii], bc);

          /*
           *  0   1   2   3   4   5   6   7   8   9
           * 300 210 201 120 111 102 030 021 012 003
           */
          tmp = bc[0] * bc[0];
          aii[0] = tmp * bc[0];
          tmp *= 3.0;
          aii[1] = tmp * bc[1];
          aii[2] = tmp * bc[2];
          tmp = bc[1] * bc[1];
          aii[6] = tmp * bc[1];
          tmp *= 3.0;
          aii[3] = tmp * bc[0];
          aii[7] = tmp * bc[2];
          tmp = bc[2] * bc[2];
          aii[9] = tmp * bc[2];
          tmp *= 3.0;
          aii[5] = tmp * bc[0];
          aii[8] = tmp * bc[1];
          aii[4] = bc[0] * bc[1] * bc[2] * 6.0;
        }

        svd_lsq(A, 10, npoints, z, std, w, b);

        wmin = w[0];
        wmax = w[0];
        for (ii = 1; ii < 10; ++ii)
        {
          if (w[ii] < wmin)
            wmin = w[ii];
          else if (w[ii] > wmax)
            wmax = w[ii];
        }
        if (wmin < wmax / a->k)
          ok = 0;

        free2d(A);

      }
      else if (q == 2)
      {
        double** A = (double **)alloc2d(6, npoints, sizeof(double));
        double w[6];

        for (ii = 0; ii < npoints; ++ii)
        {
          double* aii = A[ii];

          triangle_calculatebc(s, 0, points[ii], bc);

          /*
           *  0   1   2   3   4   5
           * 200 110 101 020 011 002
           */

          aii[0] = bc[0] * bc[0];
          aii[1] = bc[0] * bc[1] * 2.0;
          aii[2] = bc[0] * bc[2] * 2.0;
          aii[3] = bc[1] * bc[1];
          aii[4] = bc[1] * bc[2] * 2.0;
          aii[5] = bc[2] * bc[2];
        }

        svd_lsq(A, 6, npoints, z, std, w, b1);

        wmin = w[0];
        wmax = w[0];
        for (ii = 1; ii < 6; ++ii)
        {
          if (w[ii] < wmin)
            wmin = w[ii];
          else if (w[ii] > wmax)
            wmax = w[ii];
        }
        if (wmin < wmax / a->k)
          ok = 0;
        else
        {              /* degree raising */
          ok = 1;
          b[0] = b1[0];
          b[1] = (b1[0] + 2.0 * b1[1]) / 3.0;
          b[2] = (b1[0] + 2.0 * b1[2]) / 3.0;
          b[3] = (b1[3] + 2.0 * b1[1]) / 3.0;
          b[4] = (b1[1] + b1[2] + b1[4]) / 3.0;
          b[5] = (b1[5] + 2.0 * b1[2]) / 3.0;
          b[6] = b1[3];
          b[7] = (b1[3] + 2.0 * b1[4]) / 3.0;
          b[8] = (b1[5] + 2.0 * b1[4]) / 3.0;
          b[9] = b1[5];
        }

        free2d(A);

      }
      else if (q == 1)
      {
        double** A = (double **)alloc2d(3, npoints, sizeof(double));
        double w[3];

        for (ii = 0; ii < npoints; ++ii)
        {
          double* aii = A[ii];

          triangle_calculatebc(s, 0, points[ii], bc);

          aii[0] = bc[0];
          aii[1] = bc[1];
          aii[2] = bc[2];
        }

        svd_lsq(A, 3, npoints, z, std, w, b1);

        wmin = w[0];
        wmax = w[0];
        for (ii = 1; ii < 3; ++ii)
        {
          if (w[ii] < wmin)
            wmin = w[ii];
          else if (w[ii] > wmax)
            wmax = w[ii];
        }
        if (wmin < wmax / a->k)
          ok = 0;
        else
        {              /* degree raising */
          ok = 1;
          b[0] = b1[0];
          b[1] = (2.0 * b1[0] + b1[1]) / 3.0;
          b[2] = (2.0 * b1[0] + b1[2]) / 3.0;
          b[3] = (2.0 * b1[1] + b1[0]) / 3.0;
          b[4] = (b1[0] + b1[1] + b1[2]) / 3.0;
          b[5] = (2.0 * b1[2] + b1[0]) / 3.0;
          b[6] = b1[1];
          b[7] = (2.0 * b1[1] + b1[2]) / 3.0;
          b[8] = (2.0 * b1[2] + b1[1]) / 3.0;
          b[9] = b1[2];
        }

        free2d(A);
      }
      else if (q == 0)
      {
        double** A = (double **)alloc2d(1, npoints, sizeof(double));
        double w[1];

        for (ii = 0; ii < npoints; ++ii)
          A[ii][0] = 1.0;

        svd_lsq(A, 1, npoints, z, std, w, b1);

        ok = 1;
        b[0] = b1[0];
        b[1] = b1[0];
        b[2] = b1[0];
        b[3] = b1[0];
        b[4] = b1[0];
        b[5] = b1[0];
        b[6] = b1[0];
        b[7] = b1[0];
        b[8] = b1[0];
        b[9] = b1[0];

        free2d(A);
      }
    }
    while (!ok);

    a->norder[q]++;
    s->order = q;

    {
      square* s = t->parent;
      double* coeffs = s->coeffs;

      coeffs[12] = b[0];
      coeffs[9] = b[1];
      coeffs[6] = b[3];
      coeffs[3] = b[6];
      coeffs[2] = b[7];
      coeffs[1] = b[8];
      coeffs[0] = b[9];
      coeffs[4] = b[5];
      coeffs[8] = b[2];
      coeffs[5] = b[4];
    }

    free(z);
    if (std != NULL)
      free(std);

    if (t->points != NULL)
    {
      free(t->points);
      t->points = NULL;
      t->npoints = 0;
      t->nallocated = 0;
    }
    if (t->std != NULL)
    {
      free(t->std);
      t->std = NULL;
    }
  }

  /* Calculates spline coefficients in each primary triangle by least squares
   * fitting to data attached by csa_attachpointstriangle().
   */
  static void csa_findprimarycoeffs(csa* a)
  {
    int i;

    if (csa_verbose)
      fprintf(stderr, "calculating spline coefficients for primary triangles:\n  ");

    for (i = 0; i < a->npt; ++i)
    {
      triangle* t = a->pt[i];

      csa_attachpointstriangle(a, t);
      csa_findprimarycoeffstriangle(a, t);
    }

    if (csa_verbose)
    {
      fprintf(stderr, "\n  3rd order -- %d sets\n", a->norder[3]);
      fprintf(stderr, "  2nd order -- %d sets\n", a->norder[2]);
      fprintf(stderr, "  1st order -- %d sets\n", a->norder[1]);
      fprintf(stderr, "  0th order -- %d sets\n", a->norder[0]);
      fflush(stderr);
    }

    if (csa_verbose == 2)
    {
      int j;

      fprintf(stderr, " j\\i");
      for (i = 0; i < a->ni; ++i)
        fprintf(stderr, "%2d ", i);
      fprintf(stderr, "\n");
      for (j = a->nj - 1; j >= 0; --j)
      {
        fprintf(stderr, "%2d  ", j);
        for (i = 0; i < a->ni; ++i)
        {
          square* s = a->squares[j][i];

          if (s->primary)
            fprintf(stderr, "%2d ", s->order);
          else
            fprintf(stderr, " . ");
        }
        fprintf(stderr, "\n");
      }
    }
  }

  /* Finds spline coefficients in (adjacent to primary triangles) secondary
   * triangles from C1 smoothness conditions.
   */
  static void csa_findsecondarycoeffs(csa* a)
  {
    square*** squares = a->squares;
    int ni = a->ni;
    int nj = a->nj;
    int ii;

    if (csa_verbose)
    {
      fprintf(stderr, "propagating spline coefficients to the remaining triangles:\n");
      fflush(stderr);
    }

    /*
     * red
     */
    for (ii = 0; ii < a->npt; ++ii)
    {
      triangle* t = a->pt[ii];
      square* s = t->parent;
      int i = s->i;
      int j = s->j;
      double* c = s->coeffs;
      double* cl = (i > 0) ? squares[j][i - 1]->coeffs : NULL;
      double* cb = (j > 0) ? squares[j - 1][i]->coeffs : NULL;
      double* cbl = (i > 0 && j > 0) ? squares[j - 1][i - 1]->coeffs : NULL;
      double* ca = (j < nj - 1) ? squares[j + 1][i]->coeffs : NULL;
      double* cal = (j < nj - 1 && i > 0) ? squares[j + 1][i - 1]->coeffs : NULL;

      c[7] = 2.0 * c[4] - c[1];
      c[11] = 2.0 * c[8] - c[5];
      c[15] = 2.0 * c[12] - c[9];

      c[10] = 2.0 * c[6] - c[2];
      c[13] = 2.0 * c[9] - c[5];
      c[16] = 2.0 * c[12] - c[8];

      c[19] = 2.0 * c[15] - c[11];

      if (cl != NULL)
      {
        cl[21] = c[0];
        cl[22] = c[1];
        cl[23] = c[2];
        cl[24] = c[3];

        cl[18] = c[0] + c[1] - c[4];
        cl[19] = c[1] + c[2] - c[5];
        cl[20] = c[2] + c[3] - c[6];

        cl[17] = 2.0 * cl[20] - cl[23];
        cl[14] = 2.0 * cl[18] - cl[22];
      }

      if (cb != NULL)
      {
        cb[3] = c[0];
        cb[10] = c[7];

        cb[6] = c[0] + c[7] - c[4];
        cb[2] = 2.0 * cb[6] - cb[10];
      }

      if (cbl != NULL)
      {
        cbl[23] = cb[2];
        cbl[24] = cb[3];

        cbl[20] = cb[2] + cb[3] - cb[6];
        cbl[17] = cl[14];
      }

      if (ca != NULL)
      {
        ca[0] = c[3];
        ca[7] = c[10];

        ca[4] = c[3] + c[10] - c[6];
        ca[1] = 2.0 * ca[4] - ca[7];
      }

      if (cal != NULL)
      {
        cal[21] = c[3];
        cal[22] = ca[1];

        cal[18] = ca[0] + ca[1] - ca[4];
        cal[14] = cl[17];
      }
    }

    /*
     * blue
     */
    for (ii = 0; ii < a->npt; ++ii)
    {
      triangle* t = a->pt[ii];
      square* s = t->parent;
      int i = s->i;
      int j = s->j;
      double* c = s->coeffs;
      double* cr = (i < ni - 1) ? squares[j][i + 1]->coeffs : NULL;
      double* car = (i < ni - 1 && j < nj - 1) ? squares[j + 1][i + 1]->coeffs : NULL;
      double* cbr = (i < ni - 1 && j > 0) ? squares[j - 1][i + 1]->coeffs : NULL;

      if (car != NULL)
        cr[13] = car[7] + car[14] - car[11];

      if (cbr != NULL)
        cr[11] = cbr[10] + cbr[17] - cbr[13];

      if (cr != NULL)
        cr[5] = c[22] + c[23] - c[19];
    }

    /*
     * green & yellow
     */
    for (ii = 0; ii < a->npt; ++ii)
    {
      triangle* t = a->pt[ii];
      square* s = t->parent;
      int i = s->i;
      int j = s->j;
      double* cr = (i < ni - 1) ? squares[j][i + 1]->coeffs : NULL;

      if (cr != NULL)
      {
        cr[9] = (cr[5] + cr[13]) / 2.0;
        cr[8] = (cr[5] + cr[11]) / 2.0;
        cr[15] = (cr[11] + cr[19]) / 2.0;
        cr[16] = (cr[13] + cr[19]) / 2.0;
        cr[12] = (cr[8] + cr[16]) / 2.0;
      }
    }

    if (csa_verbose)
    {
      fprintf(stderr, "checking that all coefficients have been set:\n");
      fflush(stderr);
    }

    for (ii = 0; ii < ni * nj; ++ii)
    {
      square* s = squares[0][ii];
      double* c = s->coeffs;
      int i;

      if (s->npoints == 0)
        continue;
      for (i = 0; i < 25; ++i)
        if (__isNaN(c[i]))
          fprintf(stderr, "  squares[%d][%d]->coeffs[%d] = __makeNaN()\n", s->j, s->i, i);
    }
  }

  static int i300[] = { 12, 12, 12, 12 };
  static int i030[] = { 3, 24, 21, 0 };
  static int i003[] = { 0, 3, 24, 21 };
  static int i210[] = { 9, 16, 15, 8 };
  static int i021[] = { 2, 17, 22, 7 };
  static int i102[] = { 4, 6, 20, 18 };
  static int i120[] = { 6, 20, 18, 4 };
  static int i012[] = { 1, 10, 23, 14 };
  static int i201[] = { 8, 9, 16, 15 };
  static int i111[] = { 5, 13, 19, 11 };

  static int* iall[] = { i300, i030, i003, i210, i021, i102, i120, i012, i201, i111 };

  static void csa_sethascoeffsflag(csa* a)
  {
    int i, j;

    for (j = 0; j < a->nj; ++j)
    {
      for (i = 0; i < a->ni; ++i)
      {
        square* s = a->squares[j][i];
        double* coeffs = s->coeffs;
        int ii;

        for (ii = 0; ii < 4; ++ii)
        {
          int cc;

          for (cc = 0; cc < 10; ++cc)
            if (__isNaN(coeffs[iall[cc][ii]]))
              break;
          if (cc == 10)
            s->hascoeffs[ii] = 1;
        }
      }
    }
  }

  void csa_calculatespline(csa* a)
  {
    if (a->std != NULL)
      assert(a->nstd == a->npoints);
    csa_squarize(a);
    csa_findprimarycoeffs(a);
    csa_findsecondarycoeffs(a);
    csa_sethascoeffsflag(a);
  }

  void csa_approximatepoint(csa* a, point* p)
  {
    double h = a->h;
    double ii = (p->x - a->xmin) / h + 1.0;
    double jj = (p->y - a->ymin) / h + 1.0;
    int i, j;
    square* s;
    double fi, fj;
    int ti;
    double bc[3];

    if (fabs(floor(ii) - ii) / h < EPS)
      ii = floor(ii);
    if (fabs(floor(jj) - jj) / h < EPS)
      jj = floor(jj);

    if (ii < 0.0 || jj < 0.0 || ii > (double) a->ni - 1.0 || jj > (double) a->nj - 1.0)
    {
      p->z = 0.0;
      return;
    }

    i = (int) floor(ii);
    j = (int) floor(jj);
    s = a->squares[j][i];
    fi = ii - i;
    fj = jj - j;

    if (fj < fi)
    {
      if (fi + fj < 1.0)
        ti = 3;
      else
        ti = 2;
    }
    else
    {
      if (fi + fj < 1.0)
        ti = 0;
      else
        ti = 1;
    }

    if (!s->hascoeffs[ti])
    {
      p->z = 0.0;
      return;
    }

    triangle_calculatebc(s, ti, p, bc);

    {
      double* c = s->coeffs;
      double bc1 = bc[0];
      double bc2 = bc[1];
      double bc3 = bc[2];
      double tmp1 = bc1 * bc1;
      double tmp2 = bc2 * bc2;
      double tmp3 = bc3 * bc3;

      switch (ti)
      {
      case 0:
        p->z = c[12] * bc1 * tmp1 + c[3] * bc2 * tmp2 + c[0] * bc3 * tmp3 + 3.0 * (c[9] * tmp1 * bc2 + c[2] * tmp2 * bc3 + c[4] * tmp3 * bc1 + c[6] * bc1 * tmp2 + c[1] * bc2 * tmp3 + c[8] * tmp1 * bc3) + 6.0 * c[5] * bc1 * bc2 * bc3;
        break;
      case 1:
        p->z = c[12] * bc1 * tmp1 + c[24] * bc2 * tmp2 + c[3] * bc3 * tmp3 + 3.0 * (c[16] * tmp1 * bc2 + c[17] * tmp2 * bc3 + c[6] * tmp3 * bc1 + c[20] * bc1 * tmp2 + c[10] * bc2 * tmp3 + c[9] * tmp1 * bc3) + 6.0 * c[13] * bc1 * bc2 * bc3;
        break;
      case 2:
        p->z = c[12] * bc1 * tmp1 + c[21] * bc2 * tmp2 + c[24] * bc3 * tmp3 + 3.0 * (c[15] * tmp1 * bc2 + c[22] * tmp2 * bc3 + c[20] * tmp3 * bc1 + c[18] * bc1 * tmp2 + c[23] * bc2 * tmp3 + c[16] * tmp1 * bc3) + 6.0 * c[19] * bc1 * bc2 * bc3;
        break;
      default:               /* 3 */
        p->z = c[12] * bc1 * tmp1 + c[0] * bc2 * tmp2 + c[21] * bc3 * tmp3 + 3.0 * (c[8] * tmp1 * bc2 + c[7] * tmp2 * bc3 + c[18] * tmp3 * bc1 + c[4] * bc1 * tmp2 + c[14] * bc2 * tmp3 + c[15] * tmp1 * bc3) + 6.0 * c[11] * bc1 * bc2 * bc3;
      }
    }
  }

  void csa_approximatepoints(csa* a, int n, point* points)
  {
    int ii;

    for (ii = 0; ii < n; ++ii)
      csa_approximatepoint(a, &points[ii]);
  }

  void csa_setnpmin(csa* a, int npmin)
  {
    a->npmin = npmin;
  }

  void csa_setnpmax(csa* a, int npmax)
  {
    a->npmax = npmax;
  }

  void csa_setk(csa* a, int k)
  {
    a->k = k;
  }

  void csa_setnppc(csa* a, int nppc)
  {
    a->nppc = nppc;
  }

#define BUFSIZE 10240
#define STRBUFSIZE 64
#define NALLOCATED_START 1024
  /* Reads array of points from a columnar file.
   *
   * @param fname File name (can be "stdin" or "-" for standard input)
   * @param dim Number of dimensions (must be 2 or 3)
   * @param n Pointer to number of points (output)
   * @param points Pointer to array of points [*n] (output) (to be freed)
   * @param std Pointer to array of data std (to be freed)
   */
  void points_read(char* fname, int dim, int* n, point** points, double** std)
  {
    FILE* f = NULL;
    int nallocated = NALLOCATED_START;
    char buf[BUFSIZE];
    char seps[] = " ,;\t";
    char* token;
    double stdval = __makeNaN();

    if (dim < 2 || dim > 3)
    {
      *n = 0;
      *points = NULL;
      return;
    }

    if (fname == NULL)
      f = stdin;
    else
    {
      if (strcmp(fname, "stdin") == 0 || strcmp(fname, "-") == 0)
        f = stdin;
      else
      {
        f = fopen(fname, "r");
        if (f == NULL)
          quit("%s: %s\n", fname, strerror(errno));
      }
    }

    *points = (point *)malloc(nallocated * sizeof(point));
    *n = 0;
    while (fgets(buf, BUFSIZE, f) != NULL)
    {
      point* p;

      if (*n == nallocated)
      {
        nallocated *= 2;
        *points = (point *)realloc(*points, nallocated * sizeof(point));
        if (std != NULL && *std != NULL)
          *std = (double *)realloc(*std, nallocated * sizeof(double));
      }

      p = &(*points)[*n];

      if (buf[0] == '#')
        continue;
      if ((token = strtok(buf, seps)) == NULL)
        continue;
      if (!str2double(token, &p->x))
        continue;
      if ((token = strtok(NULL, seps)) == NULL)
        continue;
      if (!str2double(token, &p->y))
        continue;
      if (dim == 2)
        p->z = __makeNaN();
      else
      {
        /*
         * z
         */
        if ((token = strtok(NULL, seps)) == NULL)
          continue;
        if (!str2double(token, &p->z))
          continue;
        /*
         * std
         */
        if (std == NULL)
          continue;

        if (*n == 0)
        {
          if ((token = strtok(NULL, seps)) != NULL)
            *std = (double *)malloc(nallocated * sizeof(double));
          if (csa_verbose)
            fprintf(stderr, "%s std data\n", (token != NULL) ? "found" : "no");
        }

        if (*std != NULL)
        {
          if (*n != 0)
            token = strtok(NULL, seps);
          if (token != NULL && !str2double(token, &stdval))
            quit("%s: could not convert \"%s\" to std\n", fname, token);
          (*std)[*n] = stdval;
        }
      }
      (*n)++;
    }

    if (*n == 0)
    {
      free(*points);
      *points = NULL;
    }
    else
      *points = (point *)realloc(*points, *n * sizeof(point));
    if (std != NULL && *std != NULL)
      *std = (double *)realloc(*std, *n * sizeof(point));

    if (f != stdin)
      if (fclose(f) != 0)
        quit("%s: %s\n", fname, strerror(errno));
  }

  void points_write(int n, point* points)
  {
    int i;

    for (i = 0; i < n; ++i)
    {
      point* p = &points[i];

      printf("%.15g %.15g %.15g\n", p->x, p->y, p->z);
    }
  }

  void parse_commandline(char* arg, int *invY, int *logz, int* invariant, int* square, int* nppc, int* k)
  {
    int i=0,j;
    char workstring[2048];
    for (i = 0; arg[i]!= '\0'; i++)
    {
      /* -- Y- force to inverse Y coordinate */
      if (arg[i] == 'Y')
      {
        if(arg[i+1] == '+')
          *invY = 0;
        if(arg[i+1] == '-')
          *invY = 1;
        i++;
      }
      /* -- e- use log preprocess to z */
      if (arg[i] == 'e')
      {
        *logz=1;
      }
      /* -- l- use linear preprocess to z, against to e */
      if (arg[i] == 'l')
      {
        *logz=0;
      }
      /* -- scale internally so that the minimal ellipse turns into a circle */
      if (arg[i] == 'c')
      {
        *invariant = 1;
        *square = 0;
      }
      /*  -- scale internally so that Xmax - Xmin = Ymax - Ymin */
      if (arg[i] == 's')
      {
        *square = 1;
        *invariant = 0;
      }
      /* -- set the spline sensitivity (default = 140, reduce to get smoother results */
      if (arg[i] == 'k')
      {
        if (((arg[i + 1] >= '0') && (arg[i + 1] <= '9')) )
        {
          j = 0;
          while (((arg[i + 1] >= '0') && (arg[i + 1] <= '9')) )
          {
            i++;
            workstring[j] = arg[i];
            j++;
          }
          workstring[j] = '\0';
          *k = strtol(workstring, (char **) NULL,0);
        }
      }
      /* -- set the average number of points per cell (default = 5)
                   works best for uniform data. Decrease to get smaller
                   cells or increase to get larger cells
             */
      if (arg[i] == 'n')
      {
        if (((arg[i + 1] >= '0') && (arg[i + 1] <= '9')) )
        {
          j = 0;
          while (((arg[i + 1] >= '0') && (arg[i + 1] <= '9')) )
          {
            i++;
            workstring[j] = arg[i];
            j++;
          }
          workstring[j] = '\0';
          *nppc = strtol(workstring, (char **) NULL,0);
        }
      }

    }

  }

}
