#ifndef __qshep_h__
#define __qshep_h__

#ifdef __cplusplus
extern "C"
{
#endif

  void   QSHEP2   ( int &n, double *x, double *y, double *f, int &nq, int &nw, int &nr,
                    int *lcell, int *lnext, double &xmin, double &ymin, double &dx, double &dy, double &rmax, double *rsq, double *a, int &ier );
  
  double GETNP2   ( double &px, double &py, double *x, double *y, int &nr, int *lcell, int *lnext,
                    double &xmin, double &ymin, double &dx, double &dy, int &np, double &dsq );
  
  void   QSHEP3   ( int &n, double *x, double *y, double *z, double *f, int &nq, int &nw, int &nr,
                    int *lcell, int *lnext, double *xyzmin, double *xyzdel, double &rmax, double *rsq, double *a, int &ier );

  double QS3VAL   ( double &px, double &py, double &pz, int &n, double *x, double *y, double *z, double *f,
                    int &nr, int *lcell, int *lnext, double *xyzmin, double *xyzdel, double &rmax, double *rsq, double *a );

#ifdef __cplusplus
}
#endif


#define  QSHEP2  QSHEP2
#define  GETNP2  GETNP2
#define  QSHEP3  QSHEP3
#define  QS3VAL  QS3VAL

#endif // #ifndef __qshep_h__
