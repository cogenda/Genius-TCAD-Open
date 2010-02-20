#ifndef __qshep_h__
#define __qshep_h__

#ifdef __cplusplus
extern "C"
{
#endif

  void   qshep2_   ( int &n, double *x, double *y, double *f, int &nq, int &nw, int &nr,
                     int *lcell, int *lnext, double &xmin, double &ymin, double &dx, double &dy, double &rmax, double *rsq, double *a, int &ier );
  
  double getnp2_   ( double &px, double &py, double *x, double *y, int &nr, int *lcell, int *lnext,
                     double &xmin, double &ymin, double &dx, double &dy, int &np, double &dsq );
  
  void   qshep3_   ( int &n, double *x, double *y, double *z, double *f, int &nq, int &nw, int &nr,
                     int *lcell, int *lnext, double *xyzmin, double *xyzdel, double &rmax, double *rsq, double *a, int &ier );

  double qs3val_   ( double &px, double &py, double &pz, int &n, double *x, double *y, double *z, double *f,
                     int &nr, int *lcell, int *lnext, double *xyzmin, double *xyzdel, double &rmax, double *rsq, double *a );

#ifdef __cplusplus
}
#endif


#define  QSHEP2  qshep2_
#define  GETNP2  getnp2_
#define  QSHEP3  qshep3_
#define  QS3VAL  qs3val_

#endif // #ifndef __qshep_h__
