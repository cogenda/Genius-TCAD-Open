/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/

#ifndef __flux3_h__
#define __flux3_h__

// inline function for discrete energy balance equation with SG scheme



#include "mathfunc.h"


/* ----------------------------------------------------------------------------
 * Theta function
 *
 *                                           T2-T1
 *                         Theta(T1,T2) = -------------
 *                                         log(T2/T1)
 */
inline Real Theta(Real T1, Real T2)
{
        Real x = T2/T1-1;
        if(fabs(x)>1e-6)
                return (T2-T1)/log(fabs(T2/T1));
        else
                return T1/(1-0.5*x);
}
inline AutoDScalar Theta(const AutoDScalar &T1, const AutoDScalar &T2)
{
        AutoDScalar x = T2/T1-1;
        if(fabs(x)>1e-6)
                return (T2-T1)/log(fabs(T2/T1));
        else
                return T1/(1-0.5*x);
}

//-----------------------------------------------------------------------------
// FIXME I am very afraid about float exception of exp operator here.
inline Real nmid_eb(Real kb, Real e,  Real V1,  Real V2,
                           Real n1, Real n2, Real Tn1, Real Tn2)
{
  Real Tn = 0.5*(Tn1+Tn2);
  Real Q;
  if( fabs(Tn2-Tn1)/(Tn2+Tn1) > 1e-6 )
    Q = (1-exp((2-e/kb*(V2-V1)/(Tn2-Tn1))*log(Tn1/Tn)))/(1-exp((2-e/kb*(V2-V1)/(Tn2-Tn1))*log(Tn1/Tn2)));
  else
    Q = (1-exp(e*(V2-V1)/(2*kb*Tn1)))/(1-exp(e*(V2-V1)/(kb*Tn1)));
  return n1/Tn1*Tn*(1-Q) + n2/Tn2*Tn*Q;
}

inline AutoDScalar nmid_eb(Real kb, Real e, const AutoDScalar &V1, const AutoDScalar &V2,
                           const AutoDScalar &n1,const AutoDScalar &n2, const AutoDScalar &Tn1,const AutoDScalar &Tn2)
{
  AutoDScalar Tn = 0.5*(Tn1+Tn2);
  AutoDScalar alpha = 2-e/kb*(V2-V1)/(Tn2-Tn1);
  AutoDScalar Q;
  if( fabs(Tn2-Tn1)/(Tn2+Tn1) > 1e-6 )
    Q = (1-exp((2-e/kb*(V2-V1)/(Tn2-Tn1))*log(Tn1/Tn)))/(1-exp((2-e/kb*(V2-V1)/(Tn2-Tn1))*log(Tn1/Tn2)));
  else
    Q = (1-exp(e*(V2-V1)/(2*kb*Tn1)))/(1-exp(e*(V2-V1)/(kb*Tn1)));
  return n1/Tn1*Tn*(1-Q) + n2/Tn2*Tn*Q;
}


//-----------------------------------------------------------------------------

inline Real pmid_eb(Real kb, Real e,  Real V1,  Real V2,
                           Real p1, Real p2, Real Tp1, Real Tp2)
{
  Real Tp = 0.5*(Tp1+Tp2);
  Real Q;
  if( fabs(Tp2-Tp1)/(Tp2+Tp1) > 1e-6 )
    Q = (1-exp((2+e/kb*(V2-V1)/(Tp2-Tp1))*log(Tp1/Tp)))/(1-exp((2+e/kb*(V2-V1)/(Tp2-Tp1))*log(Tp1/Tp2)));
  else
    Q = (1-exp(e*(V1-V2)/(2*kb*Tp1)))/(1-exp(e*(V1-V2)/(kb*Tp1)));
  return p1/Tp1*Tp*(1-Q) + p2/Tp2*Tp*Q;
}

inline AutoDScalar pmid_eb(Real kb, Real e, const AutoDScalar &V1, const AutoDScalar &V2,
                           const AutoDScalar &p1,const AutoDScalar &p2, const AutoDScalar &Tp1,const AutoDScalar &Tp2)
{
  AutoDScalar Tp = 0.5*(Tp1+Tp2);
  AutoDScalar Q;
  if( fabs(Tp2-Tp1)/(Tp2+Tp1) > 1e-6 )
    Q = (1-exp((2+e/kb*(V2-V1)/(Tp2-Tp1))*log(Tp1/Tp)))/(1-exp((2+e/kb*(V2-V1)/(Tp2-Tp1))*log(Tp1/Tp2)));
  else
    Q = (1-exp(e*(V1-V2)/(2*kb*Tp1)))/(1-exp(e*(V1-V2)/(kb*Tp1)));
  return p1/Tp1*Tp*(1-Q) + p2/Tp2*Tp*Q;
}


//-----------------------------------------------------------------------------

inline Real In_eb(Real kb, Real e,  Real V1,  Real V2,
                         Real n1, Real n2, Real Tn1, Real Tn2, Real h)
{
  Real theta = Theta(Tn1,Tn2);
  Real alpha = (e/kb*(V2-V1)-2*(Tn2-Tn1))/theta;
  return kb*0.5*(Tn1+Tn2)*theta*(bern(alpha)*n2/Tn2 - bern(-alpha)*n1/Tn1)/h;
}

inline AutoDScalar In_eb(Real kb, Real e, const AutoDScalar &V1, const AutoDScalar &V2,
                         const AutoDScalar &n1,const AutoDScalar &n2, const AutoDScalar &Tn1,const AutoDScalar &Tn2, Real h)
{
  AutoDScalar theta = Theta(Tn1,Tn2);
  AutoDScalar alpha = (e/kb*(V2-V1)-2*(Tn2-Tn1))/theta;
  return kb*0.5*(Tn1+Tn2)*theta*(bern(alpha)*n2/Tn2 - bern(-alpha)*n1/Tn1)/h;
}



//-----------------------------------------------------------------------------

inline Real Ip_eb(Real kb, Real e,  Real V1,  Real V2,
                         Real p1, Real p2, Real Tp1, Real Tp2, Real h)
{
  Real theta = Theta(Tp1,Tp2);
  Real alpha = (e/kb*(V2-V1)+2*(Tp2-Tp1))/theta;
  return kb*0.5*(Tp1+Tp2)*theta*(bern(alpha)*p1/Tp1 - bern(-alpha)*p2/Tp2)/h;
}

inline AutoDScalar Ip_eb(Real kb, Real e, const AutoDScalar &V1, const AutoDScalar &V2,
                         const AutoDScalar &p1,const AutoDScalar &p2, const AutoDScalar &Tp1,const AutoDScalar &Tp2, Real h)
{
  AutoDScalar theta = Theta(Tp1,Tp2);
  AutoDScalar alpha = (e/kb*(V2-V1)+2*(Tp2-Tp1))/theta;
  return kb*0.5*(Tp1+Tp2)*theta*(bern(alpha)*p1/Tp1 - bern(-alpha)*p2/Tp2)/h;
}




//-----------------------------------------------------------------------------

inline Real Sn_eb(Real kb, Real e,  Real V1,  Real V2,
                         Real n1, Real n2, Real Tn1, Real Tn2, Real h)
{
  Real theta = Theta(Tn1,Tn2);
  Real alpha = (e/kb*(V2-V1)-2*(Tn2-Tn1))/theta;
  Real phi   = (e/kb*(V2-V1)-(Tn2-Tn1))/theta-log(fabs(n2/n1));
  Real Dn    = kb*0.5*(Tn1+Tn2)/e;
  if(alpha > BP4_BERN || 1.25*phi > BP4_BERN)
    return -2.0*kb*Dn/h*theta*( - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*n1);
  return -2.0*kb*Dn/h*theta*(bern(alpha)*bern(1.25*phi)/bern(phi)*n2 - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*n1);
}

inline AutoDScalar Sn_eb(Real kb, Real e, const  AutoDScalar &V1,const  AutoDScalar &V2,
                         const AutoDScalar &n1, const AutoDScalar &n2, const AutoDScalar &Tn1,const AutoDScalar &Tn2, Real h)
{
  AutoDScalar theta = Theta(Tn1,Tn2);
  AutoDScalar alpha = (e/kb*(V2-V1)-2*(Tn2-Tn1))/theta;
  AutoDScalar phi   = (e/kb*(V2-V1)-(Tn2-Tn1))/theta-log(fabs(n2/n1));
  AutoDScalar Dn    = kb*0.5*(Tn1+Tn2)/e;
  if(alpha > BP4_BERN || 1.25*phi > BP4_BERN)
    return -2.0*kb*Dn/h*theta*( - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*n1);
  return -2.0*kb*Dn/h*theta*(bern(alpha)*bern(1.25*phi)/bern(phi)*n2 - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*n1);
}



//-----------------------------------------------------------------------------

inline Real Sp_eb(Real kb, Real e,  Real V1,  Real V2,
                         Real p1, Real p2, Real Tp1, Real Tp2, Real h)
{
  Real theta = Theta(Tp1,Tp2);
  Real alpha = (-e/kb*(V2-V1)-2*(Tp2-Tp1))/theta;
  Real phi   = (-e/kb*(V2-V1)-(Tp2-Tp1))/theta-log(fabs(p2/p1));
  Real Dp    = kb*0.5*(Tp1+Tp2)/e;
  if(alpha > BP4_BERN || 1.25*phi > BP4_BERN)
    return -2.0*kb*Dp/h*theta*(- bern(-alpha)*bern(-1.25*phi)/bern(-phi)*p1);
  return  -2.0*kb*Dp/h*theta*(bern(alpha)*bern(1.25*phi)/bern(phi)*p2 - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*p1);
}

inline AutoDScalar Sp_eb(Real kb, Real e, const  AutoDScalar &V1,const  AutoDScalar &V2,
                         const AutoDScalar &p1, const AutoDScalar &p2, const AutoDScalar &Tp1,const AutoDScalar &Tp2, Real h)
{
  AutoDScalar theta = Theta(Tp1,Tp2);
  AutoDScalar alpha = (-e/kb*(V2-V1)-2*(Tp2-Tp1))/theta;
  AutoDScalar phi   = (-e/kb*(V2-V1)-(Tp2-Tp1))/theta-log(fabs(p2/p1));
  AutoDScalar Dp    = kb*0.5*(Tp1+Tp2)/e;
  if(alpha > BP4_BERN || 1.25*phi > BP4_BERN)
    return -2.0*kb*Dp/h*theta*(- bern(-alpha)*bern(-1.25*phi)/bern(-phi)*p1);
  return   -2.0*kb*Dp/h*theta*(bern(alpha)*bern(1.25*phi)/bern(phi)*p2 - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*p1);
}



#endif // #define __flux3_h__
