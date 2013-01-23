/*****************************************************************************/
/*              8888888         88888888         88888888                    */
/*            8                8                8                            */
/*           8                 8                8                            */
/*           8                  88888888         88888888                    */
/*           8      8888                8                8                   */
/*            8       8                 8                8                   */
/*              888888         888888888        888888888                    */
/*                                                                           */
/*       A Two-Dimensional General Purpose Semiconductor Simulator.          */
/*                                                                           */
/*  GSS 0.4x                                                                 */
/*  Last update: July 18, 2007                                               */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _jflux1q_h_
#define _jflux1q_h_
#include "mathfunc.h"


inline PetscScalar nmid(PetscScalar Vt, PetscScalar Vc1, PetscScalar Vc2, PetscScalar n1, PetscScalar n2)
{
  return n1*aux2((Vc2-Vc1)/(2*Vt)) + n2*aux2((Vc1-Vc2)/(2*Vt));
}

inline PetscScalar pmid(PetscScalar Vt, PetscScalar Vv1, PetscScalar Vv2, PetscScalar p1, PetscScalar p2)
{
  return p1*aux2((Vv1-Vv2)/(2*Vt)) + p2*aux2((Vv2-Vv1)/(2*Vt));
}

inline AutoDScalar nmid(PetscScalar Vt, const AutoDScalar &Vc1, const AutoDScalar &Vc2, const AutoDScalar &n1, const AutoDScalar &n2)
{
  return n1*aux2((Vc2-Vc1)/(2*Vt)) + n2*aux2((Vc1-Vc2)/(2*Vt));
}

inline AutoDScalar pmid(PetscScalar Vt, const AutoDScalar &Vv1, const AutoDScalar &Vv2, const AutoDScalar &p1, const AutoDScalar &p2)
{
  return p1*aux2((Vv1-Vv2)/(2*Vt)) + p2*aux2((Vv2-Vv1)/(2*Vt));
}


//-----------------------------------------------------------------------------
inline PetscScalar In(PetscScalar Vt,PetscScalar Vc1,PetscScalar Vc2,PetscScalar n1,PetscScalar n2,PetscScalar h)
{
  PetscScalar E    = (Vc2-Vc1)/h;
  PetscScalar n    = n1*aux2((Vc2-Vc1)/(2*Vt)) + n2*aux2((Vc1-Vc2)/(2*Vt));
  PetscScalar dndx = aux1((Vc2-Vc1)/(2*Vt))*(n2-n1)/h;
  return (E*n+Vt*dndx);
}

inline PetscScalar Ip(PetscScalar Vt,PetscScalar Vv1,PetscScalar Vv2,PetscScalar p1,PetscScalar p2,PetscScalar h)
{
  PetscScalar E    = (Vv2-Vv1)/h;
  PetscScalar p    = p1*aux2((Vv1-Vv2)/(2*Vt)) + p2*aux2((Vv2-Vv1)/(2*Vt));
  PetscScalar dpdx = aux1((Vv2-Vv1)/(2*Vt))*(p2-p1)/h;
  return (E*p-Vt*dpdx);
}

inline AutoDScalar In(PetscScalar Vt, const AutoDScalar &Vc1, const AutoDScalar &Vc2, const AutoDScalar &n1, const AutoDScalar &n2, PetscScalar h)
{
  AutoDScalar E    = (Vc2-Vc1)/h;
  AutoDScalar n    = n1*aux2((Vc2-Vc1)/(2*Vt)) + n2*aux2((Vc1-Vc2)/(2*Vt));
  AutoDScalar dndx = aux1((Vc2-Vc1)/(2*Vt))*(n2-n1)/h;
  return (E*n+Vt*dndx);
}

inline AutoDScalar Ip(PetscScalar Vt, const AutoDScalar &Vv1, const AutoDScalar &Vv2, const AutoDScalar &p1, const AutoDScalar &p2, PetscScalar h)
{
  AutoDScalar E    = (Vv2-Vv1)/h;
  AutoDScalar p    = p1*aux2((Vv1-Vv2)/(2*Vt)) + p2*aux2((Vv2-Vv1)/(2*Vt));
  AutoDScalar dpdx = aux1((Vv2-Vv1)/(2*Vt))*(p2-p1)/h;
  return (E*p-Vt*dpdx);
}





//-----------------------------------------------------------------------------
inline PetscScalar In(PetscScalar Vt,PetscScalar dVc,PetscScalar n1,PetscScalar n2,PetscScalar h)
{
  return Vt/h*(n2*bern(-dVc/Vt)-n1*bern(dVc/Vt));
}

inline PetscScalar Ip(PetscScalar Vt,PetscScalar dVv,PetscScalar p1,PetscScalar p2,PetscScalar h)
{
  return Vt/h*(p1*bern(-dVv/Vt)-p2*bern(dVv/Vt));
}

inline AutoDScalar In(PetscScalar Vt,const AutoDScalar &dVc,const AutoDScalar &n1,const AutoDScalar &n2, PetscScalar h)
{
  return Vt/h*(n2*bern(-dVc/Vt)-n1*bern(dVc/Vt));
}

inline AutoDScalar Ip(PetscScalar Vt,const AutoDScalar &dVv,const AutoDScalar &p1,const AutoDScalar &p2, PetscScalar h)
{
  return Vt/h*(p1*bern(-dVv/Vt)-p2*bern(dVv/Vt));
}




//-----------------------------------------------------------------------------
inline PetscScalar qV(PetscScalar q1,PetscScalar q2,PetscScalar h)
{
  PetscScalar x = 0.5*(q1-q2);
  if(x<0)
   return (1-exp(x))/h;
  else
   return (-x-x*x/2)/h;  
}

inline AutoDScalar qV(const AutoDScalar &q1, const AutoDScalar &q2, PetscScalar h)
{
  AutoDScalar x = 0.5*(q1-q2);
  if(x<0)
   return (1-exp(x))/h;
  else
   return (-x-x*x/2)/h;  
}


#endif
