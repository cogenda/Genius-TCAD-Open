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

#ifndef __jflux1e_h__
#define __jflux1e_h__

// inline function for discrete drift-diffusion equation with SG scheme



#include "mathfunc.h"


inline PetscScalar nmid_dd(PetscScalar Vt, PetscScalar Vc1, PetscScalar Vc2, PetscScalar n1, PetscScalar n2)
{
  return n1*aux2((Vc2-Vc1)/(2*Vt)) + n2*aux2((Vc1-Vc2)/(2*Vt));
}

inline PetscScalar pmid_dd(PetscScalar Vt, PetscScalar Vv1, PetscScalar Vv2, PetscScalar p1, PetscScalar p2)
{
  return p1*aux2((Vv1-Vv2)/(2*Vt)) + p2*aux2((Vv2-Vv1)/(2*Vt));
}

inline AutoDScalar nmid_dd(PetscScalar Vt, const AutoDScalar &Vc1, const AutoDScalar &Vc2, const AutoDScalar &n1, const AutoDScalar &n2)
{
  return n1*aux2((Vc2-Vc1)/(2*Vt)) + n2*aux2((Vc1-Vc2)/(2*Vt));
}

inline AutoDScalar pmid_dd(PetscScalar Vt, const AutoDScalar &Vv1, const AutoDScalar &Vv2, const AutoDScalar &p1, const AutoDScalar &p2)
{
  return p1*aux2((Vv1-Vv2)/(2*Vt)) + p2*aux2((Vv2-Vv1)/(2*Vt));
}


//-----------------------------------------------------------------------------
inline PetscScalar In_dd(PetscScalar Vt,PetscScalar Vc1,PetscScalar Vc2,PetscScalar n1,PetscScalar n2,PetscScalar h)
{
  PetscScalar E    = (Vc2-Vc1)/h;
  PetscScalar n    = n1*aux2((Vc2-Vc1)/(2*Vt)) + n2*aux2((Vc1-Vc2)/(2*Vt));
  PetscScalar dndx = aux1((Vc2-Vc1)/(2*Vt))*(n2-n1)/h;
  return (E*n+Vt*dndx);
}

inline PetscScalar Ip_dd(PetscScalar Vt,PetscScalar Vv1,PetscScalar Vv2,PetscScalar p1,PetscScalar p2,PetscScalar h)
{
  PetscScalar E    = (Vv2-Vv1)/h;
  PetscScalar p    = p1*aux2((Vv1-Vv2)/(2*Vt)) + p2*aux2((Vv2-Vv1)/(2*Vt));
  PetscScalar dpdx = aux1((Vv2-Vv1)/(2*Vt))*(p2-p1)/h;
  return (E*p-Vt*dpdx);
}

inline AutoDScalar In_dd(PetscScalar Vt, const AutoDScalar &Vc1, const AutoDScalar &Vc2, const AutoDScalar &n1, const AutoDScalar &n2, PetscScalar h)
{
  AutoDScalar E    = (Vc2-Vc1)/h;
  AutoDScalar n    = n1*aux2((Vc2-Vc1)/(2*Vt)) + n2*aux2((Vc1-Vc2)/(2*Vt));
  AutoDScalar dndx = aux1((Vc2-Vc1)/(2*Vt))*(n2-n1)/h;
  return (E*n+Vt*dndx);
}

inline AutoDScalar Ip_dd(PetscScalar Vt, const AutoDScalar &Vv1, const AutoDScalar &Vv2, const AutoDScalar &p1, const AutoDScalar &p2, PetscScalar h)
{
  AutoDScalar E    = (Vv2-Vv1)/h;
  AutoDScalar p    = p1*aux2((Vv1-Vv2)/(2*Vt)) + p2*aux2((Vv2-Vv1)/(2*Vt));
  AutoDScalar dpdx = aux1((Vv2-Vv1)/(2*Vt))*(p2-p1)/h;
  return (E*p-Vt*dpdx);
}


//-----------------------------------------------------------------------------

inline PetscScalar In_drift(PetscScalar Vt,PetscScalar Vc1,PetscScalar Vc2,PetscScalar n1,PetscScalar n2,PetscScalar h)
{
  PetscScalar E    = (Vc2-Vc1)/h;
  PetscScalar n    = n1*aux2((Vc2-Vc1)/(2*Vt)) + n2*aux2((Vc1-Vc2)/(2*Vt));
  return E*n;
}

inline PetscScalar Ip_drift(PetscScalar Vt,PetscScalar Vv1,PetscScalar Vv2,PetscScalar p1,PetscScalar p2,PetscScalar h)
{
  PetscScalar E    = (Vv2-Vv1)/h;
  PetscScalar p    = p1*aux2((Vv1-Vv2)/(2*Vt)) + p2*aux2((Vv2-Vv1)/(2*Vt));
  return E*p;
}

inline AutoDScalar In_drift(PetscScalar Vt, const AutoDScalar &Vc1, const AutoDScalar &Vc2, const AutoDScalar &n1, const AutoDScalar &n2, PetscScalar h)
{
  AutoDScalar E    = (Vc2-Vc1)/h;
  AutoDScalar n    = n1*aux2((Vc2-Vc1)/(2*Vt)) + n2*aux2((Vc1-Vc2)/(2*Vt));
  return E*n;
}

inline AutoDScalar Ip_drift(PetscScalar Vt, const AutoDScalar &Vv1, const AutoDScalar &Vv2, const AutoDScalar &p1, const AutoDScalar &p2, PetscScalar h)
{
  AutoDScalar E    = (Vv2-Vv1)/h;
  AutoDScalar p    = p1*aux2((Vv1-Vv2)/(2*Vt)) + p2*aux2((Vv2-Vv1)/(2*Vt));
  return E*p;
}


inline PetscScalar In_diffusion(PetscScalar Vt,PetscScalar Vc1,PetscScalar Vc2,PetscScalar n1,PetscScalar n2,PetscScalar h)
{
  PetscScalar dndx = aux1((Vc2-Vc1)/(2*Vt))*(n2-n1)/h;
  return Vt*dndx;
}

inline PetscScalar Ip_diffusion(PetscScalar Vt,PetscScalar Vv1,PetscScalar Vv2,PetscScalar p1,PetscScalar p2,PetscScalar h)
{
  PetscScalar dpdx = aux1((Vv2-Vv1)/(2*Vt))*(p2-p1)/h;
  return -Vt*dpdx;
}

inline AutoDScalar In_diffusion(PetscScalar Vt, const AutoDScalar &Vc1, const AutoDScalar &Vc2, const AutoDScalar &n1, const AutoDScalar &n2, PetscScalar h)
{
  AutoDScalar dndx = aux1((Vc2-Vc1)/(2*Vt))*(n2-n1)/h;
  return Vt*dndx;
}

inline AutoDScalar Ip_diffusion(PetscScalar Vt, const AutoDScalar &Vv1, const AutoDScalar &Vv2, const AutoDScalar &p1, const AutoDScalar &p2, PetscScalar h)
{
  AutoDScalar dpdx = aux1((Vv2-Vv1)/(2*Vt))*(p2-p1)/h;
  return Vt*dpdx;
}



//-----------------------------------------------------------------------------
inline PetscScalar In_dd(PetscScalar Vt,PetscScalar dVc,PetscScalar n1,PetscScalar n2,PetscScalar h)
{
  return Vt*(n2*bern(-dVc/Vt)-n1*bern(dVc/Vt))/h;
}

inline PetscScalar dVIn_dd(PetscScalar Vt,PetscScalar dVc,PetscScalar n1,PetscScalar n2,PetscScalar h)
{
  return (-n2*pd1bern(-dVc/Vt)-n1*pd1bern(dVc/Vt))/h;
}

inline PetscScalar Ip_dd(PetscScalar Vt,PetscScalar dVv,PetscScalar p1,PetscScalar p2,PetscScalar h)
{
  return Vt*(p1*bern(-dVv/Vt)-p2*bern(dVv/Vt))/h;
}

inline PetscScalar dVIp_dd(PetscScalar Vt,PetscScalar dVv,PetscScalar p1,PetscScalar p2,PetscScalar h)
{
  return (-p1*pd1bern(-dVv/Vt)-p2*pd1bern(dVv/Vt))/h;
}

inline AutoDScalar In_dd(PetscScalar Vt,const AutoDScalar &dVc,const AutoDScalar &n1,const AutoDScalar &n2, PetscScalar h)
{
  return Vt*(n2*bern(-dVc/Vt)-n1*bern(dVc/Vt))/h;
}

inline AutoDScalar Ip_dd(PetscScalar Vt,const AutoDScalar &dVv,const AutoDScalar &p1,const AutoDScalar &p2, PetscScalar h)
{
  return Vt*(p1*bern(-dVv/Vt)-p2*bern(dVv/Vt))/h;
}


inline PetscScalar In_uw(PetscScalar ,PetscScalar dVc,PetscScalar n1,PetscScalar n2,PetscScalar h)
{
  if(dVc >0)
    return n2*dVc/h;
  else
    return n1*dVc/h;
}

inline AutoDScalar In_uw(PetscScalar , const AutoDScalar &dVc,const AutoDScalar &n1,const AutoDScalar &n2, PetscScalar h)
{
  if(dVc >0)
    return n2*dVc/h;
  else
    return n1*dVc/h;
}

inline PetscScalar Ip_uw(PetscScalar ,PetscScalar dVv,PetscScalar p1,PetscScalar p2,PetscScalar h)
{
  if(dVv >0)
    return p1*dVv/h;
  else
    return p2*dVv/h;
}

inline AutoDScalar Ip_uw(PetscScalar , const AutoDScalar &dVv,const AutoDScalar &p1,const AutoDScalar &p2, PetscScalar h)
{
  if(dVv >0)
    return p1*dVv/h;
  else
    return p2*dVv/h;
}


#endif // #define __flux1_h__
