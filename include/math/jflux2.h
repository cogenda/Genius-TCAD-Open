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

#ifndef __jflux2_h__
#define __jflux2_h__

// inline function for discrete lattice temperature corrected drift-diffusion equation with SG scheme



#include "mathfunc.h"


inline Real nmid_lt(Real kb,Real e,Real dV,Real n1,Real n2,Real T,Real dT)
{
  Real Vt = kb*T/e;
  Real alpha = -dV/(2*Vt)+ dT/(2*T);
  return n1*aux2(alpha) + n2*aux2(-alpha);;
}

inline AutoDScalar nmid_lt(Real kb,Real e, const AutoDScalar &dV, const AutoDScalar &n1, const AutoDScalar &n2,
                        const AutoDScalar &T, const AutoDScalar &dT)
{
  AutoDScalar Vt = kb*T/e;
  AutoDScalar alpha = -dV/(2*Vt)+ dT/(2*T);
  return n1*aux2(alpha) + n2*aux2(-alpha);;
}



//-----------------------------------------------------------------------------
inline Real pmid_lt(Real kb,Real e,Real dV,Real p1,Real p2,Real T,Real dT)
{
  Real Vt = kb*T/e;
  Real alpha = -dV/(2*Vt)- dT/(2*T);
  return p1*aux2(-alpha) + p2*aux2(alpha);
}

inline AutoDScalar pmid_lt(Real kb,Real e, const AutoDScalar &dV, const AutoDScalar &p1, const AutoDScalar &p2,
                        const AutoDScalar &T, const AutoDScalar &dT)
{
  AutoDScalar Vt = kb*T/e;
  AutoDScalar alpha = -dV/(2*Vt)- dT/(2*T);
  return p1*aux2(-alpha) + p2*aux2(alpha);
}



//-----------------------------------------------------------------------------

inline Real In_lt(Real kb,Real e,Real dV,Real n1,Real n2,Real T, Real dT,Real h)
{
  Real E  = -dV/h;
  Real Vt = kb*T/e;
  Real alpha = -dV/(2*Vt)+ dT/(2*T);
  Real n  = n1*aux2(alpha) + n2*aux2(-alpha);
  Real dndx = aux1(alpha)*(n2-n1)/h;
  return (E*n + Vt*dndx + kb*n/e*dT/h);
}

inline AutoDScalar In_lt(Real kb,Real e, const AutoDScalar &dV, const AutoDScalar &n1, const AutoDScalar &n2,
                      const AutoDScalar &T, const AutoDScalar &dT, Real h)
{
  AutoDScalar E  = -dV/h;
  AutoDScalar Vt = kb*T/e;
  AutoDScalar alpha = -dV/(2*Vt)+ dT/(2*T);
  AutoDScalar n  = n1*aux2(alpha) + n2*aux2(-alpha);
  AutoDScalar dndx = aux1(alpha)*(n2-n1)/h;
  return (E*n + Vt*dndx + kb*n/e*dT/h);
}



//-----------------------------------------------------------------------------

inline Real Ip_lt(Real kb,Real e,Real dV,Real p1,Real p2,Real T, Real dT,Real h)
{
  Real E  = -dV/h;
  Real Vt = kb*T/e;
  Real alpha = -dV/(2*Vt)- dT/(2*T);
  Real p  = p1*aux2(-alpha) + p2*aux2(alpha);
  Real dpdx = aux1(alpha)*(p2-p1)/h;
  return (E*p-Vt*dpdx - kb*p/e*dT/h);
}

inline AutoDScalar Ip_lt(Real kb,Real e, const AutoDScalar &dV, const AutoDScalar &p1, const AutoDScalar &p2,
                      const AutoDScalar &T, const AutoDScalar &dT,Real h)
{
  AutoDScalar E  = -dV/h;
  AutoDScalar Vt = kb*T/e;
  AutoDScalar alpha = -dV/(2*Vt)- dT/(2*T);
  AutoDScalar p  = p1*aux2(-alpha) + p2*aux2(alpha);
  AutoDScalar dpdx = aux1(alpha)*(p2-p1)/h;
  return (E*p-Vt*dpdx - kb*p/e*dT/h);
}


#endif // #define __flux2_h__
