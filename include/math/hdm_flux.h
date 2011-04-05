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

#ifndef __hdm_flux_h__
#define __hdm_flux_h__

#include "vector_value.h"


class HDMVector
{
  public:
    /**
     * empty constructor
     */
    HDMVector() {}

    /**
     * construct from solution array
     */
    HDMVector(const Real *array, unsigned int size=4)
    {
      for(unsigned int n=0; n<size; ++n)
        _v[n] = array[n];
    }

    /**
     * access density
     */
    Real   density () const                         { return _v[0]; }
    Real & density()                                { return _v[0]; }
    void   set_density(Real d)                      { _v[0] = d; }

    /**
     * access velocity
     */
    Real   speed()   const                          { return sqrt(_v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3]); }
    VectorValue<Real> velocity() const              { return VectorValue<Real>(_v[1], _v[2], _v[3]); }
    void set_velocity(const VectorValue<Real> &ve)  { _v[1]=ve[0]; _v[2]=ve[1]; _v[3]=ve[2]; }
    void add_velocity(const VectorValue<Real> &ve)  { _v[1]+=ve[0]; _v[2]+=ve[1]; _v[3]+=ve[2]; }

    /**
     * access energy
     */
    Real   energy () const                          { return _v[4]; }
    Real & energy()                                 { return _v[4]; }
    void   set_energy(Real e)                       { _v[4] = e; }

    /**
     * access it as an array
     */
    Real & operator [] (unsigned int n)             {return _v[n];}
    const Real & operator [] (unsigned int n) const {return _v[n];}

    /**
     * scale each component
     */
    HDMVector & operator *= (Real scale)
    {
      for(unsigned int n=0; n<4; ++n)
        _v[n] *= scale;
      return *this;
    }

    /**
     * mirror the HDM vector, Vg=V-2(V*n)n,
     * useful for processing solid wall boundary by ghost cell method
     */
    void mirror(const VectorValue<Real> & norm )
    {
      VectorValue<Real> v(_v[1], _v[2], _v[3]);
      v -= 2*(v*norm)*norm;
      set_velocity(v);
    }

    void print(std::ostream& os) const
    {
      os << "n=" << _v[0] << '\n';
      os << "nv=(" << _v[1] << ',' << _v[2] << ',' << _v[3] << ')' << '\n';
      //os << "nw=" << _v[4] << '\n';
    }

    friend std::ostream& operator << (std::ostream& os, const HDMVector& t)
    {
      t.print(os);
      return os;
    }

  private:

    Real   _v[5];
};


// inline function for discrete hydrodynamic equation

/**
 * flux function for AUSM+ scheme
 */

inline Real MPolynomial2plus(Real M)
{    return 0.25*(M+1)*(M+1);  }

inline Real MPolynomial2neg(Real M)
{    return -0.25*(M-1)*(M-1);  }

inline Real MPolynomial4plus(Real M)
{
  if(fabs(M)>=1)
    return 0.5*(M+fabs(M));
  else
    return MPolynomial2plus(M)*(1-2*MPolynomial2neg(M));
}

inline Real MPolynomial4neg(Real M)
{
  if(fabs(M)>=1)
    return 0.5*(M-fabs(M));
  else
    return MPolynomial2neg(M)*(1+2*MPolynomial2plus(M));
}

inline Real PPolynomial5plus(Real M)
{
  if(fabs(M)>=1)
    return 1/M*0.5*(M+fabs(M));
  else
    return MPolynomial2plus(M)*(2-M-3*M*MPolynomial2neg(M));
}

inline Real PPolynomial5neg(Real M)
{
  if(fabs(M)>=1)
    return 1/M*0.5*(M-fabs(M));
  else
    return MPolynomial2neg(M)*(-2-M+3*M*MPolynomial2plus(M));
}


/**
 * AUSM+ flux at interface
 * @param U1 hydrodynamic variable at left
 * @param U2 hydrodynamic variable at right
 * @param m  mass of carrier
 * @param kt Kb*T, we don't consider carrier tempeature
 * @param norm the norm of interface from left to right
 * @param d  distance from left to right
 * @param S  control volumn surface area
 * @param dt max time step
 *
 */
inline HDMVector AUSM_if_flux(const HDMVector &U1, const HDMVector &U2,
                              Real m, Real kT, Real d, Real S, const Point &norm, Real &dt)
{
  HDMVector flux;

  Real  n1 = U1.density();  //
  Real  v1 = U1.speed()/(m*n1);
  Real  p1 = kT*n1;

  Real  n2 = U2.density();  //
  Real  v2 = U2.speed()/(m*n2);
  Real  p2 = kT*n2;

  Real  a  = sqrt(5.0/3.0*kT/m);

  dt=d/(a+std::max(v1,v2));

  Real M1 = (U1.velocity()*norm)/n1/m/a;
  Real M2 = (U2.velocity()*norm)/n2/m/a;

  Real Mavg = 0.5*(M1*M1+M2*M2);
  Real M  = MPolynomial4plus(M1) + MPolynomial4neg(M2)-std::max(1.0-0.8*Mavg,0.0)*(p2-p1)/(p2+p1);
  Real massflow  = M>0? a*M*n1: a*M*n2;

  if(massflow > 0)
  {
    flux.set_density(massflow);
    flux.set_velocity(U1.velocity()/n1*massflow);
  }
  else
  {
    flux.set_density(massflow);
    flux.set_velocity(U2.velocity()/n2*massflow);
  }

  Real p =  PPolynomial5plus(M1)*p1 + PPolynomial5neg(M2)*p2
           -PPolynomial5plus(M1)*PPolynomial5neg(M2)*m*sqrt(n1*n2)*a*a*(M2-M1);

  flux.add_velocity(p*norm);
  flux *= S;

  return flux;
}


#if 0
inline HDMVector AUSM_if_flux(const HDMVector &U1, const HDMVector &U2,
                              Real m, Real k, Real d, Real S, const Point &norm, Real &dt)
{
  HDMVector flux;

  Real  n1 = U1.density();  //
  Real  v1 = U1.speed()/(m*n1);
  Real  p1 = (2.0/3.0*(U1.energy() - 0.5*m*n1*v1*v1));
  Real  a1 = sqrt((5.0/3.0*p1/(n1*m)));
  Real  H1 = U1.energy()/n1+p1/n1;

  Real  n2 = U2.density();  //
  Real  v2 = U2.speed()/(m*n2);
  Real  p2 = (2.0/3.0*(U2.energy() - 0.5*m*n2*v2*v2));
  Real  a2 = sqrt((5.0/3.0*p2/(n2*m)));
  Real  H2 = U2.energy()/n2+p2/n2;

  Real  a  = sqrt(a1*a2);

  dt=d/(a+std::max(v1,v2));

  Real M1 = (U1.velocity()*norm)/n1/m/a;
  Real M2 = (U2.velocity()*norm)/n2/m/a;

  Real Mavg = 0.5*(M1*M1+M2*M2);
  Real M  = MPolynomial4plus(M1) + MPolynomial4neg(M2)-std::max(1.0-0.8*Mavg,0.0)*(p2-p1)/(p2+p1);
  Real massflow  = M>0? a*M*n1: a*M*n2;

  if(massflow > 0)
  {
    flux.set_density(massflow);
    flux.set_velocity(U1.velocity()/n1*massflow);
    flux.set_energy(H1*massflow);
  }
  else
  {
    flux.set_density(massflow);
    flux.set_velocity(U2.velocity()/n2*massflow);
    flux.set_energy(H2*massflow);
  }

  Real p =  PPolynomial5plus(M1)*p1 + PPolynomial5neg(M2)*p2
      -PPolynomial5plus(M1)*PPolynomial5neg(M2)*m*sqrt(n1*n2)*a*a*(M2-M1);

  flux.add_velocity(p*norm);
  flux *= S;

  return flux;
}
#endif

#endif
