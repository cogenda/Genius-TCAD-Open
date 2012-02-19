/*----------------------------------------------------------------------------
 This file is borrowed from ADOL-C tapless version and slightly modified.
 The original information is:

 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     adouble.h
 Revision: $Id: adolc.h,v 1.5 2008/07/07 08:44:19 gdiso Exp $
 Contents: adouble.h contains the basis for the class of adouble
           included here are all the possible functions defined on
           the adouble class.  Notice that, as opposed to ealier versions,
           both the class adub and the class adouble are derived from a base
           class (badouble).  See below for further explanation.

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing

 The author promises that I can use adouble.h in gss code covered by BSD license.
 > There is a license problem I should write to you.
 > My code is covered by BSD license, but ADOLC is CPL. I don't like to
 > change BSD license.
 > As a resuit, I'd like to get your written permission of using adouble.h
 > in my code.

 Is it OK for you, if I just confirm it in this email? Then:
 Hereby I confirm that you can use adouble.h for your project.
 -- Andrea Walther

 History:
 ----->   20041110 kowarz: tapeless (scalar/vector) forward version added
          20040423 kowarz: adapted to configure - make - make install
          20000107 olvo:   iostream.h instaed of stream.h
          19991210 olvo:   checking the changes
          19991122 olvo:   new op_codes eq_plus_prod eq_min_prod
                           for  y += x1 * x2
                           and  y -= x1 * x2
          19981201 olvo:   last check:
                           - taputil things changed, includes
          19980820 olvo:   new comparison strategy & some inlines
          19980709 olvo:   modified sign operators
----------------------------------------------------------------------------*/


#ifndef _adolc_h_
#define _adolc_h_


#include <cassert>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

#include "asinh.hpp" // for asinh
#include "acosh.hpp" // for acosh
#include "atanh.hpp" // for atanh

#include "genius_dll.h"

// the max number of independent variable may be used by Genius
// for EBM solver on Hex8 element, there are 6*8=48 independent variables
// with 8 more independent variables for reserve purpose
// this is not a flexible method,
// However, using std::vector (or even new adval array) instead of fxied length array makes system performance greatly slow done.
#define ADTL_NUMBER_DIRECTIONS 56


extern "C"
{
  /**
   * function for set the static adtl::AutoDScalar::numdir
   */
  DLL_EXPORT_DECLARE  void  set_ad_number(const unsigned int p);
}

typedef double PetscScalar;


namespace adtl
{

  class AutoDScalar
  {
  public:
    // ctors
    inline AutoDScalar();
    inline AutoDScalar(const PetscScalar v);
    inline AutoDScalar(const PetscScalar v, const PetscScalar * adv);
    inline AutoDScalar(const PetscScalar v, const PetscScalar * adv, unsigned int n);
    inline AutoDScalar(const AutoDScalar& a);
    inline AutoDScalar(const AutoDScalar& a, unsigned int *, unsigned int n);

    /*******************  temporary results  ******************************/
    // sign
    inline const AutoDScalar operator - () const;
    inline const AutoDScalar operator + () const;

    // addition
    inline const AutoDScalar operator + (const PetscScalar v) const;
    inline const AutoDScalar operator + (const AutoDScalar& a) const;
    inline friend
    const AutoDScalar operator + (const PetscScalar v, const AutoDScalar& a);

    // substraction
    inline const AutoDScalar operator - (const PetscScalar v) const;
    inline const AutoDScalar operator - (const AutoDScalar& a) const;
    inline friend
    const AutoDScalar operator - (const PetscScalar v, const AutoDScalar& a);

    // multiplication
    inline const AutoDScalar operator * (const PetscScalar v) const;
    inline const AutoDScalar operator * (const AutoDScalar& a) const;
    inline friend
    const AutoDScalar operator * (const PetscScalar v, const AutoDScalar& a);

    // division
    inline const AutoDScalar operator / (const PetscScalar v) const;
    inline const AutoDScalar operator / (const AutoDScalar& a) const;
    inline friend
    const AutoDScalar operator / (const PetscScalar v, const AutoDScalar& a);

    // inc/dec
    inline const AutoDScalar operator ++ ();
    inline const AutoDScalar operator ++ (int);
    inline const AutoDScalar operator -- ();
    inline const AutoDScalar operator -- (int);

    // functions
    inline friend const AutoDScalar tan(const AutoDScalar &a);
    inline friend const AutoDScalar exp(const AutoDScalar &a);
    inline friend const AutoDScalar log(const AutoDScalar &a);
    inline friend const AutoDScalar sqrt(const AutoDScalar &a);
    inline friend const AutoDScalar sin(const AutoDScalar &a);
    inline friend const AutoDScalar cos(const AutoDScalar &a);
    inline friend const AutoDScalar asin(const AutoDScalar &a);
    inline friend const AutoDScalar acos(const AutoDScalar &a);
    inline friend const AutoDScalar atan(const AutoDScalar &a);

    inline friend const AutoDScalar atan2(const AutoDScalar &a, const AutoDScalar &b);
    inline friend const AutoDScalar pow(const AutoDScalar &a, PetscScalar v);
    inline friend const AutoDScalar pow(const AutoDScalar &a, const AutoDScalar &b);
    inline friend const AutoDScalar pow(PetscScalar v, const AutoDScalar &a);
    inline friend const AutoDScalar log10(const AutoDScalar &a);

    inline friend const AutoDScalar sinh (const AutoDScalar &a);
    inline friend const AutoDScalar cosh (const AutoDScalar &a);
    inline friend const AutoDScalar tanh (const AutoDScalar &a);

    inline friend const AutoDScalar asinh (const AutoDScalar &a);
    inline friend const AutoDScalar acosh (const AutoDScalar &a);
    inline friend const AutoDScalar atanh (const AutoDScalar &a);

    inline friend const AutoDScalar fabs (const AutoDScalar &a);
    inline friend const AutoDScalar ceil (const AutoDScalar &a);
    inline friend const AutoDScalar floor (const AutoDScalar &a);
    inline friend const AutoDScalar fmax (const AutoDScalar &a, const AutoDScalar &b);
    inline friend const AutoDScalar fmax (PetscScalar v, const AutoDScalar &a);
    inline friend const AutoDScalar fmax (const AutoDScalar &a, PetscScalar v);
    inline friend const AutoDScalar fmin (const AutoDScalar &a, const AutoDScalar &b);
    inline friend const AutoDScalar fmin (PetscScalar v, const AutoDScalar &a);
    inline friend const AutoDScalar fmin (const AutoDScalar &a, PetscScalar v);
    inline friend const AutoDScalar ldexp (const AutoDScalar &a, const AutoDScalar &b);
    inline friend const AutoDScalar ldexp (const AutoDScalar &a, const PetscScalar v);
    inline friend const AutoDScalar ldexp (const PetscScalar v, const AutoDScalar &a);
    inline friend       PetscScalar frexp (const AutoDScalar &a, int* v);
#ifndef WINDOWS
    inline friend const AutoDScalar erf (const AutoDScalar &a);
#endif


    /*******************  nontemporary results  ***************************/
    // assignment
    inline void operator = (const PetscScalar v);
    inline void operator = (const AutoDScalar& a);

    // addition
    inline void operator += (const PetscScalar v);
    inline void operator += (const AutoDScalar& a);

    // substraction
    inline void operator -= (const PetscScalar v);
    inline void operator -= (const AutoDScalar& a);

    // multiplication
    inline void operator *= (const PetscScalar v);
    inline void operator *= (const AutoDScalar& a);

    // division
    inline void operator /= (const PetscScalar v);
    inline void operator /= (const AutoDScalar& a);

    // not
    inline int operator ! () const;

    // comparision
    inline int operator != (const AutoDScalar&) const;
    inline int operator != (const PetscScalar) const;
    inline friend int operator != (const PetscScalar, const AutoDScalar&);

    inline int operator == (const AutoDScalar&) const;
    inline int operator == (const PetscScalar) const;
    inline friend int operator == (const PetscScalar, const AutoDScalar&);

    inline int operator <= (const AutoDScalar&) const;
    inline int operator <= (const PetscScalar) const;
    inline friend int operator <= (const PetscScalar, const AutoDScalar&);

    inline int operator >= (const AutoDScalar&) const;
    inline int operator >= (const PetscScalar) const;
    inline friend int operator >= (const PetscScalar, const AutoDScalar&);

    inline int operator >  (const AutoDScalar&) const;
    inline int operator >  (const PetscScalar) const;
    inline friend int operator >  (const PetscScalar, const AutoDScalar&);

    inline int operator <  (const AutoDScalar&) const;
    inline int operator <  (const PetscScalar) const;
    inline friend int operator <  (const PetscScalar, const AutoDScalar&);

    /*******************  getter / setter  ********************************/
    inline PetscScalar getValue() const;
    inline void setValue(const PetscScalar v);
    inline const PetscScalar * getADValue() const;
    inline std::vector<PetscScalar> getADValueVector() const;
    inline void setADValue(const PetscScalar * v);
#if defined(ADTL_NUMBER_DIRECTIONS)
    inline PetscScalar getADValue(const unsigned int p) const;
    inline void setADValue(const unsigned int p, const PetscScalar v);
#endif
    /*******************  i/o operations  *********************************/
    inline friend std::ostream& operator << ( std::ostream&, const AutoDScalar& );
    inline friend std::istream& operator >> ( std::istream&, AutoDScalar& );

    static unsigned int numdir;
    static void setNumDir(const unsigned int p)
    {
      if (p>ADTL_NUMBER_DIRECTIONS) numdir=ADTL_NUMBER_DIRECTIONS;
      else numdir=p;
    }

  private:
    // internal variables

    PetscScalar val;
    PetscScalar adval[ADTL_NUMBER_DIRECTIONS];
  };

  /*******************************  ctors  ************************************/
  AutoDScalar::AutoDScalar(): val(0)
  {
    memset(adval, 0, sizeof(PetscScalar)*ADTL_NUMBER_DIRECTIONS);
  }

  AutoDScalar::AutoDScalar(const PetscScalar v) : val(v)
  {
    memset(adval, 0, sizeof(PetscScalar)*ADTL_NUMBER_DIRECTIONS);
  }

  AutoDScalar::AutoDScalar(const PetscScalar v, const PetscScalar * adv) : val(v)
  {
    memset(adval, 0, sizeof(PetscScalar)*ADTL_NUMBER_DIRECTIONS);

    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]=adv[_i];
  }

  AutoDScalar::AutoDScalar(const PetscScalar v, const PetscScalar * adv, unsigned int n) : val(v)
  {
    memset(adval, 0, sizeof(PetscScalar)*ADTL_NUMBER_DIRECTIONS);

    for (unsigned int _i=0; _i<n; ++_i)
      adval[_i]=adv[_i];
  }

  AutoDScalar::AutoDScalar(const AutoDScalar& a) : val(a.val)
  {
    memset(adval, 0, sizeof(PetscScalar)*ADTL_NUMBER_DIRECTIONS);

    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]=a.adval[_i];
  }

  AutoDScalar::AutoDScalar(const AutoDScalar& a, unsigned int *order, unsigned int n) : val(a.val)
  {
    memset(adval, 0, sizeof(PetscScalar)*ADTL_NUMBER_DIRECTIONS);

    for (unsigned int _i=0; _i<n; ++_i)
      adval[order[_i]]=a.adval[_i];
  }


  /*************************  temporary results  ******************************/
  // sign
  const AutoDScalar AutoDScalar::operator - () const
  {
    AutoDScalar tmp;
    tmp.val=-val;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=-adval[_i];
    return tmp;
  }

  const AutoDScalar AutoDScalar::operator + () const
  {
    return *this;
  }

  // addition
  const AutoDScalar AutoDScalar::operator + (const PetscScalar v) const
  {
    return AutoDScalar(val+v, adval);
  }

  const AutoDScalar AutoDScalar::operator + (const AutoDScalar& a) const
  {
    AutoDScalar tmp;
    tmp.val=val+a.val;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=adval[_i]+a.adval[_i];
    return tmp;
  }

  const AutoDScalar operator + (const PetscScalar v, const AutoDScalar& a)
  {
    return AutoDScalar(v+a.val, a.adval);
  }

  // subtraction
  const AutoDScalar AutoDScalar::operator - (const PetscScalar v) const
  {
    return AutoDScalar(val-v, adval);
  }

  const AutoDScalar AutoDScalar::operator - (const AutoDScalar& a) const
  {
    AutoDScalar tmp;
    tmp.val=val-a.val;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=adval[_i]-a.adval[_i];
    return tmp;
  }

  const AutoDScalar operator - (const PetscScalar v, const AutoDScalar& a)
  {
    AutoDScalar tmp;
    tmp.val=v-a.val;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=-a.adval[_i];
    return tmp;
  }

  // multiplication
  const AutoDScalar AutoDScalar::operator * (const PetscScalar v) const
  {
    AutoDScalar tmp;
    tmp.val=val*v;
    for (unsigned int _i=0; _i<numdir; ++_i)
      tmp.adval[_i]=adval[_i]*v;
    return tmp;
  }

  const AutoDScalar AutoDScalar::operator * (const AutoDScalar& a) const
  {
    AutoDScalar tmp;
    tmp.val=val*a.val;
    for (unsigned int _i=0; _i<numdir; ++_i)
      tmp.adval[_i]=adval[_i]*a.val+val*a.adval[_i];
    return tmp;
  }

  const AutoDScalar operator * (const PetscScalar v, const AutoDScalar& a)
  {
    AutoDScalar tmp;
    tmp.val=v*a.val;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=v*a.adval[_i];
    return tmp;
  }

  // division
  const AutoDScalar AutoDScalar::operator / (const PetscScalar v) const
  {
    AutoDScalar tmp;
    PetscScalar t=1.0/v;
    tmp.val=val*t;
    for (unsigned int _i=0; _i<numdir; ++_i)
      tmp.adval[_i]=adval[_i]*t;
    return tmp;
  }

  const AutoDScalar AutoDScalar::operator / (const AutoDScalar& a) const
  {
    AutoDScalar tmp;
    tmp.val=val/a.val;
    for (unsigned int _i=0; _i<numdir; ++_i)
      tmp.adval[_i]=(adval[_i]*a.val-val*a.adval[_i])/a.val/a.val;
    return tmp;
  }

  const AutoDScalar operator / (const PetscScalar v, const AutoDScalar& a)
  {
    AutoDScalar tmp;
    tmp.val=v/a.val;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=(-v*a.adval[_i])/a.val/a.val;
    return tmp;
  }

  // inc/dec
  const AutoDScalar AutoDScalar::operator ++ ()
  {
    ++val;
    return *this;
  }

  const AutoDScalar AutoDScalar::operator ++ (int)
  {
    AutoDScalar tmp;
    tmp.val=val++;
    for (unsigned int _i=0; _i<numdir; ++_i)
      tmp.adval[_i]=adval[_i];
    return tmp;
  }

  const AutoDScalar AutoDScalar::operator -- ()
  {
    --val;
    return *this;
  }

  const AutoDScalar AutoDScalar::operator -- (int)
  {
    AutoDScalar tmp;
    tmp.val=val--;
    for (unsigned int _i=0; _i<numdir; ++_i)
      tmp.adval[_i]=adval[_i];
    return tmp;
  }

  // functions
  const AutoDScalar tan(const AutoDScalar& a)
  {
    AutoDScalar tmp;
    PetscScalar tmp2;
    tmp.val=::tan(a.val);
    tmp2=::cos(a.val);
    tmp2*=tmp2;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=a.adval[_i]/tmp2;
    return tmp;
  }

  const AutoDScalar exp(const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::exp(a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=tmp.val*a.adval[_i];
    return tmp;
  }

  const AutoDScalar log(const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::log(a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      if (a.val>0 || (a.val==0 && a.adval[_i]>=0)) tmp.adval[_i]=a.adval[_i]/a.val;
      else tmp.adval[_i]=std::numeric_limits<PetscScalar>::quiet_NaN();
    return tmp;
  }

  const AutoDScalar sqrt(const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::sqrt(a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
    {
      if (a.val>0)
        tmp.adval[_i]=0.5*a.adval[_i]/tmp.val;
      else if (a.val==0 && a.adval[_i]==0)
        tmp.adval[_i]=0;
      else
        tmp.adval[_i]=std::numeric_limits<PetscScalar>::quiet_NaN();
    }
    return tmp;
  }

  const AutoDScalar sin(const AutoDScalar &a)
  {
    AutoDScalar tmp;
    PetscScalar tmp2;
    tmp.val=::sin(a.val);
    tmp2=::cos(a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=tmp2*a.adval[_i];
    return tmp;
  }

  const AutoDScalar cos(const AutoDScalar &a)
  {
    AutoDScalar tmp;
    PetscScalar tmp2;
    tmp.val=::cos(a.val);
    tmp2=-::sin(a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=tmp2*a.adval[_i];
    return tmp;
  }

  const AutoDScalar asin(const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::asin(a.val);
    PetscScalar tmp2=::sqrt(1-a.val*a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=a.adval[_i]/tmp2;
    return tmp;
  }

  const AutoDScalar acos(const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::acos(a.val);
    PetscScalar tmp2=-::sqrt(1-a.val*a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=a.adval[_i]/tmp2;
    return tmp;
  }

  const AutoDScalar atan(const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::atan(a.val);
    PetscScalar tmp2=1+a.val*a.val;
    tmp2=1/tmp2;
    if (tmp2!=0)
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=a.adval[_i]*tmp2;
    else
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=0.0;
    return tmp;
  }

  const AutoDScalar atan2(const AutoDScalar &a, const AutoDScalar &b)
  {
    AutoDScalar tmp;
    tmp.val=::atan2(a.val, b.val);
    PetscScalar tmp2=a.val*a.val;
    PetscScalar tmp3=b.val*b.val;
    PetscScalar tmp4=tmp3/(tmp2+tmp3);
    if (tmp4!=0)
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=(a.adval[_i]*b.val-a.val*b.adval[_i])/tmp3*tmp4;
    else
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=0.0;
    return tmp;
  }

  const AutoDScalar pow(const AutoDScalar &a, PetscScalar v)
  {
    AutoDScalar tmp;
    tmp.val=::pow(a.val, v);
    PetscScalar tmp2=v*::pow(a.val, v-1);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=tmp2*a.adval[_i];
    return tmp;
  }

  const AutoDScalar pow(const AutoDScalar &a, const AutoDScalar &b)
  {
    AutoDScalar tmp;
    tmp.val=::pow(a.val, b.val);
    PetscScalar tmp2=b.val*::pow(a.val, b.val-1);
    PetscScalar tmp3=::log(a.val)*tmp.val;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=tmp2*a.adval[_i]+tmp3*b.adval[_i];
    return tmp;
  }

  const AutoDScalar pow(PetscScalar v, const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::pow(v, a.val);
    PetscScalar tmp2=tmp.val*::log(v);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=tmp2*a.adval[_i];
    return tmp;
  }

  const AutoDScalar log10(const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::log10(a.val);
    PetscScalar tmp2=::log((PetscScalar)10)*a.val;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=a.adval[_i]/tmp2;
    return tmp;
  }

  const AutoDScalar sinh (const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::sinh(a.val);
    PetscScalar tmp2=::cosh(a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=a.adval[_i]*tmp2;
    return tmp;
  }

  const AutoDScalar cosh (const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::cosh(a.val);
    PetscScalar tmp2=::sinh(a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=a.adval[_i]*tmp2;
    return tmp;
  }

  const AutoDScalar tanh (const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::tanh(a.val);
    PetscScalar tmp2=::cosh(a.val);
    tmp2*=tmp2;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=a.adval[_i]/tmp2;
    return tmp;
  }


  const AutoDScalar asinh (const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=boost::math::asinh(a.val);
    PetscScalar tmp2=::sqrt(a.val*a.val+1);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=a.adval[_i]/tmp2;
    return tmp;
  }

  const AutoDScalar acosh (const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=boost::math::acosh(a.val);
    PetscScalar tmp2=::sqrt(a.val*a.val-1);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=a.adval[_i]/tmp2;
    return tmp;
  }

  const AutoDScalar atanh (const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=boost::math::atanh(a.val);
    PetscScalar tmp2=1-a.val*a.val;
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=a.adval[_i]/tmp2;
    return tmp;
  }


  const AutoDScalar fabs (const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::fabs(a.val);
    int as=0;
    if (a.val>0) as=1;
    if (a.val<0) as=-1;
    if (as!=0)
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=a.adval[_i]*as;
    else
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      {
        as=0;
        if (a.adval[_i]>0) as=1;
        if (a.adval[_i]<0) as=-1;
        tmp.adval[_i]=a.adval[_i]*as;
      }
    return tmp;
  }

  const AutoDScalar ceil (const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::ceil(a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=0.0;
    return tmp;
  }

  const AutoDScalar floor (const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::floor(a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=0.0;
    return tmp;
  }

  const AutoDScalar fmax (const AutoDScalar &a, const AutoDScalar &b)
  {
    AutoDScalar tmp;
    PetscScalar tmp2=a.val-b.val;
    if (tmp2<0)
    {
      tmp.val=b.val;
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=b.adval[_i];
    }
    else
    {
      tmp.val=a.val;
      if (tmp2>0)
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
          tmp.adval[_i]=a.adval[_i];
      }
      else
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        {
          if (a.adval[_i]<b.adval[_i]) tmp.adval[_i]=b.adval[_i];
          else tmp.adval[_i]=a.adval[_i];
        }
      }
    }
    return tmp;
  }

  const AutoDScalar fmax (PetscScalar v, const AutoDScalar &a)
  {
    AutoDScalar tmp;
    PetscScalar tmp2=v-a.val;
    if (tmp2<0)
    {
      tmp.val=a.val;
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=a.adval[_i];
    }
    else
    {
      tmp.val=v;
      if (tmp2>0)
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
          tmp.adval[_i]=0.0;
      }
      else
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        {
          if (a.adval[_i]>0) tmp.adval[_i]=a.adval[_i];
          else tmp.adval[_i]=0.0;
        }
      }
    }
    return tmp;
  }

  const AutoDScalar fmax (const AutoDScalar &a, PetscScalar v)
  {
    AutoDScalar tmp;
    PetscScalar tmp2=a.val-v;
    if (tmp2<0)
    {
      tmp.val=v;
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=0.0;
    }
    else
    {
      tmp.val=a.val;
      if (tmp2>0)
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
          tmp.adval[_i]=a.adval[_i];
      }
      else
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        {
          if (a.adval[_i]>0) tmp.adval[_i]=a.adval[_i];
          else tmp.adval[_i]=0.0;
        }
      }
    }
    return tmp;
  }

  const AutoDScalar fmin (const AutoDScalar &a, const AutoDScalar &b)
  {
    AutoDScalar tmp;
    PetscScalar tmp2=a.val-b.val;
    if (tmp2<0)
    {
      tmp.val=a.val;
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=a.adval[_i];
    }
    else
    {
      tmp.val=b.val;
      if (tmp2>0)
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
          tmp.adval[_i]=b.adval[_i];
      }
      else
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        {
          if (a.adval[_i]<b.adval[_i]) tmp.adval[_i]=a.adval[_i];
          else tmp.adval[_i]=b.adval[_i];
        }
      }
    }
    return tmp;
  }

  const AutoDScalar fmin (PetscScalar v, const AutoDScalar &a)
  {
    AutoDScalar tmp;
    PetscScalar tmp2=v-a.val;
    if (tmp2<0)
    {
      tmp.val=v;
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=0.0;
    }
    else
    {
      tmp.val=a.val;
      if (tmp2>0)
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
          tmp.adval[_i]=a.adval[_i];
      }
      else
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        {
          if (a.adval[_i]<0) tmp.adval[_i]=a.adval[_i];
          else tmp.adval[_i]=0.0;
        }
      }
    }
    return tmp;
  }

  const AutoDScalar fmin (const AutoDScalar &a, PetscScalar v)
  {
    AutoDScalar tmp;
    PetscScalar tmp2=a.val-v;
    if (tmp2<0)
    {
      tmp.val=a.val;
      for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        tmp.adval[_i]=a.adval[_i];
    }
    else
    {
      tmp.val=v;
      if (tmp2>0)
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
          tmp.adval[_i]=0.0;
      }
      else
      {
        for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
        {
          if (a.adval[_i]<0) tmp.adval[_i]=a.adval[_i];
          else tmp.adval[_i]=0.0;
        }
      }
    }
    return tmp;
  }

  const AutoDScalar ldexp (const AutoDScalar &a, const AutoDScalar &b)
  {
    return a*pow(2.,b);
  }

  const AutoDScalar ldexp (const AutoDScalar &a, const PetscScalar v)
  {
    return a*::pow(2.,v);
  }

  const AutoDScalar ldexp (const PetscScalar v, const AutoDScalar &a)
  {
    return v*pow(2.,a);
  }

  PetscScalar frexp (const AutoDScalar &a, int* v)
  {
    return ::frexp(a.val, v);
  }

#ifndef WINDOWS
  const AutoDScalar erf (const AutoDScalar &a)
  {
    AutoDScalar tmp;
    tmp.val=::erf(a.val);
    PetscScalar tmp2=2.0/::sqrt(::acos(-1.0))*::exp(-a.val*a.val);
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      tmp.adval[_i]=tmp2*a.adval[_i];
    return tmp;
  }
#endif


  /*******************  nontemporary results  *********************************/
  void AutoDScalar::operator = (const PetscScalar v)
  {
    val=v;
    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]=0.0;
  }

  void AutoDScalar::operator = (const AutoDScalar& a)
  {
    val=a.val;
    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]=a.adval[_i];
  }

  void AutoDScalar::operator += (const PetscScalar v)
  {
    val+=v;
  }

  void AutoDScalar::operator += (const AutoDScalar& a)
  {
    val=val+a.val;
    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]+=a.adval[_i];
  }

  void AutoDScalar::operator -= (const PetscScalar v)
  {
    val-=v;
  }

  void AutoDScalar::operator -= (const AutoDScalar& a)
  {
    val=val-a.val;
    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]-=a.adval[_i];
  }

  void AutoDScalar::operator *= (const PetscScalar v)
  {
    val=val*v;
    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]*=v;
  }

  void AutoDScalar::operator *= (const AutoDScalar& a)
  {
    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]=adval[_i]*a.val+val*a.adval[_i];
    val*=a.val;
  }

  void AutoDScalar::operator /= (const PetscScalar v)
  {
    val/=v;
    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]/=v;
  }

  void AutoDScalar::operator /= (const AutoDScalar& a)
  {
    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]=(adval[_i]*a.val-val*a.adval[_i])/a.val/a.val;
    val=val/a.val;
  }

  // not
  int AutoDScalar::operator ! () const
  {
    return val==0.0;
  }

  // comparision
  int AutoDScalar::operator != (const AutoDScalar &a) const
  {
    return val!=a.val;
  }

  int AutoDScalar::operator != (const PetscScalar v) const
  {
    return val!=v;
  }

  int operator != (const PetscScalar v, const AutoDScalar &a)
  {
    return v!=a.val;
  }

  int AutoDScalar::operator == (const AutoDScalar &a) const
  {
    return val==a.val;
  }

  int AutoDScalar::operator == (const PetscScalar v) const
  {
    return val==v;
  }

  int operator == (const PetscScalar v, const AutoDScalar &a)
  {
    return v==a.val;
  }

  int AutoDScalar::operator <= (const AutoDScalar &a) const
  {
    return val<=a.val;
  }

  int AutoDScalar::operator <= (const PetscScalar v) const
  {
    return val<=v;
  }

  int operator <= (const PetscScalar v, const AutoDScalar &a)
  {
    return v<=a.val;
  }

  int AutoDScalar::operator >= (const AutoDScalar &a) const
  {
    return val>=a.val;
  }

  int AutoDScalar::operator >= (const PetscScalar v) const
  {
    return val>=v;
  }

  int operator >= (const PetscScalar v, const AutoDScalar &a)
  {
    return v>=a.val;
  }

  int AutoDScalar::operator >  (const AutoDScalar &a) const
  {
    return val>a.val;
  }

  int AutoDScalar::operator >  (const PetscScalar v) const
  {
    return val>v;
  }

  int operator >  (const PetscScalar v, const AutoDScalar &a)
  {
    return v>a.val;
  }

  int AutoDScalar::operator <  (const AutoDScalar &a) const
  {
    return val<a.val;
  }

  int AutoDScalar::operator <  (const PetscScalar v) const
  {
    return val<v;
  }

  int operator <  (const PetscScalar v, const AutoDScalar &a)
  {
    return v<a.val;
  }

  /*******************  getter / setter  **************************************/
  PetscScalar AutoDScalar::getValue() const
  {
    return val;
  }

  void AutoDScalar::setValue(const PetscScalar v)
  {
    val=v;
  }

  std::vector<PetscScalar> AutoDScalar::getADValueVector() const
  {
    std::vector<PetscScalar> ad_vec(numdir);
    for (unsigned int _i=0; _i<numdir; ++_i)
      ad_vec[_i]=adval[_i];
    return ad_vec;
  }

  const PetscScalar * AutoDScalar::getADValue() const
  {
    return adval;
  }

  void AutoDScalar::setADValue(const PetscScalar * v)
  {
    for (unsigned int _i=0; _i<numdir; ++_i)
      adval[_i]=v[_i];
  }

#  if defined(ADTL_NUMBER_DIRECTIONS)
  PetscScalar AutoDScalar::getADValue(const unsigned int p) const
  {
    return adval[p];
  }

  void AutoDScalar::setADValue(const unsigned int p, const PetscScalar v)
  {
    adval[p]=v;
  }
#  endif

  /*******************  i/o operations  ***************************************/
  std::ostream& operator << ( std::ostream& out, const AutoDScalar& a)
  {
    out << "Value: " << a.val;
#if !defined(ADTL_NUMBER_DIRECTIONS)
    out << " ADValue: ";
#else
    out << " ADValues (" << AutoDScalar::numdir << "): ";
#endif
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      out << a.adval[_i] << " ";
    out << "(a)";
    return out;
  }

  std::istream& operator >> ( std::istream& in, AutoDScalar& a)
  {
    char c;
    do
    {
      in >> c;
    }
    while (c!=':' && !in.eof());
    in >> a.val;
#if !defined(ADTL_NUMBER_DIRECTIONS)
    do in >> c;
    while (c!=':' && !in.eof());
#else
    unsigned int num;
    do in >> c;
    while (c!='(' && !in.eof());
    in >> num;
    if (num>ADTL_NUMBER_DIRECTIONS)
    {
      std::cout << "ADOL-C error: to many directions in input\n";
      exit(-1);
    }
    do in >> c;
    while (c!=')' && !in.eof());
#endif
    for (unsigned int _i=0; _i<AutoDScalar::numdir; ++_i)
      in >> a.adval[_i];
    do in >> c;
    while (c!=')' && !in.eof());
    return in;
  }

}

#endif
