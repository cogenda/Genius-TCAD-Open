#ifndef __adsmt_h__
#define __adsmt_h__

#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <limits>

#include "asinh.hpp" // for asinh
#include "acosh.hpp" // for acosh
#include "atanh.hpp" // for atanh

#define ADSMT_NUMBER_DIRECTIONS 32



typedef double PetscScalar;
typedef int    PetscInt;

#if 0
/**
 *
 */
class ADPattern
{
public:
  /// build a new AD pattern, return the pattern id
  static size_t build_pattern( const std::vector<PetscInt> & pattern);

  /// the index of pattern p
  static PetscInt * pattern_index(const size_t p);

  /// the size of pattern p
  static PetscInt   pattern_size ( const size_t p);

  /// for two different pattern, pad them to the same size and return the id of new pattern
  static size_t merge_patterns(const size_t p1, const size_t p2);

  /// clear all the patterns, free memory
  static void   clear_patterns();

  /// number of patterns
  static size_t num_patterns();

private:

  //static std:vector< std::vector<PetscInt> > ad_patterns;

};
#endif


class AsmtDScalar
{
public:

  /*******************  ctors / dtors  ********************************/
  inline AsmtDScalar(PetscScalar val=0.0);
  //inline AsmtDScalar(PetscScalar val, const std::vector<PetscInt> & ad_index, const PetscInt ad_index, const PetscScalar ad_value);
  inline AsmtDScalar(PetscScalar val, const size_t size, const PetscInt * ad_index, const PetscScalar * ad_value);
  inline AsmtDScalar(PetscScalar val, const std::vector<PetscInt> & ad_index, const std::vector<PetscScalar> & ad_value);
  inline AsmtDScalar(const AsmtDScalar& a);
  inline ~AsmtDScalar();

  /*******************  getter / setter  ********************************/
  inline PetscScalar getValue() const { return _val; }
  inline void setValue(PetscScalar v) { _val = v; }

  inline PetscScalar getADValue(const PetscInt p) const;
  inline void setADValue(const PetscInt p, const PetscScalar v);
  inline void setADValue(const std::vector<PetscInt> & ad_index, const PetscInt p, const PetscScalar v);
  inline PetscInt getADSize() const { return (PetscInt)_occupied; }
  inline const PetscInt * getADIndex() const { return _deriv_index; }
  inline const PetscScalar * getADValue() const { return _deriv_value; }

  /*******************  temporary results  ******************************/
  // sign
  inline const AsmtDScalar operator - () const;
  inline const AsmtDScalar operator + () const;

  // addition
  inline const AsmtDScalar operator + (const PetscScalar v) const;
  inline const AsmtDScalar operator + (const AsmtDScalar& a) const;
  inline friend const AsmtDScalar operator + (const PetscScalar v, const AsmtDScalar& a);

  // substraction
  inline const AsmtDScalar operator - (const PetscScalar v) const;
  inline const AsmtDScalar operator - (const AsmtDScalar& a) const;
  inline friend const AsmtDScalar operator - (const PetscScalar v, const AsmtDScalar& a);

  // multiplication
  inline const AsmtDScalar operator * (const PetscScalar v) const;
  inline const AsmtDScalar operator * (const AsmtDScalar& a) const;
  inline friend const AsmtDScalar operator * (const PetscScalar v, const AsmtDScalar& a);

  // division
  inline const AsmtDScalar operator / (const PetscScalar v) const;
  inline const AsmtDScalar operator / (const AsmtDScalar& a) const;
  inline friend const AsmtDScalar operator / (const PetscScalar v, const AsmtDScalar& a);

  // functions
  inline friend const AsmtDScalar exp(const AsmtDScalar &a);
  inline friend const AsmtDScalar log(const AsmtDScalar &a);
  inline friend const AsmtDScalar sqrt(const AsmtDScalar &a);

  inline friend const AsmtDScalar sin(const AsmtDScalar &a);
  inline friend const AsmtDScalar cos(const AsmtDScalar &a);
  inline friend const AsmtDScalar tan(const AsmtDScalar &a);

  inline friend const AsmtDScalar asin(const AsmtDScalar &a);
  inline friend const AsmtDScalar acos(const AsmtDScalar &a);
  inline friend const AsmtDScalar atan(const AsmtDScalar &a);

  inline friend const AsmtDScalar pow(const AsmtDScalar &a, PetscScalar v);
  inline friend const AsmtDScalar pow(const AsmtDScalar &a, const AsmtDScalar &b);
  inline friend const AsmtDScalar pow(PetscScalar v, const AsmtDScalar &a);
  inline friend const AsmtDScalar log10(const AsmtDScalar &a);

  inline friend const AsmtDScalar sinh (const AsmtDScalar &a);
  inline friend const AsmtDScalar cosh (const AsmtDScalar &a);
  inline friend const AsmtDScalar tanh (const AsmtDScalar &a);

  inline friend const AsmtDScalar asinh (const AsmtDScalar &a);
  inline friend const AsmtDScalar acosh (const AsmtDScalar &a);
  inline friend const AsmtDScalar atanh (const AsmtDScalar &a);

  inline friend const AsmtDScalar fabs (const AsmtDScalar &a);
  inline friend const AsmtDScalar ceil (const AsmtDScalar &a);
  inline friend const AsmtDScalar floor (const AsmtDScalar &a);
  inline friend const AsmtDScalar fmax (const AsmtDScalar &a, const AsmtDScalar &b);
  inline friend const AsmtDScalar fmax (PetscScalar v, const AsmtDScalar &a);
  inline friend const AsmtDScalar fmax (const AsmtDScalar &a, PetscScalar v);
  inline friend const AsmtDScalar fmin (const AsmtDScalar &a, const AsmtDScalar &b);
  inline friend const AsmtDScalar fmin (PetscScalar v, const AsmtDScalar &a);
  inline friend const AsmtDScalar fmin (const AsmtDScalar &a, PetscScalar v);

  inline friend const AsmtDScalar fermi_half (const AsmtDScalar &a);

  /*******************  nontemporary results  ***************************/
  // assignment
  inline void operator = (const PetscScalar v);
  inline void operator = (const AsmtDScalar& a);

  // addition
  inline void operator += (const PetscScalar v);
  inline void operator += (const AsmtDScalar& a);

  // substraction
  inline void operator -= (const PetscScalar v);
  inline void operator -= (const AsmtDScalar& a);

  // multiplication
  inline void operator *= (const PetscScalar v);
  inline void operator *= (const AsmtDScalar& a);

  // division
  inline void operator /= (const PetscScalar v);
  inline void operator /= (const AsmtDScalar& a);

  // not
  inline bool operator ! () const;

  // comparision
  inline bool operator != (const AsmtDScalar&) const;
  inline bool operator != (const PetscScalar) const;
  inline friend bool operator != (const PetscScalar, const AsmtDScalar&);

  inline bool operator == (const AsmtDScalar&) const;
  inline bool operator == (const PetscScalar) const;
  inline friend bool operator == (const PetscScalar, const AsmtDScalar&);

  inline bool operator <= (const AsmtDScalar&) const;
  inline bool operator <= (const PetscScalar) const;
  inline friend bool operator <= (const PetscScalar, const AsmtDScalar&);

  inline bool operator >= (const AsmtDScalar&) const;
  inline bool operator >= (const PetscScalar) const;
  inline friend bool operator >= (const PetscScalar, const AsmtDScalar&);

  inline bool operator >  (const AsmtDScalar&) const;
  inline bool operator >  (const PetscScalar) const;
  inline friend bool operator >  (const PetscScalar, const AsmtDScalar&);

  inline bool operator <  (const AsmtDScalar&) const;
  inline bool operator <  (const PetscScalar) const;
  inline friend bool operator <  (const PetscScalar, const AsmtDScalar&);

  /*******************  i/o operations  *********************************/
  inline friend std::ostream& operator << ( std::ostream&, const AsmtDScalar& );

private:

  /// value
  PetscScalar _val;

  // when AD independent variable less than 4x8=32, use fixed buffer for cache friendly access
  // or we have to use dynamic array
  size_t _occupied;
  PetscInt _deriv_index_buf[ADSMT_NUMBER_DIRECTIONS];
  PetscScalar _deriv_value_buf[ADSMT_NUMBER_DIRECTIONS];

  // when _occupied < ADSMT_NUMBER_DIRECTIONS, point to buf, or point to dynamic allocate memory
  PetscInt *_deriv_index;
  PetscScalar * _deriv_value;

  //some aux functions
  inline PetscScalar * _findADValue (const PetscInt p) const;
  inline bool          _withSamePattern(const AsmtDScalar & a) const;
  inline size_t        _patternSize(const AsmtDScalar & a) const;
  inline void          _padPattern(const AsmtDScalar & a, AsmtDScalar & pa);
};

//---------------------------------------------------------------------------

AsmtDScalar::AsmtDScalar(PetscScalar val)
    :_val(val), _occupied(0)
{
  _deriv_index = (PetscInt *)_deriv_index_buf;
  _deriv_value = (PetscScalar *)_deriv_value_buf;
}


AsmtDScalar::AsmtDScalar(PetscScalar val, const size_t size, const PetscInt * ad_index, const PetscScalar * ad_value)
    :_val(val), _occupied(size)
{
  if( _occupied <= ADSMT_NUMBER_DIRECTIONS )
  {
    _deriv_index = (PetscInt *)_deriv_index_buf;
    _deriv_value = (PetscScalar *)_deriv_value_buf;
  }
  else
  {
    _deriv_index = new PetscInt[_occupied];
    _deriv_value = new PetscScalar[_occupied];
  }
  memcpy(_deriv_index, (void*)ad_index, _occupied*sizeof(PetscInt));
  memcpy(_deriv_value, (void*)ad_value, _occupied*sizeof(PetscScalar));

}


AsmtDScalar::AsmtDScalar(PetscScalar val, const std::vector<PetscInt> & ad_index, const std::vector<PetscScalar> & ad_value)
    :_val(val)
{
  _occupied = ad_index.size();
  if( _occupied <= ADSMT_NUMBER_DIRECTIONS )
  {
    _deriv_index = (PetscInt *)_deriv_index_buf;
    _deriv_value = (PetscScalar *)_deriv_value_buf;
  }
  else
  {
    _deriv_index = new PetscInt[_occupied];
    _deriv_value = new PetscScalar[_occupied];
  }
  memcpy(_deriv_index, (void*)&ad_index[0], _occupied*sizeof(PetscInt));
  memcpy(_deriv_value, (void*)&ad_value[0], _occupied*sizeof(PetscScalar));
}


AsmtDScalar::AsmtDScalar(const AsmtDScalar& a)
    :_val(a._val), _occupied(a._occupied)
{
  if( _occupied <= ADSMT_NUMBER_DIRECTIONS )
  {
    _deriv_index = (PetscInt *)_deriv_index_buf;
    _deriv_value = (PetscScalar *)_deriv_value_buf;
  }
  else
  {
    _deriv_index = new PetscInt[_occupied];
    _deriv_value = new PetscScalar[_occupied];
  }
  memcpy(_deriv_index, (void*)a._deriv_index, _occupied*sizeof(PetscInt));
  memcpy(_deriv_value, (void*)a._deriv_value, _occupied*sizeof(PetscScalar));
}


AsmtDScalar::~AsmtDScalar()
{
  if(_occupied > ADSMT_NUMBER_DIRECTIONS)
  {
    delete [] _deriv_index;
    delete [] _deriv_value;
  }
}


/*************************  temporary results  ******************************/
// sign
const AsmtDScalar AsmtDScalar::operator - () const
{
  AsmtDScalar tmp(-_val, _occupied, _deriv_index, _deriv_value);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i] = -tmp._deriv_value[_i];
  return tmp;
}

const AsmtDScalar AsmtDScalar::operator + () const
{
  return *this;
}

// addition
const AsmtDScalar AsmtDScalar::operator + (const PetscScalar v) const
{
  return AsmtDScalar(_val+v, _occupied, _deriv_index, _deriv_value);
}

const AsmtDScalar AsmtDScalar::operator + (const AsmtDScalar& a) const
{
  AsmtDScalar tmp(*this);
  tmp._val += a._val;
  if(_withSamePattern(a))
  {
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i] += a._deriv_value[_i];
  }
  else
  {
    AsmtDScalar pa;
    tmp._padPattern(a, pa);
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i] += pa._deriv_value[_i];
  }
  return tmp;
}

const AsmtDScalar operator + (const PetscScalar v, const AsmtDScalar& a)
{
  return AsmtDScalar(v+a._val, a._occupied, a._deriv_index, a._deriv_value);
}


// subtraction
const AsmtDScalar AsmtDScalar::operator - (const PetscScalar v) const
{
  return AsmtDScalar(_val-v, _occupied, _deriv_index, _deriv_value);
}

const AsmtDScalar AsmtDScalar::operator - (const AsmtDScalar& a) const
{
  AsmtDScalar tmp(*this);
  tmp._val -= a._val;
  if(_withSamePattern(a))
  {
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i] -= a._deriv_value[_i];
  }
  else
  {
    AsmtDScalar pa;
    tmp._padPattern(a, pa);
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i] -= pa._deriv_value[_i];
  }
  return tmp;
}

const AsmtDScalar operator - (const PetscScalar v, const AsmtDScalar& a)
{
  AsmtDScalar tmp(-a._val+v, a._occupied, a._deriv_index, a._deriv_value);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i] = -tmp._deriv_value[_i];
  return tmp;
}


// multiplication
const AsmtDScalar AsmtDScalar::operator * (const PetscScalar v) const
{
  AsmtDScalar tmp(_val*v, _occupied, _deriv_index, _deriv_value);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i] *= v;
  return tmp;
}


const AsmtDScalar AsmtDScalar::operator * (const AsmtDScalar& a) const
{
  AsmtDScalar tmp(*this);
  tmp._val *= a._val;
  if(_withSamePattern(a))
  {
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i] = _deriv_value[_i]*a._val + _val*a._deriv_value[_i];
  }
  else
  {
    AsmtDScalar pa;
    tmp._padPattern(a, pa);
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i] = tmp._deriv_value[_i]*pa._val + _val*pa._deriv_value[_i];
  }
  return tmp;
}

const AsmtDScalar operator * (const PetscScalar v, const AsmtDScalar& a)
{
  AsmtDScalar tmp(a._val*v, a._occupied, a._deriv_index, a._deriv_value);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i] *= v;
  return tmp;
}

// division
const AsmtDScalar AsmtDScalar::operator / (const PetscScalar v) const
{
  AsmtDScalar tmp(_val/v, _occupied, _deriv_index, _deriv_value);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i] /= v;
  return tmp;
}

const AsmtDScalar AsmtDScalar::operator / (const AsmtDScalar& a) const
{
  AsmtDScalar tmp(*this);
  tmp._val /= a._val;
  if(_withSamePattern(a))
  {
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i] = (_deriv_value[_i]*a._val-_val*a._deriv_value[_i])/(a._val*a._val);
  }
  else
  {
    AsmtDScalar pa;
    tmp._padPattern(a, pa);
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i] = (tmp._deriv_value[_i]*pa._val-_val*pa._deriv_value[_i])/(pa._val*pa._val);
  }
  return tmp;
}

const AsmtDScalar operator / (const PetscScalar v, const AsmtDScalar& a)
{
  AsmtDScalar tmp(v/a._val, a._occupied, a._deriv_index, a._deriv_value);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=(-v*a._deriv_value[_i])/(a._val*a._val);
  return tmp;
}

// functions
const AsmtDScalar exp(const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::exp(a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=tmp._val*a._deriv_value[_i];
  return tmp;
}

const AsmtDScalar log(const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::log(a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    if (a._val>0 || (a._val==0 && a._deriv_value[_i]>=0))
      tmp._deriv_value[_i] = a._deriv_value[_i]/a._val;
    else
      tmp._deriv_value[_i] = std::numeric_limits<PetscScalar>::quiet_NaN();
  return tmp;
}

const AsmtDScalar sqrt(const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::sqrt(a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
  {
    if (a._val>0)
      tmp._deriv_value[_i]=0.5*a._deriv_value[_i]/tmp._val;
    else if (a._val==0 && a._deriv_value[_i]==0)
      tmp._deriv_value[_i]=0;
    else
      tmp._deriv_value[_i]=std::numeric_limits<PetscScalar>::quiet_NaN();
  }
  return tmp;
}

const AsmtDScalar sin(const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  PetscScalar tmp2;
  tmp._val=std::sin(a._val);
  tmp2=std::cos(a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i] = tmp2*a._deriv_value[_i];
  return tmp;
}

const AsmtDScalar cos(const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  PetscScalar tmp2;
  tmp._val=std::cos(a._val);
  tmp2=-std::sin(a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=tmp2*a._deriv_value[_i];
  return tmp;
}

const AsmtDScalar tan(const AsmtDScalar& a)
{
  AsmtDScalar tmp(a);
  PetscScalar tmp2;
  tmp._val=std::tan(a._val);
  tmp2=cos(a._val);
  tmp2*=tmp2;
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=a._deriv_value[_i]/tmp2;
  return tmp;
}

const AsmtDScalar asin(const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::asin(a._val);
  PetscScalar tmp2=std::sqrt(1-a._val*a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=a._deriv_value[_i]/tmp2;
  return tmp;
}

const AsmtDScalar acos(const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::acos(a._val);
  PetscScalar tmp2=-std::sqrt(1-a._val*a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=a._deriv_value[_i]/tmp2;
  return tmp;
}

const AsmtDScalar atan(const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::atan(a._val);
  PetscScalar tmp2=1+a._val*a._val;
  tmp2=1/tmp2;
  if (tmp2!=0)
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i]=a._deriv_value[_i]*tmp2;
  else
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i]=0.0;
  return tmp;
}

const AsmtDScalar pow(const AsmtDScalar &a, PetscScalar v)
{
  AsmtDScalar tmp(a);
  tmp._val=std::pow(a._val, v);
  PetscScalar tmp2=v*std::pow(a._val, v-1);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=tmp2*a._deriv_value[_i];
  return tmp;
}

const AsmtDScalar pow(const AsmtDScalar &a, const AsmtDScalar &b)
{
  if(a._withSamePattern(b))
  {
    AsmtDScalar tmp(a);
    tmp._val=std::pow(a._val, b._val);
    PetscScalar tmp2=b._val*::pow(a._val, b._val-1);
    PetscScalar tmp3=std::log(a._val)*tmp._val;
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i]=tmp2*a._deriv_value[_i]+tmp3*b._deriv_value[_i];
    return tmp;
  }
  else
  {
    AsmtDScalar tmp(a);
    AsmtDScalar pb;
    tmp._padPattern(b, pb);
    tmp._val=std::pow(a._val, b._val);
    PetscScalar tmp2=b._val*::pow(a._val, b._val-1);
    PetscScalar tmp3=std::log(a._val)*tmp._val;
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i]=tmp2*tmp._deriv_value[_i]+tmp3*pb._deriv_value[_i];
    return tmp;
  }
}

const AsmtDScalar pow(PetscScalar v, const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::pow(v, a._val);
  PetscScalar tmp2=tmp._val*::log(v);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=tmp2*a._deriv_value[_i];
  return tmp;
}

const AsmtDScalar log10(const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::log10(a._val);
  PetscScalar tmp2=std::log((PetscScalar)10)*a._val;
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=a._deriv_value[_i]/tmp2;
  return tmp;
}


const AsmtDScalar sinh (const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::sinh(a._val);
  PetscScalar tmp2=std::cosh(a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=a._deriv_value[_i]*tmp2;
  return tmp;
}

const AsmtDScalar cosh (const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::cosh(a._val);
  PetscScalar tmp2=std::sinh(a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=a._deriv_value[_i]*tmp2;
  return tmp;
}

const AsmtDScalar tanh (const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::tanh(a._val);
  PetscScalar tmp2=std::cosh(a._val);
  tmp2*=tmp2;
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=a._deriv_value[_i]/tmp2;
  return tmp;
}

const AsmtDScalar asinh (const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=boost::math::asinh(a._val);
  PetscScalar tmp2=std::sqrt(a._val*a._val+1);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=a._deriv_value[_i]/tmp2;
  return tmp;
}

const AsmtDScalar acosh (const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=boost::math::acosh(a._val);
  PetscScalar tmp2=std::sqrt(a._val*a._val-1);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=a._deriv_value[_i]/tmp2;
  return tmp;
}

const AsmtDScalar atanh (const AsmtDScalar &a)
{
  AsmtDScalar tmp;
  tmp._val=boost::math::atanh(a._val);
  PetscScalar tmp2=1-a._val*a._val;
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=a._deriv_value[_i]/tmp2;
  return tmp;
}


const AsmtDScalar fabs (const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=std::abs(a._val);
  int as=0;
  if (a._val>0) as=1;
  if (a._val<0) as=-1;
  if (as!=0)
    for (size_t _i=0; _i<tmp._occupied; ++_i)
      tmp._deriv_value[_i]=a._deriv_value[_i]*as;
  else
    for (size_t _i=0; _i<tmp._occupied; ++_i)
    {
      as=0;
      if (a._deriv_value[_i]>0) as=1;
      if (a._deriv_value[_i]<0) as=-1;
      tmp._deriv_value[_i]=a._deriv_value[_i]*as;
    }
  return tmp;
}

const AsmtDScalar ceil (const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=ceil(a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=0.0;
  return tmp;
}

const AsmtDScalar floor (const AsmtDScalar &a)
{
  AsmtDScalar tmp(a);
  tmp._val=floor(a._val);
  for (size_t _i=0; _i<tmp._occupied; ++_i)
    tmp._deriv_value[_i]=0.0;
  return tmp;
}


const AsmtDScalar fmax (const AsmtDScalar &a, const AsmtDScalar &b)
{
  PetscScalar tmp2=a._val-b._val;
  if (tmp2<0)
    return b;

  if (tmp2>0)
    return a;

  AsmtDScalar tmp(a);
  if( tmp._withSamePattern(b))
  {
    for (size_t _i=0; _i<tmp._occupied; ++_i)
    {
      if (a._deriv_value[_i]<b._deriv_value[_i])
        tmp._deriv_value[_i]=b._deriv_value[_i];
      //else
      //  tmp._deriv_value[_i]=a._deriv_value[_i];
    }
  }
  else
  {
    AsmtDScalar pb;
    tmp._padPattern(b, pb);
    for (size_t _i=0; _i<tmp._occupied; ++_i)
    {
      if (tmp._deriv_value[_i]<pb._deriv_value[_i])
        tmp._deriv_value[_i]=pb._deriv_value[_i];
      //else
      //  tmp._deriv_value[_i]=tmp._deriv_value[_i];
    }
  }
  return tmp;
}

const AsmtDScalar fmax (PetscScalar v, const AsmtDScalar &a)
{
  PetscScalar tmp2=v-a._val;
  if (tmp2<0)
  {
    return a;
  }
  else
  {
    AsmtDScalar tmp(a);
    tmp._val=v;
    if (tmp2>0)
    {
      for (size_t _i=0; _i<tmp._occupied; ++_i)
        tmp._deriv_value[_i]=0.0;
    }
    else
    {
      for (size_t _i=0; _i<tmp._occupied; ++_i)
      {
        if (tmp._deriv_value[_i]<0)
         tmp._deriv_value[_i]=0.0;
      }
    }
    return tmp;
  }
}

const AsmtDScalar fmax (const AsmtDScalar &a, PetscScalar v)
{
  PetscScalar tmp2=a._val-v;
  if (tmp2>0)
  {
    return a;
  }
  else
  {
    AsmtDScalar tmp(a);
    if (tmp2<0)
    {
      tmp._val=v;
      for (size_t _i=0; _i<tmp._occupied; ++_i)
        tmp._deriv_value[_i]=0.0;
    }
    else
    {
      for (size_t _i=0; _i<tmp._occupied; ++_i)
      {
        if (tmp._deriv_value[_i]<0)
          tmp._deriv_value[_i]=0.0;
      }
    }
    return tmp;
  }
}

const AsmtDScalar fmin (const AsmtDScalar &a, const AsmtDScalar &b)
{
  PetscScalar tmp2=a._val-b._val;
  if (tmp2<0)
    return a;
  if (tmp2>0)
    return b;

  AsmtDScalar tmp(a);
  if( tmp._withSamePattern(b))
  {
    for (size_t _i=0; _i<tmp._occupied; ++_i)
    {
      if (a._deriv_value[_i]>b._deriv_value[_i])
        tmp._deriv_value[_i]=b._deriv_value[_i];
      //else
      //  tmp._deriv_value[_i]=a._deriv_value[_i];
    }
  }
  else
  {
    AsmtDScalar pb;
    tmp._padPattern(b, pb);
    for (size_t _i=0; _i<tmp._occupied; ++_i)
    {
      if (tmp._deriv_value[_i]>pb._deriv_value[_i])
        tmp._deriv_value[_i]=pb._deriv_value[_i];
      //else
      //  tmp._deriv_value[_i]=tmp._deriv_value[_i];
    }
  }
  return tmp;
}

const AsmtDScalar fmin (PetscScalar v, const AsmtDScalar &a)
{
  PetscScalar tmp2=v-a._val;
  if (tmp2>0)
  {
    return a;
  }
  else
  {
    AsmtDScalar tmp(a);
    tmp._val=v;
    if (tmp2<0)
    {
      for (size_t _i=0; _i<tmp._occupied; ++_i)
        tmp._deriv_value[_i]=0.0;
    }
    else
    {
      for (size_t _i=0; _i<tmp._occupied; ++_i)
      {
        if (tmp._deriv_value[_i]>0)
          tmp._deriv_value[_i]=0.0;
      }
    }
    return tmp;
  }
}

const AsmtDScalar fmin (const AsmtDScalar &a, PetscScalar v)
{
  PetscScalar tmp2=a._val-v;
  if (tmp2<0)
  {
    return a;
  }
  else
  {
    AsmtDScalar tmp(a);
    if (tmp2>0)
    {
      tmp._val=v;
      for (size_t _i=0; _i<tmp._occupied; ++_i)
        tmp._deriv_value[_i]=0.0;
    }
    else
    {
      for (size_t _i=0; _i<tmp._occupied; ++_i)
      {
        if (tmp._deriv_value[_i]>0)
          tmp._deriv_value[_i]=0.0;
      }
    }
    return tmp;
  }
}

const AsmtDScalar fermi_half (const AsmtDScalar &a)
{
  /* use an analytic expression. The result achieves within 0.4% error in all ranges.*/
  if(a<-4.5)
  {
    //for small arguments, fhfp and exp are almost identical.
    return 1.0/exp(-a);
  }
  else if(a<0.0)
  {
    AsmtDScalar v = pow(a,4) + 50 + 33.6*a*(1-0.68*exp(-0.17*(a+1)*(a+1)));
    AsmtDScalar p = 1.329340388179*pow(v,PetscScalar(-0.375));
    return 1.0/(exp(-a) + p);
  }
  else
  {
    AsmtDScalar v = pow(a,4) + 50 + 33.6*a*(1-0.68*exp(-0.17*(a+1)*(a+1)));
    AsmtDScalar p = 1.329340388179*pow(v,PetscScalar(-0.375));
    return 1.0/(1.0/exp(a) + p);
  }
}


/*******************  nontemporary results  *********************************/
void AsmtDScalar::operator = (const PetscScalar v)
{
  _val = v;
  if(_occupied > ADSMT_NUMBER_DIRECTIONS)
  {
    delete [] _deriv_index;
    delete [] _deriv_value;
  }
  _occupied = 0;
  _deriv_index = (PetscInt *)_deriv_index_buf;
  _deriv_value = (PetscScalar *)_deriv_value_buf;
}

void AsmtDScalar::operator = (const AsmtDScalar& a)
{
  _val=a._val;

  if(_occupied > ADSMT_NUMBER_DIRECTIONS)
  {
    delete [] _deriv_index;
    delete [] _deriv_value;
  }
  _occupied = a._occupied;
  if( _occupied <= ADSMT_NUMBER_DIRECTIONS )
  {
    _deriv_index = (PetscInt *)_deriv_index_buf;
    _deriv_value = (PetscScalar *)_deriv_value_buf;
  }
  else
  {
    _deriv_index = new PetscInt[_occupied];
    _deriv_value = new PetscScalar[_occupied];
  }
  memcpy(_deriv_index, (void*)a._deriv_index, _occupied*sizeof(PetscInt));
  memcpy(_deriv_value, (void*)a._deriv_value, _occupied*sizeof(PetscScalar));
}

void AsmtDScalar::operator += (const PetscScalar v)
{
  _val += v;
}

void AsmtDScalar::operator += (const AsmtDScalar& a)
{
  _val += a._val;
  if(_withSamePattern(a))
  {
    for (size_t _i=0; _i<_occupied; ++_i)
      _deriv_value[_i] += a._deriv_value[_i];
  }
  else
  {
    AsmtDScalar pa;
    _padPattern(a, pa);
    for (size_t _i=0; _i<_occupied; ++_i)
      _deriv_value[_i] += pa._deriv_value[_i];
  }
}

void AsmtDScalar::operator -= (const PetscScalar v)
{
  _val -= v;
}

void AsmtDScalar::operator -= (const AsmtDScalar& a)
{
  _val -= a._val;
  if(_withSamePattern(a))
  {
    for (size_t _i=0; _i<_occupied; ++_i)
      _deriv_value[_i] -= a._deriv_value[_i];
  }
  else
  {
    AsmtDScalar pa;
    _padPattern(a, pa);
    for (size_t _i=0; _i<_occupied; ++_i)
      _deriv_value[_i] -= pa._deriv_value[_i];
  }
}

void AsmtDScalar::operator *= (const PetscScalar v)
{
  _val *= v;
  for (size_t _i=0; _i<_occupied; ++_i)
    _deriv_value[_i] *= v;
}

void AsmtDScalar::operator *= (const AsmtDScalar& a)
{
  if(_withSamePattern(a))
  {
    for (size_t _i=0; _i<_occupied; ++_i)
      _deriv_value[_i] = _deriv_value[_i]*a._val+_val*a._deriv_value[_i];
  }
  else
  {
    AsmtDScalar pa;
    _padPattern(a, pa);
    for (size_t _i=0; _i<_occupied; ++_i)
      _deriv_value[_i] = _deriv_value[_i]*pa._val+_val*pa._deriv_value[_i];
  }
  _val *= a._val;
}

void AsmtDScalar::operator /= (const PetscScalar v)
{
  _val /= v;
  for (size_t _i=0; _i<_occupied; ++_i)
    _deriv_value[_i] /= v;
}

void AsmtDScalar::operator /= (const AsmtDScalar& a)
{
  if(_withSamePattern(a))
  {
    for (size_t _i=0; _i<_occupied; ++_i)
      _deriv_value[_i] = (_deriv_value[_i]*a._val-_val*a._deriv_value[_i])/(a._val*a._val);
  }
  else
  {
    AsmtDScalar pa;
    _padPattern(a, pa);
    for (size_t _i=0; _i<_occupied; ++_i)
      _deriv_value[_i] = (_deriv_value[_i]*pa._val-_val*pa._deriv_value[_i])/(pa._val*pa._val);
  }
  _val /= a._val;
}

// not
bool AsmtDScalar::operator ! () const
{
  return _val == 0.0;
}

// comparision
bool AsmtDScalar::operator != (const AsmtDScalar &a) const
{
  return _val != a._val;
}

bool AsmtDScalar::operator != (const PetscScalar v) const
{
  return _val != v;
}

bool operator != (const PetscScalar v, const AsmtDScalar &a)
{
  return v != a._val;
}

bool AsmtDScalar::operator == (const AsmtDScalar &a) const
{
  return _val == a._val;
}

bool AsmtDScalar::operator == (const PetscScalar v) const
{
  return _val == v;
}

bool operator == (const PetscScalar v, const AsmtDScalar &a)
{
  return  v == a._val;
}

bool AsmtDScalar::operator <= (const AsmtDScalar &a) const
{
  return _val <= a._val;
}

bool AsmtDScalar::operator <= (const PetscScalar v) const
{
  return _val <= v;
}

bool operator <= (const PetscScalar v, const AsmtDScalar &a)
{
  return v <= a._val;
}

bool AsmtDScalar::operator >= (const AsmtDScalar &a) const
{
  return _val >= a._val;
}

bool AsmtDScalar::operator >= (const PetscScalar v) const
{
  return _val>=v;
}

bool operator >= (const PetscScalar v, const AsmtDScalar &a)
{
  return v>=a._val;
}

bool AsmtDScalar::operator >  (const AsmtDScalar &a) const
{
  return _val>a._val;
}

bool AsmtDScalar::operator >  (const PetscScalar v) const
{
  return _val>v;
}

bool operator >  (const PetscScalar v, const AsmtDScalar &a)
{
  return v>a._val;
}

bool AsmtDScalar::operator <  (const AsmtDScalar &a) const
{
  return _val<a._val;
}

bool AsmtDScalar::operator <  (const PetscScalar v) const
{
  return _val<v;
}

bool operator <  (const PetscScalar v, const AsmtDScalar &a)
{
  return v<a._val;
}

//---------------------------------------------------------------------------

PetscScalar AsmtDScalar::getADValue(const PetscInt p) const
  { return *_findADValue(p); }


void AsmtDScalar::setADValue(const PetscInt p, const PetscScalar v)
{
  assert(_occupied==0);
  _deriv_index[0] = p;
  _deriv_value[0] = v;
  ++_occupied;
}


void AsmtDScalar::setADValue(const std::vector<PetscInt> & ad_index, const PetscInt p, const PetscScalar v)
{
  if( _occupied  > ADSMT_NUMBER_DIRECTIONS)
  {
    delete [] _deriv_index;
    delete [] _deriv_value;
  }

  _occupied=ad_index.size();
  if( _occupied <= ADSMT_NUMBER_DIRECTIONS )
  {
    _deriv_index = (PetscInt *)_deriv_index_buf;
    _deriv_value = (PetscScalar *)_deriv_value_buf;
  }
  else
  {
    _deriv_index = new PetscInt[_occupied];
    _deriv_value = new PetscScalar[_occupied];
  }

  memcpy(_deriv_index, (void*)&ad_index[0], _occupied*sizeof(PetscInt));
  for (size_t _i=0; _i<_occupied; ++_i)
  {
    if(_deriv_index[_i] == p)
      _deriv_value[_i] = v;
    else
      _deriv_value[_i] = 0.0;
  }
}



/*******************  i/o operations  ***************************************/
std::ostream& operator << ( std::ostream& out, const AsmtDScalar& a)
{
  out << "Value: " << a._val;
  out << " ADValues (" << a._occupied << "): ";
  for (size_t _i=0; _i<a._occupied; ++_i)
    out << '(' << a._deriv_index[_i] << ',' << a._deriv_value[_i] << ')' << ' ';
  out << "(a)";
  return out;
}


/*******************  aux functions  ***************************************/

PetscScalar * AsmtDScalar::_findADValue (const PetscInt p) const
{
  for(size_t i=0; i<_occupied; i++)
    if( _deriv_index[i] == p ) return (PetscScalar *)&(_deriv_value[i]);
  return 0;
}

bool AsmtDScalar::_withSamePattern(const AsmtDScalar & a) const
{
  if( _occupied != a._occupied ) return false;
  return (memcmp ((void*)a._deriv_index, _deriv_index, _occupied*sizeof(PetscInt)) == 0);
}


size_t  AsmtDScalar:: _patternSize(const AsmtDScalar & a) const
{
  size_t pattern_size = 0;

  size_t i = 0;
  size_t j = 0;
  while(true)
  {
    PetscInt ix,iy;
    bool end=true;

    if (i!=_occupied)
    {
      ix = _deriv_index[i];
      end = false;
    }
    else
      ix = std::numeric_limits<int>::max();

    if (j!=a._occupied)
    {
      iy = a._deriv_index[j];
      end = false;
    }
    else
      iy = std::numeric_limits<int>::max();

    if (end) break;

    PetscInt ii = std::min(ix,iy);
    if (ix<iy)
    {
      pattern_size++;
      i++;
    }
    else if (ix>iy)
    {
      pattern_size++;
      j++;
    }
    else
    {
      pattern_size++;
      i++;
      j++;
    }
  }

  return pattern_size;
}



void  AsmtDScalar:: _padPattern(const AsmtDScalar & a, AsmtDScalar & pa)
{
  size_t pattern_size = _patternSize(a);

  PetscInt index_buf[ADSMT_NUMBER_DIRECTIONS];
  PetscScalar value_buf[ADSMT_NUMBER_DIRECTIONS];

  PetscInt * index_p;
  PetscScalar * value_p;

  if(pattern_size <= ADSMT_NUMBER_DIRECTIONS)
  {
    index_p = (PetscInt *)index_buf;
    value_p = (PetscScalar *)value_buf;
  }
  else
  {
    index_p = new PetscInt[pattern_size];
    value_p = new PetscScalar[pattern_size];
  }

  pa._val = a._val;
  pa._occupied = pattern_size;
  if(pattern_size <= ADSMT_NUMBER_DIRECTIONS)
  {
    pa._deriv_index = (PetscInt *)pa._deriv_index_buf;
    pa._deriv_value = (PetscScalar *)pa._deriv_value_buf;
  }
  else
  {
    pa._deriv_index = new PetscInt[pattern_size];
    pa._deriv_value = new PetscScalar[pattern_size];
  }

  size_t i = 0;
  size_t j = 0;
  size_t k = 0;
  while(true)
  {
    PetscInt ix,iy;
    bool end=true;

    if (i!=_occupied)
    {
      ix = _deriv_index[i];
      end = false;
    }
    else
      ix = std::numeric_limits<int>::max();

    if (j!=a._occupied)
    {
      iy = a._deriv_index[j];
      end = false;
    }
    else
      iy = std::numeric_limits<int>::max();

    if (end) break;

    PetscInt ii = std::min(ix,iy);
    if (ix<iy)
    {
      index_p[k] = ii;
      value_p[k] = _deriv_value[i];
      pa._deriv_index[k] = ii;
      pa._deriv_value[k] = 0.0;
      i++;
    }
    else if (ix>iy)
    {
      index_p[k] = ii;
      value_p[k] = 0.0;
      pa._deriv_index[k] = ii;
      pa._deriv_value[k] = a._deriv_value[j];
      j++;
    }
    else
    {
      index_p[k] = ii;
      value_p[k] = _deriv_value[i];
      pa._deriv_index[k] = ii;
      pa._deriv_value[k] = a._deriv_value[j];
      i++;
      j++;
    }

    k++;
  }

  if(pattern_size <= ADSMT_NUMBER_DIRECTIONS)
  {
    _occupied = pattern_size;
    memcpy(_deriv_index, (void*)index_p, _occupied*sizeof(PetscInt));
    memcpy(_deriv_value, (void*)value_p, _occupied*sizeof(PetscScalar));
  }
  else
  {
    if( _occupied  > ADSMT_NUMBER_DIRECTIONS)
    {
      delete [] _deriv_index;
      delete [] _deriv_value;
    }
    _occupied = pattern_size;
    _deriv_index = index_p;
    _deriv_value = value_p;
  }
}

#endif
