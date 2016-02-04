#ifndef __std_ext_h__
#define __std_ext_h__

#include <iostream>
#include "config.h"

// VS2013  will support asinh etc, waiting...
#ifdef WINDOWS
#include "asinh.hpp" // for asinh
#include "acosh.hpp" // for acosh
#include "atanh.hpp" // for atanh
using namespace boost::math;
#endif



#if (WITH_PETSCSCALAR_FLOAT128 && __GNUC__ >= 4 && __GNUC_MINOR__ >=6) 
#include <quadmath.h>  
#endif

// For some reason the real std::max, std::min
// don't handle mixed compatible types
namespace std {
#if (WITH_PETSCSCALAR_FLOAT128 && __GNUC__ >= 4 && __GNUC_MINOR__ >=6) 
  inline __float128 max(__float128 a, double b)
  { return (a>b?a:b); }
  inline __float128 min(__float128 a, double b)
  { return (a<b?a:b); }

  inline __float128 max(double a, __float128 b)
  { return (a>b?a:b); }
  inline __float128 min(double a, __float128 b)
  { return (a<b?a:b); }
  
  inline __float128 abs(__float128 a)
  { return fabsq(a); }
  
  inline __float128 sign(__float128 a)
  { return (a>0 ? 1.0 : -1.0); }
  
  inline __float128 pow(double a, __float128 x)
  { return powq((__float128)a, x); }
  inline __float128 pow(__float128 a, __float128 x)
  { return powq(a, x); }
  inline __float128 pow(__float128 a, double x)
  { return powq(a, (__float128)x); }
  inline __float128 pow(__float128 a, int x)
  { return powq(a, (__float128)x); }
  
  inline __float128 log10(__float128 a)
  { return log10q(a); }
  
  inline __float128 tanh(__float128 a)
  { return tanhq(a); }
  
  inline std::ostream& operator << (std::ostream& os, const __float128 & t)
  {
    double tmp = t;
    os << tmp;
    return os;
  }
  
  inline std::istream& operator >> (std::istream& is, __float128 & t)
  {
    double tmp;
    is >> tmp;
    t = tmp;
    return is;
  }
#endif  
  
  inline long double max(long double a, double b)
  { return (a>b?a:b); }
  inline long double min(long double a, double b)
  { return (a<b?a:b); }

  inline long double max(double a, long double b)
  { return (a>b?a:b); }
  inline long double min(double a, long double b)
  { return (a<b?a:b); }

  inline double max(double a, float b)
  { return (a>b?a:b); }
  inline double min(double a, float b)
  { return (a<b?a:b); }

  inline double max(float a, double b)
  { return (a>b?a:b); }
  inline double min(float a, double b)
  { return (a<b?a:b); }

  inline long double max(long double a, float b)
  { return (a>b?a:b); }
  inline long double min(long double a, float b)
  { return (a<b?a:b); }

  inline long double max(float a, long double b)
  { return (a>b?a:b); }
  inline long double min(float a, long double b)
  { return (a<b?a:b); }

  inline unsigned int max(unsigned int a, unsigned int b)
  { return (a>b?a:b); }
  inline unsigned int min(unsigned int a, unsigned int b)
  { return (a<b?a:b); }


  inline double sign(float a)
  { return (a>0 ? 1.0 : -1.0); }

  inline double sign(double a)
  { return (a>0 ? 1.0 : -1.0); }

  inline long double sign(long double a)
  { return (a>0 ? 1.0 : -1.0); }

}




#endif
