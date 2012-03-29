#ifndef __anti_reflection_coating_h__
#define __anti_reflection_coating_h__

#include <vector>
#include <complex>

#include "point.h"

struct ARCoatings
{
  unsigned int inner_region;

  /// total anti-reflection coating layers
  unsigned int layers;
  /// complex refractive index
  std::vector< std::complex<double> > layer_refractive_index;
  /// thickness
  std::vector<double> layer_thickness;
};

#endif

