/*****************************************************************************/
/*                                                                           */
/*              8888888         88888888         88888888                    */
/*            8                8                8                            */
/*           8                 8                8                            */
/*           8                  88888888         88888888                    */
/*           8      8888                8                8                   */
/*            8       8                 8                8                   */
/*              888888         888888888        888888888                    */
/*                                                                           */
/*       A Two-Dimensional General Purpose Semiconductor Gemulator.          */
/*                                                                           */
/*  GSS material database Version 0.4                                        */
/*  Last update: Feb 17, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: Germanium


#include "PMI.h"




class GSS_Ge_Optical : public PMIS_Optical
{
private:

  void _init_default_wave_table()
  {
    if(!_wave_table.empty()) _wave_table.clear();

    _wave_table.push_back( RefractionItem(0.1240,  0.930,  0.860) );
    _wave_table.push_back( RefractionItem(0.1350,  0.925,  1.000) );
    _wave_table.push_back( RefractionItem(0.1378,  0.920,  1.140) );
    _wave_table.push_back( RefractionItem(0.1459,  0.920,  1.200) );
    _wave_table.push_back( RefractionItem(0.1550,  0.920,  1.400) );
    _wave_table.push_back( RefractionItem(0.1653,  0.960,  1.600) );
    _wave_table.push_back( RefractionItem(0.1771,  1.000,  1.800) );
    _wave_table.push_back( RefractionItem(0.1907,  1.100,  2.050) );
    _wave_table.push_back( RefractionItem(0.2066,  1.300,  2.340) );
    _wave_table.push_back( RefractionItem(0.2254,  1.380,  2.842) );
    _wave_table.push_back( RefractionItem(0.2480,  1.394,  3.197) );
    _wave_table.push_back( RefractionItem(0.2755,  1.953,  4.297) );
    _wave_table.push_back( RefractionItem(0.3100,  3.905,  3.336) );
    _wave_table.push_back( RefractionItem(0.3306,  3.947,  2.922) );
    _wave_table.push_back( RefractionItem(0.3542,  4.020,  2.667) );
    _wave_table.push_back( RefractionItem(0.3815,  4.144,  2.405) );
    _wave_table.push_back( RefractionItem(0.4133,  4.082,  2.145) );
    _wave_table.push_back( RefractionItem(0.4275,  4.037,  2.140) );
    _wave_table.push_back( RefractionItem(0.4428,  4.035,  2.181) );
    _wave_table.push_back( RefractionItem(0.4592,  4.082,  2.240) );
    _wave_table.push_back( RefractionItem(0.4769,  4.180,  2.309) );
    _wave_table.push_back( RefractionItem(0.4959,  4.340,  2.384) );
    _wave_table.push_back( RefractionItem(0.5166,  4.610,  2.455) );
    _wave_table.push_back( RefractionItem(0.5391,  5.062,  2.318) );
    _wave_table.push_back( RefractionItem(0.5636,  5.283,  2.049) );
    _wave_table.push_back( RefractionItem(0.5904,  5.748,  1.634) );
    _wave_table.push_back( RefractionItem(0.6199,  5.588,  0.933) );
    _wave_table.push_back( RefractionItem(0.6526,  5.380,  0.638) );
    _wave_table.push_back( RefractionItem(0.6888,  5.067,  0.500) );
    _wave_table.push_back( RefractionItem(0.7293,  4.897,  0.401) );
    _wave_table.push_back( RefractionItem(0.7749,  4.763,  0.345) );
    _wave_table.push_back( RefractionItem(0.8266,  4.653,  0.298) );
    _wave_table.push_back( RefractionItem(1.2400,  4.325,  0.081) );
    _wave_table.push_back( RefractionItem(1.3780,  4.285,  0.075) );
    _wave_table.push_back( RefractionItem(1.5500,  4.275,  0.057) );
    _wave_table.push_back( RefractionItem(1.7710,  4.180,  0.028) );
    _wave_table.push_back( RefractionItem(2.0000,  4.180,  0.000) );
  }

public:

  std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Tl, PetscScalar Eg=0) const
  {
    std::complex<PetscScalar> n(1.0,0.0);
    unsigned int table_size = _wave_table.size();

    if( lamda < _wave_table[0].wavelength*um )
      return std::complex<PetscScalar> (_wave_table[0].RefractionIndexRe, _wave_table[0].RefractionIndexIm);

    if( lamda > _wave_table[table_size-1].wavelength*um )
      return std::complex<PetscScalar> (_wave_table[table_size-1].RefractionIndexRe, _wave_table[table_size-1].RefractionIndexIm);

    for(unsigned int i=0; i<table_size-1; i++)
    {
      // do a linear interpolation
      if(lamda>=_wave_table[i].wavelength*um && lamda<=_wave_table[i+1].wavelength*um)
      {
        std::complex<PetscScalar> n1(_wave_table[i].RefractionIndexRe, _wave_table[i].RefractionIndexIm);
        std::complex<PetscScalar> n2(_wave_table[i+1].RefractionIndexRe, _wave_table[i+1].RefractionIndexIm);
        PetscScalar d1 = lamda - _wave_table[i].wavelength*um;
        PetscScalar d2 = _wave_table[i+1].wavelength*um - lamda;
        n = (n1*d2 + n2*d1)/(d1+d2);
        break;
      }
    }
    return n;
  }

  // constructions
public:
  GSS_Ge_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
    _init_default_wave_table();
  }

  ~GSS_Ge_Optical() {}
}
;



extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Optical*  PMIS_Ge_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_Ge_Optical(env);
  }
}
