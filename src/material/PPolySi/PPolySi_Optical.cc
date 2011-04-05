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
/*       A Two-Dimensional General Purpose Semiconductor Simulator.          */
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
// Material Type: Amorphous Silicon


#include "PMI.h"



class GSS_PPolySi_Optical : public PMIC_Optical
{
private:

  void _init_default_wave_table()
  {
    if(!_wave_table.empty()) _wave_table.clear();

    _wave_table.push_back( RefractionItem(0.1240,  0.466,  1.090) );
    _wave_table.push_back( RefractionItem(0.1378,  0.499,  1.330) );
    _wave_table.push_back( RefractionItem(0.1550,  0.554,  1.660) );
    _wave_table.push_back( RefractionItem(0.1771,  0.670,  2.080) );
    _wave_table.push_back( RefractionItem(0.2066,  0.961,  2.650) );
    _wave_table.push_back( RefractionItem(0.2480,  1.660,  3.380) );
    _wave_table.push_back( RefractionItem(0.2755,  2.310,  3.710) );
    _wave_table.push_back( RefractionItem(0.3100,  3.360,  3.920) );
    _wave_table.push_back( RefractionItem(0.3543,  4.590,  3.380) );
    _wave_table.push_back( RefractionItem(0.4133,  5.430,  2.190) );
    _wave_table.push_back( RefractionItem(0.4960,  5.250,  0.992) );
    _wave_table.push_back( RefractionItem(0.6199,  4.710,  0.217) );
    _wave_table.push_back( RefractionItem(0.8190,  3.900,  2.23e-2) );
    _wave_table.push_back( RefractionItem(0.9613,  3.660,  1.1e-2) );
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
  GSS_PPolySi_Optical(const PMIC_Environment &env):PMIC_Optical(env)
  {
    _init_default_wave_table();
  }

  ~GSS_PPolySi_Optical() {}
}
;



extern "C"
{
  DLL_EXPORT_DECLARE  PMIC_Optical*  PMIC_PPolySi_Optical_Default (const PMIC_Environment& env)
  {
    return new GSS_PPolySi_Optical(env);
  }
}
