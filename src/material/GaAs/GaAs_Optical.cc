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
/*       A Two-Dimensional General Purpose Semiconductor GaAsmulator.        */
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
// Material Type: GaAs


#include "PMI.h"


class GSS_GaAs_Optical : public PMIS_Optical
{
private:

  void _init_default_wave_table()
  {
    if(!_wave_table.empty()) _wave_table.clear();

    _wave_table.push_back( RefractionItem(0.0620,  1.025,  0.144) );
    _wave_table.push_back( RefractionItem(0.0827,  0.836,  0.282) );
    _wave_table.push_back( RefractionItem(0.1240,  0.913,  0.974) );
    _wave_table.push_back( RefractionItem(0.1378,  0.901,  1.136) );
    _wave_table.push_back( RefractionItem(0.1550,  0.899,  1.435) );
    _wave_table.push_back( RefractionItem(0.1771,  1.063,  1.838) );
    _wave_table.push_back( RefractionItem(0.2066,  1.264,  2.472) );
    _wave_table.push_back( RefractionItem(0.2254,  1.383,  2.936) );
    _wave_table.push_back( RefractionItem(0.2480,  2.273,  4.084) );
    _wave_table.push_back( RefractionItem(0.2755,  3.913,  2.919) );
    _wave_table.push_back( RefractionItem(0.3100,  3.601,  1.920) );
    _wave_table.push_back( RefractionItem(0.3542,  3.531,  2.013) );
    _wave_table.push_back( RefractionItem(0.4133,  4.509,  1.948) );
    _wave_table.push_back( RefractionItem(0.4275,  5.052,  1.721) );
    _wave_table.push_back( RefractionItem(0.4428,  4.959,  0.991) );
    _wave_table.push_back( RefractionItem(0.4592,  4.694,  0.696) );
    _wave_table.push_back( RefractionItem(0.4769,  4.492,  0.539) );
    _wave_table.push_back( RefractionItem(0.4959,  4.333,  0.441) );
    _wave_table.push_back( RefractionItem(0.5166,  4.205,  0.371) );
    _wave_table.push_back( RefractionItem(0.5391,  4.100,  0.320) );
    _wave_table.push_back( RefractionItem(0.5636,  4.013,  0.276) );
    _wave_table.push_back( RefractionItem(0.5904,  3.904,  0.240) );
    _wave_table.push_back( RefractionItem(0.6199,  3.878,  0.211) );
    _wave_table.push_back( RefractionItem(0.6526,  3.826,  0.179) );
    _wave_table.push_back( RefractionItem(0.6888,  3.785,  0.151) );
    _wave_table.push_back( RefractionItem(0.7293,  3.742,  0.112) );
    _wave_table.push_back( RefractionItem(0.7749,  3.700,  0.091) );
    _wave_table.push_back( RefractionItem(0.8266,  3.666,  0.080) );
    _wave_table.push_back( RefractionItem(0.8856,  3.614,  0.002) );
    _wave_table.push_back( RefractionItem(1.0000,  3.509,  8.5e-5) );
    _wave_table.push_back( RefractionItem(2.0000,  3.509,  0.000) );
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
  GSS_GaAs_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
    _init_default_wave_table();
  }

  ~GSS_GaAs_Optical()
  {}
}
;


extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Optical*  PMIS_GaAs_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_GaAs_Optical(env);
  }
}

