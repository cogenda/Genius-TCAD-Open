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
// Material Type: Silicon


#include "PMI.h"



class GSS_Si_Optical : public PMIS_Optical
{
private:

  /**
   * build in refraction data in wavelength(um), n, k
   * NOTE: um = 1.2398/eV
   */
  void _init_default_wave_table()
  {
    // from CRC Handbook of Chemistry and Physics
    if(!_wave_table.empty()) _wave_table.clear();

    _wave_table.push_back( RefractionItem(0.1240,  0.306 , 1.38) );
    _wave_table.push_back( RefractionItem(0.1305,  0.332 , 1.51) );
    _wave_table.push_back( RefractionItem(0.1378,  0.367 , 1.66) );
    _wave_table.push_back( RefractionItem(0.1459,  0.414 , 1.82) );
    _wave_table.push_back( RefractionItem(0.1550,  0.478 , 2.00) );
    _wave_table.push_back( RefractionItem(0.1653,  0.563 , 2.21) );
    _wave_table.push_back( RefractionItem(0.1771,  0.682 , 2.45) );
    _wave_table.push_back( RefractionItem(0.1907,  0.847 , 2.73) );
    _wave_table.push_back( RefractionItem(0.2066,  1.11  , 3.05) );
    _wave_table.push_back( RefractionItem(0.2254,  1.340 , 3.302) );
    _wave_table.push_back( RefractionItem(0.2480,  1.570 , 3.565) );
    _wave_table.push_back( RefractionItem(0.2755,  2.451 , 5.082) );
    _wave_table.push_back( RefractionItem(0.3100,  5.010 , 3.650) );
    _wave_table.push_back( RefractionItem(0.3306,  5.105 , 3.111) );
    _wave_table.push_back( RefractionItem(0.3542,  5.610 , 3.014) );
    _wave_table.push_back( RefractionItem(0.3815,  6.389 , 0.880) );
    _wave_table.push_back( RefractionItem(0.4133,  5.222 , 0.269) );
    _wave_table.push_back( RefractionItem(0.4275,  4.961 , 0.203) );
    _wave_table.push_back( RefractionItem(0.4428,  4.753 , 0.163) );
    _wave_table.push_back( RefractionItem(0.4592,  4.583 , 0.130) );
    _wave_table.push_back( RefractionItem(0.4769,  4.442 , 0.090) );
    _wave_table.push_back( RefractionItem(0.4959,  4.320 , 0.073) );
    _wave_table.push_back( RefractionItem(0.5166,  4.215 , 0.060) );
    _wave_table.push_back( RefractionItem(0.5391,  4.123 , 0.048) );
    _wave_table.push_back( RefractionItem(0.5636,  4.042 , 0.032) );
    _wave_table.push_back( RefractionItem(0.5904,  3.969 , 0.030) );
    _wave_table.push_back( RefractionItem(0.6199,  3.906 , 0.022) );
    _wave_table.push_back( RefractionItem(0.6526,  3.847 , 0.016) );
    _wave_table.push_back( RefractionItem(0.6888,  3.796 , 0.013) );
    _wave_table.push_back( RefractionItem(0.7293,  3.752 , 0.010) );
    _wave_table.push_back( RefractionItem(0.7749,  3.714 , 0.008) );
    _wave_table.push_back( RefractionItem(0.8266,  3.673 , 0.005) );   // 1.5eV
    _wave_table.push_back( RefractionItem(0.8856,  3.673 , 7.75e-3) ); // 1.4
    _wave_table.push_back( RefractionItem(0.9537,  3.673 , 2.26e-3) ); // 1.3
    _wave_table.push_back( RefractionItem(1.0332,  3.673 , 1.80e-4) ); // 1.2
    _wave_table.push_back( RefractionItem(1.1271,  3.5341, 1.30e-5) ); // 1.1
    _wave_table.push_back( RefractionItem(1.200 ,  3.5193, 2.50e-9));  // 1.033
  }

  PetscScalar  _concentration;

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

  PetscScalar FreeCarrierAbsorption(PetscScalar lamda, PetscScalar n, PetscScalar p, PetscScalar ) const
  {
    // source:
    // Free Carrier Absorption in Silicon
    // IEEE JOURNAL OF SOLID-STATE CIRCUITS, VOL. SC-13, NO. 1,FEBRUARY 1978
    // use experiment value

    // should use um and cm^-3
    lamda /= um;
    n /= _concentration;
    p /= _concentration;

    return (1e-18*lamda*lamda*n + 2.7e-18*lamda*lamda*p)/cm;
  }



  // constructions
public:
  GSS_Si_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
    PMI_Info = "This is the optical model of Silicon";
    _init_default_wave_table();

    _concentration = std::pow(cm, -3);
  }

  ~GSS_Si_Optical() {}

};




extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Optical*  PMIS_Si_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_Si_Optical(env);
  }
}
