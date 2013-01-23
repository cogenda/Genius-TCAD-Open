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
// Material Type: Cu


#include "PMI.h"



class GSS_Cu_Optical : public PMIC_Optical
{
private:

  void _init_default_wave_table()
  {
    if(!_wave_table.empty()) _wave_table.clear();

    //Source: FreeSnell
    _wave_table.push_back( RefractionItem(0.0200,  0.97,    0.11) );
    _wave_table.push_back( RefractionItem(0.0400,  0.8799,  0.2401) );
    _wave_table.push_back( RefractionItem(0.0600,  0.8935,  0.423) );
    _wave_table.push_back( RefractionItem(0.0800,  0.9801,  0.6901) );
    _wave_table.push_back( RefractionItem(0.0900,  1.064,   0.72) );
    _wave_table.push_back( RefractionItem(0.1000,  1.086,   0.726) );
    _wave_table.push_back( RefractionItem(0.1500,  1.03,    1.00) );
    _wave_table.push_back( RefractionItem(0.2000,  0.988,   1.504) );
    _wave_table.push_back( RefractionItem(0.3000,  1.393,   1.670) );
    _wave_table.push_back( RefractionItem(0.4000,  1.175,   2.130) );
    _wave_table.push_back( RefractionItem(0.5000,  1.134,   2.570 ) );
    _wave_table.push_back( RefractionItem(0.6000,  0.747,   3.362 ) );
    _wave_table.push_back( RefractionItem(0.7000,  0.412,   4.203 ) );
    _wave_table.push_back( RefractionItem(0.8000,  0.454,   4.978 ) );
    _wave_table.push_back( RefractionItem(0.9000,  0.496,   5.754 ) );
    _wave_table.push_back( RefractionItem(1.0000,  0.538,   6.530 ) );
    _wave_table.push_back( RefractionItem(1.5000,  0.733,   10.017) );
    _wave_table.push_back( RefractionItem(2.0000,  0.879,   13.400) );
    _wave_table.push_back( RefractionItem(3.0000,  1.562,   20.055) );
    _wave_table.push_back( RefractionItem(4.0000,  2.250,   26.500) );
    _wave_table.push_back( RefractionItem(5.0000,  3.260,   33.000) );
    _wave_table.push_back( RefractionItem(6.0000,  3.944,   39.078) );
    _wave_table.push_back( RefractionItem(7.0000,  4.942,   45.336) );
    _wave_table.push_back( RefractionItem(8.0000,  6.074,   51.292) );
    _wave_table.push_back( RefractionItem(9.0000,  7.118,   57.167) );
    _wave_table.push_back( RefractionItem(10.000,  8.310,   63.000) );
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
  GSS_Cu_Optical(const PMIC_Environment &env):PMIC_Optical(env)
  {
    _init_default_wave_table();
  }

  ~GSS_Cu_Optical(){}
}
;


extern "C"
{
  DLL_EXPORT_DECLARE  PMIC_Optical*  PMIC_Cu_Optical_Default (const PMIC_Environment& env)
  {
    return new GSS_Cu_Optical(env);
  }
}
