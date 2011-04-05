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
// Material Type: Metal (Al)


#include "PMI.h"


class GSS_Elec_Optical : public PMIC_Optical
{
private:

  void _init_default_wave_table()
  {
    if(!_wave_table.empty()) _wave_table.clear();

    //Source: MEDICI data base
    _wave_table.push_back( RefractionItem(0.0165,  1.010,  0.024) );
    _wave_table.push_back( RefractionItem(0.0172,  1.020,  0.0035) );
    _wave_table.push_back( RefractionItem(0.0376,  0.902,  0.0102) );
    _wave_table.push_back( RefractionItem(0.0620,  0.668,  0.0268) );
    _wave_table.push_back( RefractionItem(0.0800,  0.258,  0.0777) );
    _wave_table.push_back( RefractionItem(0.1033,  0.0328, 0.791) );
    _wave_table.push_back( RefractionItem(0.1600,  0.0765, 1.740) );
    _wave_table.push_back( RefractionItem(0.2066,  0.130,  2.39) );
    _wave_table.push_back( RefractionItem(0.3000,  0.276,  3.61) );
    _wave_table.push_back( RefractionItem(0.4000,  0.490,  4.86) );
    _wave_table.push_back( RefractionItem(0.5000,  0.769,  6.08) );
    _wave_table.push_back( RefractionItem(0.6000,  1.20,   7.26) );
    _wave_table.push_back( RefractionItem(0.7000,  1.83,   8.31) );
    _wave_table.push_back( RefractionItem(0.8000,  2.80,   8.45) );
    _wave_table.push_back( RefractionItem(0.9000,  2.06,   8.30) );
    _wave_table.push_back( RefractionItem(1.0000,  1.35,   9.58) );
    _wave_table.push_back( RefractionItem(1.2000,  1.21,   12.0) );
    _wave_table.push_back( RefractionItem(1.5000,  1.38,   15.4) );
    _wave_table.push_back( RefractionItem(2.0000,  2.15,   20.7) );
    _wave_table.push_back( RefractionItem(4.0000,  6.43,   39.8) );
    _wave_table.push_back( RefractionItem(10.000,  25.3,   89.8) );
    _wave_table.push_back( RefractionItem(15.000,  44.01,  1.2e2) );
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
  GSS_Elec_Optical(const PMIC_Environment &env):PMIC_Optical(env)
  {
    _init_default_wave_table();
  }

  ~GSS_Elec_Optical(){}
}
;


extern "C"
{
  DLL_EXPORT_DECLARE  PMIC_Optical*  PMIC_Elec_Optical_Default (const PMIC_Environment& env)
  {
    return new GSS_Elec_Optical(env);
  }
}
