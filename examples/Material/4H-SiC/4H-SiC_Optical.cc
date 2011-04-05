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
// Material Type: 4H-SiC.


#include "PMI.h"



class GSS_SiC4H_Optical : public PMIS_Optical
{
private:

  void _init_default_wave_table()
  {
    if(!_wave_table.empty()) _wave_table.clear();

    //Source: FreeSnell
    _wave_table.push_back( RefractionItem(0.1500,  4.021 , 1.720));
    _wave_table.push_back( RefractionItem(0.2000,  3.960 , 1.06));
    _wave_table.push_back( RefractionItem(0.2100,  3.759 , 0.702));
    _wave_table.push_back( RefractionItem(0.2200,  3.565 , 0.487));
    _wave_table.push_back( RefractionItem(0.2400,  3.252 , 0.307));
    _wave_table.push_back( RefractionItem(0.2600,  3.044 , 0.206));
    _wave_table.push_back( RefractionItem(0.2800,  2.990 , 0.195));
    _wave_table.push_back( RefractionItem(0.3000,  2.953 , 0.18));
    _wave_table.push_back( RefractionItem(0.3300,  2.896 , 0.156));
    _wave_table.push_back( RefractionItem(0.3500,  2.844 , 0.119));
    _wave_table.push_back( RefractionItem(0.4000,  2.764 , 0.1));
    _wave_table.push_back( RefractionItem(0.4200,  2.744 , 0.0443));
    _wave_table.push_back( RefractionItem(0.4500,  2.720 , 7.59e-4));
    _wave_table.push_back( RefractionItem(0.5000,  2.686 , 0.000));
    _wave_table.push_back( RefractionItem(0.6000,  2.644 , 0.000));
    _wave_table.push_back( RefractionItem(0.7000,  2.620 , 0.000));
    _wave_table.push_back( RefractionItem(0.8000,  2.601 , 0.000));
    _wave_table.push_back( RefractionItem(0.9000,  2.589 , 0.000));
    _wave_table.push_back( RefractionItem(1.0000 , 2.573,  0.000));
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


public:

  // constructions
  GSS_SiC4H_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
    _init_default_wave_table();
  }

  ~GSS_SiC4H_Optical()
  {}

};


extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Optical*  PMIS_SiC4H_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_SiC4H_Optical(env);
  }
}
