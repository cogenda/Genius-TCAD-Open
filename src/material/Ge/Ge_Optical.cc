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
// Material Type: Ge


#include "PMI.h"




class GSS_Ge_Optical : public PMIS_Optical
{
private:

  static const RefractionItem WaveTable[];

  static const unsigned int table_size;

public:

  std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=1.12, PetscScalar Tl=1.0) const
  {
    std::complex<PetscScalar> n(1.0,0.0);

    if( lamda < WaveTable[0].wavelength*um )
      return std::complex<PetscScalar> (WaveTable[0].RefractionIndexRe, WaveTable[0].RefractionIndexIm);

    if( lamda > WaveTable[table_size-1].wavelength*um )
      return std::complex<PetscScalar> (WaveTable[table_size-1].RefractionIndexRe, WaveTable[table_size-1].RefractionIndexIm);

    for(unsigned int i=0; i<table_size-1; i++)
    {
      // do a linear interpolation
      if(lamda>=WaveTable[i].wavelength*um && lamda<=WaveTable[i+1].wavelength*um)
      {
        std::complex<PetscScalar> n1(WaveTable[i].RefractionIndexRe, WaveTable[i].RefractionIndexIm);
        std::complex<PetscScalar> n2(WaveTable[i+1].RefractionIndexRe, WaveTable[i+1].RefractionIndexIm);
        PetscScalar d1 = lamda - WaveTable[i].wavelength*um;
        PetscScalar d2 = WaveTable[i+1].wavelength*um - lamda;
        n = (n1*d2 + n2*d1)/(d1+d2);
        break;
      }
    }
    return n;
  }

  // constructions
public:
  GSS_Ge_Optical(const PMIS_Environment &env):PMIS_Optical(env) {}

  ~GSS_Ge_Optical() {}
}
;


const RefractionItem GSS_Ge_Optical::WaveTable[] =
  {
    {0.1240,  0.930,  0.860},
    {0.1350,  0.925,  1.000},
    {0.1378,  0.920,  1.140},
    {0.1459,  0.920,  1.200},
    {0.1550,  0.920,  1.400},
    {0.1653,  0.960,  1.600},
    {0.1771,  1.000,  1.800},
    {0.1907,  1.100,  2.050},
    {0.2066,  1.300,  2.340},
    {0.2254,  1.380,  2.842},
    {0.2480,  1.394,  3.197},
    {0.2755,  1.953,  4.297},
    {0.3100,  3.905,  3.336},
    {0.3306,  3.947,  2.922},
    {0.3542,  4.020,  2.667},
    {0.3815,  4.144,  2.405},
    {0.4133,  4.082,  2.145},
    {0.4275,  4.037,  2.140},
    {0.4428,  4.035,  2.181},
    {0.4592,  4.082,  2.240},
    {0.4769,  4.180,  2.309},
    {0.4959,  4.340,  2.384},
    {0.5166,  4.610,  2.455},
    {0.5391,  5.062,  2.318},
    {0.5636,  5.283,  2.049},
    {0.5904,  5.748,  1.634},
    {0.6199,  5.588,  0.933},
    {0.6526,  5.380,  0.638},
    {0.6888,  5.067,  0.500},
    {0.7293,  4.897,  0.401},
    {0.7749,  4.763,  0.345},
    {0.8266,  4.653,  0.298},
    {1.2400,  4.325,  0.081},
    {1.3780,  4.285,  0.075},
    {1.5500,  4.275,  0.057},
    {1.7710,  4.180,  0.028},
  };

const unsigned int GSS_Ge_Optical::table_size = sizeof(GSS_Ge_Optical::WaveTable)/sizeof(RefractionItem);


extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Optical*  PMIS_Ge_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_Ge_Optical(env);
  }
}
