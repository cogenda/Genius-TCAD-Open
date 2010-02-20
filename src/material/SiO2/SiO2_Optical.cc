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
// Material Type: SiO2


#include "PMI.h"



class GSS_SiO2_Optical : public PMII_Optical
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
  GSS_SiO2_Optical(const PMII_Environment &env):PMII_Optical(env) {}

  ~GSS_SiO2_Optical() {}
}
;


const RefractionItem GSS_SiO2_Optical::WaveTable[] =
  {
    {0.0496,  0.733,  0.325},
    {0.0550,  0.822,  0.408},
    {0.0620,  0.859,  0.585},
    {0.0689,  0.957,  0.717},
    {0.0775,  1.172,  0.808},
    {0.0886,  1.265,  0.861},
    {0.1033,  1.475,  0.323},
    {0.1240,  2.330,  0.056},
    {0.1305,  2.092,  1.89e-2},
    {0.1378,  1.904,  5.57e-3},
    {0.1459,  1.702,  3.2e-5},
    {0.1550,  1.600,  0.e0},
    {0.1771,  1.567,  0.e0},
    {0.1907,  1.543,  0.e0},
    {0.2066,  1.520,  0.e0},
    {0.2302,  1.500,  0.e0},
    {0.2652,  1.490,  0.e0},
    {0.2894,  1.480,  0.e0},
    {0.3303,  1.475,  0.e0},
    {0.3610,  1.469,  0.e0},
    {0.4046,  1.460,  0.e0},
    {0.5460,  1.455,  0.e0},
    {1.0139,  1.450,  0.e0},
    {1.5295,  1.444,  0.e0},
    {1.8131,  1.440,  0.e0},
    {2.4374,  1.430,  0.e0},
    {3.3026,  1.411,  0.e0},
    {3.5564,  1.404,  0.e0},
    {5.0000,  1.342,  0.e0},
    {6.0000,  1.266,  0.e0},
    {7.0000,  1.093,  0.e0},
    {7.5000,  0.9025, 0.e0},
    {7.5000,  0.9025, 0.e0},
  };

const unsigned int GSS_SiO2_Optical::table_size = sizeof(GSS_SiO2_Optical::WaveTable)/sizeof(RefractionItem);


extern "C"
{
  DLL_EXPORT_DECLARE  PMII_Optical*  PMII_SiO2_Optical_Default (const PMII_Environment& env)
  {
    return new GSS_SiO2_Optical(env);
  }
}
