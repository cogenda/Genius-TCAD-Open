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
  GSS_PPolySi_Optical(const PMIC_Environment &env):PMIC_Optical(env)  {}

  ~GSS_PPolySi_Optical() {}
}
;


const RefractionItem GSS_PPolySi_Optical::WaveTable[] =
  {
    {0.1240,  0.466,  1.090},
    {0.1378,  0.499,  1.330},
    {0.1550,  0.554,  1.660},
    {0.1771,  0.670,  2.080},
    {0.2066,  0.961,  2.650},
    {0.2480,  1.660,  3.380},
    {0.2755,  2.310,  3.710},
    {0.3100,  3.360,  3.920},
    {0.3543,  4.590,  3.380},
    {0.4133,  5.430,  2.190},
    {0.4960,  5.250,  0.992},
    {0.6199,  4.710,  0.217},
    {0.8190,  3.900,  2.23e-2},
    {0.9613,  3.660,  1.1e-2},
  };


const unsigned int GSS_PPolySi_Optical::table_size = sizeof(GSS_PPolySi_Optical::WaveTable)/sizeof(RefractionItem);


extern "C"
{
  DLL_EXPORT_DECLARE  PMIC_Optical*  PMIC_PPolySi_Optical_Default (const PMIC_Environment& env)
  {
    return new GSS_PPolySi_Optical(env);
  }
}
