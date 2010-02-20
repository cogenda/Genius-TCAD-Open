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
/*                                                                           */
/*  Gong Ding                                                                */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: 4H-SiC.


#include "PMI.h"



class GSS_SiC4H_Optical : public PMIS_Optical
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


public:

  // constructions
  GSS_SiC4H_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
  }

  ~GSS_SiC4H_Optical()
  {}
}
;

//Source: FreeSnell
const RefractionItem GSS_SiC4H_Optical::WaveTable[] =
  {
    {0.1500,  4.021 , 1.720},
    {0.2000,  3.960 , 1.06},
    {0.2100,  3.759 , 0.702},
    {0.2200,  3.565 , 0.487},
    {0.2400,  3.252 , 0.307},
    {0.2600,  3.044 , 0.206},
    {0.2800,  2.990 , 0.195},
    {0.3000,  2.953 , 0.18},
    {0.3300,  2.896 , 0.156},
    {0.3500,  2.844 , 0.119},
    {0.4000,  2.764 , 0.1},
    {0.4200,  2.744 , 0.0443},
    {0.4500,  2.720 , 7.59e-4},
    {0.5000,  2.686 , 0.000},
    {0.6000,  2.644 , 0.000},
    {0.7000,  2.620 , 0.000},
    {0.8000,  2.601 , 0.000},
    {0.9000,  2.589 , 0.000},
    {1.0000 , 2.573,  0.000},
  };

const unsigned int GSS_SiC4H_Optical::table_size = sizeof(GSS_SiC4H_Optical::WaveTable)/sizeof(RefractionItem);

extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Optical*  PMIS_SiC4H_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_SiC4H_Optical(env);
  }
}
