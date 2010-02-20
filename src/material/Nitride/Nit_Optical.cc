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
// Material Type: Nitride


#include "PMI.h"





class GSS_Nitride_Optical : public PMII_Optical
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
  GSS_Nitride_Optical(const PMII_Environment &env):PMII_Optical(env)  {}

  ~GSS_Nitride_Optical() {}
}
;



const RefractionItem GSS_Nitride_Optical::WaveTable[] =
  {
    {0.0512,  0.655,  0.420   },
    {0.0539,  0.625,  0.481   },
    {0.0564,  0.611,  0.560   },
    {0.0590,  0.617,  0.647   },
    {0.0620,  0.635,  0.743   },
    {0.0653,  0.676,  0.841   },
    {0.0689,  0.735,  0.936   },
    {0.0729,  0.810,  1.03    },
    {0.0775,  0.902,  1.11    },
    {0.0827,  1.001,  1.18    },
    {0.0886,  1.111,  1.26    },
    {0.0953,  1.247,  1.35    },
    {0.1033,  1.417,  1.43    },
    {0.1127,  1.657,  1.52    },
    {0.1181,  1.827,  1.53    },
    {0.1240,  2.000,  1.49    },
    {0.1305,  2.162,  1.44    },
    {0.1378,  2.326,  1.32    },
    {0.1459,  2.492,  1.16    },
    {0.1550,  2.651,  0.962   },
    {0.1600,  2.711,  0.866   },//7.75eV
    {0.1653,  2.753,  0.750   },
    {0.1710,  2.766,  0.612   },
    {0.1771,  2.752,  0.493   },
    {0.1837,  2.724,  0.380   },
    {0.1907,  2.682,  0.273   },
    {0.1984,  2.620,  0.174   },
    {0.2066,  2.541,  0.102   },
    {0.2156,  2.464,  5.7e-2  },
    {0.2254,  2.393,  2.9e-2  },
    {0.2342,  2.331,  1.1e-2  },
    {0.2380,  2.278,  4.9e-2  },
    {0.2610,  2.234,  1.2e-3  },
    {0.2755,  2.198,  2.2e-4  },
    {0.2917,  2.167,  0.0     },
    {0.3100,  2.141,  0.0     },
    {0.3542,  2.099,  0.0     },
    {0.4133,  2.066,  0.0     },
    {0.4959,  2.041,  0.0     },
    {0.6199,  2.022,  0.0     },
    {0.8266,  2.008,  0.0     },
    {1.240 ,  1.998,  0.0     },//1.0eV
  };


const unsigned int GSS_Nitride_Optical::table_size = sizeof(GSS_Nitride_Optical::WaveTable)/sizeof(RefractionItem);


extern "C"
{
  DLL_EXPORT_DECLARE  PMII_Optical*  PMII_Nitride_Optical_Default (const PMII_Environment& env)
  {
    return new GSS_Nitride_Optical(env);
  }
}
