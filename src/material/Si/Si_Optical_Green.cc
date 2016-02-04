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
// srouce:
// Green M A. Self-consistent optical parameters of intrinsic silicon at 300K including temperature coefficients[J].
// Solar Energy Materials and Solar Cells, 2008, 92(11): 1305-1310.

#include "PMI.h"



class GSS_Si_Optical_Green : public PMIS_Optical
{
private:


  struct GreenOpticalProperties
  {
     ///empty constructor
    GreenOpticalProperties() {}

    ///constructor
    GreenOpticalProperties(PetscScalar lambda, PetscScalar n, PetscScalar k, PetscScalar cn, PetscScalar ck)
      :wavelength(lambda), refraction(n, k), Tn(cn), Tk(ck)  {}

    ///wave length
    PetscScalar  wavelength;

    ///refraction index (n, k)
    std::complex<PetscScalar>  refraction;

    /// temperature coefficient for n
    PetscScalar Tn;

    /// temperature coefficient for k
    PetscScalar Tk;
  };

  std::vector<GreenOpticalProperties> optical_data_array;

  PetscScalar  _concentration;

  /**
   * build in refraction data in wavelength(um), n, alpha, k, c_n, c_alpha
   * NOTE: um = 1.23984190/eV
   */
  void _init_green_nk_wave_table()
  {
    double green_nk_table[][6] = {
      {    0.25, 1.84E6, 1.665, 3.665,  2.9, -0.9 },
      {    0.26, 1.97E6, 1.757, 4.084,  2.0, -1.5 },
      {    0.27, 2.18E6, 2.068, 4.680,   0 , -3.1 },
      {    0.28, 2.37E6, 2.959, 5.287, -4.8, -3.3 },
      {    0.29, 2.29E6, 4.356, 5.286, -9.0,  0.8 },
      {    0.30, 1.77E6, 4.976, 4.234, -3.8,  2.5 },
      {    0.31, 1.46E6, 5.121, 3.598, -1.6,  3.2 },
      {    0.32, 1.30E6, 5.112, 3.303, -1.3,  1.5 },
      {    0.33, 1.18E6, 5.195, 3.100, -1.2,  0.7 },
      {    0.34, 1.10E6, 5.301, 2.977, -1.0,  0.3 },
      {    0.35, 1.06E6, 5.494, 2.938, -1.8,  0.0 },
      {    0.36, 1.04E6, 6.026, 2.966, -4.1, -1.4 },
      {    0.37, 7.37E5, 6.891, 2.171, -4.4,  4.2 },
      {    0.38, 3.13E5, 6.616, 0.946, -2.3,  9.1 },
      {    0.39, 1.43E5, 6.039, 0.445, 1   , 26   },
      {    0.40, 9.30E4, 5.613, 0.296, 2.1 , 33   },
      {    0.41, 6.95E4, 5.330, 0.227, 2.1 , 31   },
      {    0.42, 5.27E4, 5.119, 0.176, 1.9 , 29   },
      {    0.43, 4.02E4, 4.949, 0.138, 1.8 , 29   },
      {    0.44, 3.07E4, 4.812, 0.107, 1.7 , 28   },
      {    0.45, 2.41E4, 4.691, 0.086, 1.6 , 28   },
      {    0.46, 1.95E4, 4.587, 0.071, 1.6 , 29   },
      {    0.47, 1.66E4, 4.497, 0.062, 1.5 , 29   },
      {    0.48, 1.44E4, 4.419, 0.055, 1.4 , 30   },
      {    0.49, 1.26E4, 4.350, 0.049, 1.4 , 30   },
      {    0.50, 1.11E4, 4.294, 0.044, 1.3 , 31   },
      {    0.51, 9700  , 4.241, 0.039, 1.3 , 31   },
      {    0.52, 8800  , 4.193, 0.036, 1.2 , 32   },
      {    0.53, 7850  , 4.151, 0.033, 1.2 , 33   },
      {    0.54, 7050  , 4.112, 0.030, 1.2 , 33   },
      {    0.55, 6390  , 4.077, 0.028, 1.1 , 33   },
      {    0.56, 5780  , 4.045, 0.026, 1.1 , 34   },
      {    0.57, 5320  , 4.015, 0.024, 1.1 , 34   },
      {    0.58, 4880  , 3.988, 0.023, 1.1 , 34   },
      {    0.59, 4490  , 3.963, 0.021, 1 ,   34   },
      {    0.60, 4175  , 3.940, 0.020, 1 ,   34   },
      {    0.61, 3800  , 3.918, 0.018, 1 ,   35   },
      {    0.62, 3520  , 3.898, 0.017, 1 ,   35   },
      {    0.63, 3280  , 3.879, 0.016, 1 ,   35   },
      {    0.64, 3030  , 3.861, 0.015, 1 ,   35   },
      {    0.65, 2790  , 3.844, 0.014, 0.9 , 35   },
      {    0.66, 2570  , 3.828, 0.013, 0.9 , 35   },
      {    0.67, 2390  , 3.813, 0.013, 0.9 , 36   },
      {    0.68, 2200  , 3.798, 0.012, 0.9 , 36   },
      {    0.69, 2040  , 3.784, 0.011, 0.9 , 36   },
      {    0.70, 1890  , 3.772, 0.011, 0.9 , 37   },
      {    0.71, 1780  , 3.759, 0.010, 0.9 , 37   },
      {    0.72, 1680  , 3.748, 0.010, 0.9 , 37   },
      {    0.73, 1540  , 3.737, 0.009, 0.8 , 37   },
      {    0.74, 1420  , 3.727, 0.008, 0.8 , 37   },
      {    0.75, 1310  , 3.717, 0.008, 0.8 , 37   },
      {    0.76, 1190  , 3.708, 0.007, 0.8 , 37   },
      {    0.77, 1100  , 3.699, 0.007, 0.8 , 37   },
      {    0.78, 1030  , 3.691, 0.006, 0.8 , 37   },
      {    0.79, 928   , 3.683, 0.006, 0.8 , 38   },
      {    0.8 , 850   , 3.675, 0.005, 0.8 , 40   },
      {    0.81, 775   , 3.668, 0.005, 0.8 , 41   },
      {    0.82, 707   , 3.661, 0.005, 0.8 , 42   },
      {    0.83, 647   , 3.654, 0.004, 0.7 , 44   },
      {    0.84, 590   , 3.647, 0.004, 0.7 , 45   },
      {    0.85, 534   , 3.641, 0.004, 0.7 , 46   },
      {    0.86, 479   , 3.635, 0.003, 0.7 , 47   },
      {    0.87, 431   , 3.630, 0.003, 0.7 , 49   },
      {    0.88, 383   , 3.624, 0.003, 0.7 , 51   },
      {    0.89, 343   , 3.619, 0.002, 0.7 , 52   },
      {    0.9 , 303   , 3.614, 0.002, 0.7 , 54   },
      {    0.91, 271   , 3.609, 0.002, 0.7 , 56   },
      {    0.92, 240   , 3.604, 0.002, 0.7 , 57   },
      {    0.93, 209   , 3.600, 0.002, 0.7 , 59   },
      {    0.94, 183   , 3.595, 0.001, 0.7 , 62   },
      {    0.95, 156   , 3.591, 0.001, 0.7 , 65   },
      {    0.96, 134   , 3.587, 0.001, 0.7 , 69   },
      {    0.97, 113   , 3.583, 0.001, 0.7 , 73   },
      {    0.98, 96.0  , 3.579, 0.001, 0.7 , 78   },
      {    0.99, 79.0  , 3.575, 0.001, 0.7 , 83   },
      {    1   , 64.0  , 3.572, 0.001, 0.7 , 90   },
      {    1.01, 51.1  , 3.568, 0.0  , 0.7 , 97   },
      {    1.02, 39.9  , 3.565, 0.0  , 0.7 , 105  },
      {    1.03, 30.2  , 3.562, 0.0  , 0.7 , 112  },
      {    1.04, 22.6  , 3.559, 0.0  , 0.6 , 120  },
      {    1.05, 16.3  , 3.556, 0.0  , 0.6 , 135  },
      {    1.06, 11.1  , 3.553, 0.0  , 0.6 , 145  },
      {    1.07, 8.0   , 3.550, 0.0  , 0.6 , 155  },
      {    1.08, 6.2   , 3.547, 0.0  , 0.6 , 160  },
      {    1.09, 4.7   , 3.545, 0.0  , 0.6 , 165  },
      {    1.1 , 3.5   , 3.542, 0.0  , 0.6 , 175  },
      {    1.11, 2.7   , 3.540, 0.0  , 0.6 , 180  },
      {    1.12, 2.0   , 3.537, 0.0  , 0.6 , 185  },
      {    1.13, 1.5   , 3.535, 0.0  , 0.6 , 190  },
      {    1.14, 1.0   , 3.532, 0.0  , 0.6 , 200  },
      {    1.15, 0.68  , 3.530, 0.0  , 0.6 , 210  },
      {    1.16, 0.42  , 3.528, 0.0  , 0.6 , 230  },
      {    1.17, 0.22  , 3.526, 0.0  , 0.6 , 260  },
      {    1.18, 6.5E-2, 3.524, 0.0  , 0.6 , 320  },
      {    1.19, 3.6E-2, 3.522, 0.0  , 0.6 , 345  },
      {    1.2 , 2.2E-2, 3.520, 0.0  , 0.6 , 355  },
      {    1.21, 1.3E-2, 3.518, 0.0  , 0.6 , 380  },
      {    1.22, 8.2E-3, 3.517, 0.0  , 0.6 , 390  },
      {    1.23, 4.7E-3, 3.515, 0.0  , 0.6 , 405  },
      {    1.24, 2.4E-3, 3.513, 0.0  , 0.6 , 410  },
      {    1.25, 1.0E-3, 3.511, 0.0  , 0.6 , 430  },
      {    1.26, 3.6E-4, 3.509, 0.0  , 0.6 , 440  },
      {    1.27, 2.0E-4, 3.508, 0.0  , 0.6 , 455  },
      {    1.28, 1.2E-4, 3.506, 0.0  , 0.6 , 470  },
      {    1.29, 7.1E-5, 3.505, 0.0  , 0.6 , 500  },
      {    1.3 , 4.5E-5, 3.503, 0.0  , 0.6 , 525  },
      {    1.31, 2.7E-5, 3.502, 0.0  , 0.6 , 550  },
      {    1.32, 1.6E-5, 3.500, 0.0  , 0.6 , 580  },
      {    1.33, 8.0E-6, 3.499, 0.0  , 0.6 , 610  },
      {    1.34, 3.5E-6, 3.497, 0.0  , 0.6 , 650  },
      {    1.35, 1.7E-6, 3.496, 0.0  , 0.6 , 670  },
      {    1.36, 9.5E-7, 3.495, 0.0  , 0.6 , 675  },
      {    1.37, 6.0E-7, 3.494, 0.0  , 0.6 , 680  },
      {    1.38, 3.8E-7, 3.492, 0.0  , 0.6 , 685  },
      {    1.39, 2.3E-7, 3.491, 0.0  , 0.6 , 690  },
      {    1.4 , 1.4E-7, 3.490, 0.0  , 0.6 , 700  },
      {    1.41, 8.5E-8, 3.489, 0.0  , 0.6 , 710  },
      {    1.42, 5.0E-8, 3.488, 0.0  , 0.6 , 720  },
      {    1.43, 2.5E-8, 3.487, 0.0  , 0.6 , 730  },
      {    1.44, 1.8E-8, 3.486, 0.0  , 0.6 , 740  },
      {    1.45, 1.2E-8, 3.485, 0.0  , 0.6 , 750  },
    };

    unsigned int dim = sizeof(green_nk_table)/(sizeof(green_nk_table[0]));
    for(unsigned int i=0; i<dim; i++)
    {
      PetscScalar lambda = green_nk_table[i][0]*um;
      PetscScalar alpha  = green_nk_table[i][1]/cm;
      PetscScalar n      = green_nk_table[i][2];
      PetscScalar k      = alpha*lambda/12.5663706144;
      PetscScalar cn     = green_nk_table[i][4]*1e-4/K;
      PetscScalar ck     = green_nk_table[i][5]*1e-4/K;

      optical_data_array.push_back(GreenOpticalProperties(lambda, n, k, cn, ck));
    }
  }


public:

  std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Tl, PetscScalar Eg=0) const
  {
    std::complex<PetscScalar> n(1.0,0.0);
    PetscScalar  cn=1, ck=1;

    unsigned int table_size = optical_data_array.size();

    if( lamda < optical_data_array[0].wavelength )
    {
      n = optical_data_array[0].refraction;
      cn = optical_data_array[0].Tn;
      ck = optical_data_array[0].Tk;
    }


    if( lamda > optical_data_array[table_size-1].wavelength )
    {
      n = optical_data_array[table_size-1].refraction;
      cn = optical_data_array[table_size-1].Tn;
      ck = optical_data_array[table_size-1].Tk;
    }


    for(unsigned int i=0; i<table_size-1; i++)
    {
      // do a linear interpolation
      if(lamda>=optical_data_array[i].wavelength && lamda<=optical_data_array[i+1].wavelength)
      {
        std::complex<PetscScalar> n1 = optical_data_array[i].refraction;
        PetscScalar cn1 = optical_data_array[i].Tn;
        PetscScalar ck1 = optical_data_array[i].Tk;

        std::complex<PetscScalar> n2 = optical_data_array[i+1].refraction;
        PetscScalar cn2 = optical_data_array[i+1].Tn;
        PetscScalar ck2 = optical_data_array[i+1].Tk;

        PetscScalar d1 = lamda - optical_data_array[i].wavelength;
        PetscScalar d2 = optical_data_array[i+1].wavelength - lamda;

        n = (n1*d2 + n2*d1)/(d1+d2);
        cn = (cn1*d2 + cn2*d1)/(d1+d2);
        ck = (ck1*d2 + ck2*d1)/(d1+d2);

        break;
      }
    }

    PetscScalar bn = cn*300*K;
    PetscScalar bk = ck*300*K;

    std::complex<PetscScalar> nn(n.real()*std::pow(Tl/(300*K), bn), n.imag()*std::pow(Tl/(300*K), bk));

    return nn;
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
  GSS_Si_Optical_Green(const PMIS_Environment &env):PMIS_Optical(env)
  {
    PMI_Info = "This is the Green optical model of Silicon";
    _init_green_nk_wave_table();

    _concentration = std::pow(cm, -3);
  }

  ~GSS_Si_Optical_Green() {}

};




extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Optical*  PMIS_Si_Optical_Green (const PMIS_Environment& env)
  {
    return new GSS_Si_Optical_Green(env);
  }
}
