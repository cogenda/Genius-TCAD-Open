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

class GSS_SiO2_BandStructure : public PMII_BandStructure
{
private:
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar BANDGAP;   // The bandgap for the material.

  PetscScalar HCI_ECN;   // critical electric field for electron scattering in the insulator
  PetscScalar HCI_ECP;   // critical electric field for hole scattering in the insulator
  PetscScalar HCI_BARLN; // represents barrier lowering effects to electron due to the image field
  PetscScalar HCI_TUNLN; // accounts for the possibility of tunneling to electron
  PetscScalar HCI_BARLP; // represents barrier lowering effects to hole due to the image field
  PetscScalar HCI_TUNLP; // accounts for the possibility of tunneling to hole
  PetscScalar HCI_MFP;   // mean-free-path for electrons in silicon dioxide

  PetscScalar FN_A;      // Coefficient of the pre-exponential term for the Fowler-Nordheim  tunneling model
  PetscScalar FN_B;      // Coefficient of the exponential term for the Fowler-Nordheim tunneling model.

  void   Band_Init()
  {
    AFFINITY = 1.070000e+00*eV;
    BANDGAP  = 9.0*eV;

    HCI_ECN  = 1.650000e+05*V/cm;
    HCI_ECP  = 1.650000e+05*V/cm;
    HCI_BARLN = 2.590000E-04*pow(V*cm, 0.5);
    HCI_TUNLN = 3.000000E-05*pow(V*cm*cm, 1.0/3.0);
    HCI_BARLP = 2.590000E-04*pow(V*cm, 0.5);
    HCI_TUNLP = 3.000000E-05*pow(V*cm*cm, 1.0/3.0);
    HCI_MFP   = 3.200000E-07*cm;

    FN_A      = 6.320000E-07*A/(V*V);
    FN_B      = 2.210000E+08*V/cm;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("AFFINITY",  PARA("AFFINITY", "The electron affinity for the material", "eV", eV, &AFFINITY)) );
    parameter_map.insert(para_item("BANDGAP",   PARA("BANDGAP",  "The bandgap for the material", "eV", eV, &BANDGAP)) );
#endif
  }


public:

  PetscScalar Eg            (const PetscScalar &Tl) const { return BANDGAP;  }

  PetscScalar HCI_Probability_Insulator_n(const PetscScalar &t_ins, const PetscScalar &E_ins) const
  {
    if(E_ins > 0)
    {
      PetscScalar alpha = std::min(t_ins/HCI_MFP, sqrt(HCI_ECN/E_ins) );
      if( alpha > 30.0 ) return 0;
      return exp( -alpha  );
    }
    return exp( -t_ins/HCI_MFP  );
  }

  PetscScalar HCI_Probability_Insulator_p(const PetscScalar &t_ins, const PetscScalar &E_ins) const
  {
    if(E_ins < 0)
    {
      PetscScalar alpha = std::min(t_ins/HCI_MFP, sqrt(HCI_ECN/(-E_ins)) );
      if( alpha > 30.0 ) return 0;
      return exp( -alpha  );
    }
    return exp( -t_ins/HCI_MFP  );
  }

  PetscScalar HCI_Barrier_n(const PetscScalar &affinity_semi, const PetscScalar &,
                            const PetscScalar &t_ins, const PetscScalar &E_ins) const
  {
    if(E_ins > 0)
      return (affinity_semi - AFFINITY) - HCI_BARLN*pow(E_ins, 0.5) - HCI_TUNLN*pow(E_ins, 2.0/3.0) ;
    return (affinity_semi - AFFINITY) - E_ins*t_ins;
  }

  PetscScalar HCI_Barrier_p(const PetscScalar &affinity_semi, const PetscScalar & Eg_semi,
                            const PetscScalar &t_ins, const PetscScalar &E_ins) const
  {
    if(E_ins < 0)
      return (AFFINITY + BANDGAP - affinity_semi -  Eg_semi) - HCI_BARLP*pow(fabs(E_ins), 0.5) - HCI_TUNLP*pow(fabs(E_ins), 2.0/3.0) ;
    return (AFFINITY + BANDGAP - affinity_semi -  Eg_semi) + E_ins*t_ins;
  }


  PetscScalar J_FN_Tunneling(const PetscScalar &E_ins, const PetscScalar &alpha) const
  {
    PetscScalar E = fabs(E_ins) + 1*V/cm;
    if( FN_B/E > 30.0 ) return 0.0;

    if(alpha == 1.0)
      return FN_A*E_ins*E_ins*exp( - FN_B/E );
    else
      return FN_A*E_ins*E_ins/std::pow(1.0-sqrt(1.0-alpha), 2.0)*exp( - FN_B/E*(1-std::pow(1.0-alpha, 1.5)) );
  }

public:
  GSS_SiO2_BandStructure(const PMII_Environment &env):PMII_BandStructure(env)
  {
    Band_Init();
  }
  ~GSS_SiO2_BandStructure()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMII_BandStructure* PMII_SiO2_BandStructure_Default (const PMII_Environment& env)
  {
    return new GSS_SiO2_BandStructure(env);
  }
}
