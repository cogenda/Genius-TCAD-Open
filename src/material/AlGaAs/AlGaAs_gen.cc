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
// Material Type: AlGaAs

#include "PMI.h"

class GSS_AlGaAs_Avalanche : public PMIS_Avalanche
{
private:
  PetscScalar N_IONIZA; // The constant term in the multiplicative prefactor of the electron ionization coefficient.
  PetscScalar EXN_II;   // The exponent of the ratio of the critical electrical field to the local electric field.
  PetscScalar P_IONIZA; // The constant term in the multiplicative prefactor of the hole ionization coefficient.
  PetscScalar EXP_II;   // The exponent of the ratio of the critical electrical field to the local electric field.
  PetscScalar N_ION_1;  // The coefficient multiplying T in the multiplicative prefactor of the electron ionization coefficient.
  PetscScalar N_ION_2;  // The coefficient multiplying T^2 in the multiplicative prefactor of the electron ionization coefficient.
  PetscScalar P_ION_1;  // The coefficient multiplying T in the multiplicative prefactor of the hole ionization coefficient.
  PetscScalar P_ION_2;  // The coefficient multiplying T^2 in the multiplicative prefactor of the hole ionization coefficient.
  //// Impact Ionization Model Depending on Lattice Temperature.
  PetscScalar LAN300;   // Energy free path for electrons at 300 K, used for the impact ionization model depending on lattice temperature.
  PetscScalar LAP300;   // Energy free path for holes at 300 K, used for the impact ionization model depending on lattice temperature.
  PetscScalar OP_PH_EN; // Mean optical phonon energy used for the impact ionization model depending on lattice temperature.
  //
  PetscScalar ElecTauw;
  PetscScalar HoleTauw;
  PetscScalar T300    ;

  void 	Avalanche_Init()
  {
    N_IONIZA  =  2.994000e+05/cm;
    EXN_II    =  1.600000e+00;
    P_IONIZA  =  2.215000e+05/cm;
    EXP_II    =  1.750000e+00;
    N_ION_1   =  0.000000e+00/cm/K;
    N_ION_2   =  0.000000e+00/cm/K/K;
    P_ION_1   =  0.000000e+00/cm/K;
    P_ION_2   =  0.000000e+00/cm/K/K;
    LAN300    =  3.527240e-06*cm;
    LAP300    =  3.676490e-06*cm;
    OP_PH_EN  =  3.500000e-02*eV;

    ElecTauw  = 1e-12*s;
    HoleTauw  = 1e-12*s;
    T300      = 300.0*K;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("N.IONIZA",    PARA("N.IONIZA",    "The constant term in the multiplicative prefactor of the electron ionization coefficient", "/cm", 1.0/cm, &N_IONIZA)) );
    parameter_map.insert(para_item("EXN.II",  PARA("EXN.II",  "The exponent of the ratio of the critical electrical field to the local electric field", "-", 1.0, &EXN_II)) );
    parameter_map.insert(para_item("P.IONIZA",  PARA("P.IONIZA",  "The constant term in the multiplicative prefactor of the hole ionization coefficient", "/cm",1.0/cm , &P_IONIZA)) );
    parameter_map.insert(para_item("EXP.II",  PARA("EXP.II",  "The exponent of the ratio of the critical electrical field to the local electric field", "-",1.0 , &EXP_II)) );
    parameter_map.insert(para_item("N.ION.1",  PARA("N.ION.1",  "The coefficient multiplying T in the multiplicative prefactor of the electron ionization coefficient", "/cm/K",1.0/cm/K , &N_ION_1)) );
    parameter_map.insert(para_item("N.ION.2",  PARA("N.ION.2",  "The coefficient multiplying T^2 in the multiplicative prefactor of the electron ionization coefficient", "/cm/K/K",1.0/cm/K/K , &N_ION_2)) );
    parameter_map.insert(para_item("P.ION.1",  PARA("P.ION.1",  "The coefficient multiplying T in the multiplicative prefactor of the hole ionization coefficient", "/cm/K", 1.0/cm/K, &P_ION_1)) );
    parameter_map.insert(para_item("P.ION.2",  PARA("P.ION.2",  "The coefficient multiplying T^2 in the multiplicative prefactor of the hole ionization coefficient", "/cm/K/K",1.0/cm/K/K , &P_ION_2)) );
    parameter_map.insert(para_item("LAN300",  PARA("LAN300",  "Energy free path for electrons at 300 K, used for the impact ionization model depending on lattice temperature", "cm", cm, &LAN300)) );
    parameter_map.insert(para_item("LAP300",  PARA("LAP300",  "Energy free path for holes at 300 K, used for the impact ionization model depending on lattice temperature", "cm", cm, &LAP300)) );
    parameter_map.insert(para_item("OP.PH.EN",  PARA("OP.PH.EN",  "Mean optical phonon energy used for the impact ionization model depending on lattice temperature", "eV", eV, &OP_PH_EN)) );
    parameter_map.insert(para_item("ElecTauw",  PARA("ElecTauw",  "", "s", s, &ElecTauw)) );
    parameter_map.insert(para_item("HoleTauw",  PARA("HoleTauw",  "", "s", s, &HoleTauw)) );
#endif
  }
public:
  //---------------------------------------------------------------------------
  // Electron Impact Ionization rate for DDM
  PetscScalar ElecGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const
  {
    if (Ep < 1e3*V/cm)
    {
      return 0;
    }
    else
    {
      PetscScalar alpha = N_IONIZA + N_ION_1*Tl + N_ION_2*Tl*Tl;
      PetscScalar L=LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      return alpha*exp(-std::pow(Eg/(e*L)/Ep,EXN_II));
    }
  }
  AutoDScalar ElecGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const
  {
    if (Ep < 1e3*V/cm)
    {
      return 0;
    }
    else
    {
      AutoDScalar alpha = N_IONIZA + N_ION_1*Tl + N_ION_2*Tl*Tl;
      AutoDScalar L = LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      AutoDScalar Ecrit = Eg/(e*L);
      return alpha*exp(-adtl::pow(Ecrit/Ep,EXN_II));
    }
  }

  //---------------------------------------------------------------------------
  // Hole Impact Ionization rate for DDM
  PetscScalar HoleGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const
  {
    if (Ep < 1e3*V/cm)
    {
      return 0;
    }
    else
    {
      PetscScalar alpha = P_IONIZA+P_ION_1*Tl+P_ION_2*Tl*Tl;
      PetscScalar L = LAP300*tanh(OP_PH_EN/(2*kb*Tl));
      return alpha*exp(-std::pow(Eg/(e*L)/Ep,EXP_II));
    }
  }
  AutoDScalar HoleGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const
  {
    if (Ep < 1e3*V/cm)
    {
      return 0;
    }
    else
    {
      AutoDScalar alpha = P_IONIZA+P_ION_1*Tl+P_ION_2*Tl*Tl;
      AutoDScalar L = LAP300*tanh(OP_PH_EN/(2*kb*Tl));
      AutoDScalar Ecrit = Eg/(e*L);
      return alpha*exp(-adtl::pow(Ecrit/Ep,EXP_II));
    }
  }


  //---------------------------------------------------------------------------
  // Electron Impact Ionization rate for EBM
  PetscScalar ElecGenRateEBM (const PetscScalar &Tn,const PetscScalar &Tl,const PetscScalar &Eg) const
  {
    if (fabs(Tn - Tl)<1*K)
    {
      return 0;
    }
    else
    {
      PetscScalar vsat = 11.3E6*cm/s - 1.2E4*cm/s/K*Tl;
      PetscScalar L = LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      PetscScalar Ecrit = Eg/(e*L);
      PetscScalar uc   = 2*vsat*ElecTauw/3*Ecrit;
      PetscScalar ut   = kb/e*(Tn-Tl);
      return N_IONIZA/e*exp(-std::pow(uc/ut,EXN_II));
    }
  }
  AutoDScalar ElecGenRateEBM (const AutoDScalar &Tn,const AutoDScalar &Tl,const AutoDScalar &Eg) const
  {
    if (fabs(Tn - Tl)<1*K)
    {
      return 0;
    }
    else
    {
      AutoDScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
      AutoDScalar L = LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      AutoDScalar Ecrit = Eg/(e*L);
      AutoDScalar uc   = 2*vsat*ElecTauw/3*Ecrit;
      AutoDScalar ut   = kb/e*(Tn-Tl);
      return N_IONIZA/e*exp(-adtl::pow(uc/ut,EXN_II));
    }
  }

  //---------------------------------------------------------------------------
  // Hole Impact Ionization rate for EBM
  PetscScalar HoleGenRateEBM (const PetscScalar &Tp,const PetscScalar &Tl,const PetscScalar &Eg) const
  {
    if (fabs(Tp - Tl)<1*K)
    {
      return 0;
    }
    else
    {
      PetscScalar vsat = 11.3E6*cm/s - 1.2E4*cm/s/K*Tl;
      PetscScalar L = LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      PetscScalar Ecrit = Eg/(e*L);
      PetscScalar uc   = 2*vsat*HoleTauw/3*Ecrit;
      PetscScalar ut   = kb/e*(Tp-Tl);
      return P_IONIZA/e*exp(-std::pow(uc/ut,EXP_II));
    }
  }
  AutoDScalar HoleGenRateEBM (const AutoDScalar &Tp,const AutoDScalar &Tl,const AutoDScalar &Eg) const
  {
    if (fabs(Tp - Tl)<1*K)
    {
      return 0;
    }
    else
    {
      AutoDScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
      AutoDScalar L = LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      AutoDScalar Ecrit = Eg/(e*L);
      AutoDScalar uc   = 2*vsat*HoleTauw/3*Ecrit;
      AutoDScalar ut   = kb/e*(Tp-Tl);
      return P_IONIZA/e*exp(-adtl::pow(uc/ut,EXP_II));
    }
  }

//----------------------------------------------------------------
// constructor and destructor
public:
  GSS_AlGaAs_Avalanche(const PMIS_Environment &env):PMIS_Avalanche(env)
  {
    Avalanche_Init();
  }
  ~GSS_AlGaAs_Avalanche()
  {
  }

}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Avalanche* PMIS_AlGaAs_Avalanche_Selberherr (const PMIS_Environment& env)
  {
    return new GSS_AlGaAs_Avalanche(env);
  }
  DLL_EXPORT_DECLARE  PMIS_Avalanche* PMIS_AlGaAs_Avalanche_Default (const PMIS_Environment& env)
  {
    return new GSS_AlGaAs_Avalanche(env);
  }
}
