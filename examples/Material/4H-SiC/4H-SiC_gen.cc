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
// Material Type: 4H-SiC

#include "PMI.h"

class GSS_SiC4H_Avalanche : public PMIS_Avalanche
{
private:
  PetscScalar    AN_II_LO;
  PetscScalar    BN_II_LO;
  PetscScalar    AN_II_HI;
  PetscScalar    BN_II_HI;
  PetscScalar    E0N_II;
  PetscScalar    EN_OP_PH;

  PetscScalar    AP_II_LO;
  PetscScalar    BP_II_LO;
  PetscScalar    AP_II_HI;
  PetscScalar    BP_II_HI;
  PetscScalar    E0P_II;
  PetscScalar    EP_OP_PH;
  //
  PetscScalar ElecTauw;
  PetscScalar HoleTauw;
  PetscScalar T300    ;

  void 	Avalanche_Init()
  {
    AN_II_LO  =  4.200000E+05/cm;
    BN_II_LO  =  1.670000E+07*V/cm;
    AN_II_HI  =  4.200000E+05/cm;
    BN_II_HI  =  1.670000E+07*V/cm;
    E0N_II    =  4.000000E+05*V/cm;
    EN_OP_PH  =  1.000000E+00*eV;

    AP_II_LO  =  4.200000E+05/cm;
    BP_II_LO  =  1.670000E+07*V/cm;
    AP_II_HI  =  4.200000E+05/cm;
    BP_II_HI  =  1.670000E+07*V/cm;
    E0P_II    =  4.000000E+05*V/cm;
    EP_OP_PH  =  1.000000E+00*eV;
    //
    ElecTauw  = 2e-13*s;
    HoleTauw  = 2e-13*s;
    T300      = 300.0*K;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("AN.II.LO",    PARA("AN.II.LO", "", "cm^-1", 1.0/cm, &AN_II_LO)) );
    parameter_map.insert(para_item("BN.II.LO",    PARA("BN.II.LO", "", "V/cm", V/cm, &BN_II_LO)) );
    parameter_map.insert(para_item("AN.II.HI",    PARA("AN.II.HI", "", "cm^-1", 1.0/cm, &AN_II_HI)) );
    parameter_map.insert(para_item("BN.II.HI",    PARA("BN.II.HI", "", "V/cm", V/cm, &BN_II_HI)) );
    parameter_map.insert(para_item("E0N.II",      PARA("E0N.II", "", "V/cm", V/cm, &E0N_II)) );
    parameter_map.insert(para_item("EN.OP.PH",    PARA("EN.OP.PH", "", "eV", eV, &EN_OP_PH)) );

    parameter_map.insert(para_item("AP.II.LO",    PARA("AP.II.LO", "", "cm^-1", 1.0/cm, &AP_II_LO)) );
    parameter_map.insert(para_item("BP.II.LO",    PARA("BP.II.LO", "", "V/cm", V/cm, &BP_II_LO)) );
    parameter_map.insert(para_item("AP.II.HI",    PARA("AP.II.HI", "", "cm^-1", 1.0/cm, &AP_II_HI)) );
    parameter_map.insert(para_item("BP.II.HI",    PARA("BP.II.HI", "", "V/cm", V/cm, &BP_II_HI)) );
    parameter_map.insert(para_item("E0P.II",      PARA("E0P.II", "", "V/cm", V/cm, &E0P_II)) );
    parameter_map.insert(para_item("EP.OP.PH",    PARA("EP.OP.PH", "", "eV", eV, &EP_OP_PH)) );
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
      PetscScalar gamma = tanh(EN_OP_PH/(2*kb*T300))/tanh(EN_OP_PH/(2*kb*Tl));
      if (Ep < E0N_II)
        return gamma * AN_II_LO * exp (-BN_II_LO * gamma /Ep);
      else
        return gamma * AN_II_HI * exp (-BN_II_HI * gamma /Ep);
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
      AutoDScalar gamma = tanh(EN_OP_PH/(2*kb*T300))/tanh(EN_OP_PH/(2*kb*Tl));
      if (Ep < E0N_II)
        return gamma * AN_II_LO * exp (-BN_II_LO * gamma /Ep);
      else
        return gamma * AN_II_HI * exp (-BN_II_HI * gamma /Ep);
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
      PetscScalar gamma = tanh(EP_OP_PH/(2*kb*T300))/tanh(EP_OP_PH/(2*kb*Tl));
      if (Ep < E0P_II)
        return gamma * AP_II_LO * exp (-BP_II_LO * gamma /Ep);
      else
        return gamma * AP_II_HI * exp (-BP_II_HI * gamma /Ep);
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
      AutoDScalar gamma = tanh(EP_OP_PH/(2*kb*T300))/tanh(EP_OP_PH/(2*kb*Tl));
      if (Ep < E0P_II)
        return gamma * AP_II_LO * exp (-BP_II_LO * gamma /Ep);
      else
        return gamma * AP_II_HI * exp (-BP_II_HI * gamma /Ep);
    }
  }



  //---------------------------------------------------------------------------
  // Electron Impact Ionization rate for EBM
  PetscScalar ElecGenRateEBM (const PetscScalar &Tn,const PetscScalar &Tl,const PetscScalar &Eg) const
  {
    return 0;
  }
  AutoDScalar ElecGenRateEBM (const AutoDScalar &Tn,const AutoDScalar &Tl,const AutoDScalar &Eg) const
  {
    return 0;
  }

  //---------------------------------------------------------------------------
  // Hole Impact Ionization rate for EBM
  PetscScalar HoleGenRateEBM (const PetscScalar &Tp,const PetscScalar &Tl,const PetscScalar &Eg) const
  {
    return 0;
  }
  AutoDScalar HoleGenRateEBM (const AutoDScalar &Tp,const AutoDScalar &Tl,const AutoDScalar &Eg) const
  {
    return 0;
  }

//----------------------------------------------------------------
// constructor and destructor
public:

  GSS_SiC4H_Avalanche(const PMIS_Environment &env):PMIS_Avalanche(env)
  {
    Avalanche_Init();
  }
  ~GSS_SiC4H_Avalanche()
  {
  }

}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Avalanche* PMIS_SiC4H_Avalanche_vanOverstraetendeMan (const PMIS_Environment& env)
  {
    return new GSS_SiC4H_Avalanche(env);
  }
  DLL_EXPORT_DECLARE  PMIS_Avalanche* PMIS_SiC4H_Avalanche_Default (const PMIS_Environment& env)
  {
    return new GSS_SiC4H_Avalanche(env);
  }
}
