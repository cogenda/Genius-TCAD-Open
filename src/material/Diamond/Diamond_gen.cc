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
// Material Type: Diamond

#include "PMI.h"

class GSS_Diamond_Avalanche_Default : public PMIS_Avalanche
{
private:
  PetscScalar   N_IONIZA;
  PetscScalar   EN_II;
  PetscScalar   P_IONIZA;
  PetscScalar   EP_II;

  PetscScalar cut_low;
  PetscScalar cut_end;

  void 	Avalanche_Init()
  {
    N_IONIZA  =  1.935E8/cm;
    EN_II     =  7.749E6*V/cm;

    P_IONIZA  =  1.935E8/cm;
    EP_II     =  7.749E6*V/cm;

    cut_low   = 0.1;
    cut_end   = 0.2;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("N.IONIZA",PARA("N.IONIZA",    "The constant term in the multiplicative prefactor of the electron ionization coefficient", "/cm", 1.0/cm, &N_IONIZA)) );
    parameter_map.insert(para_item("EN_II",   PARA("EN_II",  "The critical electrical field to the local electric field", "V/cm", V/cm, &EN_II)) );
    parameter_map.insert(para_item("P.IONIZA",PARA("P.IONIZA",  "The constant term in the multiplicative prefactor of the hole ionization coefficient", "/cm",1.0/cm , &P_IONIZA)) );
    parameter_map.insert(para_item("EP.II",   PARA("EP.II",  "The exponent of the ratio of the critical electrical field to the local electric field", "-",1.0 , &EP_II)) );
    parameter_map.insert(para_item("CutLow",  PARA("CutLow",  "Disable impact ionization when E field lower than CutLow*Ecrit", "-", 1.0, &cut_low)) );
    parameter_map.insert(para_item("CutEnd",  PARA("CutEnd",  "Damp impact ionization when E field between CutLow*Ecrit and CutEnd*Ecrit", "-", 1.0, &cut_end)) );
#endif
  }


public:
  //---------------------------------------------------------------------------
  // Electron Impact Ionization rate for DDM
  PetscScalar ElecGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const
  {
    if (Ep < cut_low*EN_II)
    {
      return 0;
    }
    else if (Ep < cut_end*EN_II)
    {
      PetscScalar smooth_rate = 1/(cut_end-cut_low)*(Ep/EN_II-cut_low);
      return smooth_rate*N_IONIZA*exp(-EN_II/Ep);
    }
    else
      return N_IONIZA*exp(-EN_II/Ep);
  }

  AutoDScalar ElecGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const
  {
    if (Ep < cut_low*EN_II)
    {
      return 0;
    }
    else if (Ep < cut_end*EN_II)
    {
      AutoDScalar smooth_rate = 1/(cut_end-cut_low)*(Ep/EN_II-cut_low);
      return smooth_rate*N_IONIZA*exp(-EN_II/Ep);
    }
    else
      return N_IONIZA*exp(-EN_II/Ep);
  }

  //---------------------------------------------------------------------------
  // Hole Impact Ionization rate for DDM
  PetscScalar HoleGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const
  {
    if (Ep < cut_low*EP_II)
    {
      return 0;
    }
    else if (Ep < cut_end*EP_II)
    {
      PetscScalar smooth_rate = 1/(cut_end-cut_low)*(Ep/EP_II-cut_low);
      return smooth_rate*P_IONIZA*exp(-EP_II/Ep);
    }
    else
      return P_IONIZA*exp(-EP_II/Ep);
  }
  AutoDScalar HoleGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const
  {
    if (Ep < cut_low*EP_II)
    {
      return 0;
    }
    else if (Ep < cut_end*EP_II)
    {
      AutoDScalar smooth_rate = 1/(cut_end-cut_low)*(Ep/EP_II-cut_low);
      return smooth_rate*P_IONIZA*exp(-EP_II/Ep);
    }
    else
      return P_IONIZA*exp(-EP_II/Ep);
  }



  //---------------------------------------------------------------------------
  // Electron Impact Ionization rate for EBM
  PetscScalar ElecGenRateEBM (const PetscScalar &Tn,const PetscScalar &Tl,const PetscScalar &Eg) const
  {
    return 0.0;
  }
  AutoDScalar ElecGenRateEBM (const AutoDScalar &Tn,const AutoDScalar &Tl,const AutoDScalar &Eg) const
  {
    return 0.0;
  }

  //---------------------------------------------------------------------------
  // Hole Impact Ionization rate for EBM
  PetscScalar HoleGenRateEBM (const PetscScalar &Tp,const PetscScalar &Tl,const PetscScalar &Eg) const
  {
    return 0.0;
  }
  AutoDScalar HoleGenRateEBM (const AutoDScalar &Tp,const AutoDScalar &Tl,const AutoDScalar &Eg) const
  {
    return 0.0;
  }

//----------------------------------------------------------------
// constructor and destructor
public:

  GSS_Diamond_Avalanche_Default(const PMIS_Environment &env):PMIS_Avalanche(env)
  {
    Avalanche_Init();
  }
  ~GSS_Diamond_Avalanche_Default()
  {
  }

}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Avalanche* PMIS_Diamond_Avalanche_Selberherr (const PMIS_Environment& env)
  {
    return new GSS_Diamond_Avalanche_Default(env);
  }
  DLL_EXPORT_DECLARE  PMIS_Avalanche* PMIS_Diamond_Avalanche_Default (const PMIS_Environment& env)
  {
    return new GSS_Diamond_Avalanche_Default(env);
  }
}
