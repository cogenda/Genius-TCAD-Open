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
// Material Type: 4H-SiC


#include "PMI.h"

class GSS_SiC4H_Thermal : public PMIS_Thermal
{
private:
  PetscScalar A_SP_HEA;	// First parameter for the specific heat model of the material.
  PetscScalar B_SP_HEA;	// Second parameter for the specific heat model of the material.
  PetscScalar C_SP_HEA;	// Third parameter for the specific heat model of the material.
  PetscScalar D_SP_HEA;	// Fourth parameter for the specific heat model of the material.
  PetscScalar F_SP_HEA;	// Fifth parameter for the specific heat model of the material.
  PetscScalar G_SP_HEA;	// Sixth parameter for the specific heat model of the material.	 
  PetscScalar T300;
  void   Thermal_Init()
  {
     A_SP_HEA =  6.900000e+02*J/kg/K;
     B_SP_HEA =  0.000000e+00*J/kg/std::pow(K,2);
     C_SP_HEA =  0.000000e+00*J/kg/std::pow(K,3);
     D_SP_HEA =  0.000000e+00*J/kg*K;
     F_SP_HEA =  0.000000e+00*J/kg/std::pow(K,4);
     G_SP_HEA =  0.000000e+00*J/kg/std::pow(K,5);
     T300     =  300.0*K;
     
#ifdef __CALIBRATE__
	 parameter_map.insert(para_item("A.SP.HEA",    PARA("A.SP.HEA",    "First parameter for the specific heat model of the material", "J/kg/K", J/kg/K, &A_SP_HEA)) );
	 parameter_map.insert(para_item("B.SP.HEA",    PARA("B.SP.HEA",    "Second parameter for the specific heat model of the material", "J/kg/std::pow(K,2)",J/kg/std::pow(K,2) , &B_SP_HEA)) );
	 parameter_map.insert(para_item("C.SP.HEA",    PARA("C.SP.HEA",    "Third parameter for the specific heat model of the material", "J/kg/std::pow(K,3)", J/kg/std::pow(K,3), &C_SP_HEA)) );
	 parameter_map.insert(para_item("D.SP.HEA",    PARA("D.SP.HEA",    "Fourth parameter for the specific heat model of the material", "J/kg*K",J/kg*K , &D_SP_HEA)) );
	 parameter_map.insert(para_item("F.SP.HEA",    PARA("F.SP.HEA",    "Fifth parameter for the specific heat model of the material", "J/kg/std::pow(K,4)",J/kg/std::pow(K,4) , &F_SP_HEA)) );
	 parameter_map.insert(para_item("G.SP.HEA",    PARA("G.SP.HEA",    "Sixth parameter for the specific heat model of the material", "J/kg/std::pow(K,5)",J/kg/std::pow(K,5) , &G_SP_HEA)) );
#endif       
  }
public:
  //---------------------------------------------------------------------------
  // Heat Capacity, no exact values, use paramters for Si.
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const
  {
    return A_SP_HEA + B_SP_HEA*Tl + C_SP_HEA*Tl*Tl + D_SP_HEA/Tl/Tl
           + F_SP_HEA*Tl*Tl*Tl + G_SP_HEA*Tl*Tl*Tl*Tl;
  }
  AutoDScalar HeatCapacity  (const AutoDScalar &Tl) const
  {
    return A_SP_HEA + B_SP_HEA*Tl + C_SP_HEA*Tl*Tl + D_SP_HEA/Tl/Tl
           + F_SP_HEA*Tl*Tl*Tl + G_SP_HEA*Tl*Tl*Tl*Tl;
  }

  //---------------------------------------------------------------------------
  // Heat Conduction, source: Semiconductors on NSM
  PetscScalar HeatConduction(const PetscScalar &Tl) const
  {
    return 3.7*W/cm/K;
  }
  AutoDScalar HeatConduction(const AutoDScalar &Tl) const
  {
    return 3.7*W/cm/K;
  }

// constructor and destructor  
public:     
  GSS_SiC4H_Thermal(const PMIS_Environment &env):PMIS_Thermal(env)
  {
    Thermal_Init();
  }
  ~GSS_SiC4H_Thermal()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Thermal* PMIS_SiC4H_Thermal_Default (const PMIS_Environment& env)
  {
    return new GSS_SiC4H_Thermal(env);
  }
}
