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


#include "PMI.h"

class GSS_Si_Thermal : public PMIS_Thermal
{
private:
  PetscScalar A_SP_HEA;	// First parameter for the specific heat model of the material.
  PetscScalar B_SP_HEA;	// Second parameter for the specific heat model of the material.
  PetscScalar C_SP_HEA;	// Third parameter for the specific heat model of the material.
  PetscScalar D_SP_HEA;	// Fourth parameter for the specific heat model of the material.
  PetscScalar F_SP_HEA;	// Fifth parameter for the specific heat model of the material.
  PetscScalar G_SP_HEA;	// Sixth parameter for the specific heat model of the material.
  PetscScalar H_SP_HEA;	// Seventh parameter for the specific heat model of the material.
  PetscScalar A_TH_CON;	// First parameter for the thermal conductivity model of the material.
  PetscScalar B_TH_CON;	// Second parameter for the thermal conductivity model of the material.
  PetscScalar C_TH_CON;	// Third parameter for the thermal conductivity model of the material.
  PetscScalar E_TH_CON;	// Fifth parameter for the thermal conductivity model of the material.
  PetscScalar D_TH_CON;	// Fourth parameter for the thermal conductivity model of the material.
  PetscScalar T300;
  void   Thermal_Init()
  {
    //f(x)=a+b*x+c*x^2+d*x^3+e*x^4+f*x^5
    //a=-40.783173616135490*J/kg/K;
    //b=4.326055304269036*J/kg/pow(K,2);
    //c=-0.008164731922246*J/kg/pow(K,3);
    //d=7.631085506760700e-06*J/kg/pow(K,4);
    //e=-3.394906879769543e-09*J/kg/pow(K,5);
    //f=5.692999601098961e-13*J/kg/pow(K,6);

    A_SP_HEA = -4.783173616135490e+01*J/kg/K;
    B_SP_HEA =  4.326055304269036*J/kg/std::pow(K,2);
    C_SP_HEA = -0.008164731922246*J/kg/std::pow(K,3);
    D_SP_HEA =  0.000000e+00*J/kg*K;
    F_SP_HEA =  7.631085506760700e-06*J/kg/std::pow(K,4);
    G_SP_HEA = -3.394906879769543e-09*J/kg/std::pow(K,5);
    H_SP_HEA =  5.692999601098961e-13*J/kg/std::pow(K,6);
    A_TH_CON =  3.000000e-02*cm*K/J*s;
    B_TH_CON =  1.560000e-03*cm/J*s;
    C_TH_CON =  1.650000e-06*cm/J*s/K;
    E_TH_CON =  0.000000e+00;
    D_TH_CON =  0.000000e+00*cm/J*s*std::pow(K,1-E_TH_CON);
    T300     =  300.0*K;

#ifdef __CALIBRATE__
	 parameter_map.insert(para_item("A.SP.HEA",    PARA("A.SP.HEA",    "First parameter for the specific heat model of the material", "J/kg/K", J/kg/K, &A_SP_HEA)) );
	 parameter_map.insert(para_item("B.SP.HEA",    PARA("B.SP.HEA",    "Second parameter for the specific heat model of the material", "J/kg/std::pow(K,2)",J/kg/std::pow(K,2) , &B_SP_HEA)) );
	 parameter_map.insert(para_item("C.SP.HEA",    PARA("C.SP.HEA",    "Third parameter for the specific heat model of the material", "J/kg/std::pow(K,3)", J/kg/std::pow(K,3), &C_SP_HEA)) );
	 parameter_map.insert(para_item("D.SP.HEA",    PARA("D.SP.HEA",    "Fourth parameter for the specific heat model of the material", "J/kg*K",J/kg*K , &D_SP_HEA)) );
	 parameter_map.insert(para_item("F.SP.HEA",    PARA("F.SP.HEA",    "Fifth parameter for the specific heat model of the material", "J/kg/std::pow(K,4)",J/kg/std::pow(K,4) , &F_SP_HEA)) );
	 parameter_map.insert(para_item("G.SP.HEA",    PARA("G.SP.HEA",    "Sixth parameter for the specific heat model of the material", "J/kg/std::pow(K,5)",J/kg/std::pow(K,5) , &G_SP_HEA)) );
         parameter_map.insert(para_item("H.SP.HEA",    PARA("H.SP.HEA",    "Seventh parameter for the specific heat model of the material", "J/kg/std::pow(K,6)",J/kg/std::pow(K,6) , &H_SP_HEA)) );
	 parameter_map.insert(para_item("A.TH.CON",    PARA("A.TH.CON",    "First parameter for the thermal conductivity model of the material", "cm*K/J*s", cm*K/J*s, &A_TH_CON)) );
	 parameter_map.insert(para_item("B.TH.CON",    PARA("B.TH.CON",    "Second parameter for the thermal conductivity model of the material", "cm/J*s", cm/J*s, &B_TH_CON)) );
	 parameter_map.insert(para_item("C.TH.CON",    PARA("C.TH.CON",    "Third parameter for the thermal conductivity model of the material", "cm/J*s/K",cm/J*s/K , &C_TH_CON)) );
	 parameter_map.insert(para_item("E.TH.CON",    PARA("E.TH.CON",    "Fifth parameter for the thermal conductivity model of the material", "-", 1.0, &E_TH_CON)) );
	 parameter_map.insert(para_item("D.TH.CON",    PARA("D.TH.CON",    "Fourth parameter for the thermal conductivity model of the material", "cm/J*s*std::pow(K,1-E.TH.CON)", cm/J*s*std::pow(K,1-E_TH_CON), &D_TH_CON)) );
#endif        	
  }
public:
  //---------------------------------------------------------------------------
  // Heat Capacity, cover the T range of [20,1500]
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const
  {
    if(Tl > 20*K)
      return A_SP_HEA + B_SP_HEA*Tl + C_SP_HEA*Tl*Tl + F_SP_HEA*Tl*Tl*Tl + G_SP_HEA*Tl*Tl*Tl*Tl + H_SP_HEA*Tl*Tl*Tl*Tl*Tl;
    else
    {
      PetscScalar T  = 20*K;
      PetscScalar Cp = A_SP_HEA + B_SP_HEA*T + C_SP_HEA*Tl*T + F_SP_HEA*T*T*T + G_SP_HEA*T*T*T*T + H_SP_HEA*T*T*T*T*T;
      return Cp*Tl/T;
    }
    return 0.0;
  }
  AutoDScalar HeatCapacity  (const AutoDScalar &Tl) const
  {
    if(Tl > 20*K)
      return A_SP_HEA + B_SP_HEA*Tl + C_SP_HEA*Tl*Tl + F_SP_HEA*Tl*Tl*Tl + G_SP_HEA*Tl*Tl*Tl*Tl + H_SP_HEA*Tl*Tl*Tl*Tl*Tl;
    else
    {
      PetscScalar T  = 20*K;
      PetscScalar Cp = A_SP_HEA + B_SP_HEA*T + C_SP_HEA*T*T + F_SP_HEA*T*T*T + G_SP_HEA*T*T*T*T + H_SP_HEA*T*T*T*T*T;
      return Cp*Tl/T;
    }
    return 0.0;
  }

  //---------------------------------------------------------------------------
  // Heat Conduction
  PetscScalar HeatConduction(const PetscScalar &Tl) const
  {
    return 1.0/(A_TH_CON + B_TH_CON*Tl + C_TH_CON*Tl*Tl);
  }
  AutoDScalar HeatConduction(const AutoDScalar &Tl) const
  {
    return 1.0/(A_TH_CON + B_TH_CON*Tl + C_TH_CON*Tl*Tl);
  }

// constructor and destructor  
public:     
  GSS_Si_Thermal(const PMIS_Environment &env):PMIS_Thermal(env)
  {
    PMI_Info = "This is the Default thermal model of Silicon"; 
    Thermal_Init();
  }
  ~GSS_Si_Thermal()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Thermal* PMIS_Si_Thermal_Default (const PMIS_Environment& env)
  {
    return new GSS_Si_Thermal(env);
  }
}
