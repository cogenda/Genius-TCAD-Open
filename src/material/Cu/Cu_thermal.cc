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
/*  Last update: May 27, 2007                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: Cu


#include "PMI.h"

class GSS_Cu_Thermal : public PMIC_Thermal
{
public:
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const 
  {
    return 3.42*J/(K*std::pow(cm,3));
  }
  PetscScalar HeatConduction(const PetscScalar &Tl) const 
  {
    return 3.85*W/(K*cm);
  }
  GSS_Cu_Thermal(const PMIC_Environment &env):PMIC_Thermal(env)
  {
    
  }
  ~GSS_Cu_Thermal()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIC_Thermal* PMIC_Cu_Thermal_Default (const PMIC_Environment& env)
  {
    return new GSS_Cu_Thermal(env);
  }
}
