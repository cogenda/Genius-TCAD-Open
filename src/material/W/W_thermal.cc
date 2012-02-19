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
// Material Type: W


#include "PMI.h"

class GSS_W_Thermal : public PMIC_Thermal
{
public:
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const
  {
    return 24.27*J/(K*(183.84*g));
  }
  PetscScalar HeatConduction(const PetscScalar &Tl) const
  {
    return 173*W/(K*m);
  }
  GSS_W_Thermal(const PMIC_Environment &env):PMIC_Thermal(env)
  {

  }
  ~GSS_W_Thermal()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIC_Thermal* PMIC_W_Thermal_Default (const PMIC_Environment& env)
  {
    return new GSS_W_Thermal(env);
  }
}
