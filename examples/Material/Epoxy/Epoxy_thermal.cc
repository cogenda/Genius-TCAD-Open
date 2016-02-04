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
// Material Type: Epoxy


#include "PMI.h"

class GSS_Epoxy_Thermal : public PMII_Thermal
{
private:

  void   Thermal_Init()
  {

  }

public:
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const 
  {
    return 1.00*J/(g*K);
  }
  PetscScalar HeatConduction(const PetscScalar &Tl) const 
  {
    //S. Okamoto and H. Ishida, “Thermal Conductivity of PET/(LDPE/AI) Composites Determined by MDSC”,
    // Macromolecules, 34(2001), p. 7392.
    return 0.4*W/(m*K);
  }

  GSS_Epoxy_Thermal(const PMII_Environment &env):PMII_Thermal(env)
  {
    Thermal_Init();
  }
  ~GSS_Epoxy_Thermal()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMII_Thermal* PMII_Epoxy_Thermal_Default (const PMII_Environment& env)
  {
    return new GSS_Epoxy_Thermal(env);
  }
}
