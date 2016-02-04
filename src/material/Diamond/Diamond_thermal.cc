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

class GSS_Diamond_Thermal : public PMIS_Thermal
{
private:

  void   Thermal_Init()
  {

  }
public:
  //---------------------------------------------------------------------------
  // Heat Capacity, no exact values, use paramters for Si.
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const
  {
    return 0.52*J/g/K;
  }
  AutoDScalar HeatCapacity  (const AutoDScalar &Tl) const
  {
    return 0.52*J/g/K;
  }

  //---------------------------------------------------------------------------
  // Heat Conduction, source: Semiconductors on NSM
  PetscScalar HeatConduction(const PetscScalar &Tl) const
  {
    return 10*W/cm/K;
  }
  AutoDScalar HeatConduction(const AutoDScalar &Tl) const
  {
    return 10*W/cm/K;
  }

// constructor and destructor
public:
  GSS_Diamond_Thermal(const PMIS_Environment &env):PMIS_Thermal(env)
  {
    Thermal_Init();
  }
  ~GSS_Diamond_Thermal()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Thermal* PMIS_Diamond_Thermal_Default (const PMIS_Environment& env)
  {
    return new GSS_Diamond_Thermal(env);
  }
}
