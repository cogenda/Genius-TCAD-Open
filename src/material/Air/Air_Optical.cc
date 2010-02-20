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
// Material Type: Air


#include "PMI.h"


class GSS_Air_Optical : public PMII_Optical
{
public:
  std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=9.0, PetscScalar Tl=1.0) const 
  {
      return std::complex<PetscScalar> (1.0,0.0);
  }                                            
  
  // constructions
public:
  GSS_Air_Optical(const PMII_Environment &env):PMII_Optical(env) {}

  ~GSS_Air_Optical(){}
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMII_Optical*  PMII_Air_Optical_Default (const PMII_Environment& env)
  {
    return new GSS_Air_Optical(env);
  }
}
