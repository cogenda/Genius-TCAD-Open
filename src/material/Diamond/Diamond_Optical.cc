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
// Material Type: Diamond.


#include "PMI.h"



class GSS_Diamond_Optical : public PMIS_Optical
{
private:

  void _init_default_wave_table()
  {
    if(!_wave_table.empty()) _wave_table.clear();
  }

public:
  std::complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Tl, PetscScalar Eg=0) const
  {
    std::complex<PetscScalar> n(2.4,0.0);
    return n;
  }


public:

  // constructions
  GSS_Diamond_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
    _init_default_wave_table();
  }

  ~GSS_Diamond_Optical()
  {}

};


extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Optical*  PMIS_Diamond_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_Diamond_Optical(env);
  }
}
