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
// Material Type: GaAs


#include "PMI.h"

class GSS_GaAs_BasicParameter : public PMIS_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of GaAs.
  PetscScalar PERMEABI;  // The relative megnetic permeability of GaAs.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.
  void   Basic_Init()
  {
    PERMITTI = 1.310000e+01;
    PERMEABI = 1.0;
    AFFINITY = 4.070000e+00*eV;
    DENSITY  = 5.317600e-03*kg*std::pow(cm, -3);
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("PERMITTI", PARA("PERMITTI", "Relative permittivity", "-", 1.0, &PERMITTI)) );
    parameter_map.insert(para_item("PERMEABI", PARA("PERMEABI", "Relative permeability", "-", 1.0, &PERMEABI)) );
    parameter_map.insert(para_item("AFFINITY", PARA("AFFINITY", "Electron affinity", "eV", eV, &AFFINITY)) );
    parameter_map.insert(para_item("DENSITY",  PARA("DENSITY",  "DENSITY", "kg cm^-3", 1.0, &DENSITY)) );
#endif
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }

  void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction) const
  {
    atoms.push_back("Gallium");
    atoms.push_back("Arsenic");
    fraction.push_back(1.0);
    fraction.push_back(1.0);
  }

  GSS_GaAs_BasicParameter(const PMIS_Environment &env):PMIS_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_GaAs_BasicParameter()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_BasicParameter* PMIS_GaAs_BasicParameter_Default (const PMIS_Environment& env)
  {
    return new GSS_GaAs_BasicParameter(env);
  }
}
