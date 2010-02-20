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
/*  Last update: July 25, 2007                                               */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: Nitride


#include "PMI.h"

class GSS_Nitride_BasicParameter : public PMII_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of Nitride.
  PetscScalar PERMEABI;  // The relative megnetic permeability of Nitride.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar BANDGAP;   // The banggap for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.
  void   Basic_Init()
  {
    PERMITTI = 7.500000E+00;
    PERMEABI = 1.0;
    AFFINITY = 2.5*eV;
    BANDGAP  = 4.7*eV;
    DENSITY  = 3.440000E-03*kg*std::pow(cm,-3);

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("PERMITTI",    PARA("PERMITTI",    "The relative dielectric permittivity", "-", 1.0, &PERMITTI)) );
    parameter_map.insert(para_item("PERMEABI",  PARA("PERMEABI",  "The relative megnetic permeability ", "-", 1.0, &PERMEABI)) );
    parameter_map.insert(para_item("AFFINITY", PARA("AFFINITY", "The electron affinity for the material", "eV", eV, &AFFINITY)) );
    parameter_map.insert(para_item("DENSITY", PARA("DENSITY", "Specific mass density for the material", "kg*cm^-3", kg*std::pow(cm,-3), &DENSITY)) );
#endif
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }
  PetscScalar Eg            (const PetscScalar &Tl) const { return BANDGAP;  }

  void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction) const
  {
    atoms.push_back("Silicon");
    atoms.push_back("Nitrogen");

    fraction.push_back(3.0);
    fraction.push_back(4.0);
  }

  GSS_Nitride_BasicParameter(const PMII_Environment &env):PMII_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_Nitride_BasicParameter()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMII_BasicParameter* PMII_Nitride_BasicParameter_Default (const PMII_Environment& env)
  {
    return new GSS_Nitride_BasicParameter(env);
  }
}
