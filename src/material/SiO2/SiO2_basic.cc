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
// Material Type: SiO2


#include "PMI.h"

class GSS_SiO2_BasicParameter : public PMII_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of SiO2.
  PetscScalar PERMEABI;  // The relative megnetic permeability of SiO2.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.

  void   Basic_Init()
  {
    PERMITTI = 3.900000e+00;
    PERMEABI = 1.0;
    AFFINITY = 1.070000e+00*eV;
    DENSITY  = 2.260000e-03*kg*std::pow(cm,-3);

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("PERMITTI",  PARA("PERMITTI", "The relative dielectric permittivity of SiO2", "-", 1.0, &PERMITTI)) );
    parameter_map.insert(para_item("PERMEABI",  PARA("PERMEABI", "The relative megnetic permeability of SiO2", "-", 1.0, &PERMEABI)) );
    parameter_map.insert(para_item("AFFINITY",  PARA("AFFINITY", "The electron affinity for the material", "eV", eV, &AFFINITY)) );
    parameter_map.insert(para_item("DENSITY",   PARA("DENSITY",  "Specific mass density for the material", "kg*cm^-3", kg*std::pow(cm,-3), &DENSITY)) );
#endif
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }

  PetscScalar Conductance   (const PetscScalar &Tl, const PetscScalar &E) const
  { return 1.0/(1e+15*V/A*m); }
  AutoDScalar Conductance   (const AutoDScalar &Tl, const AutoDScalar &E) const
  { return 1.0/(1e+15*V/A*m); }

  PetscScalar RadConductance(const PetscScalar &DRate)const
  {
    return 0.0;
  }

  PetscScalar Mobility      (const PetscScalar &Tl) const
  { return 5e-11*cm*cm/V/s; }

  void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction) const
  {
    atoms.push_back("Si");//Silicon
    atoms.push_back("O");//Oxygen

    fraction.push_back(1.0);
    fraction.push_back(2.0);
  }

public:
  GSS_SiO2_BasicParameter(const PMII_Environment &env):PMII_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_SiO2_BasicParameter()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMII_BasicParameter* PMII_SiO2_BasicParameter_Default (const PMII_Environment& env)
  {
    return new GSS_SiO2_BasicParameter(env);
  }
}
