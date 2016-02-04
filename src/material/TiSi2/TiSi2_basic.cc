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
// Material Type: TiSi2


#include "PMI.h"

class GSS_TiSi2_BasicParameter : public PMIC_BasicParameter
{
private:

  PetscScalar PERMITTI;      // The relative dielectric permittivity.
  PetscScalar PERMEABI;      // The relative megnetic permeability.
  PetscScalar AFFINITY;      // The electron affinity for the material.
  PetscScalar DENSITY;       // Specific mass density for the material.
  PetscScalar IONDENSITY;    // Specific ion (free electron) density for the material.
  PetscScalar CONDUCTANCE;   // Specific conductance for the material.

  void   Basic_Init()
  {
    PERMITTI = 10.0;//approx?
    PERMEABI = 1.0;
    AFFINITY = 0.970000e+00*eV;//use value for Al, should be corrected later.
    DENSITY  = 4.040000e-03*kg*std::pow(cm,-3);//TiSi2
    IONDENSITY  = 1.000000e+23*std::pow(cm,-3);
    CONDUCTANCE = 1.0/(37e-6*V/A*cm); //TiSi2

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("PERMITTI",    PARA("PERMITTI",    "The relative dielectric permittivity", "-", 1.0, &PERMITTI)) );
    parameter_map.insert(para_item("PERMEABI",  PARA("PERMEABI",  "The relative megnetic permeability", "-", 1.0, &PERMEABI)) );
    parameter_map.insert(para_item("AFFINITY", PARA("AFFINITY", "The electron affinity for the material", "eV", eV, &AFFINITY)) );
    parameter_map.insert(para_item("DENSITY", PARA("DENSITY", "Specific mass density for the material", "kg*cm^-3", kg*std::pow(cm,-3), &DENSITY)) );
    parameter_map.insert(para_item("IONDENSITY", PARA("IONDENSITY", "Specific ion density for the material", "cm^-3", std::pow(cm,-3), &IONDENSITY)) );
    parameter_map.insert(para_item("CONDUCTANCE", PARA("CONDUCTANCE", "Specific conductance for the material", "(ohmic*m)^-1", A/V/m, &CONDUCTANCE)) );
#endif
  }
public:

  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;     }
  PetscScalar IonDensity    (const PetscScalar &Tl) const { return IONDENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI;    }
  PetscScalar Permeability  ()                      const { return PERMEABI;    }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY;    }
  PetscScalar Conductance   ()                      const { return CONDUCTANCE; }
  PetscScalar ThermalVn     (const PetscScalar &Tl) const { return 1e6*cm/s;    }
    
  PetscScalar CurrentDensity(const PetscScalar &E, const PetscScalar &Tl) const 
  { return E*CONDUCTANCE; }
  
  AutoDScalar CurrentDensity(const AutoDScalar &E, const AutoDScalar &Tl) const 
  { return E*CONDUCTANCE; }
  
  void G4Material(std::vector<Atom> &atoms, std::vector<double> & fraction) const
  {
    atoms.push_back(Atom("Titanium",  "Ti", 22, 47.9));//Titanium
    atoms.push_back(Atom("Silicon",   "Si", 14, 28.086));//Silicon

    fraction.push_back(1.0);
    fraction.push_back(2.0);
  }

public:
  GSS_TiSi2_BasicParameter(const PMIC_Environment &env):PMIC_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_TiSi2_BasicParameter()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIC_BasicParameter* PMIC_TiSi2_BasicParameter_Default (const PMIC_Environment& env)
  {
    return new GSS_TiSi2_BasicParameter(env);
  }
}
