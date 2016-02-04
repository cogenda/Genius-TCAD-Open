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

  PetscScalar rad;
  PetscScalar RadGen;   // radiative generated elec-hole pairs per rad
  PetscScalar RadYield_E;
  PetscScalar RadYield_P;
  PetscScalar ElecMob;  // Electron mobility
  PetscScalar HoleMob;  // Hole mobility
  PetscScalar HIonMob;  // H+ mobility

  void   Basic_Init()
  {
    PERMITTI = 3.900000e+00;
    PERMEABI = 1.0;
    AFFINITY = 1.070000e+00*eV;
    DENSITY  = 2.260000e-03*kg*std::pow(cm,-3);

    rad      = 0.01*J/kg;
    RadGen   = 7.6e12/rad*std::pow(cm,-3);
    
    // there are two set of parameters
    //Fernández-Martínez P, Cortés I, Hidalgo S, et al. Simulation of total ionising dose in MOS capacitors[C]
    //Electron Devices (CDE), 2011 Spanish Conference on. IEEE, 2011: 1-4.
    //   for Co60 gamma      for X-ray
    //    0.55e6*V/cm       1.35e6*V/cm
    //      0.7                0.9
    
    //defaut use parameter set for Co60
    RadYield_E =  0.55e6*V/cm;       //1.35e6*V/cm;
    RadYield_P =  0.7;               //0.9;
    
    
    ElecMob  = 20*cm*cm/V/s;
    HoleMob  = 1e-5*cm*cm/V/s;
    HIonMob  = 1.14e-11*cm*cm/V/s;
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("PERMITTI",  PARA("PERMITTI", "The relative dielectric permittivity of SiO2", "-", 1.0, &PERMITTI)) );
    parameter_map.insert(para_item("PERMEABI",  PARA("PERMEABI", "The relative megnetic permeability of SiO2", "-", 1.0, &PERMEABI)) );
    parameter_map.insert(para_item("AFFINITY",  PARA("AFFINITY", "The electron affinity for the material", "eV", eV, &AFFINITY)) );
    parameter_map.insert(para_item("DENSITY",   PARA("DENSITY",  "Specific mass density for the material", "kg*cm^-3", kg*std::pow(cm,-3), &DENSITY)) );

    parameter_map.insert(para_item("RadGen",    PARA("RadGen",    "Radiative generated elec-hole pairs per rad", "1/rad/cm^-3", 1.0/rad/std::pow(cm,-3), &RadGen)) );
    parameter_map.insert(para_item("RadYield.E",PARA("RadYield.E","E parameter for eh escaped yield", "V/cm", V/cm, &RadYield_E)) );
    parameter_map.insert(para_item("RadYield.P",PARA("RadYield.P","Power parameter for eh escaped yield", "-", 1.0, &RadYield_P)) );
    parameter_map.insert(para_item("ElecMob",   PARA("ElecMob",   "The electron mobility", "cm*cm/V/s", cm*cm/V/s, &ElecMob)) );
    parameter_map.insert(para_item("HoleMob",   PARA("HoleMob",   "The hole mobility",     "cm*cm/V/s", cm*cm/V/s, &HoleMob)) );
    parameter_map.insert(para_item("HIonMob",   PARA("HIonMob",   "The H+ mobility",       "cm*cm/V/s", cm*cm/V/s, &HIonMob)) );
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
  { return 0.0; }

  PetscScalar RadGenRate   (const PetscScalar &E) const
  {
    PetscScalar Y = std::pow((E+1*V/cm)/(E+RadYield_E), RadYield_P);
    return RadGen * Y;
  }

  PetscScalar ElecMobility      (const PetscScalar &Tl) const
  { return ElecMob; }

  PetscScalar HoleMobility      (const PetscScalar &Tl) const
  { return HoleMob; }
  
  PetscScalar HIonMobility      (const PetscScalar &Tl) const
  { return HIonMob; }

  void G4Material(std::vector<Atom> &atoms, std::vector<double> & fraction) const
  {
    atoms.push_back(Atom("Silicon",    "Si", 14, 28.086));//Silicon
    atoms.push_back(Atom("Oxygen",   "O", 8, 16.0));//Oxygen

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
