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

class GSS_Epoxy_BasicParameter : public PMII_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of Kapton.
  PetscScalar PERMEABI;  // The relative megnetic permeability of Kapton.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.

  PetscScalar CONDUCTANCE;     // Intrinsic conductance
  PetscScalar MOBILITY;  // Electron mobility
  PetscScalar Kp;        // RIC dose rate factor
  PetscScalar Delta;     // power parameter of RIC dose rate factor

  void   Basic_Init()
  {
    PERMITTI = 3.6;
    PERMEABI = 1.0;
    AFFINITY = 2.5*eV;
    DENSITY  = 1.5*g*std::pow(cm,-3);

    PetscScalar rad = 0.01*J/kg;

    CONDUCTANCE = 0.5e-15/(Ohm*m);
    MOBILITY    = 5e-13*cm*cm/V/s;
    Delta       = 1.0;
    Kp          = 6e-14/(Ohm*m)/std::pow(rad/s, Delta);
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("PERMITTI",    PARA("PERMITTI",    "The relative dielectric permittivity", "-", 1.0, &PERMITTI)) );
    parameter_map.insert(para_item("PERMEABI",  PARA("PERMEABI",  "The relative megnetic permeability ", "-", 1.0, &PERMEABI)) );
    parameter_map.insert(para_item("AFFINITY", PARA("AFFINITY", "The electron affinity for the material", "eV", eV, &AFFINITY)) );
    parameter_map.insert(para_item("DENSITY", PARA("DENSITY", "Specific mass density for the material", "kg*cm^-3", kg*std::pow(cm,-3), &DENSITY)) );
    parameter_map.insert(para_item("CONDUCTANCE", PARA("CONDUCTANCE", "Specific intrinsic conductance", "Ohm^-1*m^-1", 1.0/(Ohm*m), &CONDUCTANCE)) );
    parameter_map.insert(para_item("DELTA", PARA("DELTA", "Power parameter of RIC dose rate factor", "-", 1.0, &Delta)) );
    parameter_map.insert(para_item("KP", PARA("KP", "RIC dose rate factor", "Ohm^-1*m^-1/(rad/s)^DELTA", 1.0/(Ohm*m)/std::pow(rad/s, Delta), &Kp)) );
#endif
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }

  PetscScalar Conductance   (const PetscScalar &Tl, const PetscScalar &E) const
  {
    PetscScalar Ep = E + 1*V/cm;
    if(Ep > 1e5*V/cm) Ep=1e5*V/cm;
    PetscScalar alpha = 2*kb*Tl/(e*Ep*1e-3*um);
    PetscScalar beta  = sqrt(Ep*e*e*e/(3.14159*PERMITTI*eps0))/(2*kb*Tl);
    PetscScalar eh = (2+cosh(beta))/3.0*(alpha*sinh(1/alpha));
    return CONDUCTANCE  * eh ;//* exp(-1*eV/(kb*Tl));
  }

  AutoDScalar Conductance   (const AutoDScalar &Tl, const AutoDScalar &E) const
  {
    AutoDScalar Ep = E + 1*V/cm;
    if(Ep > 1e5*V/cm) Ep=1e5*V/cm;
    AutoDScalar alpha = 2*kb*Tl/(e*Ep*1e-3*um);
    AutoDScalar beta  = sqrt(Ep*e*e*e/(3.14159*PERMITTI*eps0))/(2*kb*Tl);
    AutoDScalar eh = (2+cosh(beta))/3.0*(alpha*sinh(1/alpha));
    return CONDUCTANCE  * eh ;//* exp(-1*eV/(kb*Tl));
  }

  PetscScalar RadConductance(const PetscScalar &DRate)const
  { return Kp*std::pow(DRate, Delta); }

  PetscScalar RadGenRate   (const PetscScalar &E) const
  { return 0.0; }

  PetscScalar ElecMobility      (const PetscScalar &Tl) const
  { return MOBILITY; }

  PetscScalar HoleMobility      (const PetscScalar &Tl) const
  { return MOBILITY; }

  void G4Material(std::vector<Atom> &atoms, std::vector<double> & fraction) const
  {
    // http://images2.freshpatents.com/pdf/US20060204760A1.pdf
    atoms.push_back(Atom("Hydrogen", "H", 1, 1.00797)); //Hydrogen
    atoms.push_back(Atom("Carbon",   "C", 6, 12.01115)); //Carbon
    atoms.push_back(Atom("Nitrogen", "N", 7, 14.01)); //Oxygen
    atoms.push_back(Atom("Oxygen",   "O", 8, 16.0)); //Nitrogen

    fraction.push_back(15.0);
    fraction.push_back(12.0);
    fraction.push_back(3.0);
    fraction.push_back(6.0);
  }


  GSS_Epoxy_BasicParameter(const PMII_Environment &env):PMII_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_Epoxy_BasicParameter()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMII_BasicParameter* PMII_Epoxy_BasicParameter_Default (const PMII_Environment& env)
  {
    return new GSS_Epoxy_BasicParameter(env);
  }
}
