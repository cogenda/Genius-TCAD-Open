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
// Material Type: Al(x)Ga(1-x)As


#include "PMI.h"

class GSS_AlGaAs_BasicParameter : public PMIS_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of AlGaAs.
  PetscScalar EPS_X1;
  PetscScalar EPS_X2;
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar AF_X0;
  PetscScalar AF_X1;
  PetscScalar AF_X2;
  PetscScalar AF_X3;
  PetscScalar AF_X4;
  PetscScalar AF_X5;
  PetscScalar AF_XL;
  PetscScalar PERMEABI;  // The relative megnetic permeability of AlGaAs.
  PetscScalar DENSITY;   // Specific mass density for the material.

  void   Basic_Init()
  {
    PERMITTI =  1.310000E+01;
    EPS_X1   =  0.000000E+00;
    EPS_X2   =  0.000000E+00;
    AFFINITY =  4.070000e+00*eV;
    AF_X0    =  0.000000E+00;
    AF_X1    = -1.100000E+00;
    AF_X2    =  0.000000E+00;
    AF_X3    = -4.300000E-01;
    AF_X4    = -1.400000E-01;
    AF_X5    =  0.000000E+00;
    AF_XL    =  4.500000E-01;
    PERMEABI =  1.0;
    DENSITY  =  5.317600E-03*kg*std::pow(cm,-3);

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("PERMITTI",    PARA("PERMITTI",    "The relative dielectric permittivity", "-", 1.0, &PERMITTI)) );
    parameter_map.insert(para_item("PERMEABI",  PARA("PERMEABI",  "The relative megnetic permeability ", "-", 1.0, &PERMEABI)) );
    parameter_map.insert(para_item("AFFINITY", PARA("AFFINITY", "The electron affinity for the material", "eV", eV, &AFFINITY)) );
    parameter_map.insert(para_item("DENSITY", PARA("DENSITY", "Specific mass density for the material", "kg*cm^-3", kg*std::pow(cm,-3), &DENSITY)) );
    parameter_map.insert(para_item("EPS.X1", PARA("EPS.X1", "", "-",1.0 , &EPS_X1)) );
    parameter_map.insert(para_item("EPS.X2", PARA("EPS.X2", "", "-",1.0 , &EPS_X2)) );
    parameter_map.insert(para_item("AF.X0", PARA("AF.X0", "", "-",1.0 , &AF_X0)) );
    parameter_map.insert(para_item("AF.X1", PARA("AF.X1", "", "-",1.0 , &AF_X1)) );
    parameter_map.insert(para_item("AF.X2", PARA("AF.X2", "", "-",1.0 , &AF_X2)) );
    parameter_map.insert(para_item("AF.X3", PARA("AF.X3", "", "-",1.0 , &AF_X3)) );
    parameter_map.insert(para_item("AF.X4", PARA("AF.X4", "", "-",1.0 , &AF_X4)) );
    parameter_map.insert(para_item("AF.X5", PARA("AF.X5", "", "-",1.0 , &AF_X5)) );
    parameter_map.insert(para_item("AF.XL", PARA("AF.XL", "", "-",1.0 , &AF_XL)) );
#endif
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const
  {
  	return DENSITY;
  }
  PetscScalar Permittivity() const
  {
        PetscScalar mole_x = ReadxMoleFraction();
        return PERMITTI + EPS_X1*mole_x + EPS_X2*mole_x*mole_x;
  }
  PetscScalar Permeability() const
  {
  	return PERMEABI;
  }
  PetscScalar Affinity      (const PetscScalar &Tl) const
  {
        PetscScalar mole_x = ReadxMoleFraction();
        if(mole_x<AF_XL)
                return AFFINITY + AF_X0 + AF_X1*mole_x + AF_X2*mole_x*mole_x;
        else
                return AFFINITY + AF_X3 + AF_X4*mole_x + AF_X5*mole_x*mole_x;
  }

  void G4Material(std::vector<Atom> &atoms, std::vector<double> & fraction) const
  {
    atoms.push_back(Atom("Aluminum",   "Al", 13, 26.98)); //Aluminum
    atoms.push_back(Atom("Gallium",   "Ga", 31, 69.72)); //Gallium
    atoms.push_back(Atom("Arsenic",   "As", 33, 74.922)); //Arsenic

    PetscScalar mole_x = ReadxMoleFraction();
    fraction.push_back(mole_x);
    fraction.push_back(1.0-mole_x);
    fraction.push_back(1.0);
  }

  GSS_AlGaAs_BasicParameter(const PMIS_Environment &env):PMIS_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_AlGaAs_BasicParameter()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_BasicParameter* PMIS_AlGaAs_BasicParameter_Default (const PMIS_Environment& env)
  {
    return new GSS_AlGaAs_BasicParameter(env);
  }
}
