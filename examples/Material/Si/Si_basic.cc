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
// Material Type: Silicon


#include "PMI.h"

class GSS_Si_BasicParameter : public PMIS_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of silicon.
  PetscScalar PERMEABI;  // The relative megnetic permeability of silicon.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.


  TensorValue<PetscScalar> T1;
  TensorValue<PetscScalar> T2;

  void   Basic_Init()
  {
    PERMITTI = 1.170000e+01;
    PERMEABI = 1.0;
    AFFINITY = 4.170000e+00*eV;
    DENSITY  = 2.320000e-03*kg*std::pow(cm,-3);

    //S[1][1] = 0.77  # [1e-12 cm^2/din]
    //S[1][2] = -2.1000e-01   # [1e-12 cm^2/din]
    //S[4][4] = 1.25  # [1e-12 cm^2/din]
    PetscScalar din      = 1e-5*kg*m/s/s;
    PetscScalar S11      = 0.77*1e-12*cm*cm/din;
    PetscScalar S12      = -2.1000e-01*1e-12*cm*cm/din;
    PetscScalar S44      = 1.25*1e-12*cm*cm/din;

    T1 = TensorValue<PetscScalar>(S11, S12, S12, S12, S11, S12, S12, S12, S11);
    T2 = TensorValue<PetscScalar>(S44, 0, 0, 0, S44, 0, 0, 0, S44);

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("PERMITTI",    PARA("PERMITTI",    "The relative dielectric permittivity of silicon", "-", 1.0, &PERMITTI)) );
    parameter_map.insert(para_item("PERMEABI",  PARA("PERMEABI",  "The relative megnetic permeability of silicon", "-", 1.0, &PERMEABI)) );
    parameter_map.insert(para_item("AFFINITY", PARA("AFFINITY", "The electron affinity for the material", "eV", eV, &AFFINITY)) );
    parameter_map.insert(para_item("DENSITY", PARA("DENSITY", "Specific mass density for the material", "kg*cm^-3", kg*std::pow(cm,-3), &DENSITY)) );
#endif
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }


  TensorValue<PetscScalar> Strain(const TensorValue<PetscScalar> & stress) const
  {
    //std::cout<<stress;

    VectorValue<PetscScalar> V1(stress[0], stress[4], stress[8]);
    VectorValue<PetscScalar> V2(stress[1], stress[2], stress[5]);

    //std::cout<<V1;
    //std::cout<<V2;

    VectorValue<PetscScalar> S1 = T1*V1;
    VectorValue<PetscScalar> S2 = T2*V2;

    //std::cout<<S1;
    //std::cout<<S2;

    TensorValue<PetscScalar> strain=TensorValue<PetscScalar>(S1[0], S2[0], S2[1], S2[0], S1[1], S2[2], S2[1], S2[2], S1[2]);
    //std::cout<<strain;
    return strain;
  }


  void G4Material(std::vector<Atom> &atoms, std::vector<double> & fraction) const
  {
    atoms.push_back(Atom("Silicon",    "Si", 14, 28.086));//Silicon
    fraction.push_back(1.0);
  }

  GSS_Si_BasicParameter(const PMIS_Environment &env):PMIS_BasicParameter(env)
  {
    PMI_Info = "This is the Default model for basic physical parameters of Silicon";
    Basic_Init();
  }
  ~GSS_Si_BasicParameter()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE PMIS_BasicParameter* PMIS_Si_BasicParameter_Default (const PMIS_Environment& env)
  {
    return new GSS_Si_BasicParameter(env);
  }
}
