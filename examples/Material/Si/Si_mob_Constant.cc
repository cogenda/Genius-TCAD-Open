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

class GSS_Si_Mob_Constant : public PMIS_Mobility
{
private:
  // parameters for constant mobility
  PetscScalar MUN ; //Constant mobility for electron
  PetscScalar MUP ; //Constant mobility for hole

  void Mob_Constant_Init()
  {
    MUN = 1.400E+03*cm*cm/V/s;
    MUP = 4.500E+02*cm*cm/V/s;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("MUN", PARA("MUN", "Constant mobility for electron", "cm*cm/V/s", cm*cm/V/s, &MUN)) );
    parameter_map.insert(para_item("MUP", PARA("MUP", "Constant mobility for hole",     "cm*cm/V/s", cm*cm/V/s, &MUP)) );
#endif
  }


public:

  // Hall mobility factor  for electrons
  PetscScalar RH_ELEC()  { return 1.1; }

  // Hall mobility factor  for holes
  PetscScalar RH_HOLE()  { return 0.7; }

  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
     return MUN;
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    return MUN;
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
     return MUP;
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    return MUP;
  }

// constructor
public:
  GSS_Si_Mob_Constant(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    PMI_Info = "This is the Constant low-field mobility model of Silicon";
    Mob_Constant_Init();
  }


  ~GSS_Si_Mob_Constant()
  {
  }

}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 */
extern "C"
{
  DLL_EXPORT_DECLARE PMIS_Mobility* PMIS_Si_Mob_Constant (const PMIS_Environment& env)
  {
    return new GSS_Si_Mob_Constant(env);
  }
}
