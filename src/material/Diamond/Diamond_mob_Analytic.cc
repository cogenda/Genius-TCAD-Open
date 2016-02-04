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
// Material Type: Diamond


#include "PMI.h"

class GSS_Diamond_Mob_Constant : public PMIS_Mobility
{
private:

  PetscScalar T300;
  // parameters for constant mobility
  PetscScalar MUN ; //Constant mobility for electron
  PetscScalar MUP ; //Constant mobility for hole

  // parameters for high field modification
  PetscScalar BETAN;
  PetscScalar BETAP;
  PetscScalar VSATN0;
  PetscScalar VSATP0;
  PetscScalar VSATN_A;
  PetscScalar VSATP_A;

  void Mob_Constant_Init()
  {
    T300 = 300*K;
    MUN = 2200*cm*cm/V/s;
    MUP = 1600*cm*cm/V/s;

    BETAN   =  2.000000E+00;
    BETAP   =  1.000000E+00;
    VSATN0  =  2.400000E7*cm/s;
    VSATN_A =  0.8;
    VSATP0  =  1.200000E7*cm/s;
    VSATP_A =  0.8;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("MUN", PARA("MUN", "Constant mobility for electron", "cm*cm/V/s", cm*cm/V/s, &MUN)) );
    parameter_map.insert(para_item("MUP", PARA("MUP", "Constant mobility for hole",     "cm*cm/V/s", cm*cm/V/s, &MUP)) );
    parameter_map.insert(para_item("BETAN",    PARA("BETAN",    "", "-", 1.0, &BETAN)) );
    parameter_map.insert(para_item("BETAP",    PARA("BETAP",    "", "-",1.0 , &BETAP)) );
    parameter_map.insert(para_item("VSATN0",    PARA("VSATN0",    "", "cm/s", cm/s, &VSATN0)) );
    parameter_map.insert(para_item("VSATP0",    PARA("VSATP0",    "", "cm/s", cm/s , &VSATP0)) );
    parameter_map.insert(para_item("VSATN.A",    PARA("VSATN.A",    "", "-", 1.0, &VSATN_A)) );
    parameter_map.insert(para_item("VSATP.A",    PARA("VSATP.A",    "", "-",1.0 , &VSATP_A)) );
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
    PetscScalar mu0  = MUN;
    PetscScalar vsat = VSATN0/(1+VSATN_A*exp(Tl/(2*T300)));
    return mu0/std::pow(1+std::pow(mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN);
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    PetscScalar mu0  = MUN;
    AutoDScalar vsat = VSATN0/(1+VSATN_A*exp(Tl/(2*T300)));
    return mu0/adtl::pow(1+adtl::pow(mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN);
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar mu0  = MUP;
    PetscScalar vsat = VSATP0/(1+VSATP_A*exp(Tl/(2*T300)));
    return mu0/std::pow(1+std::pow(mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    PetscScalar mu0  = MUP;
    AutoDScalar vsat = VSATP0/(1+VSATP_A*exp(Tl/(2*T300)));
    return mu0/adtl::pow(1+adtl::pow(mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP);
  }

  // constructor
public:
  GSS_Diamond_Mob_Constant(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    PMI_Info = "This is the Constant low-field mobility model of Diamond";
    Mob_Constant_Init();
  }


  ~GSS_Diamond_Mob_Constant()
  {}

}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  use Analytic model as default mobility model
 */
extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_Diamond_Mob_Default (const PMIS_Environment& env)
  {
    return new GSS_Diamond_Mob_Constant(env);
  }
}
/* alias */
extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_Diamond_Mob_Analytic (const PMIS_Environment& env)
  {
    return new GSS_Diamond_Mob_Constant(env);
  }
}
