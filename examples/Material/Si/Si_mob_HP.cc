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
/*  Last update: July 23, 2007                                               */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: Silicon


#include "PMI.h"

class GSS_Si_Mob_HP : public PMIS_Mobility
{
private:
  // parameters for Analytic mobility
  PetscScalar MUN_MIN ;
  PetscScalar MUN_MAX ;
  PetscScalar NREFN   ;
  PetscScalar NUN     ;
  PetscScalar XIN     ;
  PetscScalar ALPHAN  ;
  PetscScalar MUP_MIN ;
  PetscScalar MUP_MAX ;
  PetscScalar NREFP   ;
  PetscScalar NUP     ;
  PetscScalar XIP     ;
  PetscScalar ALPHAP  ;
  PetscScalar T300    ;
  // parameters for HP high field modification
  PetscScalar   MUN0_HP;
  PetscScalar   ECN_HP;
  PetscScalar   VSN_HP;
  PetscScalar   VCN_HP;
  PetscScalar   GN_HP;
  PetscScalar   MUP0_HP;
  PetscScalar   ECP_HP;
  PetscScalar   VSP_HP;
  PetscScalar   VCP_HP;
  PetscScalar   GP_HP;
  PetscScalar   NREF_HP;

  void Mob_HP_Init()
  {
    MUN_MIN = 5.524000E+01*cm*cm/V/s;
    MUN_MAX = 1.429230E+03*cm*cm/V/s;
    NREFN   = 1.072000E+17*std::pow(cm,-3);
    NUN     = -2.300000E+00;
    XIN     = -3.800000E+00;
    ALPHAN  = 7.300000E-01;
    MUP_MIN = 4.970000E+01*cm*cm/V/s;
    MUP_MAX = 4.793700E+02*cm*cm/V/s;
    NREFP   = 1.606000E+17*std::pow(cm,-3);
    NUP     = -2.200000E+00;
    XIP     = -3.700000E+00;
    ALPHAP  = 7.000000E-01;
    T300    = 300.0*K;

    MUN0_HP = 7.740000E+02*cm*cm/V/s;
    ECN_HP  = 5.500000E+05*V/cm;
    VSN_HP  = 1.036000E+07*cm/s;
    VCN_HP  = 4.900000E+06*cm/s;
    GN_HP   = 8.800000E+00;
    MUP0_HP = 2.500000E+02*cm*cm/V/s;
    ECP_HP  = 2.780000E+05*V/cm;
    VSP_HP  = 1.200000E+07*cm/s;
    VCP_HP  = 2.928000E+06*cm/s;
    GP_HP   = 1.600000E+00;
    NREF_HP = 5.000000E+17*std::pow(cm,-3);

#ifdef __CALIBRATE__
         parameter_map.insert(para_item("MUN.MIN",    PARA("MUN.MIN",    "", "cm*cm/V/s", cm*cm/V/s, &MUN_MIN)) );
         parameter_map.insert(para_item("MUN.MAX",    PARA("MUN.MAX",    "", "cm*cm/V/s",cm*cm/V/s , &MUN_MAX)) );
         parameter_map.insert(para_item("NREFN",    PARA("NREFN",    "", "cm^-3", std::pow(cm,-3), &NREFN)) );
         parameter_map.insert(para_item("NUN",    PARA("NUN",    "", "-", 1.0, &NUN)) );
         parameter_map.insert(para_item("XIN",    PARA("XIN",    "", "-", 1.0, &XIN)) );
         parameter_map.insert(para_item("ALPHAN",    PARA("ALPHAN",    "", "-", 1.0, &ALPHAN)) );
         parameter_map.insert(para_item("MUP.MIN",    PARA("MUP.MIN",    "", "cm*cm/V/s", cm*cm/V/s, &MUP_MIN)) );
         parameter_map.insert(para_item("MUP.MAX",    PARA("MUP.MAX",    "", "cm*cm/V/s",cm*cm/V/s , &MUP_MAX)) );
         parameter_map.insert(para_item("NREFP",    PARA("NREFP",    "", "cm^-3",std::pow(cm,-3) , &NREFP)) );
         parameter_map.insert(para_item("NUP",    PARA("NUP",    "", "-",1.0 , &NUP)) );
         parameter_map.insert(para_item("XIP",    PARA("XIP",    "", "-",1.0 , &XIP)) );
         parameter_map.insert(para_item("ALPHAP",    PARA("ALPHAP",    "", "-", 1.0, &ALPHAP)) );
         parameter_map.insert(para_item("MUN0.HP",    PARA("MUN0.HP",    "", "cm*cm/V/s",cm*cm/V/s , &MUN0_HP)) );
         parameter_map.insert(para_item("ECN.HP",    PARA("ECN.HP",    "", "V/cm",V/cm , &ECN_HP)) );
         parameter_map.insert(para_item("VSN.HP",    PARA("VSN.HP",    "", "cm/s",cm/s , &VSN_HP)) );
         parameter_map.insert(para_item("VCN.HP ",    PARA("VCN.HP ",    "", "cm/s",cm/s , &VCN_HP )) );
         parameter_map.insert(para_item("GN.HP",    PARA("GN.HP",    "", "-", 1.0, &GN_HP)) );
         parameter_map.insert(para_item("MUP0.HP",    PARA("MUP0.HP",    "", "cm*cm/V/s",cm*cm/V/s , &MUP0_HP)) );
         parameter_map.insert(para_item("ECP.HP",    PARA("ECP.HP",    "", "V/cm",V/cm , &ECP_HP)) );
         parameter_map.insert(para_item("VSP.HP",    PARA("VSP.HP",    "", "cm/s",cm/s , &VSP_HP)) );
         parameter_map.insert(para_item("VCP.HP",    PARA("VCP.HP",    "", "cm/s", cm/s, &VCP_HP)) );
         parameter_map.insert(para_item("GP.HP",    PARA("GP.HP",    "", "-",1.0 , &GP_HP)) );
         parameter_map.insert(para_item("NREF.HP",    PARA("NREF.HP",    "", "cm^-3",std::pow(cm,-3) , &NREF_HP)) );
#endif
  }

private:
  //---------------------------------------------------------------------------
  // Electron low field mobility, Analytic model
  PetscScalar ElecMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUN_MIN+(MUN_MAX*std::pow(Tl/T300,NUN)-MUN_MIN)/ \
           (1+std::pow(Tl/T300,XIN)*std::pow((Na+Nd)/NREFN,ALPHAN));
  }
  AutoDScalar ElecMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUN_MIN+(MUN_MAX*adtl::pow(Tl/T300,NUN)-MUN_MIN)/ \
           (1+adtl::pow(Tl/T300,XIN)*std::pow((Na+Nd)/NREFN,ALPHAN));
  }

  //---------------------------------------------------------------------------
  // Hole low field mobility, Analytic model
  PetscScalar HoleMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUP_MIN+(MUP_MAX*std::pow(Tl/T300,NUP)-MUP_MIN)/ \
           (1+std::pow(Tl/T300,XIP)*std::pow((Na+Nd)/NREFP,ALPHAP));
  }
  AutoDScalar HoleMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUP_MIN+(MUP_MAX*adtl::pow(Tl/T300,NUP)-MUP_MIN)/ \
           (1+adtl::pow(Tl/T300,XIP)*std::pow((Na+Nd)/NREFP,ALPHAP));
  }


public:
  /**
   * A factor used in determining the effective electric field at interfaces
   * used in the field-dependent mobility models for electrons.
   */
  PetscScalar ZETAN() { return 1.0; }

  /**
   * A factor used in determining the effective electric field at interfaces
   * used in the field-dependent mobility models for electrons.
   */
  PetscScalar ETAN()  { return 0.5; }

  /**
   * A factor used in determining the effective electric field at interfaces
   * used in field-dependent mobility models for holes.
   */
  PetscScalar ZETAP() { return 1.0; }

  /**
   * A factor used in determining the effective electric field at interfaces
   * used in field-dependent mobility models for holes.
   */
  PetscScalar ETAP()  { return 0.333; }


  // Hall mobility factor  for electrons
  PetscScalar RH_ELEC()  { return 1.1; }

  // Hall mobility factor  for holes
  PetscScalar RH_HOLE()  { return 0.7; }


  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
    PetscScalar Ntotal = ReadDopingNa()+ReadDopingNd();
    PetscScalar mu0;
    if(Ntotal > NREF_HP) mu0=ElecMobLowField(Tl);
    else mu0 = MUN0_HP/(1+Et/ECN_HP);
    PetscScalar alpha = mu0*Ep/VCN_HP;
    PetscScalar beta  = mu0*Ep/VSN_HP;
    return mu0/sqrt(1+alpha*alpha/(alpha+GN_HP)+beta*beta);
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    PetscScalar Ntotal = ReadDopingNa()+ReadDopingNd();
    AutoDScalar mu0;
    if(Ntotal > NREF_HP) mu0=ElecMobLowField(Tl);
    else mu0 = MUN0_HP/(1+Et/ECN_HP);
    AutoDScalar alpha = mu0*Ep/VCN_HP;
    AutoDScalar beta  = mu0*Ep/VSN_HP;
    return mu0/sqrt(1+alpha*alpha/(alpha+GN_HP)+beta*beta);
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar Ntotal = ReadDopingNa()+ReadDopingNd();
    PetscScalar mu0;
    if(Ntotal > NREF_HP) mu0=HoleMobLowField(Tl);
    else mu0 = MUP0_HP/(1+Et/ECP_HP);
    PetscScalar alpha = mu0*Ep/VCP_HP;
    PetscScalar beta  = mu0*Ep/VSP_HP;
    return mu0/sqrt(1+alpha*alpha/(alpha+GP_HP)+beta*beta);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    PetscScalar Ntotal = ReadDopingNa()+ReadDopingNd();
    AutoDScalar mu0;
    if(Ntotal > NREF_HP) mu0=HoleMobLowField(Tl);
    else mu0 = MUP0_HP/(1+Et/ECP_HP);
    AutoDScalar alpha = mu0*Ep/VCP_HP;
    AutoDScalar beta  = mu0*Ep/VSP_HP;
    return mu0/sqrt(1+alpha*alpha/(alpha+GP_HP)+beta*beta);
  }

// constructor
public:
  GSS_Si_Mob_HP(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    PMI_Info = "This is the Hewlett-Packard mobility model of Silicon";
    Mob_HP_Init();
  }


  ~GSS_Si_Mob_HP()
  {
  }

}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  load HP mobility model
 */
/* alias */
extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_Si_Mob_HP (const PMIS_Environment& env)
  {
    return new GSS_Si_Mob_HP(env);
  }
}
