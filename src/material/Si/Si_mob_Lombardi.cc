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
/*  Last update: July 26, 2007                                               */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: Silicon


#include "PMI.h"

class GSS_Si_Mob_Lombardi : public PMIS_Mobility
{
private:
  PetscScalar EXN1_LSM;
  PetscScalar EXN2_LSM;
  PetscScalar EXN3_LSM;
  PetscScalar EXN4_LSM;
  PetscScalar EXN8_LSM;
  PetscScalar MUN0_LSM;
  PetscScalar MUN1_LSM;
  PetscScalar MUN2_LSM;
  PetscScalar CRN_LSM;
  PetscScalar CSN_LSM;
  PetscScalar BN_LSM;
  PetscScalar CN_LSM;
  PetscScalar DN_LSM;

  PetscScalar EXP1_LSM;
  PetscScalar EXP2_LSM;
  PetscScalar EXP3_LSM;
  PetscScalar EXP4_LSM;
  PetscScalar EXP8_LSM;
  PetscScalar MUP0_LSM;
  PetscScalar MUP1_LSM;
  PetscScalar MUP2_LSM;
  PetscScalar CRP_LSM;
  PetscScalar CSP_LSM;
  PetscScalar BP_LSM;
  PetscScalar CP_LSM;
  PetscScalar DP_LSM;
  PetscScalar PC_LSM;

  // parameters for parallel field modification
  PetscScalar BETAN;
  PetscScalar BETAP;
  PetscScalar VSATN0;
  PetscScalar VSATP0;
  PetscScalar VSATN_A;
  PetscScalar VSATP_A;
  PetscScalar T300;

  void Mob_Lombardi_Init()
  {
    EXN1_LSM=6.800000E-01;
    EXN2_LSM=2.000000E+00;
    EXN3_LSM=2.500000E+00;
    EXN4_LSM=1.250000E-01;
    EXN8_LSM=2.000000E+00;
    MUN0_LSM=5.220000E+01*cm*cm/V/s;
    MUN1_LSM=4.340000E+01*cm*cm/V/s;
    MUN2_LSM=1.417000E+03*cm*cm/V/s;
    CRN_LSM=9.680000E+16*std::pow(cm,-3);
    CSN_LSM=3.430000E+20*std::pow(cm,-3);
    BN_LSM=4.750000E+07*cm/s;
    CN_LSM=1.740000E+05*K*cm/s*std::pow(V/cm,PetscScalar(-2.0/3.0))*std::pow(cm,3*EXN4_LSM);
    DN_LSM=5.820000E+14*cm*cm/V/s*std::pow(V/cm,EXN8_LSM);

    EXP1_LSM=7.190000E-01;
    EXP2_LSM=2.000000E+00;
    EXP3_LSM=2.200000E+00;
    EXP4_LSM=3.170000E-02;
    EXP8_LSM=2.000000E+00;
    MUP0_LSM=4.490000E+01*cm*cm/V/s;
    MUP1_LSM=2.900000E+01*cm*cm/V/s;
    MUP2_LSM=4.705000E+02*cm*cm/V/s;
    CRP_LSM=2.230000E+17*std::pow(cm,-3);
    CSP_LSM=6.100000E+20*std::pow(cm,-3);
    BP_LSM=9.930000E+06*cm/s;
    CP_LSM=8.840000E+05*K*cm/s*std::pow(V/cm,PetscScalar(-2.0/3.0))*std::pow(cm,3*EXP4_LSM);
    DP_LSM=2.050000E+14*cm*cm/V/s*std::pow(V/cm,EXN8_LSM);
    PC_LSM=9.230000E+16*std::pow(cm,-3);

    BETAN = 2.000000E+00;
    BETAP = 1.000000E+00;
    VSATN0  =  2.400000E7*cm/s;
    VSATN_A =  0.8;
    VSATP0  =  2.400000E7*cm/s;
    VSATP_A =  0.8;
    T300  = 300.0*K;

#ifdef __CALIBRATE__
         parameter_map.insert(para_item("EXN1.LSM",    PARA("EXN1.LSM",    "", "-", 1.0, &EXN1_LSM)) );
         parameter_map.insert(para_item("EXN2.LSM",    PARA("EXN2.LSM",    "", "-", 1.0, &EXN2_LSM)) );
         parameter_map.insert(para_item("EXN3.LSM",    PARA("EXN3.LSM",    "", "-", 1.0, &EXN3_LSM)) );
         parameter_map.insert(para_item("EXN4.LSM",    PARA("EXN4.LSM",    "", "-", 1.0, &EXN4_LSM)) );
         parameter_map.insert(para_item("EXN8.LSM",    PARA("EXN8.LSM",    "", "-", 1.0, &EXN8_LSM)) );
         parameter_map.insert(para_item("MUN0.LSM",    PARA("MUN0.LSM",    "", "cm*cm/V/s", cm*cm/V/s, &MUN0_LSM)) );
         parameter_map.insert(para_item("MUN1.LSM",    PARA("MUN1.LSM",    "", "cm*cm/V/s",cm*cm/V/s , &MUN1_LSM)) );
         parameter_map.insert(para_item("MUN2.LSM",    PARA("MUN2.LSM",    "", "cm*cm/V/s", cm*cm/V/s, &MUN2_LSM)) );
         parameter_map.insert(para_item("CRN.LSM",    PARA("CRN.LSM",    "", "cm^-3",std::pow(cm,-3) , &CRN_LSM)) );
         parameter_map.insert(para_item("CSN.LSM",    PARA("CSN.LSM",    "", "cm^-3",std::pow(cm,-3) , &CSN_LSM)) );
         parameter_map.insert(para_item("BN.LSM",    PARA("BN.LSM",    "", "cm/s", cm/s, &BN_LSM)) );
         parameter_map.insert(para_item("CN.LSM",    PARA("CN.LSM",    "", "K*cm/s*std::pow(V/cm,PetscScalar(-2.0/3.0))*std::pow(cm,3*EXN4.LSM)", K*cm/s*std::pow(V/cm,PetscScalar(-2.0/3.0))*std::pow(cm,3*EXN4_LSM), &CN_LSM)) );
         parameter_map.insert(para_item("DN.LSM",    PARA("DN.LSM",    "", "cm*cm/V/s*std::pow(V/cm,EXN8.LSM)", cm*cm/V/s*std::pow(V/cm,EXN8_LSM), &DN_LSM)) );
         parameter_map.insert(para_item("EXP1.LSM",    PARA("EXP1.LSM",    "", "-",1.0 , &EXP1_LSM)) );
         parameter_map.insert(para_item("EXP2.LSM",    PARA("EXP2.LSM",    "", "-",1.0 , &EXP2_LSM)) );
         parameter_map.insert(para_item("EXP3.LSM",    PARA("EXP3.LSM",    "", "-",1.0 , &EXP3_LSM)) );
         parameter_map.insert(para_item("EXP4.LSM",    PARA("EXP4.LSM",    "", "-",1.0 , &EXP4_LSM)) );
         parameter_map.insert(para_item("EXP8.LSM",    PARA("EXP8.LSM",    "", "-",1.0 , &EXP8_LSM)) );
         parameter_map.insert(para_item("MUP0.LSM",    PARA("MUP0.LSM",    "", "cm*cm/V/s", cm*cm/V/s, &MUP0_LSM)) );
         parameter_map.insert(para_item("MUP1.LSM",    PARA("MUP1.LSM",    "", "cm*cm/V/s",cm*cm/V/s , &MUP1_LSM)) );
         parameter_map.insert(para_item("MUP2.LSM",    PARA("MUP2.LSM",    "", "cm*cm/V/s",cm*cm/V/s , &MUP2_LSM)) );
         parameter_map.insert(para_item("CRP.LSM",    PARA("CRP.LSM",    "", "cm^-3",std::pow(cm,-3) , &CRP_LSM)) );
         parameter_map.insert(para_item("CSP.LSM",    PARA("CSP.LSM",    "", "cm^-3",std::pow(cm,-3) , &CSP_LSM)) );
         parameter_map.insert(para_item("BP.LSM",    PARA("BP.LSM",    "", "cm/s", cm/s, &BP_LSM)) );
         parameter_map.insert(para_item("CP.LSM",    PARA("CP.LSM",    "", "K*cm/s*std::pow(V/cm,PetscScalar(-2.0/3.0))*std::pow(cm,3*EXP4.LSM)", K*cm/s*std::pow(V/cm,PetscScalar(-2.0/3.0))*std::pow(cm,3*EXP4_LSM), &CP_LSM)) );
         parameter_map.insert(para_item("DP.LSM",    PARA("DP.LSM",    "", "cm*cm/V/s*std::pow(V/cm,EXN8.LSM)",cm*cm/V/s*std::pow(V/cm,EXN8_LSM) , &DP_LSM)) );
         parameter_map.insert(para_item("PC.LSM",    PARA("PC.LSM",    "", "std::pow(cm,-3)", std::pow(cm,-3), &PC_LSM)) );
         parameter_map.insert(para_item("BETAN",    PARA("BETAN",    "", "-",1.0 , &BETAN)) );
         parameter_map.insert(para_item("BETAP",    PARA("BETAP",    "", "-",1.0 , &BETAP)) );
         parameter_map.insert(para_item("VSATN0",    PARA("VSATN0",    "", "cm/s", cm/s, &VSATN0)) );
         parameter_map.insert(para_item("VSATP0",    PARA("VSATP0",    "", "cm/s", cm/s , &VSATP0)) );
         parameter_map.insert(para_item("VSATN.A",    PARA("VSATN.A",    "", "-", 1.0, &VSATN_A)) );
         parameter_map.insert(para_item("VSATP.A",    PARA("VSATP.A",    "", "-",1.0 , &VSATP_A)) );
#endif
  }
  //---------------------------------------------------------------------------
  // Electron low field mobility
  PetscScalar ElecMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    PetscScalar mu_max = MUN2_LSM*std::pow(Tl/T300,-EXN3_LSM);
    return MUN0_LSM+(mu_max-MUN0_LSM)/(1+std::pow(N_total/CRN_LSM,EXN1_LSM))-MUN1_LSM/(1+std::pow(CSN_LSM/N_total,EXN2_LSM));
  }
  AutoDScalar ElecMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    AutoDScalar mu_max = MUN2_LSM*adtl::pow(Tl/T300,-EXN3_LSM);
    return MUN0_LSM+(mu_max-MUN0_LSM)/(1+std::pow(N_total/CRN_LSM,EXN1_LSM))-MUN1_LSM/(1+std::pow(CSN_LSM/N_total,EXN2_LSM));
  }

  //---------------------------------------------------------------------------
  // Hole low field mobility, Analytic model
  PetscScalar HoleMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    PetscScalar mu_max = MUP2_LSM*std::pow(Tl/T300,-EXP3_LSM);
    return MUP0_LSM*exp(-PC_LSM/N_total)+mu_max/(1+std::pow(N_total/CRP_LSM,EXP1_LSM))-MUP1_LSM/(1+std::pow(CSP_LSM/N_total,EXP2_LSM));
  }
  AutoDScalar HoleMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    AutoDScalar mu_max = MUP2_LSM*adtl::pow(Tl/T300,-EXP3_LSM);
    return MUP0_LSM*exp(-PC_LSM/N_total)+mu_max/(1+std::pow(N_total/CRP_LSM,EXP1_LSM))-MUP1_LSM/(1+std::pow(CSP_LSM/N_total,EXP2_LSM));
  }

  //---------------------------------------------------------------------------
  // Electron surface mobility, acoustical phono scattering and roughness scattering
  PetscScalar ElecMobSurface(const PetscScalar &Tl,const PetscScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    PetscScalar ET = Et+1.0*V/cm;
    PetscScalar mu_ac = BN_LSM/ET + CN_LSM*std::pow(N_total,EXN4_LSM)/Tl*std::pow(ET,PetscScalar(-1.0/3.0));
    PetscScalar mu_sr = DN_LSM*std::pow(ET,-EXN8_LSM);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }
  AutoDScalar ElecMobSurface(const AutoDScalar &Tl,const AutoDScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    AutoDScalar ET = Et+1.0*V/cm;
    AutoDScalar mu_ac = BN_LSM/ET + CN_LSM*std::pow(N_total,EXN4_LSM)/Tl*adtl::pow(ET,PetscScalar(-1.0/3.0));
    AutoDScalar mu_sr = DN_LSM*adtl::pow(ET,-EXN8_LSM);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }

  //---------------------------------------------------------------------------
  // Hole surface mobility, acoustical phono scattering and roughness scattering
  PetscScalar HoleMobSurface(const PetscScalar &Tl,const PetscScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    PetscScalar ET = Et+1.0*V/cm;
    PetscScalar mu_ac = BP_LSM/ET + CP_LSM*std::pow(N_total,EXP4_LSM)/Tl*std::pow(ET,PetscScalar(-1.0/3.0));
    PetscScalar mu_sr = DP_LSM*std::pow(ET,-EXP8_LSM);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }
  AutoDScalar HoleMobSurface(const AutoDScalar &Tl,const AutoDScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    AutoDScalar ET = Et+1.0*V/cm;
    AutoDScalar mu_ac = BP_LSM/ET + CP_LSM*std::pow(N_total,EXP4_LSM)/Tl*adtl::pow(ET+1.0*V/cm,PetscScalar(-1.0/3.0));
    AutoDScalar mu_sr = DP_LSM*adtl::pow(ET,-EXP8_LSM);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
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
    PetscScalar vsat = VSATN0/(1+VSATN_A*exp(Tl/(2*T300)));
    PetscScalar mu0  = 1.0/(1.0/ElecMobLowField(Tl)+1.0/ElecMobSurface(Tl,Et));
    return mu0/std::pow(1+std::pow(mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN);
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    AutoDScalar vsat = VSATN0/(1+VSATN_A*exp(Tl/(2*T300)));
    AutoDScalar mu0  = 1.0/(1.0/ElecMobLowField(Tl)+1.0/ElecMobSurface(Tl,Et));
    return mu0/adtl::pow(1+adtl::pow(mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN);
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = VSATP0/(1+VSATP_A*exp(Tl/(2*T300)));
    PetscScalar mu0  = 1.0/(1.0/HoleMobLowField(Tl)+1.0/HoleMobSurface(Tl,Et));
    return mu0/std::pow(1+std::pow(mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    AutoDScalar vsat = VSATP0/(1+VSATP_A*exp(Tl/(2*T300)));
    AutoDScalar mu0  = 1.0/(1.0/HoleMobLowField(Tl)+1.0/HoleMobSurface(Tl,Et));
    return mu0/adtl::pow(1+adtl::pow(mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP);
  }

  // constructor
public:
  GSS_Si_Mob_Lombardi(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    PMI_Info = "This is the Lombardi mobility model of Silicon";
    Mob_Lombardi_Init();
  }

  ~GSS_Si_Mob_Lombardi(){}

};

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  and it setup Lombardi mobility model
 */
extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_Si_Mob_Lombardi (const PMIS_Environment& env)
  {
    return new GSS_Si_Mob_Lombardi(env);
  }
}
