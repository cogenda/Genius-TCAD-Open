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
// Material Type: 4H-SiC, use parameters for Si now.


#include "PMI.h"

class GSS_SiC4H_Mob_Masetti : public PMIS_Mobility
{
private:
  // parameters for Analytic mobility
  PetscScalar MUN_MAX  ;
  PetscScalar MUN_ZETA ;
  PetscScalar MUN_MIN1 ;
  PetscScalar MUN_MIN2 ;
  PetscScalar MUN1     ;
  PetscScalar PCN      ;
  PetscScalar CRN      ;
  PetscScalar CSN      ;
  PetscScalar MUN_ALPHA;
  PetscScalar MUN_BETA ;

  PetscScalar MUP_MAX  ;
  PetscScalar MUP_ZETA ;
  PetscScalar MUP_MIN1 ;
  PetscScalar MUP_MIN2 ;
  PetscScalar MUP1     ;
  PetscScalar PCP      ;
  PetscScalar CRP      ;
  PetscScalar CSP      ;
  PetscScalar MUP_ALPHA;
  PetscScalar MUP_BETA ;

  PetscScalar T300     ;

  // parameters for high field modification
  PetscScalar VSATN_A;
  PetscScalar VSATN_B;
  PetscScalar VSATP_A;
  PetscScalar VSATP_B;

  PetscScalar BETAN;
  PetscScalar BETAP;

  void Mob_Analytic_Init()
  {
    MUN_MAX   = 9.470000E+02 *cm*cm/V/s;
    MUN_ZETA  = 1.962000E+00;
    MUN_MIN1  = 0 *cm*cm/V/s;
    MUN_MIN2  = 0 *cm*cm/V/s;
    MUN1      = 0 *cm*cm/V/s;
    PCN       = 0 *std::pow(cm,-3);
    CRN       = 1.940000E+17 *std::pow(cm,-3);
    CSN       = 0 *std::pow(cm,-3);
    MUN_ALPHA = 0.61;
    MUN_BETA  = 0;

    MUP_MAX   = 1.240000E+02 *cm*cm/V/s;
    MUN_ZETA  = 1.424000E+00;
    MUP_MIN1  = 1.590000E+01 *cm*cm/V/s;
    MUP_MIN2  = 1.590000E+01 *cm*cm/V/s;
    MUP1      = 0 *cm*cm/V/s;
    PCP       = 0 *std::pow(cm,-3);
    CRP       = 1.760000E+19 *std::pow(cm,-3);
    CSP       = 0 *std::pow(cm,-3);
    MUP_ALPHA = 0.34;
    MUP_BETA  = 0;

    VSATN_A   =  1.070000E+07 * cm/s;
    VSATN_B   =  0.000000E+00 * cm/s;
    VSATP_A   =  8.370000E+06 * cm/s;
    VSATP_B   =  0.000000E+00 * cm/s;

    BETAN     = 0.84;
    BETAP     = 0.84;

    T300 = 300.0*K;
#ifdef __CALIBRATE__
	 parameter_map.insert(para_item("MUN.MIN1",   PARA("MUN.MIN1",    "", "cm*cm/V/s",cm*cm/V/s , &MUN_MIN1)) );
	 parameter_map.insert(para_item("MUN.MIN2",   PARA("MUN.MIN2",    "", "cm*cm/V/s",cm*cm/V/s , &MUN_MIN2)) );
	 parameter_map.insert(para_item("MUN.MAX",    PARA("MUN.MAX",    "", "cm*cm/V/s", cm*cm/V/s, &MUN_MAX)) );
	 parameter_map.insert(para_item("MUN1",       PARA("MUN1",     "", "cm*cm/V/s",cm*cm/V/s , &MUN1)) );

	 parameter_map.insert(para_item("PCN",        PARA("PCN",    "", "cm^-3", std::pow(cm,-3), &PCN)) );
	 parameter_map.insert(para_item("CRN",        PARA("CRN",    "", "cm^-3", std::pow(cm,-3), &CRN)) );
	 parameter_map.insert(para_item("CSN",        PARA("CSN",    "", "cm^-3", std::pow(cm,-3), &CSN)) );
	 parameter_map.insert(para_item("MUN.ALPHA",  PARA("MUN.ALPHA",    "", "-", 1.0, &MUN_ALPHA)) );
	 parameter_map.insert(para_item("MUN.BETA",   PARA("MUN.BETA",    "", "-", 1.0, &MUN_BETA)) );
	 parameter_map.insert(para_item("MUN.ZETA",   PARA("MUN.ZETA",    "", "-", 1.0, &MUN_ZETA)) );

	 parameter_map.insert(para_item("MUP.MIN1",   PARA("MUP.MIN1",    "", "cm*cm/V/s",cm*cm/V/s , &MUP_MIN1)) );
	 parameter_map.insert(para_item("MUP.MIN2",   PARA("MUP.MIN2",    "", "cm*cm/V/s",cm*cm/V/s , &MUP_MIN2)) );
	 parameter_map.insert(para_item("MUP.MAX",    PARA("MUP.MAX",    "", "cm*cm/V/s", cm*cm/V/s, &MUP_MAX)) );
	 parameter_map.insert(para_item("MUP1",       PARA("MUP1",     "", "cm*cm/V/s",cm*cm/V/s , &MUP1)) );

	 parameter_map.insert(para_item("PCP",        PARA("PCP",    "", "cm^-3", std::pow(cm,-3), &PCP)) );
	 parameter_map.insert(para_item("CRP",        PARA("CRP",    "", "cm^-3", std::pow(cm,-3), &CRP)) );
	 parameter_map.insert(para_item("CSP",        PARA("CSP",    "", "cm^-3", std::pow(cm,-3), &CSP)) );
	 parameter_map.insert(para_item("MUP.ALPHA",  PARA("MUP.ALPHA",    "", "-", 1.0, &MUP_ALPHA)) );
	 parameter_map.insert(para_item("MUP.BETA",   PARA("MUP.BETA",    "", "-", 1.0, &MUP_BETA)) );
	 parameter_map.insert(para_item("MUP.ZETA",   PARA("MUP.ZETA",    "", "-", 1.0, &MUP_ZETA)) );

   parameter_map.insert(para_item("VSATN_A",    PARA("VSATN_A",    "", "cm/s", cm/s, &VSATN_A)) );
	 parameter_map.insert(para_item("VSATN_B",    PARA("VSATN_B",    "", "cm/s", cm/s, &VSATN_B)) );
   parameter_map.insert(para_item("VSATP_A",    PARA("VSATP_A",    "", "cm/s", cm/s, &VSATP_A)) );
	 parameter_map.insert(para_item("VSATP_B",    PARA("VSATP_B",    "", "cm/s", cm/s, &VSATP_B)) );

   parameter_map.insert(para_item("BETAN",      PARA("BETAN",      "", "-", 1.0, &BETAN)) );
   parameter_map.insert(para_item("BETAP",      PARA("BETAP",      "", "-", 1.0, &BETAP)) );
#endif
  }

public:
  //---------------------------------------------------------------------------
  // Electron low field mobility
  PetscScalar ElecMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar mu_const = MUN_MAX * std::pow(Tl/T300,-MUN_ZETA);
    return MUN_MIN1 * exp(-PCN/(Na+Nd)) \
           + (mu_const - MUN_MIN2)/(1+std::pow((Na+Nd)/CRN,MUN_ALPHA)) \
           - MUN1/(1+ std::pow(CSN/(Na+Nd),MUN_BETA));
  }
  AutoDScalar ElecMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    AutoDScalar mu_const = MUN_MAX * adtl::pow(Tl/T300,-MUN_ZETA);
    return MUN_MIN1 * exp(-PCN/(Na+Nd)) \
           + (mu_const - MUN_MIN2)/(1+std::pow((Na+Nd)/CRN, MUN_ALPHA)) \
           - MUN1/(1+ std::pow(CSN/(Na+Nd),MUN_BETA));
  }

  //---------------------------------------------------------------------------
  // Hole low field mobility
  PetscScalar HoleMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar mu_const = MUP_MAX * std::pow(Tl/T300,-MUP_ZETA);
    return MUP_MIN1 * exp(-PCP/(Na+Nd)) \
           + (mu_const - MUP_MIN2)/(1+std::pow((Na+Nd)/CRP,MUP_ALPHA)) \
           - MUP1/(1+ std::pow(CSP/(Na+Nd),MUP_BETA));
  }
  AutoDScalar HoleMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    AutoDScalar mu_const = MUP_MAX * adtl::pow(Tl/T300,-MUP_ZETA);
    return MUP_MIN1 * exp(-PCP/(Na+Nd)) \
           + (mu_const - MUP_MIN2)/(1+std::pow((Na+Nd)/CRP,MUP_ALPHA)) \
           - MUP1/(1+ std::pow(CSP/(Na+Nd),MUP_BETA));
  }


public:
  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
    PetscScalar vsat = VSATN_A - VSATN_B * (Tl/T300);
    PetscScalar mu0  = ElecMobLowField(Tl);
    return mu0/std::pow(1+std::pow(fabs(mu0*Ep/vsat)+1e-10,BETAN),1.0/BETAN);
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    AutoDScalar vsat = VSATN_A - VSATN_B * (Tl/T300);
    AutoDScalar mu0  = ElecMobLowField(Tl);
    return mu0/adtl::pow(1+adtl::pow(fabs(mu0*Ep/vsat)+1e-10,BETAN),1.0/BETAN);
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = VSATN_A - VSATN_B * (Tl/T300);
    PetscScalar mu0  = HoleMobLowField(Tl);
    return mu0/std::pow(1+std::pow(fabs(mu0*Ep/vsat)+1e-10,BETAP),1.0/BETAP);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    AutoDScalar vsat = VSATN_A - VSATN_B * (Tl/T300);
    AutoDScalar mu0  = HoleMobLowField(Tl);
    return mu0/adtl::pow(1+adtl::pow(fabs(mu0*Ep/vsat)+1e-10,BETAP),1.0/BETAP);
  }


// constructor and destructor
public:
  GSS_SiC4H_Mob_Masetti(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    Mob_Analytic_Init();
  }

  ~GSS_SiC4H_Mob_Masetti()
  {
  }

}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  use Analytic model as default mobility model
 */
extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_SiC4H_Mob_Default (const PMIS_Environment& env)
  {
    return new GSS_SiC4H_Mob_Masetti(env);
  }
}
/* alias */
extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_SiC4H_Mob_Masetti (const PMIS_Environment& env)
  {
    return new GSS_SiC4H_Mob_Masetti(env);
  }
}
