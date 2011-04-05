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
// Material Type: GaAs


#include "PMI.h"

class GSS_GaAs_Mob_Analytic : public PMIS_Mobility
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
  // parameters for high field modification
  PetscScalar E0N     ;
  PetscScalar E0P     ;
  PetscScalar VSATN_A;
  PetscScalar VSATP_A;
  PetscScalar VSATN_B;
  PetscScalar VSATP_B;

  void Mob_Analytic_Init()
  {
    MUN_MIN =  0.000000E+00*cm*cm/V/s;
    MUN_MAX =  8.500000E+03*cm*cm/V/s;
    NREFN   =  1.690000E+17*std::pow(cm,-3);
    NUN     = -1.000000E+00;
    XIN     =  0.000000E+00;
    ALPHAN  =  4.360000E-01;
    MUP_MIN =  0.000000E+00*cm*cm/V/s;
    MUP_MAX =  4.000000E+02*cm*cm/V/s;
    NREFP   =  2.750000E+17*std::pow(cm,-3);
    NUP     = -2.100000E+00;
    XIP     =  0.000000E+00;
    ALPHAP  =  3.950000E-01;
    T300    =  300.0*K;
    E0N     =  4.000000E+03*V/cm;
    E0P     =  1.000000E+06*V/cm;
    VSATN_A =  1.130000E+07*cm/s;
    VSATP_A =  1.130000E+07*cm/s;
    VSATN_B =  1.200000E+04*cm/s/K;
    VSATP_B =  1.200000E+04*cm/s/K;
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("MUN_MIN",  PARA("MUN_MIN",    "Minimum low-field electron mobility.", "cm^2 V^-1 s^-1", cm*cm/V/s, &MUN_MIN)) );
    parameter_map.insert(para_item("MUN_MAX",  PARA("MUN_MAX",    "Maximum low-field electron mobility.", "cm^2 V^-1 s^-1", cm*cm/V/s, &MUN_MAX)) );
    parameter_map.insert(para_item("NREFN",    PARA("NREFN",      "Reference doping concentration for low-field electron mobility.", "cm^-3", std::pow(cm,-3), &NREFN)) );
    parameter_map.insert(para_item("NUN",      PARA("NUN",        "Exponent in the direct temperature dependence of low-field electron mobility.", "-", 1.0, &NUN)) );
    parameter_map.insert(para_item("XIN",      PARA("XIN",        "Exponent in the concentration-linked temperature dependence of low-field electron mobility.", "-", 1.0, &XIN)) );
    parameter_map.insert(para_item("ALPHAN",   PARA("ALPHAN",     "Exponent in the concentration dependence of low-field electron mobility.", "-", 1.0, &ALPHAN)) );
    parameter_map.insert(para_item("E0N",      PARA("E0N",        "Reference electric field in the electron velocity saturation model.", "V/cm", V/cm, &E0N)) );
    parameter_map.insert(para_item("VSATN.A",    PARA("VSATN.A",  "", "cm/s",  cm/s,    &VSATN_A)) );
    parameter_map.insert(para_item("VSATN.B",    PARA("VSATN.B",  "", "cm/s/K",cm/s/K , &VSATN_B)) );

    parameter_map.insert(para_item("MUP_MIN",  PARA("MUP_MIN",    "Minimum low-field hole mobility.", "cm^2 V^-1 s^-1", cm*cm/V/s, &MUP_MIN)) );
    parameter_map.insert(para_item("MUP_MAX",  PARA("MUP_MAX",    "Maximum low-field hole mobility.", "cm^2 V^-1 s^-1", cm*cm/V/s, &MUP_MAX)) );
    parameter_map.insert(para_item("NREFP",    PARA("NREFP",      "Reference doping concentration for low-field hole mobility.", "cm^-3", std::pow(cm,-3), &NREFP)) );
    parameter_map.insert(para_item("NUP",      PARA("NUP",        "Exponent in the direct temperature dependence of low-field hole mobility.", "-", 1.0, &NUP)) );
    parameter_map.insert(para_item("XIP",      PARA("XIP",        "Exponent in the concentration-linked temperature dependence of low-field electron mobility.", "-", 1.0, &XIP)) );
    parameter_map.insert(para_item("ALPHAP",   PARA("ALPHAP",     "Exponent in the concentration dependence of low-field hole mobility.", "-", 1.0, &ALPHAP)) );
    parameter_map.insert(para_item("E0P",      PARA("E0P",        "Reference electric field in the hole velocity saturation model.", "V/cm", V/cm, &E0P)) );
    parameter_map.insert(para_item("VSATP.A",    PARA("VSATP.A",  "", "cm/s",  cm/s,    &VSATP_A)) );
    parameter_map.insert(para_item("VSATP.B",    PARA("VSATP.B",  "", "cm/s/K",cm/s/K , &VSATP_B)) );
#endif
  }

private:
  //---------------------------------------------------------------------------
  // Electron low field mobility
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
  // Hole low field mobility
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
  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
    PetscScalar vsat = VSATN_A - VSATN_B*Tl;
    PetscScalar mu0  = ElecMobLowField(Tl);
    PetscScalar E = Ep > 0 ? Ep : 0 ;
    return (mu0+vsat*std::pow(E,3)/std::pow(E0N,4))/(1+std::pow(E/E0N,4));
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    AutoDScalar vsat = VSATN_A - VSATN_B*Tl;
    AutoDScalar mu0  = ElecMobLowField(Tl);
    AutoDScalar E = fmax(Ep, 0.0) ;
    return (mu0+vsat*adtl::pow(E,3)/std::pow(E0N,4))/(1+adtl::pow(E/E0N,4));
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = VSATP_A - VSATP_B*Tl;
    PetscScalar mu0  = HoleMobLowField(Tl);
    PetscScalar E = Ep > 0 ? Ep : 0 ;
    return (mu0+vsat*std::pow(E,3)/std::pow(E0P,4))/(1+std::pow(E/E0P,4));
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    AutoDScalar vsat = VSATP_A - VSATP_B*Tl;
    AutoDScalar mu0  = HoleMobLowField(Tl);
    AutoDScalar E = fmax(Ep, 0.0) ;
    return (mu0+vsat*adtl::pow(E,3)/std::pow(E0P,4))/(1+adtl::pow(E/E0P,4));
  }


// constructor and destructor
public:
  GSS_GaAs_Mob_Analytic(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    Mob_Analytic_Init();
  }

  ~GSS_GaAs_Mob_Analytic()
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
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_GaAs_Mob_Default (const PMIS_Environment& env)
  {
    return new GSS_GaAs_Mob_Analytic(env);
  }
}
/* alias */
extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_GaAs_Mob_Analytic (const PMIS_Environment& env)
  {
    return new GSS_GaAs_Mob_Analytic(env);
  }
}
