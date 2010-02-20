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
// Material Type: AlGaAs


#include "PMI.h"

class GSS_AlGaAs_Mob_Analytic : public PMIS_Mobility
{
private:
  // parameters for Analytic mobility
  PetscScalar MUN_MIN;
  PetscScalar MIN_X1;
  PetscScalar MIN_X2;
  PetscScalar MUN_MAX;
  PetscScalar MAN_X1;
  PetscScalar MAN_X2;
  PetscScalar NREFN;
  PetscScalar NREFN2;
  PetscScalar NUN;
  PetscScalar XIN;
  PetscScalar ALPHAN;
  PetscScalar MUP_MIN;
  PetscScalar MIP_X1;
  PetscScalar MIP_X2;
  PetscScalar MUP_MAX;
  PetscScalar MAP_X1;
  PetscScalar MAP_X2;
  PetscScalar NREFP;
  PetscScalar NREFP2;
  PetscScalar NUP;
  PetscScalar XIP;
  PetscScalar ALPHAP;
  PetscScalar T300;
  // parameters for high field modification
  PetscScalar VSATN;
  PetscScalar VSN1;
  PetscScalar VSN2;
  PetscScalar VSATP;
  PetscScalar VSP1;
  PetscScalar VSP2;
  PetscScalar E0N;
  PetscScalar EN1;
  PetscScalar EN2;
  PetscScalar E0P;

  void Mob_Analytic_Init()
  {
    MUN_MIN =  2.366000E+03*cm*cm/V/s;
    MIN_X1  = -9.831000E-01;
    MIN_X2  =  0.000000E+00;
    MUN_MAX =  9.891700E+03*cm*cm/V/s;
    MAN_X1  = -9.495000E-01;
    MAN_X2  =  0.000000E+00;
    NREFN   =  3.629700E+17*std::pow(cm,-3);
    NREFN2  =  1.745100E+18*std::pow(cm,-3);
    NUN     =  0.000000E+00;
    XIN     =  0.000000E+00;
    ALPHAN  =  1.000000E+00;

    MUP_MIN =  0.000000E+00*cm*cm/V/s;
    MIP_X1  =  0.000000E+00;
    MIP_X2  =  0.000000E+00;
    MUP_MAX =  4.000000E+02*cm*cm/V/s;
    MAP_X1  =  0.000000E+00;
    MAP_X2  =  0.000000E+00;
    NREFP   =  2.750000E+17*std::pow(cm,-3);
    NREFP2  =  1.000000E+30*std::pow(cm,-3);
    NUP     = -2.100000E+00;
    XIP     =  0.000000E+00;
    ALPHAP  =  3.950000E-01;
    T300    =  300.0*K;

    VSATN   =  6.351800E+06*cm/s;
    VSN1    = -5.304000E-01;
    VSN2    = -7.484000E-02;
    VSATP   =  0.000000E+00*cm/s;
    VSP1    =  0.000000E+00;
    VSP2    =  0.000000E+00;
    E0N     =  5.418400E+03*V/cm;
    EN1     = -2.471000E+00;
    EN2     =  7.194200E+00;
    E0P     =  1.000000E+06*V/cm;

#ifdef __CALIBRATE__
   parameter_map.insert(para_item("MUN.MIN",    PARA("MUN.MIN",    "", "cm*cm/V/s",cm*cm/V/s , &MUN_MIN)) );
   parameter_map.insert(para_item("MIN.X1", PARA("MIN.X1", "", "-", 1.0, &MIN_X1)) );
   parameter_map.insert(para_item("MIN.X2", PARA("MIN.X2", "", "-", 1.0, &MIN_X2)) );
	 parameter_map.insert(para_item("MUN.MAX",    PARA("MUN.MAX",    "", "cm*cm/V/s", cm*cm/V/s, &MUN_MAX)) );
	 parameter_map.insert(para_item("MAN.X1", PARA("MAN.X1", "", "-", 1.0, &MAN_X1)) );
   parameter_map.insert(para_item("MAN.X2", PARA("MAN.X2", "", "-", 1.0, &MAN_X2)) );
	 parameter_map.insert(para_item("NREFN",    PARA("NREFN",    "", "cm^-3", std::pow(cm,-3), &NREFN)) );
	 parameter_map.insert(para_item("NREFN2",    PARA("NREFN2",    "", "cm^-3", std::pow(cm,-3), &NREFN2)) );
	 parameter_map.insert(para_item("NUN",    PARA("NUN",    "", "-", 1.0, &NUN)) );
	 parameter_map.insert(para_item("XIN",    PARA("XIN",    "", "-", 1.0, &XIN)) );
	 parameter_map.insert(para_item("ALPHAN",    PARA("ALPHAN",    "", "-", 1.0, &ALPHAN)) );

	 parameter_map.insert(para_item("MUP.MIN",    PARA("MUP.MIN",    "", "cm*cm/V/s", cm*cm/V/s, &MUP_MIN)) );
	 parameter_map.insert(para_item("MIP.X1", PARA("MIP.X1", "", "-", 1.0, &MIP_X1)) );
   parameter_map.insert(para_item("MIP.X2", PARA("MIP.X2", "", "-", 1.0, &MIP_X2)) );
	 parameter_map.insert(para_item("MUP.MAX",    PARA("MUP.MAX",    "", "cm*cm/V/s",cm*cm/V/s , &MUP_MAX)) );
	 parameter_map.insert(para_item("MAP.X1", PARA("MAP.X1", "", "-", 1.0, &MAP_X1)) );
   parameter_map.insert(para_item("MAP.X2", PARA("MAP.X2", "", "-", 1.0, &MAP_X2)) );
	 parameter_map.insert(para_item("NREFP",    PARA("NREFP",    "", "cm^-3",std::pow(cm,-3) , &NREFP)) );
	 parameter_map.insert(para_item("NREFP2",    PARA("NREFP2",    "", "cm^-3",std::pow(cm,-3) , &NREFP2)) );
	 parameter_map.insert(para_item("NUP",    PARA("NUP",    "", "-",1.0 , &NUP)) );
	 parameter_map.insert(para_item("XIP",    PARA("XIP",    "", "-",1.0 , &XIP)) );
	 parameter_map.insert(para_item("ALPHAP",    PARA("ALPHAP",    "", "-", 1.0, &ALPHAP)) );

	 parameter_map.insert(para_item("VSATN",    PARA("VSATN",    "", "cm/s", cm/s, &VSATN)) );
	 parameter_map.insert(para_item("VSN1",    PARA("VSN1",    "", "-",1.0 , &VSN1)) );
	 parameter_map.insert(para_item("VSN2",    PARA("VSN2",    "", "-",1.0 , &VSN2)) );
	 parameter_map.insert(para_item("VSATP",    PARA("VSATP",    "", "cm/s", cm/s, &VSATP)) );
	 parameter_map.insert(para_item("VSP1",    PARA("VSP1",    "", "-",1.0 , &VSP1)) );
	 parameter_map.insert(para_item("VSP2",    PARA("VSP2",    "", "-",1.0 , &VSP2)) );
    parameter_map.insert(para_item("E0N", PARA("E0N", "", "V/cm", V/cm, &E0N)) );
    parameter_map.insert(para_item("EN1", PARA("EN1", "", "-",1.0 , &EN1)) );
    parameter_map.insert(para_item("EN2", PARA("EN2", "", "-", 1.0, &EN2)) );
    parameter_map.insert(para_item("E0P", PARA("E0P", "", "V/cm",V/cm , &E0P)) );
#endif
  }

public:
  //---------------------------------------------------------------------------
  // Electron low field mobility
  PetscScalar ElecMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar x = ReadxMoleFraction();
    PetscScalar mu_min = MUN_MIN*(1+MIN_X1*x+MIN_X2*x*x);
    PetscScalar mu_max = MUN_MAX*(1+MAN_X1*x+MAN_X2*x*x);
    return mu_min+(mu_max*std::pow(Tl/T300,NUN)-mu_min)/ \
           (1+std::pow(Tl/T300,XIN)*(std::pow((Na+Nd)/NREFN,ALPHAN)+std::pow((Na+Nd)/NREFN2,3)));
  }
  AutoDScalar ElecMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar x = ReadxMoleFraction();
    PetscScalar mu_min = MUN_MIN*(1+MIN_X1*x+MIN_X2*x*x);
    PetscScalar mu_max = MUN_MAX*(1+MAN_X1*x+MAN_X2*x*x);
    return mu_min+(mu_max*adtl::pow(Tl/T300,NUN)-mu_min)/ \
           (1+adtl::pow(Tl/T300,XIN)*(std::pow((Na+Nd)/NREFN,ALPHAN)+std::pow((Na+Nd)/NREFN2,3)));
  }

  //---------------------------------------------------------------------------
  // Hole low field mobility
  PetscScalar HoleMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar x = ReadxMoleFraction();
    PetscScalar mu_min = MUP_MIN*(1+MIP_X1*x+MIP_X2*x*x);
    PetscScalar mu_max = MUP_MAX*(1+MAP_X1*x+MAP_X2*x*x);
    return mu_min+(mu_max*std::pow(Tl/T300,NUP)-mu_min)/ \
           (1+std::pow(Tl/T300,XIP)*(std::pow((Na+Nd)/NREFP,ALPHAP)+std::pow((Na+Nd)/NREFP2,3)));
  }
  AutoDScalar HoleMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar x = ReadxMoleFraction();
    PetscScalar mu_min = MUP_MIN*(1+MIP_X1*x+MIP_X2*x*x);
    PetscScalar mu_max = MUP_MAX*(1+MAP_X1*x+MAP_X2*x*x);
    return mu_min+(mu_max*adtl::pow(Tl/T300,NUP)-mu_min)/ \
           (1+adtl::pow(Tl/T300,XIP)*(std::pow((Na+Nd)/NREFP,ALPHAP)+std::pow((Na+Nd)/NREFP2,3)));
  }


public:
  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
    PetscScalar x = ReadxMoleFraction();
    PetscScalar vsat = VSATN*(1+ VSN1*x + VSN2*x*x);
    PetscScalar E0   = E0N*(1+EN1*x+EN2*x*x);
    PetscScalar mu0  = ElecMobLowField(Tl);
    return (mu0+vsat*std::pow(Ep,3)/std::pow(E0,4))/(1+std::pow(Ep/E0,4));
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    PetscScalar x = ReadxMoleFraction();
    PetscScalar vsat = VSATN*(1+ VSN1*x + VSN2*x*x);
    PetscScalar E0   = E0N*(1+EN1*x+EN2*x*x);
    AutoDScalar mu0  = ElecMobLowField(Tl);
    return (mu0+vsat*adtl::pow(Ep,3)/std::pow(E0,4))/(1+adtl::pow(Ep/E0,4));
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    return HoleMobLowField(Tl);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    return HoleMobLowField(Tl);
  }


// constructor and destructor
public:
  GSS_AlGaAs_Mob_Analytic(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    Mob_Analytic_Init();
  }


  ~GSS_AlGaAs_Mob_Analytic()
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
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_AlGaAs_Mob_Default (const PMIS_Environment& env)
  {
    return new GSS_AlGaAs_Mob_Analytic(env);
  }
}
/* alias */
extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_AlGaAs_Mob_Analytic (const PMIS_Environment& env)
  {
    return new GSS_AlGaAs_Mob_Analytic(env);
  }
}
