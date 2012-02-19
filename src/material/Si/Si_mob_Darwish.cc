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

//source:
//  An Improved Electron and Hole Mobility Model for General Purpose Device Simulation
//  Mohamed N. Darwish, IEEE TRANSACTIONS ON ELECTRON DEVICES, VOL. 44, NO. 9, SEPTEMBER 1997

#include "PMI.h"

class GSS_Si_Mob_Darwish : public PMIS_Mobility
{
private:

  // constant
  PetscScalar T300;
  PetscScalar E1;
  PetscScalar N1;

  // parameters for Philips low field mobility
  PetscScalar MMNN_UM;
  PetscScalar MMXN_UM;
  PetscScalar NRFN_UM;
  PetscScalar ALPN_UM;
  PetscScalar TETN_UM;
  PetscScalar NRFD_UM;
  PetscScalar CRFD_UM;

  PetscScalar MMNP_UM;
  PetscScalar MMXP_UM;
  PetscScalar NRFP_UM;
  PetscScalar ALPP_UM;
  PetscScalar TETP_UM;
  PetscScalar NRFA_UM;
  PetscScalar CRFA_UM;
  PetscScalar NSC_REF;
  PetscScalar CAR_REF;
  PetscScalar me_over_m0;
  PetscScalar mh_over_m0;
  PetscScalar me_over_mh;

  //
  std::vector<PetscScalar>  Pn_temperature_table;
  std::vector<PetscScalar>  Pp_temperature_table;

  PetscScalar Gn(PetscScalar P, PetscScalar Tl)
  {
    return 1-0.89233/std::pow(0.41372+P*std::pow(Tl/T300/me_over_m0,0.28227),0.19778)+0.005978/std::pow(P*std::pow(T300/Tl*me_over_m0,0.72169),1.80618);
  }


  PetscScalar Gp(PetscScalar P, PetscScalar Tl)
  {
    return 1-0.89233/std::pow(0.41372+P*std::pow(Tl/T300/mh_over_m0,0.28227),0.19778)+0.005978/std::pow(P*std::pow(T300/Tl*mh_over_m0,0.72169),1.80618);
  }

  PetscScalar Pn_Limiter_Init(PetscScalar Tl)
  {
    PetscScalar Pmin=0.01, Pmax=1.0;
    do
    {
      PetscScalar length = (Pmax-Pmin);
      PetscScalar P1=Pmin + 0.33*length;
      PetscScalar P2=Pmax - 0.33*length;
      if( Gn(P1, Tl) < Gn(P2, Tl) )
      {
        Pmax = P2;
      }
      else
      {
        Pmin = P1;
      }
    }
    while(Pmax-Pmin>1e-3);

    return 0.5*(Pmin + Pmax);
  }

  PetscScalar Pp_Limiter_Init(PetscScalar Tl)
  {
    PetscScalar Pmin=0.01, Pmax=1.0;
    do
    {
      PetscScalar length = (Pmax-Pmin);
      PetscScalar P1=Pmin + 0.33*length;
      PetscScalar P2=Pmax - 0.33*length;
      if( Gp(P1, Tl) < Gp(P2, Tl) )
      {
        Pmax = P2;
      }
      else
      {
        Pmin = P1;
      }
    }
    while(Pmax-Pmin>1e-3);

    return 0.5*(Pmin + Pmax);
  }


  PetscScalar Pn_Limiter(PetscScalar Tl) const
  {
    int index = static_cast<int>(Tl/(10*K))-1;
    if(index<0) index=0;
    if(index>Pn_temperature_table.size()-1) index=Pn_temperature_table.size()-1;
    return Pn_temperature_table[index];
  }

  PetscScalar Pp_Limiter(PetscScalar Tl) const
  {
    int index = static_cast<int>(Tl/(10*K))-1;
    if(index<0) index=0;
    if(index>Pp_temperature_table.size()-1) index=Pp_temperature_table.size()-1;
    return Pp_temperature_table[index];
  }

  void Mob_Philips_Init()
  {
    MMNN_UM = 5.220000E+01*cm*cm/V/s;
    MMXN_UM = 1.417000E+03*cm*cm/V/s;
    NRFN_UM = 9.680000E+16*std::pow(cm,-3);
    ALPN_UM = 6.800000E-01;
    TETN_UM = 2.285000E+00;
    NRFD_UM = 4.000000E+20*std::pow(cm,-3);
    CRFD_UM = 2.100000E-01;

    MMNP_UM = 4.490000E+01*cm*cm/V/s;
    MMXP_UM = 4.705000E+02*cm*cm/V/s;
    NRFP_UM = 2.230000E+17*std::pow(cm,-3);
    ALPP_UM = 7.190000E-01;
    TETP_UM = 2.247000E+00;
    NRFA_UM = 7.200000E+20*std::pow(cm,-3);
    CRFA_UM = 5.000000E-01;

    NSC_REF = 3.97e13*std::pow(cm,-2);
    CAR_REF = 1.36e20*std::pow(cm,-3);
    me_over_m0 = 1.0;
    mh_over_m0 = 1.258;
    me_over_mh = 1.0/1.258;

    for(PetscScalar T=10*K; T<1000*K; T+=10*K)
    {
      Pn_temperature_table.push_back(Pn_Limiter_Init(T));
      Pp_temperature_table.push_back(Pp_Limiter_Init(T));
    }

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("MMNN.UM",    PARA("MMNN.UM",    "", "cm*cm/V/s",cm*cm/V/s , &MMNN_UM)) );
    parameter_map.insert(para_item("MMXN.UM",    PARA("MMXN.UM",    "", "cm*cm/V/s", cm*cm/V/s, &MMXN_UM)) );
    parameter_map.insert(para_item("NRFN.UM",    PARA("NRFN.UM",    "", "std::pow(cm,-3)", std::pow(cm,-3), &NRFN_UM)) );
    parameter_map.insert(para_item("ALPN.UM",    PARA("ALPN.UM",    "", "-", 1.0, &ALPN_UM)) );
    parameter_map.insert(para_item("TETN.UM",    PARA("TETN.UM",    "", "-",1.0 , &TETN_UM)) );
    parameter_map.insert(para_item("NRFD.UM",    PARA("NRFD.UM",    "", "std::pow(cm,-3)",std::pow(cm,-3) , &NRFD_UM)) );
    parameter_map.insert(para_item("CRFD.UM",    PARA("CRFD.UM",    "", "-", 1.0, &CRFD_UM)) );
    parameter_map.insert(para_item("MMNP.UM",    PARA("MMNP.UM",    "", "cm*cm/V/s", cm*cm/V/s, &MMNP_UM)) );
    parameter_map.insert(para_item("MMXP.UM",    PARA("MMXP.UM",    "", "cm*cm/V/s", cm*cm/V/s, &MMXP_UM)) );
    parameter_map.insert(para_item("NRFP.UM",    PARA("NRFP.UM",    "", "std::pow(cm,-3)",std::pow(cm,-3) , &NRFP_UM)) );
    parameter_map.insert(para_item("ALPP.UM",    PARA("ALPP.UM",    "", "-",1.0 , &ALPP_UM)) );
    parameter_map.insert(para_item("TETP.UM",    PARA("TETP.UM",    "", "-", 1.0, &TETP_UM)) );
    parameter_map.insert(para_item("NRFA.UM",    PARA("NRFA.UM",    "", "std::pow(cm,-3)",std::pow(cm,-3) , &NRFA_UM)) );
    parameter_map.insert(para_item("CRFA.UM",    PARA("CRFA.UM",    "", "-",1.0 , &CRFA_UM)) );
    parameter_map.insert(para_item("NSC.REF",    PARA("NSC.REF",    "", "std::pow(cm,-2)", std::pow(cm,-2), &NSC_REF)) );
    parameter_map.insert(para_item("CAR.REF",    PARA("CAR.REF",    "", "std::pow(cm,-2)",std::pow(cm,-3) , &CAR_REF)) );
    parameter_map.insert(para_item("me_over_m0",    PARA("me_over_m0",    "", "-", 1.0, &me_over_m0)) );
    parameter_map.insert(para_item("mh_over_m0",    PARA("mh_over_m0",    "", "-",1.0 , &mh_over_m0)) );
    parameter_map.insert(para_item("me_over_mh",    PARA("me_over_mh",    "", "-", 1.0, &me_over_mh)) );
#endif

  }

  // parameter for surface mobility
  PetscScalar EXN4_LSM;
  PetscScalar EXN5_LSM;
  PetscScalar EXN8_LSM;

  PetscScalar AN_LSM;
  PetscScalar BN_LSM;
  PetscScalar CN_LSM;
  PetscScalar DN_LSM;
  PetscScalar EN_LSM;

  PetscScalar EXP4_LSM;
  PetscScalar EXP5_LSM;
  PetscScalar EXP8_LSM;

  PetscScalar AP_LSM;
  PetscScalar BP_LSM;
  PetscScalar CP_LSM;
  PetscScalar DP_LSM;
  PetscScalar EP_LSM;

  void Mob_Darwish_Init()
  {
    EXN4_LSM=2.330000E-02;
    EXN5_LSM=1.700000E+00;
    EXN8_LSM=2.580000E+00;

    AN_LSM=6.850000E-21*std::pow(cm,3);
    BN_LSM=3.610000E+07*cm/s;
    CN_LSM=1.700000E+04*cm/s*std::pow(V/cm,PetscScalar(-2.0/3.0));
    DN_LSM=3.580000E+18*cm*cm/V/s;
    EN_LSM=7.670000E-2;

    EXP4_LSM=1.190000E-02;
    EXP5_LSM=9.000000E-01;
    EXP8_LSM=2.180000E+00;

    AP_LSM=7.820000E-21*std::pow(cm,3);
    BP_LSM=1.510000E+07*cm/s;
    CP_LSM=4.180000E+03*cm/s*std::pow(V/cm,PetscScalar(-2.0/3.0));
    DP_LSM=4.100000E+15*cm*cm/V/s;
    EP_LSM=1.230000E-1;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("EXN4.LSM",    PARA("EXN4.LSM",    "", "-", 1.0, &EXN4_LSM)) );
    parameter_map.insert(para_item("EXN8.LSM",    PARA("EXN8.LSM",    "", "-", 1.0, &EXN8_LSM)) );
    parameter_map.insert(para_item("AN.LSM",    PARA("AN.LSM",    "", "cm^-3", std::pow(cm,-3), &AN_LSM)) );
    parameter_map.insert(para_item("BN.LSM",    PARA("BN.LSM",    "", "cm/s", cm/s, &BN_LSM)) );
    parameter_map.insert(para_item("CN.LSM",    PARA("CN.LSM",    "", "cm/s*pow(V/cm,(-2.0/3.0))", cm/s*std::pow(V/cm,PetscScalar(-2.0/3.0)), &CN_LSM)) );
    parameter_map.insert(para_item("DN.LSM",    PARA("DN.LSM",    "", "cm*cm/V/s", cm*cm/V/s, &DN_LSM)) );
    parameter_map.insert(para_item("EXP4.LSM",    PARA("EXP4.LSM",    "", "-",1.0 , &EXP4_LSM)) );
    parameter_map.insert(para_item("EXP8.LSM",    PARA("EXP8.LSM",    "", "-",1.0 , &EXP8_LSM)) );
    parameter_map.insert(para_item("AP.LSM",    PARA("AP.LSM",    "", "cm^-3", std::pow(cm,-3), &AP_LSM)) );
    parameter_map.insert(para_item("BP.LSM",    PARA("BP.LSM",    "", "cm/s", cm/s, &BP_LSM)) );
    parameter_map.insert(para_item("CP.LSM",    PARA("CP.LSM",    "", "cm/s*pow(V/cm,(-2.0/3.0))", cm/s*std::pow(V/cm,PetscScalar(-2.0/3.0)), &CP_LSM)) );
    parameter_map.insert(para_item("DP.LSM",    PARA("DP.LSM",    "", "cm*cm/V/s",cm*cm/V/s , &DP_LSM)) );
#endif

  }

  // parameters for parallel field modification
  PetscScalar BETAN;
  PetscScalar BETAP;
  PetscScalar VSATN0;
  PetscScalar VSATP0;
  PetscScalar VSATN_A;
  PetscScalar VSATP_A;

  void Mob_Parallel_Init()
  {
    BETAN = 2.000000E+00;
    BETAP = 1.000000E+00;
    VSATN0  =  2.400000E7*cm/s;
    VSATN_A =  0.8;
    VSATP0  =  2.400000E7*cm/s;
    VSATP_A =  0.8;
  }

  //---------------------------------------------------------------------------
  // Electron low field mobility
  PetscScalar ElecMobPhilips(const PetscScalar &p,const PetscScalar &n,const PetscScalar &Tl) const
  {
    PetscScalar mu_lattice = MMXN_UM*std::pow(Tl/T300,-TETN_UM);
    PetscScalar mu1 = MMXN_UM*MMXN_UM/(MMXN_UM-MMNN_UM)*std::pow(Tl/T300,3*ALPN_UM-1.5);
    PetscScalar mu2 = MMXN_UM*MMNN_UM/(MMXN_UM-MMNN_UM)*sqrt(T300/Tl);
    PetscScalar Na  = ReadDopingNa()+1e0*std::pow(cm,-3);
    PetscScalar Nd  = ReadDopingNd()+1e0*std::pow(cm,-3);
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    PetscScalar Nsc = Nds+Nas+fabs(p);

    PetscScalar P   = 1.0/(2.459/(NSC_REF/std::pow(Nsc,PetscScalar(2.0/3.0)))+3.828/(CAR_REF/fabs(n+p)*me_over_m0))*(Tl/T300)*(Tl/T300);
    PetscScalar pp1 = std::pow(P,PetscScalar(0.6478));
    PetscScalar F   = (0.7643*pp1+2.2999+6.5502*me_over_mh)/(pp1+2.3670-0.8552*me_over_mh);
    PetscScalar Pl  = Pn_Limiter(Tl);
    PetscScalar PG  = std::max(P, Pl);
    PetscScalar G   = 1-0.89233/std::pow(0.41372+PG*std::pow(Tl/T300/me_over_m0,0.28227),0.19778)+0.005978/std::pow(PG*std::pow(T300/Tl*me_over_m0, 0.72169),1.80618);
    //PetscScalar G   = 1-4.41804/std::pow(39.9014+P*std::pow(Tl/T300/me_over_m0,PetscScalar(0.0001)),PetscScalar(0.38297))+0.52896/std::pow(P*std::pow(T300/Tl*me_over_m0,PetscScalar(1.595787)),PetscScalar(0.25948));
    PetscScalar Nsce = Nds+Nas*G+fabs(p)/F;
    PetscScalar mu_scatt = mu1*(Nsc/Nsce)*std::pow(NRFN_UM/Nsc,ALPN_UM)+mu2*(fabs(n+p)/Nsce);
    return 1.0/(1.0/mu_lattice+1.0/mu_scatt);
  }
  AutoDScalar ElecMobPhilips(const AutoDScalar &p,const AutoDScalar &n,const AutoDScalar &Tl) const
  {
    AutoDScalar mu_lattice = MMXN_UM*adtl::pow(Tl/T300,-TETN_UM);
    AutoDScalar mu1 = MMXN_UM*MMXN_UM/(MMXN_UM-MMNN_UM)*adtl::pow(Tl/T300,3*ALPN_UM-1.5);
    AutoDScalar mu2 = MMXN_UM*MMNN_UM/(MMXN_UM-MMNN_UM)*sqrt(T300/Tl);
    PetscScalar Na  = ReadDopingNa()+1e0*std::pow(cm,-3);
    PetscScalar Nd  = ReadDopingNd()+1e0*std::pow(cm,-3);
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    AutoDScalar Nsc = Nds+Nas+fabs(p);

    AutoDScalar P   = 1.0/(2.459/(NSC_REF/adtl::pow(Nsc,PetscScalar(2.0/3.0)))+3.828/(CAR_REF/fabs(n+p)*me_over_m0))*(Tl/T300)*(Tl/T300);
    AutoDScalar pp1 = adtl::pow(P,PetscScalar(0.6478));
    AutoDScalar F   = (0.7643*pp1+2.2999+6.5502*me_over_mh)/(pp1+2.3670-0.8552*me_over_mh);
    PetscScalar Pl  = Pn_Limiter(Tl.getValue());
    AutoDScalar PG  = P.getValue() < Pl ? AutoDScalar(Pl) : P;
    AutoDScalar G   = 1-0.89233/adtl::pow(0.41372+PG*adtl::pow(Tl/T300/me_over_m0,0.28227),0.19778)+0.005978/adtl::pow(PG*adtl::pow(T300/Tl*me_over_m0,0.72169),1.80618);
    //AutoDScalar G   = 1-4.41804/adtl::pow(39.9014+P*adtl::pow(Tl/T300/me_over_m0,PetscScalar(0.0001)),PetscScalar(0.38297))+0.52896/adtl::pow(P*adtl::pow(T300/Tl*me_over_m0,PetscScalar(1.595787)),PetscScalar(0.25948));
    AutoDScalar Nsce = Nds+Nas*G+fabs(p)/F;
    AutoDScalar mu_scatt = mu1*(Nsc/Nsce)*adtl::pow(NRFN_UM/Nsc,ALPN_UM)+mu2*(fabs(n+p)/Nsce);
    return 1.0/(1.0/mu_lattice+1.0/mu_scatt);
  }

  //---------------------------------------------------------------------------
  // Hole low field mobility
  PetscScalar HoleMobPhilips(const PetscScalar &p,const PetscScalar &n,const PetscScalar &Tl) const
  {
    PetscScalar mu_lattice = MMXP_UM*std::pow(Tl/T300,-TETP_UM);
    PetscScalar mu1 = MMXP_UM*MMXP_UM/(MMXP_UM-MMNP_UM)*std::pow(Tl/T300,3*ALPP_UM-1.5);
    PetscScalar mu2 = MMXP_UM*MMNP_UM/(MMXP_UM-MMNP_UM)*sqrt(T300/Tl);
    PetscScalar Na  = ReadDopingNa()+1e0*std::pow(cm,-3);
    PetscScalar Nd  = ReadDopingNd()+1e0*std::pow(cm,-3);
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    PetscScalar Nsc = Nds+Nas+fabs(n);

    PetscScalar P   = 1.0/(2.459/(NSC_REF/std::pow(Nsc,PetscScalar(2.0/3.0)))+3.828/(CAR_REF/fabs(n+p)*mh_over_m0))*(Tl/T300)*(Tl/T300);
    PetscScalar pp1 = std::pow(P,PetscScalar(0.6478));
    PetscScalar F   = (0.7643*pp1+2.2999+6.5502/me_over_mh)/(pp1+2.3670-0.8552/me_over_mh);
    PetscScalar Pl  = Pp_Limiter(Tl);
    PetscScalar PG  = std::max(P, Pl);
    PetscScalar G = 1-0.89233/std::pow(0.41372+PG*std::pow(Tl/T300/mh_over_m0,0.28227),0.19778)+0.005978/std::pow(PG*std::pow(T300/Tl*mh_over_m0,0.72169),1.80618);
    //PetscScalar G   = 1-4.41804/std::pow(39.9014+P*std::pow(Tl/T300/mh_over_m0,PetscScalar(0.0001)),PetscScalar(0.38297))+0.52896/std::pow(P*std::pow(T300/Tl*mh_over_m0,PetscScalar(1.595787)),PetscScalar(0.25948));
    PetscScalar Nsce = Nas+Nds*G+fabs(n)/F;
    PetscScalar mu_scatt = mu1*(Nsc/Nsce)*std::pow(NRFP_UM/Nsc,ALPP_UM)+mu2*(fabs(n+p)/Nsce);
    return 1.0/(1.0/mu_lattice+1.0/mu_scatt);
  }
  AutoDScalar HoleMobPhilips(const AutoDScalar &p,const AutoDScalar &n,const AutoDScalar &Tl) const
  {
    AutoDScalar mu_lattice = MMXP_UM*adtl::pow(Tl/T300,-TETP_UM);
    AutoDScalar mu1 = MMXP_UM*MMXP_UM/(MMXP_UM-MMNP_UM)*adtl::pow(Tl/T300,3*ALPP_UM-1.5);
    AutoDScalar mu2 = MMXP_UM*MMNP_UM/(MMXP_UM-MMNP_UM)*sqrt(T300/Tl);
    PetscScalar Na  = ReadDopingNa()+1e0*std::pow(cm,-3);
    PetscScalar Nd  = ReadDopingNd()+1e0*std::pow(cm,-3);
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    AutoDScalar Nsc = Nds+Nas+fabs(n);

    AutoDScalar P   = 1.0/(2.459/(NSC_REF/adtl::pow(Nsc,PetscScalar(2.0/3.0)))+3.828/(CAR_REF/fabs(n+p)*mh_over_m0))*(Tl/T300)*(Tl/T300);
    AutoDScalar pp1 = adtl::pow(P,PetscScalar(0.6478));
    AutoDScalar F   = (0.7643*pp1+2.2999+6.5502/me_over_mh)/(pp1+2.3670-0.8552/me_over_mh);
    PetscScalar Pl  = Pp_Limiter(Tl.getValue());
    AutoDScalar PG  = P.getValue() < Pl ? AutoDScalar(Pl) : P;
    AutoDScalar G = 1-0.89233/adtl::pow(0.41372+PG*adtl::pow(Tl/T300/mh_over_m0,0.28227),0.19778)+0.005978/adtl::pow(PG*adtl::pow(T300/Tl*mh_over_m0,0.72169),1.80618);
    //AutoDScalar G   = 1-4.41804/adtl::pow(39.9014+P*adtl::pow(Tl/T300/mh_over_m0,PetscScalar(0.0001)),PetscScalar(0.38297))+0.52896/adtl::pow(P*adtl::pow(T300/Tl*mh_over_m0,PetscScalar(1.595787)),PetscScalar(0.25948));
    AutoDScalar Nsce = Nas+Nds*G+fabs(n)/F;
    AutoDScalar mu_scatt = mu1*(Nsc/Nsce)*adtl::pow(NRFP_UM/Nsc,ALPP_UM)+mu2*(fabs(n+p)/Nsce);
    return 1.0/(1.0/mu_lattice+1.0/mu_scatt);
  }

  //---------------------------------------------------------------------------
  // Electron surface mobility, acoustical phono scattering and roughness scattering
  PetscScalar ElecMobSurface(const PetscScalar &c, const PetscScalar &Tl,const PetscScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    PetscScalar ET = Et+1.0*V/cm;
    PetscScalar mu_ac = BN_LSM/ET + CN_LSM*std::pow(N_total/N1,EXN4_LSM)/std::pow(Tl/T300, EXN5_LSM)*std::pow(ET,PetscScalar(-1.0/3.0));
    PetscScalar mu_sr = DN_LSM*std::pow(ET/E1,-(EXN8_LSM + AN_LSM*c/std::pow(N_total/N1,EN_LSM)));
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }
  AutoDScalar ElecMobSurface(const AutoDScalar &c, const AutoDScalar &Tl,const AutoDScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    AutoDScalar ET = Et+1.0*V/cm;
    AutoDScalar mu_ac = BN_LSM/ET + CN_LSM*std::pow(N_total/N1,EXN4_LSM)/adtl::pow(Tl/T300, EXN5_LSM)*adtl::pow(ET,PetscScalar(-1.0/3.0));
    AutoDScalar mu_sr = DN_LSM*adtl::pow(ET/E1,-(EXN8_LSM + AN_LSM*c/std::pow(N_total/N1,EN_LSM)));
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }

  //---------------------------------------------------------------------------
  // Hole surface mobility, acoustical phono scattering and roughness scattering
  PetscScalar HoleMobSurface(const PetscScalar &c, const PetscScalar &Tl,const PetscScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    PetscScalar ET = Et+1.0*V/cm;
    PetscScalar mu_ac = BP_LSM/ET + CP_LSM*std::pow(N_total/N1,EXP4_LSM)/std::pow(Tl/T300, EXP5_LSM)*std::pow(ET,PetscScalar(-1.0/3.0));
    PetscScalar mu_sr = DP_LSM*std::pow(ET/E1,-(EXP8_LSM + AP_LSM*c/std::pow(N_total/N1,EP_LSM)));
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }
  AutoDScalar HoleMobSurface(const AutoDScalar &c, const AutoDScalar &Tl,const AutoDScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*std::pow(cm,-3);
    AutoDScalar ET = Et+1.0*V/cm;
    AutoDScalar mu_ac = BP_LSM/ET + CP_LSM*std::pow(N_total/N1,EXP4_LSM)/adtl::pow(Tl/T300, EXP5_LSM)*adtl::pow(ET+1.0*V/cm,PetscScalar(-1.0/3.0));
    AutoDScalar mu_sr = DP_LSM*adtl::pow(ET/E1,-(EXP8_LSM + AP_LSM*c/std::pow(N_total/N1,EP_LSM)));
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
    PetscScalar mu0  = 1.0/(1.0/ElecMobPhilips(p, n, Tl)+1.0/ElecMobSurface(n+p,Tl,Et));
    return mu0/std::pow(1+std::pow(mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN);
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    AutoDScalar vsat = VSATN0/(1+VSATN_A*exp(Tl/(2*T300)));
    AutoDScalar mu0  = 1.0/(1.0/ElecMobPhilips(p, n, Tl)+1.0/ElecMobSurface(n+p,Tl,Et));
    return mu0/adtl::pow(1+adtl::pow(mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN);
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = VSATP0/(1+VSATP_A*exp(Tl/(2*T300)));
    PetscScalar mu0  = 1.0/(1.0/HoleMobPhilips(p,n,Tl)+1.0/HoleMobSurface(n+p,Tl,Et));
    return mu0/std::pow(1+std::pow(mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    AutoDScalar vsat = VSATP0/(1+VSATP_A*exp(Tl/(2*T300)));
    AutoDScalar mu0  = 1.0/(1.0/HoleMobPhilips(p,n,Tl)+1.0/HoleMobSurface(n+p,Tl,Et));
    return mu0/adtl::pow(1+adtl::pow(mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP);
  }

  // constructor
public:
  GSS_Si_Mob_Darwish(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    PMI_Info = "This is the Darwish mobility model of Silicon";

    // init constant
    T300  = 300.0*K;
    E1    = 1*V/cm;
    N1    = 1.0*std::pow(cm,-3);
    Mob_Philips_Init();
    Mob_Darwish_Init();
    Mob_Parallel_Init();
  }

  ~GSS_Si_Mob_Darwish(){}
}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  and it setup Lombardi mobility model
 */
extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_Mobility* PMIS_Si_Mob_Darwish (const PMIS_Environment& env)
  {
    return new GSS_Si_Mob_Darwish(env);
  }
}
