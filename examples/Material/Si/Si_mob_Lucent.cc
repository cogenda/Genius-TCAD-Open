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

class GSS_Si_Mob_Lucent : public PMIS_Mobility
{
private:
  // parameters for Lucent mobility
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


  // temperature
  PetscScalar T300;
  // 1V/cm
  PetscScalar E1;
  // 1 cm^-3
  PetscScalar N1;

  // parameters for transvers field modification
  PetscScalar AN_LUC;
  PetscScalar AP_LUC;
  PetscScalar BN_LUC;
  PetscScalar BP_LUC;
  PetscScalar CN_LUC;
  PetscScalar CP_LUC;
  PetscScalar DN_LUC;
  PetscScalar DP_LUC;
  PetscScalar FN_LUC;
  PetscScalar FP_LUC;
  PetscScalar KN_LUC;
  PetscScalar KP_LUC;
  PetscScalar EXN4_LUC;
  PetscScalar EXP4_LUC;
  PetscScalar EXN9_LUC;
  PetscScalar EXP9_LUC;
  // parameters for parallel field modification
  PetscScalar BETAN;
  PetscScalar BETAP;
  PetscScalar VSATN0;
  PetscScalar VSATP0;
  PetscScalar VSATN_A;
  PetscScalar VSATP_A;

  void Mob_Lucent_Init()
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
    T300       = 300.0*K;
    E1         = 1.0*V/cm;
    N1         = 1e0*std::pow(cm,-3);

    for(PetscScalar T=10*K; T<1000*K; T+=10*K)
    {
      Pn_temperature_table.push_back(Pn_Limiter_Init(T));
      Pp_temperature_table.push_back(Pp_Limiter_Init(T));
    }

    EXN4_LUC = 2.330000E-02;
    EXP4_LUC = 1.190000E-02;
    EXN9_LUC = 7.670000E-02;
    EXP9_LUC = 1.230000E-01;
    AN_LUC   = 2.580000E+00;
    AP_LUC   = 2.180000E+00;
    BN_LUC   = 3.610000E+07*cm/s;
    BP_LUC   = 1.510000E+07*cm/s;
    CN_LUC   = 1.700000E+04*cm*cm/V/s*std::pow(V/cm,PetscScalar(1.0/3.0))*std::pow(cm,3*EXN4_LUC);
    CP_LUC   = 4.180000E+03*cm*cm/V/s*std::pow(V/cm,PetscScalar(1.0/3.0))*std::pow(cm,3*EXP4_LUC);
    DN_LUC   = 3.580000E+18*cm*cm/V/s;
    DP_LUC   = 4.100000E+15*cm*cm/V/s;
    FN_LUC   = 6.850000E-21*std::pow(cm,3*(1-EXN9_LUC));
    FP_LUC   = 7.820000E-21*std::pow(cm,3*(1-EXP9_LUC));
    KN_LUC   = 1.700000E+00;
    KP_LUC   = 9.000000E-01;

    BETAN = 2.000000E+00;
    BETAP = 1.000000E+00;
    VSATN0  =  2.400000E7*cm/s;
    VSATN_A =  0.8;
    VSATP0  =  2.400000E7*cm/s;
    VSATP_A =  0.8;

#ifdef __CALIBRATE__
         parameter_map.insert(para_item("MMNN.UM",    PARA("MMNN.UM",    "", "cm*cm/V/s", cm*cm/V/s, &MMNN_UM)) );
         parameter_map.insert(para_item("MMXN.UM",    PARA("MMXN.UM",    "", "cm*cm/V/s",cm*cm/V/s , &MMXN_UM)) );
         parameter_map.insert(para_item("NRFN.UM",    PARA("NRFN.UM",    "", "std::pow(cm,-3)",std::pow(cm,-3) , &NRFN_UM)) );
         parameter_map.insert(para_item("ALPN.UM",    PARA("ALPN.UM",    "", "-", 1.0, &ALPN_UM)) );
         parameter_map.insert(para_item("TETN.UM",    PARA("TETN.UM",    "", "-",1.0 , &TETN_UM)) );
         parameter_map.insert(para_item("NRFD.UM",    PARA("NRFD.UM",    "", "std::pow(cm,-3)",std::pow(cm,-3) , &NRFD_UM)) );
         parameter_map.insert(para_item("CRFD.UM",    PARA("CRFD.UM",    "", "-", 1.0, &CRFD_UM)) );
         parameter_map.insert(para_item("MMNP.UM",    PARA("MMNP.UM",    "", "cm*cm/V/s", cm*cm/V/s, &MMNP_UM)) );
         parameter_map.insert(para_item("MMXP.UM",    PARA("MMXP.UM",    "", "cm*cm/V/s", cm*cm/V/s, &MMXP_UM)) );
         parameter_map.insert(para_item("NRFP.UM",    PARA("NRFP.UM",    "", "std::pow(cm,-3)", std::pow(cm,-3), &NRFP_UM)) );
         parameter_map.insert(para_item("ALPP.UM",    PARA("ALPP.UM",    "", "-",1.0 , &ALPP_UM)) );
         parameter_map.insert(para_item("TETP.UM",    PARA("TETP.UM",    "", "-",1.0 , &TETP_UM)) );
         parameter_map.insert(para_item("NRFA.UM",    PARA("NRFA.UM",    "", "std::pow(cm,-3)",std::pow(cm,-3) , &NRFA_UM)) );
         parameter_map.insert(para_item("CRFA.UM",    PARA("CRFA.UM",    "", "-",1.0 , &CRFA_UM)) );
         parameter_map.insert(para_item("NSC.REF",    PARA("NSC.REF",    "", "std::pow(cm,-2)", std::pow(cm,-2), &NSC_REF)) );
         parameter_map.insert(para_item("CAR.REF",    PARA("CAR.REF",    "", "std::pow(cm,-3)",std::pow(cm,-3) , &CAR_REF)) );
         parameter_map.insert(para_item("me_over_m0",    PARA("me_over_m0",    "", "-",1.0 , &me_over_m0)) );
         parameter_map.insert(para_item("mh_over_m0",    PARA("mh_over_m0",    "", "-", 1.0, &mh_over_m0)) );
         parameter_map.insert(para_item("me_over_mh",    PARA("me_over_mh",    "", "-", 1.0, &me_over_mh)) );
         parameter_map.insert(para_item("EXN4.LUC",    PARA("EXN4.LUC",    "", "-", 1.0, &EXN4_LUC)) );
         parameter_map.insert(para_item("EXP4.LUC",    PARA("EXP4.LUC",    "", "-", 1.0, &EXP4_LUC)) );
         parameter_map.insert(para_item("EXN9.LUC",    PARA("EXN9.LUC",    "", "-",1.0 , &EXN9_LUC)) );
         parameter_map.insert(para_item("AN.LUC",    PARA("AN.LUC",    "", "-", 1.0, &AN_LUC)) );
         parameter_map.insert(para_item("AP.LUC",    PARA("AP.LUC",    "", "-", 1.0, &AP_LUC)) );
         parameter_map.insert(para_item("EXP9.LUC",    PARA("EXP9.LUC",    "", "-", 1.0, &EXP9_LUC)) );
         parameter_map.insert(para_item("BN.LUC",    PARA("BN.LUC",    "", "cm/s", cm/s, &BN_LUC)) );
         parameter_map.insert(para_item("BP.LUC",    PARA("BP.LUC",    "", "cm/s",cm/s , &BP_LUC)) );
         parameter_map.insert(para_item("CN.LUC",    PARA("CN.LUC",    "", "cm*cm/V/s*std::pow(V/cm,PetscScalar(1.0/3.0))*std::pow(cm,3*EXN4.LUC)",cm*cm/V/s*std::pow(V/cm,PetscScalar(1.0/3.0))*std::pow(cm,3*EXN4_LUC) , &CN_LUC)) );
         parameter_map.insert(para_item("CP.LUC",    PARA("CP.LUC",    "", "cm*cm/V/s*std::pow(V/cm,PetscScalar(1.0/3.0))*std::pow(cm,3*EXP4.LUC)", cm*cm/V/s*std::pow(V/cm,PetscScalar(1.0/3.0))*std::pow(cm,3*EXP4_LUC), &CP_LUC)) );
         parameter_map.insert(para_item("FN.LUC",    PARA("FN.LUC",    "", "std::pow(cm,3*(1-EXN9.LUC))",std::pow(cm,3*(1-EXN9_LUC)) , &FN_LUC)) );
         parameter_map.insert(para_item("FP.LUC",    PARA("FP.LUC",    "", "std::pow(cm,3*(1-EXP9.LUC))", std::pow(cm,3*(1-EXP9_LUC)), &FP_LUC)) );
         parameter_map.insert(para_item("KN.LUC",    PARA("KN.LUC",    "", "-", 1.0, &KN_LUC)) );
         parameter_map.insert(para_item("KP.LUC",    PARA("KP.LUC",    "", "-", 1.0, &KP_LUC)) );
         parameter_map.insert(para_item("BETAN",    PARA("BETAN",    "", "-", 1.0, &BETAN)) );
         parameter_map.insert(para_item("BETAP",    PARA("BETAP",    "", "-", 1.0, &BETAP)) );
         parameter_map.insert(para_item("VSATN0",    PARA("VSATN0",    "", "cm/s", cm/s, &VSATN0)) );
         parameter_map.insert(para_item("VSATP0",    PARA("VSATP0",    "", "cm/s", cm/s , &VSATP0)) );
         parameter_map.insert(para_item("VSATN.A",    PARA("VSATN.A",    "", "-", 1.0, &VSATN_A)) );
         parameter_map.insert(para_item("VSATP.A",    PARA("VSATP.A",    "", "-",1.0 , &VSATP_A)) );
#endif

  }
  //---------------------------------------------------------------------------
  // Electron low field mobility
  PetscScalar ElecMobPhilips(const PetscScalar &p,const PetscScalar &n,const PetscScalar &Tl) const
  {
    PetscScalar mu_lattice = MMXN_UM*std::pow(Tl/T300,-TETN_UM);
    PetscScalar mu1 = MMXN_UM*MMXN_UM/(MMXN_UM-MMNN_UM)*std::pow(Tl/T300,3*ALPN_UM-1.5);
    PetscScalar mu2 = MMXN_UM*MMNN_UM/(MMXN_UM-MMNN_UM)*sqrt(T300/Tl);
    PetscScalar Na  = ReadDopingNa()+N1;
    PetscScalar Nd  = ReadDopingNd()+N1;
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    PetscScalar Nsc = Nds+Nas+fabs(p);

    PetscScalar P   = 1.0/(2.459/(NSC_REF/std::pow(Nsc,PetscScalar(2.0/3.0)))+3.828*fabs(n+p)/(CAR_REF*me_over_m0))*(Tl/T300)*(Tl/T300);
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
    PetscScalar Na  = ReadDopingNa()+N1;
    PetscScalar Nd  = ReadDopingNd()+N1;
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    AutoDScalar Nsc = Nds+Nas+fabs(p);

    AutoDScalar P   = 1.0/(2.459/(NSC_REF/adtl::pow(Nsc,PetscScalar(2.0/3.0)))+3.828*fabs(n+p)/(CAR_REF*me_over_m0))*(Tl/T300)*(Tl/T300);
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
    PetscScalar Na  = ReadDopingNa()+N1;
    PetscScalar Nd  = ReadDopingNd()+N1;
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    PetscScalar Nsc = Nds+Nas+fabs(n);

    PetscScalar P   = 1.0/(2.459/(NSC_REF/std::pow(Nsc,PetscScalar(2.0/3.0)))+3.828*fabs(n+p)/(CAR_REF*mh_over_m0))*(Tl/T300)*(Tl/T300);
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
    PetscScalar Na  = ReadDopingNa()+N1;
    PetscScalar Nd  = ReadDopingNd()+N1;
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    AutoDScalar Nsc = Nds+Nas+fabs(n);

    AutoDScalar P   = 1.0/(2.459/(NSC_REF/adtl::pow(Nsc,PetscScalar(2.0/3.0)))+3.828*fabs(n+p)/(CAR_REF*mh_over_m0))*(Tl/T300)*(Tl/T300);
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
  //the surface acoustical phono scattering and roughness scattering for electron
  PetscScalar ElecMobLombardi(const PetscScalar &p,const PetscScalar &n,const PetscScalar &Tl,const PetscScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+N1;
    PetscScalar ET = fabs(Et)+E1;
    PetscScalar mu_ac = BN_LUC/ET + CN_LUC*std::pow(N_total,EXN4_LUC)*std::pow(Tl/T300,-KN_LUC)*std::pow(ET,PetscScalar(-1.0/3.0));
    PetscScalar r = AN_LUC + FN_LUC*fabs(n+p)/std::pow(N_total,EXN9_LUC);
    if(r > 2.62) r = 2.62;
    PetscScalar mu_sr = DN_LUC*std::pow(ET/E1,-r);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }
  AutoDScalar ElecMobLombardi(const AutoDScalar &p,const AutoDScalar &n,const AutoDScalar &Tl,const AutoDScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+N1;
    AutoDScalar ET = fabs(Et)+E1;
    AutoDScalar mu_ac = BN_LUC/ET + CN_LUC*std::pow(N_total,EXN4_LUC)*adtl::pow(Tl/T300,-KN_LUC)*adtl::pow(ET,PetscScalar(-1.0/3.0));
    AutoDScalar r = AN_LUC + FN_LUC*fabs(n+p)/std::pow(N_total,EXN9_LUC);
    if(r > 2.62) r = 2.62;
    AutoDScalar mu_sr = DN_LUC*adtl::pow(ET/E1,-r);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }

  //---------------------------------------------------------------------------
  //the surface acoustical phono scattering and roughness scattering for hole
  PetscScalar HoleMobLombardi(const PetscScalar &p,const PetscScalar &n,const PetscScalar &Tl,const PetscScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+N1;
    PetscScalar ET = fabs(Et)+E1;
    PetscScalar mu_ac = BP_LUC/ET + CP_LUC*std::pow(N_total,EXP4_LUC)*std::pow(Tl/T300,-KP_LUC)*std::pow(ET,PetscScalar(-1.0/3.0));
    PetscScalar r = AP_LUC + FP_LUC*fabs(n+p)/std::pow(N_total,EXP9_LUC);
    PetscScalar mu_sr = DP_LUC*std::pow(ET/E1,-r);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }
  AutoDScalar HoleMobLombardi(const AutoDScalar &p,const AutoDScalar &n,const AutoDScalar &Tl,const AutoDScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+N1;
    AutoDScalar ET = fabs(Et)+E1;
    AutoDScalar mu_ac = BP_LUC/ET + CP_LUC*std::pow(N_total,EXP4_LUC)*adtl::pow(Tl/T300,-KP_LUC)*adtl::pow(ET,PetscScalar(-1.0/3.0));
    AutoDScalar r = AP_LUC + FP_LUC*fabs(n+p)/std::pow(N_total,EXP9_LUC);
    AutoDScalar mu_sr = DP_LUC*adtl::pow(ET/E1,-r);
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
    PetscScalar mu0  = 1.0/(1.0/ElecMobPhilips(p,n,Tl)+1.0/ElecMobLombardi(p,n,Tl,Et));
    return 2*mu0/(1+std::pow(1+std::pow(2*mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN));
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    AutoDScalar vsat = VSATN0/(1+VSATN_A*exp(Tl/(2*T300)));
    AutoDScalar mu0  = 1.0/(1.0/ElecMobPhilips(p,n,Tl)+1.0/ElecMobLombardi(p,n,Tl,Et));
    return 2*mu0/(1+adtl::pow(1+adtl::pow(2*mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN));
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = VSATP0/(1+VSATP_A*exp(Tl/(2*T300)));
    PetscScalar mu0  = 1.0/(1.0/HoleMobPhilips(p,n,Tl)+1.0/HoleMobLombardi(p,n,Tl,Et));
    return 2*mu0/(1+std::pow(1+std::pow(2*mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP));
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    AutoDScalar vsat = VSATP0/(1+VSATP_A*exp(Tl/(2*T300)));
    AutoDScalar mu0  = 1.0/(1.0/HoleMobPhilips(p,n,Tl)+1.0/HoleMobLombardi(p,n,Tl,Et));
    return 2*mu0/(1+adtl::pow(1+adtl::pow(2*mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP));
  }


// constructor
public:
  GSS_Si_Mob_Lucent(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    PMI_Info = "This is the Lucent mobility model of Silicon";
    Mob_Lucent_Init();
  }

  ~GSS_Si_Mob_Lucent(){}
}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  and it setup Lucent mobility model
 */
extern "C"
{
  DLL_EXPORT_DECLARE PMIS_Mobility* PMIS_Si_Mob_Lucent (const PMIS_Environment& env)
  {
    return new GSS_Si_Mob_Lucent(env);
  }
}
