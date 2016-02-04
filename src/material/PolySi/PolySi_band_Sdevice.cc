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
// Material Type: Poly Silicon


#include "PMI.h"


class GSS_PolySi_BandStructure_Sdevice : public PMIS_BandStructure
{
private:
  PetscScalar T300;
  //[Bandgap]
    // Bandgap and Effective Density of States
  PetscScalar EG0;       // The energy bandgap of the material at 0 K.
  PetscScalar EGTP;      // The value of lattice temperature
  PetscScalar EGALPH;    // The value of alpha used in calculating the temperature depended energy bandgap.
  PetscScalar EGBETA;    // The value of beta  used in calculating the temperature depended energy bandgap.
  
  PetscScalar ELECMASS;  // The relative effective mass of electron
  PetscScalar HOLEMASS;  // The relative effective mass of hole
  PetscScalar NC300;     // The effective density of states in the conduction band at 300K.
  PetscScalar NV300;     // The effective density of states in the valence band at 300K.
  PetscScalar NC_F;      // The parameter for temperature depended effective density of states in the conduction band.
  PetscScalar NV_F;      // The parameter for temperature depended effective density of states in the valence band.
  
  // Model of Bandgap Narrowing due to Heavy Doping
  std::string BGN_Model;

  PetscScalar dEG0;      // The band gap correction term.
  PetscScalar N0_BGN;    // The concentration parameter used in band-gap narrowing model.
  PetscScalar E0_BGN;    // The voltage parameter used in band-gap narrowing model.

  // Init value
  void Eg_Init()
  {
    // Use parameters from Green (JAP 67, p.2945, 1990) for
    // silicon bandgap and densities of states
    EG0       = 1.16964*eV;
    dEG0      = 0.0*eV;
    EGTP      = 0.0*K;
    EGALPH    = 4.7300e-04*eV/K;
    EGBETA    = 6.3600e+02*K;

    ELECMASS  = 1.0903*me;
    HOLEMASS  = 1.1525*me;
    NC300     = 2.86E19*std::pow(cm,-3);
    NV300     = 3.10E19*std::pow(cm,-3);
    NC_F      = 1.58;
    NV_F      = 1.85;

    
    BGN_Model = "Bennett";
    dEG0      = 0.0*eV;
    N0_BGN    = 3.1620e+18*std::pow(cm,-3);
    E0_BGN    = 6.8400e-03*eV;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("EG0",    PARA("EG0",    "The energy bandgap of the material at 0 K", "eV", eV, &EG0)) );
    parameter_map.insert(para_item("EGALPH", PARA("EGALPH", "The value of alpha used in calculating the temperature depended energy bandgap", "eV/K", eV/K, &EGALPH)) );
    parameter_map.insert(para_item("EGBETA", PARA("EGBETA", "The value of beta used in calculating the temperature depended energy bandgap", "K", K, &EGBETA)) );

    parameter_map.insert(para_item("ELECMASS", PARA("ELECMASS", "The relative effective mass of electron", "electron mass", me, &ELECMASS)) );
    parameter_map.insert(para_item("HOLEMASS", PARA("HOLEMASS", "The relative effective mass of hole", "electron mass", me, &HOLEMASS)) );
    parameter_map.insert(para_item("NC300",    PARA("NC300",    "The effective density of states in the conduction band at 300K", "cm^-3", std::pow(cm,-3), &NC300)) );
    parameter_map.insert(para_item("NV300",    PARA("NV300",    "The effective density of states in the valence band at 300K", "cm^-3", std::pow(cm,-3), &NV300)) );
    parameter_map.insert(para_item("NC.F",     PARA("NC.F",     "The parameter for temperature depended effective density of states in the conduction band", "-", 1.0, &NC_F)) );
    parameter_map.insert(para_item("NV.F",     PARA("NV.F",     "The parameter for temperature depended effective density of states in the valence band", "-", 1.0, &NV_F)) );

    
    parameter_map.insert(para_item("BGN.MODEL",PARA("BGN.MODEL","Specific bandgap narrowing model, support Bennett(default), Slotboom, OldSlotboom, delAlamo", &BGN_Model)) );
#endif

  }
public:
  
  //---------------------------------------------------------------------------
  PetscScalar Eg (const PetscScalar &Tl)
  {
    return EG0 + dEG0 + EGALPH * EGTP*EGTP / (EGBETA + EGTP) - EGALPH * Tl*Tl / (EGBETA + Tl);
  }
  AutoDScalar Eg (const AutoDScalar &Tl)
  {
    return EG0 + dEG0 + EGALPH * EGTP*EGTP / (EGBETA + EGTP) - EGALPH * Tl*Tl / (EGBETA + Tl);
  }

  //---------------------------------------------------------------------------
  // procedure of Bandgap Narrowing due to Heavy Doping
  
  
  void Bandgap_Setup()
  {

//     dEg0(Bennett)	= 0.0000e+00	# [eV]
//     dEg0(Slotboom)	= -4.7950e-03	# [eV]
//     dEg0(OldSlotboom)	= -1.5950e-02	# [eV]
//     dEg0(delAlamo)	= -1.4070e-02	# [eV]
//         
//     OldSlotboom
//     { * deltaEg = dEg0 + Ebgn ( ln(N/Nref) + [ (ln(N/Nref))^2 + 0.5]^1/2 )
//       * dEg0 is defined in BandGap section 
//         Ebgn	= 9.0000e-03	# [eV]
//         Nref	= 1.0000e+17	# [cm^(-3)]
//     }
// 
// 
//     Slotboom
//     { * deltaEg = dEg0 + Ebgn ( ln(N/Nref) + [ (ln(N/Nref))^2 + 0.5]^1/2 )
//       * dEg0 is defined in BandGap section 
//         Ebgn	= 6.9200e-03	# [eV]
//         Nref	= 1.3000e+17	# [cm^(-3)]
//     }
// 
//     delAlamo
//     { * deltaEg = dEg0 + Ebgn  ln(N/Nref) 
//       * dEg0 is defined in BandGap section 
//         Ebgn	= 0.0187	# [eV]
//         Nref	= 7.0000e+17	# [cm^(-3)]
//     }
// 
//     Bennett
//     { * deltaEg = dEg0 + Ebgn (ln(N/Nref))^2
//       * dEg0 is defined in BandGap section 
//         Ebgn	= 6.8400e-03	# [eV]
//         Nref	= 3.1620e+18	# [cm^(-3)]
//     }


    if(BGN_Model == "Slotboom")
    {
      dEG0      = -4.7950e-03*eV;
      N0_BGN    = 1.3000e+17*std::pow(cm,-3);
      E0_BGN    = 6.9200e-03*eV;
    }
    else if(BGN_Model == "OldSlotboom")
    {
      dEG0      = -1.5950e-02*eV;
      N0_BGN    = 1.0000e+17*std::pow(cm,-3);
      E0_BGN    = 9.0000e-03*eV;
    }
    else if(BGN_Model == "delAlamo")
    {
      dEG0      = -1.4070e-02*eV;
      N0_BGN    = 7.0000e+17*std::pow(cm,-3);
      E0_BGN    = 0.0187*eV;
    }
    else if(BGN_Model == "Bennett")
    {
      dEG0      = 0.0*eV;
      N0_BGN    = 3.1620e+18*std::pow(cm,-3);
      E0_BGN    = 6.8400e-03*eV;
    }
    else
    {
      dEG0      = 0.0*eV;
      N0_BGN    = 1e+18*std::pow(cm,-3);
      E0_BGN    = 0.0*eV;
    }
  }
  
  PetscScalar EgNarrow_Slotboom(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N = Na+Nd+1.0*std::pow(cm,-3);
    return dEG0 + E0_BGN * ( log(N/N0_BGN) + sqrt( std::pow(log(N/N0_BGN), 2) + 0.5) );
  }
  AutoDScalar EgNarrow_Slotboom(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N = Na+Nd+1.0*std::pow(cm,-3);
    return dEG0 + E0_BGN * ( log(N/N0_BGN) + sqrt( std::pow(log(N/N0_BGN), 2) + 0.5) );
  }

  PetscScalar EgNarrow_Bennett(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N = Na+Nd+1.0*std::pow(cm,-3);
    if(N > N0_BGN)
      return dEG0 + E0_BGN * std::pow(log(N/N0_BGN), 2);
    return 0.0;
  }
  AutoDScalar EgNarrow_Bennett(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N = Na+Nd+1.0*std::pow(cm,-3);
    if(N > N0_BGN)
      return dEG0 + E0_BGN * std::pow(log(N/N0_BGN), 2);
    return 0.0;
  }
  
  PetscScalar EgNarrow_delAlamo(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N = Na+Nd+1.0*std::pow(cm,-3);
    if(N > N0_BGN)
      return  dEG0 + E0_BGN *  log(N/N0_BGN);
    return 0.0;
  }
  AutoDScalar EgNarrow_delAlamo(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N = Na+Nd+1.0*std::pow(cm,-3);
    if(N > N0_BGN)
      return  dEG0 + E0_BGN *  log(N/N0_BGN);
    return 0.0;
  }

  PetscScalar EgNarrow(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    if(BGN_Model == "Slotboom" || BGN_Model == "OldSlotboom")
      return EgNarrow_Slotboom(p,n,Tl);
    else if (BGN_Model == "delAlamo")
      return EgNarrow_delAlamo(p,n,Tl);
    else if(BGN_Model == "Bennett")
      return EgNarrow_Bennett(p,n,Tl);
    return 0.0;
  }
  AutoDScalar EgNarrow(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    if(BGN_Model == "Slotboom" || BGN_Model == "OldSlotboom")
      return EgNarrow_Slotboom(p,n,Tl);
    else if (BGN_Model == "delAlamo")
      return EgNarrow_delAlamo(p,n,Tl);
    else if(BGN_Model == "Bennett")
      return EgNarrow_Bennett(p,n,Tl);
    return 0.0;
  }

  PetscScalar EgNarrowToEc   (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl){return 0.5*EgNarrow(p, n, Tl);}
  PetscScalar EgNarrowToEv   (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl){return 0.5*EgNarrow(p, n, Tl);}

  AutoDScalar EgNarrowToEc   (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl){return 0.5*EgNarrow(p, n, Tl);}
  AutoDScalar EgNarrowToEv   (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl){return 0.5*EgNarrow(p, n, Tl);}


  //---------------------------------------------------------------------------
  //electron and hole effect mass
  PetscScalar EffecElecMass (const PetscScalar &Tl)
  {
    return ELECMASS;
  }
  AutoDScalar EffecElecMass (const AutoDScalar &Tl)
  {
    return ELECMASS;
  }
  PetscScalar EffecHoleMass (const PetscScalar &Tl)
  {
    return HOLEMASS;
  }
  AutoDScalar EffecHoleMass (const AutoDScalar &Tl)
  {
    return HOLEMASS;
  }

  //---------------------------------------------------------------------------
  // Nc and Nv
  PetscScalar Nc (const PetscScalar &Tl)
  {
    return NC300*std::pow(Tl/T300,NC_F);
  }
  AutoDScalar Nc (const AutoDScalar &Tl)
  {
    return NC300*adtl::pow(Tl/T300,NC_F);
  }
  PetscScalar Nv (const PetscScalar &Tl)
  {
    return NV300*std::pow(Tl/T300,NV_F);
  }
  AutoDScalar Nv (const AutoDScalar &Tl)
  {
    return NV300*adtl::pow(Tl/T300,NV_F);
  }

  //---------------------------------------------------------------------------
  PetscScalar ni (const PetscScalar &Tl)
  {
    PetscScalar bandgap = Eg(Tl);
    return sqrt(Nc(Tl)*Nv(Tl))*exp(-bandgap/(2*kb*Tl));
  }

  // nie, Eg narrow should be considered
  PetscScalar nie (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar bandgap = Eg(Tl);
    return sqrt(Nc(Tl)*Nv(Tl))*exp(-bandgap/(2*kb*Tl))*exp(EgNarrow(p, n, Tl)/(2*kb*Tl));
  }
  AutoDScalar nie (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar bandgap = Eg(Tl);
    return sqrt(Nc(Tl)*Nv(Tl))*exp(-bandgap/(2*kb*Tl))*exp(EgNarrow(p, n, Tl)/(2*kb*Tl));
  }

  //end of Bandgap

private:
  //[Lifetime]
  //Lifetimes
  PetscScalar TAUN0;         // The Shockley-Read-Hall electron lifetime.
  PetscScalar TAUP0;         // The Shockley-Read-Hall hole lifetime.
  PetscScalar STAUN;         // The electron surface recombination velocity.
  PetscScalar STAUP;         // The hole surface recombination velocity.
  //Concentration-Dependent Lifetimes
  PetscScalar NSRHN;         // The Shockley-Read-Hall concentration parameter for electrons.
  PetscScalar AN;            // The constant term in the concentration-dependent expression for electron lifetime.
  PetscScalar BN;            // The linear term coefficient in the concentration-dependent expression for electron lifetime.
  PetscScalar CN;            // The exponential term coefficient in the concentration-dependent expression for electron lifetime.
  PetscScalar EN;            // The exponent in the concentration-dependent expression for electron lifetime.
  PetscScalar NSRHP;         // The Shockley-Read-Hall concentration parameter for holes.
  PetscScalar AP;            // The constant term in the concentration-dependent expression for hole lifetime.
  PetscScalar BP;            // The linear term coefficient in the concentration-dependent expression for hole lifetime.
  PetscScalar CP;            // The exponential term coefficient in the concentration-dependent expression for hole lifetime.
  PetscScalar EP;            // The exponent in the concentration-dependent expression for hole lifetime.
  // Lattice Temperature-Dependent Lifetimes
  PetscScalar EXN_TAU;       // The exponent of lattice temperature dependent electron lifetime.
  PetscScalar EXP_TAU;       // The exponent of lattice temperature dependent hole lifetime.

  //Init value
  void Lifetime_Init()
  {
    TAUN0     = 1.000000e-07*s;
    TAUP0     = 1.000000e-07*s;
    STAUN     = 0.000000e+00*cm/s;
    STAUP     = 0.000000e+00*cm/s;
    NSRHN     = 5.000000e+16*std::pow(cm,-3);
    AN        = 1.000000e+00;
    BN        = 1.000000e+00;
    CN        = 0.000000e+00;
    EN        = 2.000000e+00;
    NSRHP     = 5.000000e+16*std::pow(cm,-3);
    AP        = 1.000000e+00;
    BP        = 1.000000e+00;
    CP        = 0.000000e+00;
    EP        = 2.000000e+00;
    EXN_TAU   = 0.000000e+00;
    EXP_TAU   = 0.000000e+00;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("TAUN0",    PARA("TAUN0",    "The Shockley-Read-Hall electron lifetime", "s", s, &TAUN0)) );
    parameter_map.insert(para_item("TAUP0",    PARA("TAUP0",    "The Shockley-Read-Hall hole lifetime", "s", s, &TAUP0)) );
    parameter_map.insert(para_item("STAUN",    PARA("STAUN",    "The electron surface recombination velocity", "cm/s", cm/s, &STAUN)) );
    parameter_map.insert(para_item("STAUP",    PARA("STAUP",    "The hole surface recombination velocity", "cm/s", cm/s, &STAUP)) );

    parameter_map.insert(para_item("NSRHN", PARA("NSRHN", "The Shockley-Read-Hall concentration parameter for electrons", "cm^-3", std::pow(cm,-3), &NSRHN)) );
    //    parameter_map.insert(para_item("AN",    PARA("AN",    "The constant term in the concentration-dependent expression for electron lifetime", "-", 1.0, &AN)) );
    //    parameter_map.insert(para_item("BN",    PARA("BN",    "The linear term coefficient in the concentration-dependent expression for electron lifetime", "-", 1.0, &BN)) );
    //    parameter_map.insert(para_item("CN",    PARA("CN",    "The exponential term coefficient in the concentration-dependent expression for electron lifetime", "-", 1.0, &CN)) );
    //    parameter_map.insert(para_item("EN",    PARA("EN",    "The exponent in the concentration-dependent expression for electron lifetime", "-", 1.0, &EN)) );

    parameter_map.insert(para_item("NSRHP", PARA("NSRHP", "The Shockley-Read-Hall concentration parameter for holes", "cm^-3", std::pow(cm,-3), &NSRHP)) );
    //    parameter_map.insert(para_item("AP",    PARA("AP",    "The constant term in the concentration-dependent expression for hole lifetime", "-", 1.0, &AP)) );
    //    parameter_map.insert(para_item("BP",    PARA("BP",    "The linear term coefficient in the concentration-dependent expression for hole lifetime", "-", 1.0, &BP)) );
    //    parameter_map.insert(para_item("CP",    PARA("CP",    "The exponential term coefficient in the concentration-dependent expression for hole lifetime", "-", 1.0, &CP)) );
    //    parameter_map.insert(para_item("EP",    PARA("EP",    "The exponent in the concentration-dependent expression for hole lifetime", "-", 1.0, &EP)) );

    parameter_map.insert(para_item("EXN_TAU",    PARA("EXN_TAU",    "The exponent of lattice temperature dependent electron lifetime", "-", 1.0, &EXN_TAU)) );
    parameter_map.insert(para_item("EXP_TAU",    PARA("EXP_TAU",    "The exponent of lattice temperature dependent hole lifetime", "-", 1.0, &EXP_TAU)) );
#endif

  }

public:
  //---------------------------------------------------------------------------
  // electron lift time for SHR Recombination
  PetscScalar TAUN (const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return TAUN0/(1+(Na+Nd)/NSRHN)*std::pow(Tl/T300,EXN_TAU);
  }
  AutoDScalar TAUN (const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return TAUN0/(1+(Na+Nd)/NSRHN)*adtl::pow(Tl/T300,EXN_TAU);
  }

  //---------------------------------------------------------------------------
  // hole lift time for SHR Recombination
  PetscScalar TAUP (const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return TAUP0/(1+(Na+Nd)/NSRHP)*std::pow(Tl/T300,EXP_TAU);
  }
  AutoDScalar TAUP (const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return TAUP0/(1+(Na+Nd)/NSRHP)*adtl::pow(Tl/T300,EXP_TAU);
  }
  // End of Lifetime

  //[the fit parameter for density-gradient solver]
  PetscScalar Gamman         () {return 3.6;}
  PetscScalar Gammap         () {return 5.6;}

private:
  //[Recombination]
  // SRH, Auger, and Direct Recombination
  PetscScalar ETRAP;         // The trap level (Et - Ei) used in determining the Shockley-Read-Hall recombination rate.
  PetscScalar AUGN;          // The Auger coefficient for electrons.
  PetscScalar AUGP;          // The Auger coefficient for holes.
  PetscScalar C_DIRECT;      // c.
  // Recombination Including Tunneling
  PetscScalar M_RTUN;        // The trap-assisted tunneling effective mass. *free electron rest mass m0
  PetscScalar S_RTUN;        // Band-to-band field power ratio.
  PetscScalar B_RTUN;        // Band-to-band tunneling rate proportionality factor.
  PetscScalar E_RTUN;        // Band-to-band reference electric field.

  // Init value
  void Recomb_Init()
  {
    ETRAP   = 0.000000e+00*eV;
    AUGN    =  2.800000e-31*std::pow(cm,6)/s;
    AUGP    =  9.900000e-32*std::pow(cm,6)/s;
    C_DIRECT = 0.000000e+00*std::pow(cm,3)/s;
    M_RTUN   = 2.500000e-01;
    S_RTUN   = 2.500000e+00;
    B_RTUN   = 4.000000e+14*std::pow(cm,S_RTUN-3)*std::pow(V,-S_RTUN)/s;
    E_RTUN   = 1.900000e+07*V/cm;
#ifdef __CALIBRATE__
    //    parameter_map.insert(para_item("ETRAP",    PARA("ETRAP",    "The trap level (Et - Ei) used in determining the Shockley-Read-Hall recombination rate", "eV", eV, &ETRAP)) );
    parameter_map.insert(para_item("AUGN",     PARA("AUGN",     "The Auger coefficient for electrons", "cm^6/s", std::pow(cm,6)/s, &AUGN)) );
    parameter_map.insert(para_item("AUGP",     PARA("AUGP",     "The Auger coefficient for holes", "cm^6/s", std::pow(cm,6)/s, &AUGP)) );
    parameter_map.insert(para_item("C.DIRECT", PARA("C.DIRECT", "The direct generation/recombination coefficient", "cm^3/s", std::pow(cm,3)/s, &C_DIRECT)) );

    //    parameter_map.insert(para_item("M_RTUN",   PARA("M_RTUN", "The trap-assisted tunneling effective mass", "-", 1.0, &M_RTUN)) );
    //    parameter_map.insert(para_item("S_RTUN",   PARA("S_RTUN", "Band-to-band field power ratio", "-", 1.0, &S_RTUN)) );
    //    parameter_map.insert(para_item("B_RTUN",   PARA("B_RTUN", "Band-to-band tunneling rate proportionality factor", "cm^(S_RTUN-3)V^(-S_RTUN)/s", std::pow(cm,S_RTUN -3)*std::pow(V,-S_RTUN)/s, &B_RTUN)) );
    //    parameter_map.insert(para_item("E_RTUN",   PARA("E_RTUN", "Band-to-band reference electric field", "V/cm", V/cm, &E_RTUN)) );
#endif

  }

public:

  /**
   * @return direct Recombination rate
   */
  PetscScalar CDIR           (const PetscScalar &Tl)  { return C_DIRECT; }

  /**
   * @return electron Auger Recombination rate
   */
  PetscScalar AUGERN           (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) { return AUGN; }

  /**
   * @return hole Auger Recombination rate
   */
  PetscScalar AUGERP           (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)  { return AUGP; }

  //---------------------------------------------------------------------------
  // Direct Recombination
  PetscScalar R_Direct     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni =   nie(p, n, Tl);
    return C_DIRECT*(n*p-ni*ni);
  }
  AutoDScalar R_Direct     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni =   nie(p, n, Tl);
    return C_DIRECT*(n*p-ni*ni);
  }

  //---------------------------------------------------------------------------
  // Total Auger Recombination
  PetscScalar R_Auger     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni =   nie(p, n, Tl);
    return AUGN*(p*n*n-n*ni*ni)+AUGP*(n*p*p-p*ni*ni);
  }
  AutoDScalar R_Auger     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni =   nie(p, n, Tl);
    return AUGN*(p*n*n-n*ni*ni)+AUGP*(n*p*p-p*ni*ni);
  }

  //---------------------------------------------------------------------------
  // Electron Auger Recombination
  PetscScalar R_Auger_N     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni =   nie(p, n, Tl);
    return AUGN*(p*n*n-n*ni*ni);
  }
  AutoDScalar R_Auger_N     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni =   nie(p, n, Tl);
    return AUGN*(p*n*n-n*ni*ni);
  }
  //---------------------------------------------------------------------------
  // Hole Auger Recombination
  PetscScalar R_Auger_P     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni =   nie(p, n, Tl);
    return AUGP*(n*p*p-p*ni*ni);
  }
  AutoDScalar R_Auger_P     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni =   nie(p, n, Tl);
    return AUGP*(n*p*p-p*ni*ni);
  }


  //---------------------------------------------------------------------------
  // SHR Recombination
  PetscScalar R_SHR     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni =   nie(p, n, Tl);
    PetscScalar taun = TAUN(Tl);
    PetscScalar taup = TAUP(Tl);
    return (p*n-ni*ni)/(taup*(n+ni)+taun*(p+ni));
  }
  AutoDScalar R_SHR     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni =   nie(p, n, Tl);
    AutoDScalar taun = TAUN(Tl);
    AutoDScalar taup = TAUP(Tl);
    return (p*n-ni*ni)/(taup*(n+ni)+taun*(p+ni));
  }

  //---------------------------------------------------------------------------
  // Surface SHR Recombination
  PetscScalar R_Surf     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni =   nie(p, n, Tl);

    PetscScalar seps = 1e-8 * cm/s; // a very small recomb velocity
    if (STAUN < seps || STAUP < seps)
      return 0;
    else
      return (p*n - ni*ni) / ((n+ni)/STAUP + (p+ni)/STAUN);
  }
  AutoDScalar R_Surf     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni =   nie(p, n, Tl);

    PetscScalar seps = 1e-8 * cm/s; // a very small recomb velocity
    if (STAUN < seps || STAUP < seps)
      return 0;
    else
      return (p*n - ni*ni) / ((n+ni)/STAUP + (p+ni)/STAUN);
  }

  //---------------------------------------------------------------------------
  // total Recombination
  PetscScalar Recomb (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni =   nie(p, n, Tl);
    PetscScalar taun = TAUN(Tl);
    PetscScalar taup = TAUP(Tl);
    PetscScalar dn   = p*n-ni*ni;
    PetscScalar Rshr = dn/(taup*(n+ni)+taun*(p+ni));
    PetscScalar Rdir = C_DIRECT*dn;
    PetscScalar Raug = (AUGN*n+AUGP*p)*dn;
    return Rshr+Rdir+Raug;
  }
  AutoDScalar Recomb (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni =   nie(p, n, Tl);
    AutoDScalar taun = TAUN(Tl);
    AutoDScalar taup = TAUP(Tl);
    AutoDScalar dn   = p*n-ni*ni;
    AutoDScalar Rshr = dn/(taup*(n+ni)+taun*(p+ni));
    AutoDScalar Rdir = C_DIRECT*dn;
    AutoDScalar Raug = (AUGN*n+AUGP*p)*dn;
    return Rshr+Rdir+Raug;
  }

  // End of Recombination

private:
  //[energy relax time]
  PetscScalar  WTN0;
  PetscScalar  WTN1;
  PetscScalar  WTN2;
  PetscScalar  WTN3;
  PetscScalar  WTN4;
  PetscScalar  WTN5;
  PetscScalar  WTNL;
  PetscScalar  TNL;
  PetscScalar  WTP0;
  PetscScalar  WTP1;
  PetscScalar  WTP2;
  PetscScalar  WTP3;
  PetscScalar  WTP4;
  PetscScalar  WTP5;
  PetscScalar  WTPL;
  PetscScalar  TPL;
  // Init value
  void RelaxTime_Init()
  {
    WTN0 =  1.685200E-13*s;
    WTN1 =  1.029900E-13*s;
    WTN2 = -5.184500E-15*s;
    WTN3 =  0.000000E+00*s;
    WTN4 =  0.000000E+00*s;
    WTN5 =  0.000000E+00*s;
    WTNL =  6.800000E-13*s;
    TNL  =  2.979800E+03*K;
    WTP0 = -1.560000E-14*s;
    WTP1 =  1.380000E-13*s;
    WTP2 = -2.500000E-14*s;
    WTP3 =  2.310000E-15*s;
    WTP4 = -1.050000E-16*s;
    WTP5 =  1.820000E-18*s;
    WTPL =  2.000000E-13*s;
    TPL  =  1.000000E+05*K;
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("WTN0",    PARA("WTN0",    "Constant term for electron energy relaxatioin time.", "s", s, &WTN0)) );
    parameter_map.insert(para_item("WTN1",    PARA("WTN1",    "Coefficient of the linear term for the temperature dependence of electron energy relaxatioin time.", "s", s, &WTN1)) );
    parameter_map.insert(para_item("WTN2",    PARA("WTN2",    "Coefficient of the quadratic term for the temperature dependence of electron energy relaxatioin time.", "s", s, &WTN2)) );
    parameter_map.insert(para_item("TNL",     PARA("TNL",     "Electron temperature upper reference.", "K", K, &TNL)) );
    parameter_map.insert(para_item("WTNL",    PARA("WTNL",    "Electron energy relaxation time for electron temperature higher than TNL.", "s", s, &WTNL)) );

    parameter_map.insert(para_item("WTP0",    PARA("WTP0",    "Constant term for hole energy relaxation time.", "s", s, &WTP0)) );
    parameter_map.insert(para_item("WTP1",    PARA("WTP1",    "Coefficient of the linear term for the temperature dependence of hole energy relaxatioin time.", "s", s, &WTP1)) );
    parameter_map.insert(para_item("WTP2",    PARA("WTP2",    "Coefficient of the quadratic term for the temperature dependence of hole energy relaxatioin time.", "s", s, &WTP2)) );
    parameter_map.insert(para_item("WTP3",    PARA("WTP3",    "Coefficient of the cubic term for the temperature dependence of hole energy relaxatioin time.", "s", s, &WTP3)) );
    parameter_map.insert(para_item("WTP4",    PARA("WTP4",    "Coefficient of the forth-order term for the temperature dependence of hole energy relaxatioin time.", "s", s, &WTP4)) );
    parameter_map.insert(para_item("WTP5",    PARA("WTP5",    "Coefficient of the fifth-order term for the temperature dependence of hole energy relaxatioin time.", "s", s, &WTP5)) );
    parameter_map.insert(para_item("TPL",     PARA("TPL",     "Hole temperature upper reference.", "K", K, &TPL)) );
    parameter_map.insert(para_item("WTPL",    PARA("WTPL",    "Hole energy relaxation time for electron temperature higher than TPL.", "s", s, &WTPL)) );
#endif

  }
public:
  //---------------------------------------------------------------------------
  // Electron relaxation time for EBM
  PetscScalar ElecEnergyRelaxTime(const PetscScalar &Tn,const PetscScalar &Tl)
  {
    if(Tn>TNL)     return WTNL;
    PetscScalar x = 1+(Tn-Tl)/T300;
    return WTN0+ WTN1*x + WTN2*x*x;
  }
  AutoDScalar ElecEnergyRelaxTime(const AutoDScalar &Tn,const AutoDScalar &Tl)
  {
    if(Tn>TNL)     return WTNL;
    AutoDScalar x = 1+(Tn-Tl)/T300;
    return WTN0+ WTN1*x + WTN2*x*x;
  }

  //---------------------------------------------------------------------------
  // Hole relaxation time for EBM
  PetscScalar HoleEnergyRelaxTime(const PetscScalar &Tp,const PetscScalar &Tl)
  {
    if(Tp>TPL)     return WTPL;
    PetscScalar x = 1+(Tp-Tl)/T300;
    return WTP0+ WTP1*x + WTP2*x*x + WTP3*x*x*x + WTP4*std::pow(x,4) + WTP5*std::pow(x,5);
  }
  AutoDScalar HoleEnergyRelaxTime(const AutoDScalar &Tp,const AutoDScalar &Tl)
  {
    if(Tp>TPL)     return WTPL;
    AutoDScalar x = 1+(Tp-Tl)/T300;
    return WTP0+ WTP1*x + WTP2*x*x + WTP3*x*x*x + WTP4*adtl::pow(x,4) + WTP5*adtl::pow(x,5);
  }
  // end of energy relax time

private:
  // [Schottky and Heterojunction]
  PetscScalar ARICHN;
  PetscScalar ARICHP;
  PetscScalar VSURFN;   // Thermionic emission velocity of electron
  PetscScalar VSURFP;

  void   Schottky_Init()
  {
    ARICHN = 1.100000e+02*A/(K*cm)/(K*cm);
    ARICHP = 3.000000e+01*A/(K*cm)/(K*cm);
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("ARICHN", PARA("ARICHN", "The effective Richardson constants for electrons", "A/(K^2*cm^2)", A/(K*cm)/(K*cm), &ARICHN)) );
    parameter_map.insert(para_item("ARICHP", PARA("ARICHP", "The effective Richardson constants for holes", "A/(K^2*cm^2)", A/(K*cm)/(K*cm), &ARICHP)) );
#endif

  }

public:
  PetscScalar ARichN()
  { return ARICHN; }

  PetscScalar ARichP()
  { return ARICHP; }

  PetscScalar SchottyJsn (PetscScalar n,PetscScalar Tl,PetscScalar Vb)
  {
    PetscScalar VSURFN = ARICHN*Tl*Tl/(e*Nc(Tl));
    PetscScalar nb = Nc(Tl)*exp(-e*Vb/(kb*Tl));
    return -e*VSURFN*(n-nb);
  }
  AutoDScalar SchottyJsn (AutoDScalar n,AutoDScalar Tl,AutoDScalar Vb)
  {
    AutoDScalar VSURFN = ARICHN*Tl*Tl/(e*Nc(Tl));
    AutoDScalar nb = Nc(Tl)*exp(-e*Vb/(kb*Tl));
    return -e*VSURFN*(n-nb);
  }

  PetscScalar SchottyJsp (PetscScalar p,PetscScalar Tl,PetscScalar Vb)
  {
    PetscScalar VSURFP = ARICHP*Tl*Tl/(e*Nv(Tl));
    PetscScalar pb = Nv(Tl)*exp((-Eg(Tl)+e*Vb)/(kb*Tl));
    return e*VSURFP*(p-pb);
  }
  AutoDScalar SchottyJsp (AutoDScalar p,AutoDScalar Tl,AutoDScalar Vb)
  {
    AutoDScalar VSURFP = ARICHP*Tl*Tl/(e*Nv(Tl));
    AutoDScalar pb = Nv(Tl)*exp((-Eg(Tl)+e*Vb)/(kb*Tl));
    return e*VSURFP*(p-pb);
  }

  PetscScalar SchottyBarrierLowerring (PetscScalar eps, PetscScalar E)
  {
    return sqrt(e/(4*3.1415926535*eps)*E);
  }
  PetscScalar pdSchottyJsn_pdn(PetscScalar n,PetscScalar Tl,PetscScalar Vb)
  {
    PetscScalar VSURFN = ARICHN*Tl*Tl/(e*Nc(Tl));
    return -e*VSURFN;
  }
  PetscScalar pdSchottyJsp_pdp(PetscScalar p,PetscScalar Tl,PetscScalar Vb)
  {
    PetscScalar VSURFP = ARICHP*Tl*Tl/(e*Nv(Tl));
    return e*VSURFP;
  }
  PetscScalar pdSchottyJsn_pdTl(PetscScalar n,PetscScalar Tl,PetscScalar Vb)
  {
    //use finite difference approximate
    PetscScalar dJ = SchottyJsn(n,Tl,Vb)-SchottyJsn(n,(1-1e-10)*Tl,Vb);
    return dJ/(1e-10*Tl);
  }
  PetscScalar pdSchottyJsp_pdTl(PetscScalar p,PetscScalar Tl,PetscScalar Vb)
  {
    //use finite difference approximate
    PetscScalar dJ = SchottyJsp(p,Tl,Vb)-SchottyJsp(p,(1-1e-10)*Tl,Vb);
    return dJ/(1e-10*Tl);
  }


  PetscScalar ThermalVn (PetscScalar Tl)
  {
    // the following two result should be equivalent in mathmatic.
    //return ARICHN*Tl*Tl/(e*Nc(Tl));
    return sqrt(kb*Tl/(2*3.14159265359*EffecElecMass(Tl)));
  }
  AutoDScalar ThermalVn (AutoDScalar Tl)
  {
    // the following two result should be equivalent in mathmatic.
    //return ARICHN*Tl*Tl/(e*Nc(Tl));
    return sqrt(kb*Tl/(2*3.14159265359*EffecElecMass(Tl)));
  }
  PetscScalar ThermalVp (PetscScalar Tl)
  {
    // the following two result should be equivalent in mathmatic.
    //return ARICHP*Tl*Tl/(e*Nv(Tl));
    return sqrt(kb*Tl/(2*3.14159265359*EffecHoleMass(Tl)));
  }
  AutoDScalar ThermalVp (AutoDScalar Tl)
  {
    // the following two result should be equivalent in mathmatic.
    //return ARICHP*Tl*Tl/(e*Nv(Tl));
    return sqrt(kb*Tl/(2*3.14159265359*EffecHoleMass(Tl)));
  }
  PetscScalar pdThermalVn_pdTl (PetscScalar Tl)
  {
    return 0;
  }
  PetscScalar pdThermalVp_pdTl (PetscScalar Tl)
  {
    return 0;
  }

private:
  // [Hot Carrier Injection]
  PetscScalar HCI_LAMHN; // hot-electron scattering mean-free-path
  PetscScalar HCI_LAMHP; // hot-hole scattering mean-free-path

  PetscScalar HCI_Fiegna_A; // Fiegna Constant
  PetscScalar HCI_Fiegna_X; // Fiegna Constant

  PetscScalar HCI_Classical_Lsem_n;   // scattering mean free path in the semiconductor
  PetscScalar HCI_Classical_Lsemr_n;  // redirection mean free path
  PetscScalar HCI_Classical_Lsem_p;   // scattering mean free path in the semiconductor
  PetscScalar HCI_Classical_Lsemr_p;  // redirection mean free path

  void HCI_Init()
  {
    HCI_LAMHN = 9.200000E-07*cm;
    HCI_LAMHP = 1.000000E-07*cm;

    HCI_Fiegna_A = 4.87E+02*m/s/std::pow(eV, 2.5);
    HCI_Fiegna_X = 1.30E+08*std::pow(V/(cm*eV*eV), 1.5);

    HCI_Classical_Lsem_n = 8.9E-07*cm;
    HCI_Classical_Lsemr_n = 6.2E-06*cm;
    HCI_Classical_Lsem_p = 1.0E-07*cm;
    HCI_Classical_Lsemr_p = 6.2E-06*cm;
  }

  PetscScalar Erfc(PetscScalar x)
  {
    // Compute the complementary error function erfc(x).
    // Erfc(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between x and infinity
    //
    //--- Nve 14-nov-1998 UU-SAP Utrecht

    // The parameters of the Chebyshev fit
    const PetscScalar  a1 = -1.26551223,   a2 = 1.00002368;
    const PetscScalar  a3 =  0.37409196,   a4 = 0.09678418;
    const PetscScalar  a5 = -0.18628806,   a6 = 0.27886807;
    const PetscScalar  a7 = -1.13520398,   a8 = 1.48851587;
    const PetscScalar  a9 = -0.82215223,   a10 = 0.17087277;

    PetscScalar v = 1; // The return value
    PetscScalar z = fabs(x);

    if (z <= 0) return v; // erfc(0)=1

    PetscScalar t = 1/(1+0.5*z);

    v = t*exp((-z*z) +a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10)))))))));

    if (x < 0) v = 2-v; // erfc(-x)=2-erfc(x)

    return v;
  }

public:

  PetscScalar HCI_Probability_Semiconductor_n(const PetscScalar &dis)
  {
    if( dis > 30*HCI_LAMHN  ) return 0;
    return exp( - dis/ HCI_LAMHN);
  }

  PetscScalar HCI_Probability_Semiconductor_p(const PetscScalar &dis)
  {
    if( dis > 30*HCI_LAMHP  ) return 0;
    return exp( - dis/ HCI_LAMHP);
  }

  PetscScalar HCI_Integral_Fiegna_n(const PetscScalar &phin, const PetscScalar &Eeff)
  {
    if( HCI_Fiegna_X > 30*Eeff  ) return 0;
    return HCI_Fiegna_A/(3*HCI_Fiegna_X)*std::pow(Eeff, 1.5)/sqrt(phin)*exp(-HCI_Fiegna_X*std::pow(phin, 3.0)/std::pow(Eeff, 1.5));
  }


  PetscScalar HCI_Integral_Fiegna_p(const PetscScalar &phip, const PetscScalar &Eeff)
  {
    if( HCI_Fiegna_X > 30*Eeff  ) return 0;
    return HCI_Fiegna_A/(3*HCI_Fiegna_X)*std::pow(Eeff, 1.5)/sqrt(phip)*exp(-HCI_Fiegna_X*std::pow(phip, 3.0)/std::pow(Eeff, 1.5));
  }


  PetscScalar HCI_Integral_Classical_n(const PetscScalar &phin, const PetscScalar &Eeff)
  {
    if( (HCI_Classical_Lsem_n*Eeff) < phin/30 ) return 0;
    PetscScalar a = phin/(HCI_Classical_Lsem_n*Eeff);
    return 1.0/(2*HCI_Classical_Lsemr_n)*(exp(-a) - sqrt(3.14159265359)*sqrt(a)*Erfc(sqrt(a)));
  }

  PetscScalar HCI_Integral_Classical_p(const PetscScalar &phip, const PetscScalar &Eeff)
  {
    if( (HCI_Classical_Lsem_p*Eeff) < phip/30 ) return 0;
    PetscScalar a = phip/(HCI_Classical_Lsem_p*Eeff);
    return 1.0/(2*HCI_Classical_Lsemr_p)*(exp(-a) - sqrt(3.14159265359)*sqrt(a)*Erfc(sqrt(a)));
  }
private:
  // [band to band Tunneling]
  PetscScalar  A_BTBT;
  PetscScalar  B_BTBT;
  void   BBTunneling_Init()
  {
    A_BTBT = 3.500000E+21*sqrt(e*V)/cm/s/V/V;
    B_BTBT = 2.250000E+07*V/cm/std::pow(e*V,PetscScalar(1.5));
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("A.BTBT",    PARA("A.BTBT",    "The prefactor in Kane's model of band-to-band tunneling .", "eV^(-1/2) cm^-1 s^-1 V^-2", sqrt(e*V)/cm/s/V/V, &A_BTBT)) );
    parameter_map.insert(para_item("B.BTBT",    PARA("B.BTBT",    "The prefactor in the exponential factor of Kane's model of band-to-band tunneling .", "V cm^-1 eV^-(2/3)", V/cm/std::pow(e*V,PetscScalar(1.5)), &B_BTBT)) );
#endif

  }
public:
  //----------------------------------------------------------------
  // band to band Tunneling
  PetscScalar BB_Tunneling(const PetscScalar &Tl, const PetscScalar &E)
  {
    return A_BTBT*E*E/sqrt(Eg(Tl))*exp(-B_BTBT*std::pow(Eg(Tl),PetscScalar(1.5))/(E+1*V/cm));
  }
  AutoDScalar BB_Tunneling(const AutoDScalar &Tl, const AutoDScalar &E)
  {
    return A_BTBT*E*E/sqrt(Eg(Tl))*exp(-B_BTBT*adtl::pow(Eg(Tl),PetscScalar(1.5))/(E+1*V/cm));
  }


  // constructor and destructor
public:
  GSS_PolySi_BandStructure_Sdevice(const PMIS_Environment &env):PMIS_BandStructure(env)
  {
    T300 = 300.0*K;
    PMI_Info = "This is the Sdevice compatable model for band structure parameters of Silicon";
    Eg_Init();
    Lifetime_Init();
    Recomb_Init();
    RelaxTime_Init();
    Schottky_Init();
    HCI_Init();
    BBTunneling_Init();
  }

  ~GSS_PolySi_BandStructure_Sdevice()
  {}
  
  
  void post_calibrate_process()
  {
    Bandgap_Setup();
  }

};


extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_BandStructure*  PMIS_PolySi_BandStructure_Sdevice (const PMIS_Environment& env)
  {
    return new GSS_PolySi_BandStructure_Sdevice(env);
  }
}
