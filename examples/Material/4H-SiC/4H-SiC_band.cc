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
/*                                                                           */
/*  Gong Ding                                                                */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: 4H-SiC


#include "PMI.h"


class GSS_SiC4H_BandStructure : public PMIS_BandStructure
{
  PetscScalar T300;
private:
  //[Bandgap]
  // Bandgap and Effective Density of States
  PetscScalar TREF;     // The energy bandgap of the material at 300 K.
  PetscScalar EGREF;     // The energy bandgap of the material at 300 K.
  PetscScalar EGALPH;    // The value of alpha used in calculating the temperature depended energy bandgap.
  PetscScalar EGBETA;    // The value of beta  used in calculating the temperature depended energy bandgap.
  PetscScalar ELECMASS;  // The relative effective mass of electron
  PetscScalar HOLEMASS;  // The relative effective mass of hole
  PetscScalar NC300;     // The effective density of states in the conduction band at 300K.
  PetscScalar NV300;     // The effective density of states in the valence band at 300K.
  PetscScalar NC_F;      // The parameter for temperature depended effective density of states in the conduction band.
  PetscScalar NV_F;      // The parameter for temperature depended effective density of states in the valence band.
  // Model of Bandgap Narrowing due to Heavy Doping
  PetscScalar N0_BGN;    // The concentration parameter used in Slotboom's band-gap narrowing model.
  PetscScalar V0_BGN;    // The voltage parameter used in Slotboom's band-gap narrowing model.
  PetscScalar CON_BGN;   // The const parameter used in Slotboom's band-gap narrowing model.

  // Init value
  void Eg_Init()
  {
    TREF      = 0.000000E+00*K;
    EGREF     = 3.260000E+00*eV;
    EGALPH    = 3.300000E-04*eV/K;
    EGBETA    = 0.000000E+02*K;

    ELECMASS  = 0.76*me;
    HOLEMASS  = 1.19*me;
    NC300     = 1.680000E+19*std::pow(cm,-3);
    NV300     = 4.580000E+18*std::pow(cm,-3);
    NC_F      = 1.500000e+00;
    NV_F      = 1.500000e+00;
    N0_BGN    = 1.000000e+17*std::pow(cm,-3);
    V0_BGN    = 9.000000E-03*V;
    CON_BGN   = 5.000000e-01*eV;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("TREF",   PARA("TREF",   "The reference temperature for bandgap model", "K", K, &TREF)) );
    parameter_map.insert(para_item("EGREF",  PARA("EGREF",  "The energy bandgap of the material at TREF", "eV", eV, &EGREF)) );
    parameter_map.insert(para_item("EGALPH", PARA("EGALPH", "The value of alpha used in calculating the temperature depended energy bandgap", "eV/K", eV/K, &EGALPH)) );
    parameter_map.insert(para_item("EGBETA", PARA("EGBETA", "The value of beta used in calculating the temperature depended energy bandgap", "K", K, &EGBETA)) );

    parameter_map.insert(para_item("ELECMASS", PARA("ELECMASS", "The relative effective mass of electron", "electron mass", me, &ELECMASS)) );
    parameter_map.insert(para_item("HOLEMASS", PARA("HOLEMASS", "The relative effective mass of hole", "electron mass", me, &HOLEMASS)) );
    parameter_map.insert(para_item("NC300",    PARA("NC300",    "The effective density of states in the conduction band at 300K", "cm^-3", std::pow(cm,-3), &NC300)) );
    parameter_map.insert(para_item("NV300",    PARA("NV300",    "The effective density of states in the valence band at 300K", "cm^-3", std::pow(cm,-3), &NV300)) );
    parameter_map.insert(para_item("NC.F",     PARA("NC.F",     "The parameter for temperature depended effective density of states in the conduction band", "-", 1.0, &NC_F)) );
    parameter_map.insert(para_item("NV.F",     PARA("NV.F",     "The parameter for temperature depended effective density of states in the valence band", "-", 1.0, &NV_F)) );

    parameter_map.insert(para_item("N0.BGN",   PARA("N0.BGN",   "The concentration parameter used in Slotboom's band-gap narrowing model", "cm^-3", std::pow(cm,-3), &N0_BGN)) );
    parameter_map.insert(para_item("V0.BGN",   PARA("V0.BGN",   "The voltage parameter used in Slotboom's band-gap narrowing model", "V", V, &V0_BGN)) );
    parameter_map.insert(para_item("CON.BGN",  PARA("CON.BGN",  "The const parameter used in Slotboom's band-gap narrowing model", "eV", eV, &CON_BGN)) );
#endif
  }
public:
  //---------------------------------------------------------------------------
  // procedure of Bandgap
  PetscScalar Eg (const PetscScalar &Tl)
  {
    return EGREF+EGALPH*(TREF*TREF/(TREF+EGBETA+1e-10) - Tl*Tl/(Tl+EGBETA+1e-10));
  }
  AutoDScalar Eg (const AutoDScalar &Tl)
  {
    return EGREF+EGALPH*(TREF*TREF/(TREF+EGBETA+1e-10) - Tl*Tl/(Tl+EGBETA+1e-10));
  }

  //---------------------------------------------------------------------------
  // procedure of Bandgap Narrowing due to Heavy Doping
  PetscScalar EgNarrow(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N = Na+Nd+1.0*std::pow(cm,-3);
    PetscScalar x = log(N/N0_BGN);
    return V0_BGN*(x+sqrt(x*x+CON_BGN));
  }
  PetscScalar EgNarrowToEc   (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl){return 0.5*EgNarrow(p, n, Tl);}
  PetscScalar EgNarrowToEv   (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl){return 0.5*EgNarrow(p, n, Tl);}

  AutoDScalar EgNarrow(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N = Na+Nd+1.0*std::pow(cm,-3);
    PetscScalar x = log(N/N0_BGN);
    return V0_BGN*(x+sqrt(x*x+CON_BGN));
  }
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
    PetscScalar Nc = NC300*std::pow(Tl/T300,NC_F);
    PetscScalar Nv = NV300*std::pow(Tl/T300,NV_F);
    return sqrt(Nc*Nv)*exp(-bandgap/(2*kb*Tl));
  }
  // nie, Eg narrow should be considered
  PetscScalar nie (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar bandgap = Eg(Tl);
    PetscScalar Nc = NC300*std::pow(Tl/T300,NC_F);
    PetscScalar Nv = NV300*std::pow(Tl/T300,NV_F);
    return sqrt(Nc*Nv)*exp(-bandgap/(2*kb*Tl))*exp(EgNarrow(p, n, Tl));
  }
  AutoDScalar nie (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar bandgap = Eg(Tl);
    AutoDScalar Nc = NC300*adtl::pow(Tl/T300,NC_F);
    AutoDScalar Nv = NV300*adtl::pow(Tl/T300,NV_F);
    AutoDScalar ret=sqrt(Nc*Nv)*exp(-bandgap/(2*kb*Tl))*exp(EgNarrow(p, n, Tl));
    return ret;
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
  PetscScalar NSRHP;            // The Shockley-Read-Hall concentration parameter for holes.
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
    TAUN0     = 5.000000e-07*s;
    TAUP0     = 5.000000e-07*s;
    STAUN     = 0.000000e+00*cm/s;
    STAUP     = 0.000000e+00*cm/s;
    NSRHN     = 3.000000e+30*std::pow(cm,-3);
    AN        = 1.000000e+00;
    BN        = 1.000000e+00;
    CN        = 0.000000e+00;
    EN        = 2.000000e+00;
    NSRHP     = 3.000000e+30*std::pow(cm,-3);
    AP        = 1.000000e+00;
    BP        = 1.000000e+00;
    CP        = 0.000000e+00;
    EP        = 2.000000e+00;
    EXN_TAU   = 3.000000e-01;
    EXP_TAU   = 3.000000e-01;

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

    parameter_map.insert(para_item("EXN.TAU",    PARA("EXN.TAU",    "The exponent of lattice temperature dependent electron lifetime", "-", 1.0, &EXN_TAU)) );
    parameter_map.insert(para_item("EXP.TAU",    PARA("EXP.TAU",    "The exponent of lattice temperature dependent hole lifetime", "-", 1.0, &EXP_TAU)) );
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
  PetscScalar C_DIRECT;      // The band-to-band recombination coefficient.
  // Recombination Including Tunneling
  PetscScalar M_RTUN;        // The trap-assisted tunneling effective mass. *free electron rest mass m0
  PetscScalar S_RTUN;        // Band-to-band field power ratio.
  PetscScalar B_RTUN;        // Band-to-band tunneling rate proportionality factor.
  PetscScalar E_RTUN;        // Band-to-band reference electric field.

  PetscScalar NI_MIN;
  // Init value
  void Recomb_Init()
  {
    ETRAP   = 0.000000e+00*eV;
    AUGN    =  5.000000E-31*std::pow(cm,6)/s;
    AUGP    =  2.000000e-31*std::pow(cm,6)/s;
    C_DIRECT = 0.000000e+00*std::pow(cm,3)/s;
    M_RTUN   = 2.500000e-01;
    S_RTUN   = 2.500000e+00;
    B_RTUN   = 4.000000e+14*std::pow(cm,S_RTUN -3)*std::pow(V,S_RTUN*-1)/s;
    E_RTUN   = 1.900000e+07*V/cm;
    NI_MIN   = 1.0*std::pow(cm,-3);
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("ETRAP",    PARA("ETRAP",    "The trap level (Et - Ei) used in determining the Shockley-Read-Hall recombination rate", "eV", eV, &ETRAP)) );
    parameter_map.insert(para_item("AUGN",     PARA("AUGN",     "The Auger coefficient for electrons", "cm^6/s", std::pow(cm,6)/s, &AUGN)) );
    parameter_map.insert(para_item("AUGP",     PARA("AUGP",     "The Auger coefficient for holes", "cm^6/s", std::pow(cm,6)/s, &AUGP)) );
    parameter_map.insert(para_item("C.DIRECT", PARA("C.DIRECT", "The direct generation/recombination coefficient", "cm^3/s", std::pow(cm,3)/s, &C_DIRECT)) );

    parameter_map.insert(para_item("NI.MIN", PARA("NI.MIN", "Minimum value for Nie", "cm^-3", std::pow(cm,-3), &NI_MIN)) );
//    parameter_map.insert(para_item("M_RTUN",   PARA("M_RTUN", "The trap-assisted tunneling effective mass", "-", 1.0, &M_RTUN)) );
//    parameter_map.insert(para_item("S_RTUN",   PARA("S_RTUN", "Band-to-band field power ratio", "-", 1.0, &S_RTUN)) );
//    parameter_map.insert(para_item("B_RTUN",   PARA("B_RTUN", "Band-to-band tunneling rate proportionality factor", "cm^(S_RTUN-3)V^(-S_RTUN)/s", std::pow(cm,S_RTUN -3)*std::pow(V,-S_RTUN)/s, &B_RTUN)) );
//    parameter_map.insert(para_item("E_RTUN",   PARA("E_RTUN", "Band-to-band reference electric field", "V/cm", V/cm, &E_RTUN)) );
#endif

  }

public:
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
    PetscScalar ni =   std::max(nie(p, n, Tl),NI_MIN);
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
    AutoDScalar ni =   adtl::fmax(nie(p, n, Tl),NI_MIN);
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
        return sqrt(kb*Tl/(2*3.14159265359*EffecElecMass(Tl)));
  }
  AutoDScalar ThermalVn (AutoDScalar Tl)
  {
        return sqrt(kb*Tl/(2*3.14159265359*EffecElecMass(Tl)));
  }
  PetscScalar ThermalVp (PetscScalar Tl)
  {
        return sqrt(kb*Tl/(2*3.14159265359*EffecHoleMass(Tl)));
  }
  AutoDScalar ThermalVp (AutoDScalar Tl)
  {
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
  // [band to band Tunneling]
  PetscScalar  A_BTBT;
  PetscScalar  B_BTBT;
  void   BBTunneling_Init()
  {
    A_BTBT = 0*sqrt(e*V)/cm/s/V/V;
    B_BTBT = 0*V/cm/std::pow(e*V,PetscScalar(1.5));
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
  GSS_SiC4H_BandStructure(const PMIS_Environment &env):PMIS_BandStructure(env)
  {
    T300 = 300.0*K;
    Eg_Init();
    Lifetime_Init();
    Recomb_Init();
    RelaxTime_Init();
    Schottky_Init();
    BBTunneling_Init();
  }

  ~GSS_SiC4H_BandStructure()
  {}
}
;


extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_BandStructure*  PMIS_SiC4H_BandStructure_Default (const PMIS_Environment& env)
  {
    return new GSS_SiC4H_BandStructure(env);
  }
}
