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


class GSS_AlGaAs_BandStructure : public PMIS_BandStructure
{
private:
  PetscScalar T300;
  //[Bandgap]
  // Bandgap and Effective Density of States
  PetscScalar EG300;     // The energy bandgap of the material at 300 K.
  PetscScalar EG_X0;
  PetscScalar EG_X1;
  PetscScalar EG_X2;
  PetscScalar EG_X3;
  PetscScalar EG_X4;
  PetscScalar EG_X5;
  PetscScalar EG_X6;
  PetscScalar EG_X7;
  PetscScalar EG_X8;
  PetscScalar EG_X9;
  PetscScalar EG_X10;
  PetscScalar EG_X11;
  PetscScalar EG_X12;
  PetscScalar EG_X13;
  PetscScalar EG_X14;

  PetscScalar EGGAMM;
  PetscScalar EGGAX;
  PetscScalar EGGAL;

  PetscScalar EGBETA;
  PetscScalar EGBEX;
  PetscScalar EGBEL;

  PetscScalar EGALPH;
  PetscScalar EGALX;
  PetscScalar EGALL;
  PetscScalar pm;

  PetscScalar NC300;     // The effective density of states in the conduction band at 300K.
  PetscScalar NV300;     // The effective density of states in the valence band at 300K.
  PetscScalar NC_F;      // The parameter for temperature depended effective density of states in the conduction band.
  PetscScalar NV_F;      // The parameter for temperature depended effective density of states in the valence band.
  PetscScalar MEG;
  PetscScalar MEG_X1;
  PetscScalar MEX;
  PetscScalar MEX_X1;
  PetscScalar MEL;
  PetscScalar MEL_X1;
  PetscScalar MH0;
  PetscScalar MH0_X1;
  PetscScalar ML0;
  PetscScalar ML0_X1;
  // Model of Bandgap Narrowing due to Heavy Doping
  PetscScalar N0_BGN;    // The concentration parameter used in Slotboom's band-gap narrowing model.
  PetscScalar V0_BGN;    // The voltage parameter used in Slotboom's band-gap narrowing model.
  PetscScalar CON_BGN;   // The const parameter used in Slotboom's band-gap narrowing model.

  // Init value
  void Eg_Init()
  {
    EG300     = 1.424000e+00*eV;
    EG_X0     = 0.000000E+00*eV;
    EG_X1     = 1.247000E+00*eV;
    EG_X2     = 0.000000E+00*eV;
    EG_X3     = 0.000000E+00*eV;
    EG_X4     = 0.000000E+00*eV;
    EG_X5     = 4.760000E-01*eV;
    EG_X6     = 1.250000E-01*eV;
    EG_X7     = 1.430000E-01*eV;
    EG_X8     = 0.000000E+00*eV;
    EG_X9     = 0.000000E+00*eV;
    EG_X10    = 2.860000E-01*eV;
    EG_X11    = 9.600000E-01*eV;
    EG_X12    = 0.000000E+00*eV;
    EG_X13    = 0.000000E+00*eV;
    EG_X14    = 0.000000E+00*eV;
    pm        = 1.5;
    EGGAMM    = 0.000000E+00*eV/K;
    EGGAX     = 0.000000E+00*eV/K;
    EGGAL     = 0.000000E+00*eV/K;

    EGBETA    = 2.040000E+02*K;
    EGBEX     = 2.040000E+02*K;
    EGBEL     = 2.040000E+02*K;

    EGALPH    = 5.405000E-04*eV/K;
    EGALX     = 4.600000E-04*eV/K;
    EGALL     = 6.050000E-04*eV/K;

    NC300     = 4.700000e+17*std::pow(cm,-3);
    NV300     = 7.000000e+18*std::pow(cm,-3);
    NC_F      = 1.500000e+00;
    NV_F      = 1.500000e+00;

    MEG       = 6.700000E-02*me;
    MEG_X1    = 8.300000E-02*me;
    MEX       = 8.500000E-01*me;
    MEX_X1    = -1.400000E-01*me;
    MEL       = 5.600000E-01*me;
    MEL_X1    = 1.000000E-01*me;
    MH0       = 6.200000E-01*me;
    MH0_X1    = 1.400000E-01*me;
    ML0       = 8.700000E-02*me;
    ML0_X1    = 6.300000E-02*me;

    N0_BGN    = 1.000000e+17*std::pow(cm,-3);
    V0_BGN    = 0.000000e+00*V;
    CON_BGN   = 0.000000e+00*eV;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("EG300",  PARA("EG300",  "The energy bandgap of the material at 300 K", "eV", eV, &EG300)) );
    parameter_map.insert(para_item("EG.X0", PARA("EG.X0", "", "eV",eV , &EG_X0)) );
    parameter_map.insert(para_item("EG.X1", PARA("EG.X1", "", "eV", eV, &EG_X1)) );
    parameter_map.insert(para_item("EG.X2", PARA("EG.X2", "", "eV", eV, &EG_X2)) );
    parameter_map.insert(para_item("EG.X3", PARA("EG.X3", "", "eV",eV , &EG_X3)) );
    parameter_map.insert(para_item("EG.X4", PARA("EG.X4", "", "eV",eV , &EG_X4)) );
    parameter_map.insert(para_item("EG.X5", PARA("EG.X5", "", "eV", eV, &EG_X5)) );
    parameter_map.insert(para_item("EG.X6", PARA("EG.X6", "", "eV", eV, &EG_X6)) );
    parameter_map.insert(para_item("EG.X7", PARA("EG.X7", "", "eV", eV, &EG_X7)) );
    parameter_map.insert(para_item("EG.X8", PARA("EG.X8", "", "eV", eV, &EG_X8)) );
    parameter_map.insert(para_item("EG.X9", PARA("EG.X9", "", "eV", eV, &EG_X9)) );
    parameter_map.insert(para_item("EG.X10", PARA("EG.X10", "", "eV", eV, &EG_X10)) );
    parameter_map.insert(para_item("EG.X11", PARA("EG.X11", "", "eV", eV, &EG_X11)) );
    parameter_map.insert(para_item("EG.X12", PARA("EG.X12", "", "eV", eV, &EG_X12)) );
    parameter_map.insert(para_item("EG.X13", PARA("EG.X13", "", "eV", eV, &EG_X13)) );
    parameter_map.insert(para_item("EG.X14", PARA("EG.X14", "", "eV", eV, &EG_X14)) );
    parameter_map.insert(para_item("pm", PARA("pm", "", "-",1.0 , &pm)) );

    parameter_map.insert(para_item("EGGAMM", PARA("EGGAMM", "", "eV/K", eV/K, &EGGAMM)) );
    parameter_map.insert(para_item("EGGAX", PARA("EGGAX", "", "eV/K",eV/K , &EGGAX)) );
    parameter_map.insert(para_item("EGGAL", PARA("EGGAL", "", "eV/K",eV/K , &EGGAL)) );

    parameter_map.insert(para_item("EGBETA", PARA("EGBETA", "", "K", K, &EGBETA)) );
    parameter_map.insert(para_item("EGBEX", PARA("EGBEX", "", "K",K , &EGBEX)) );
    parameter_map.insert(para_item("EGBEL", PARA("EGBEL", "", "K", K, &EGBEL)) );

    parameter_map.insert(para_item("EGALPH", PARA("EGALPH", "", "eV/K",eV/K , &EGALPH)) );
    parameter_map.insert(para_item("EGALX", PARA("EGALX", "", "eV/K", eV/K, &EGALX)) );
    parameter_map.insert(para_item("EGALL", PARA("EGALL", "", "eV/K",eV/K , &EGALL)) );

    parameter_map.insert(para_item("NC300",    PARA("NC300",    "The effective density of states in the conduction band at 300K", "cm^-3", std::pow(cm,-3), &NC300)) );
    parameter_map.insert(para_item("NV300",    PARA("NV300",    "The effective density of states in the valence band at 300K", "cm^-3", std::pow(cm,-3), &NV300)) );
    parameter_map.insert(para_item("NC.F",     PARA("NC.F",     "The parameter for temperature depended effective density of states in the conduction band", "-", 1.0, &NC_F)) );
    parameter_map.insert(para_item("NV.F",     PARA("NV.F",     "The parameter for temperature depended effective density of states in the valence band", "-", 1.0, &NV_F)) );

    parameter_map.insert(para_item("MEG", PARA("MEG", "", "electron mass", me, &MEG)) );
    parameter_map.insert(para_item("MEG.X1", PARA("MEG.X1", "", "electron mass", me, &MEG_X1)) );
    parameter_map.insert(para_item("MEX", PARA("MEX", "", "electron mass", me, &MEX)) );
    parameter_map.insert(para_item("MEX.X1", PARA("MEX.X1", "", "electron mass", me, &MEX_X1)) );
    parameter_map.insert(para_item("MEL", PARA("MEL", "", "electron mass", me, &MEL)) );
    parameter_map.insert(para_item("MEL.X1", PARA("MEL.X1", "", "electron mass", me, &MEL_X1)) );
    parameter_map.insert(para_item("MH0", PARA("MH0", "", "electron mass", me, &MH0)) );
    parameter_map.insert(para_item("MH0.X1", PARA("MH0.X1", "", "electron mass", me, &MH0_X1)) );
    parameter_map.insert(para_item("ML0", PARA("ML0", "", "electron mass", me, &ML0)) );
    parameter_map.insert(para_item("ML0.X1", PARA("ML0.X1", "", "electron mass", me, &ML0_X1)) );

    parameter_map.insert(para_item("N0.BGN",   PARA("N0.BGN",   "The concentration parameter used in Slotboom's band-gap narrowing model", "cm^-3", std::pow(cm,-3), &N0_BGN)) );
    parameter_map.insert(para_item("V0.BGN",   PARA("V0.BGN",   "The voltage parameter used in Slotboom's band-gap narrowing model", "V", V, &V0_BGN)) );
    parameter_map.insert(para_item("CON.BGN",  PARA("CON.BGN",  "The const parameter used in Slotboom's band-gap narrowing model", "eV", eV, &CON_BGN)) );
#endif
  }
public:
  //---------------------------------------------------------------------------
  // we need calculate all the bandgap valley and choose lowest
  PetscScalar E_Gamma(const PetscScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X0 + EG_X1*x + EG_X2*x*x + EG_X3*x*x*x + EG_X4*x*x*x*x
               + (T300*T300/(T300+EGBETA)-Tl*Tl/(Tl+EGBETA))*(EGALPH+EGGAMM*x);
  }
  AutoDScalar E_Gamma(const AutoDScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X0 + EG_X1*x + EG_X2*x*x + EG_X3*x*x*x + EG_X4*x*x*x*x
               + (T300*T300/(T300+EGBETA)-Tl*Tl/(Tl+EGBETA))*(EGALPH+EGGAMM*x);
  }

  PetscScalar E_X(const PetscScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X5 + EG_X6*x + EG_X7*x*x + EG_X8*x*x*x + EG_X9*x*x*x*x
               + (T300*T300/(T300+EGBEX)-Tl*Tl/(Tl+EGBEX))*(EGALX+EGGAX*x);
  }
  AutoDScalar E_X(const AutoDScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X5 + EG_X6*x + EG_X7*x*x + EG_X8*x*x*x + EG_X9*x*x*x*x
               + (T300*T300/(T300+EGBEX)-Tl*Tl/(Tl+EGBEX))*(EGALX+EGGAX*x);
  }

  PetscScalar E_L(const PetscScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X10 + EG_X11*x + EG_X12*x*x + EG_X13*x*x*x + EG_X14*x*x*x*x
               + (T300*T300/(T300+EGBEL)-Tl*Tl/(Tl+EGBEL))*(EGALL+EGGAL*x);
  }
  AutoDScalar E_L(const AutoDScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X10 + EG_X11*x + EG_X12*x*x + EG_X13*x*x*x + EG_X14*x*x*x*x
               + (T300*T300/(T300+EGBEL)-Tl*Tl/(Tl+EGBEL))*(EGALL+EGGAL*x);
  }

  //---------------------------------------------------------------------------
  // procedure of Bandgap, return the lowest valley
  PetscScalar Eg (const PetscScalar &Tl)
  {
    PetscScalar Eg1 = E_Gamma(Tl);
    PetscScalar Eg2 = E_X(Tl);
    return Eg1 < Eg2 ? Eg1 : Eg2;
  }
  AutoDScalar Eg (const AutoDScalar &Tl)
  {
    AutoDScalar Eg1 = E_Gamma(Tl);
    AutoDScalar Eg2 = E_X(Tl);
    return fmin( Eg1 , Eg2);
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
  PetscScalar EffecElecMass(const PetscScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        PetscScalar bandgap = Eg(Tl);

        PetscScalar m_Gamma = std::pow(std::pow(MEG+MEG_X1*x,pm)*exp((bandgap-E_Gamma(Tl))/(kb*Tl)),1.0/pm);
        PetscScalar m_X = std::pow(std::pow(MEX+MEX_X1*x,pm)*exp((bandgap-E_X(Tl))/(kb*Tl)),1.0/pm);
        PetscScalar m_L = std::pow(std::pow(MEL+MEL_X1*x,pm)*exp((bandgap-E_L(Tl))/(kb*Tl)),1.0/pm);
        return std::pow(std::pow(m_Gamma,pm)+std::pow(m_X,pm)+std::pow(m_L,pm),1.0/pm);
  }
  AutoDScalar EffecElecMass(const AutoDScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        AutoDScalar bandgap = Eg(Tl);

        AutoDScalar m_Gamma = adtl::pow(std::pow(MEG+MEG_X1*x,pm)*exp((bandgap-E_Gamma(Tl))/(kb*Tl)),1.0/pm);
        AutoDScalar m_X = adtl::pow(std::pow(MEX+MEX_X1*x,pm)*exp((bandgap-E_X(Tl))/(kb*Tl)),1.0/pm);
        AutoDScalar m_L = adtl::pow(std::pow(MEL+MEL_X1*x,pm)*exp((bandgap-E_L(Tl))/(kb*Tl)),1.0/pm);
        return adtl::pow(adtl::pow(m_Gamma,pm)+adtl::pow(m_X,pm)+adtl::pow(m_L,pm),1.0/pm);
  }

  PetscScalar EffecHoleMass(const PetscScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return std::pow(std::pow(MH0+MH0_X1*x,pm)+std::pow(ML0+ML0_X1*x,pm),1.0/pm);
  }
  AutoDScalar EffecHoleMass(const AutoDScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return std::pow(std::pow(MH0+MH0_X1*x,pm)+std::pow(ML0+ML0_X1*x,pm),1.0/pm);
  }


  //---------------------------------------------------------------------------
  // Nc and Nv calculated from effective mass
  PetscScalar Nc (const PetscScalar &Tl)
  {
    //return NC300*std::pow(Tl/T300,NC_F);
    return 2*std::pow(2*3.14159265359*EffecElecMass(Tl)*kb*Tl/(h*h),pm);
  }
  AutoDScalar Nc (const AutoDScalar &Tl)
  {
    //return NC300*std::pow(Tl/T300,NC_F);
    return 2*adtl::pow(2*3.14159265359*EffecElecMass(Tl)*kb*Tl/(h*h),pm);
  }

  PetscScalar Nv (const PetscScalar &Tl)
  {
    //return NV300*std::pow(Tl/T300,NV_F);
    return 2*std::pow(2*3.14159265359*EffecHoleMass(Tl)*kb*Tl/(h*h),pm);
  }
  AutoDScalar Nv (const AutoDScalar &Tl)
  {
    //return NV300*std::pow(Tl/T300,NV_F);
    return 2*adtl::pow(2*3.14159265359*EffecHoleMass(Tl)*kb*Tl/(h*h),pm);
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
    return sqrt(Nc*Nv)*exp(-bandgap/(2*kb*Tl))*exp(EgNarrow(p, n, Tl));
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
    TAUN0     = 1.000000e-09*s;
    TAUP0     = 1.000000e-09*s;
    STAUN     = 0.000000e+00*cm/s;
    STAUP     = 0.000000e+00*cm/s;
    NSRHN     = 5.000000e+16*std::pow(cm,-3);
    AN        = 1.000000e+00;
    BN        = 0.000000e+00;
    CN        = 0.000000e+00;
    EN        = 2.000000e+00;
    NSRHP     = 5.000000e+16*std::pow(cm,-3);
    AP        = 1.000000e+00;
    BP        = 0.000000e+00;
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

  // Init value
  void Recomb_Init()
  {
    //Source: Semiconductors on NSM
    ETRAP   =  0.000000e+00*eV;
    AUGN    =  1.000000e-31*std::pow(cm,6)/s;
    AUGP    =  1.000000e-30*std::pow(cm,6)/s;
    C_DIRECT = 1.800000e-10*std::pow(cm,3)/s;
    M_RTUN   = 2.500000e-01;
    S_RTUN   = 2.000000e+00;
    B_RTUN   = 0.000000e+00*std::pow(cm,S_RTUN -3)*std::pow(V,S_RTUN*-1)/s;
    E_RTUN   = 0.000000e+00*V/cm;
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
   WTN0 =  2.400000E-12*s;
   WTN1 =  4.000000E-13*s;
   WTN2 =  0.000000E+00*s;
   WTN3 =  0.000000E+00*s;
   WTN4 =  0.000000E+00*s;
   WTN5 =  0.000000E+00*s;
   WTNL =  6.800000E-13*s;
   TNL  =  1.866270E+03*K;
   WTP0 =  0.000000E+00*s;
   WTP1 =  0.000000E+00*s;
   WTP2 =  0.000000E+00*s;
   WTP3 =  0.000000E+00*s;
   WTP4 =  0.000000E+00*s;
   WTP5 =  0.000000E+00*s;
   WTPL =  1.000000E-12*s;
   TPL  =  0.000000E+00*K;
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("WTN0",    PARA("WTN0",    "Electron energy relaxatioin time when electron temperature is higher than lattice temperature by TNL.", "s", s, &WTN0)) );
    parameter_map.insert(para_item("WTN1",    PARA("WTN1",    "Electron energy relaxatioin time when electrons are nearly in equilibrium with lattice.", "s", s, &WTN1)) );
    parameter_map.insert(para_item("TNL",     PARA("TNL",     "Electron temperature reference scale.", "K", K, &TNL)) );
    parameter_map.insert(para_item("WTPL",    PARA("WTPL",    "Hole energy relaxation time.", "s", s, &WTPL)) );
#endif
  }
public:
  //---------------------------------------------------------------------------
  // Electron relaxation time for EBM
  PetscScalar ElecEnergyRelaxTime(const PetscScalar &Tn,const PetscScalar &Tl)
  {
    PetscScalar r = (Tn-Tl)/TNL;
    return WTN1+(WTN0-WTN1)*r*r*exp(2-2*r);
  }
  AutoDScalar ElecEnergyRelaxTime(const AutoDScalar &Tn,const AutoDScalar &Tl)
  {
    AutoDScalar r = (Tn-Tl)/TNL;
    return WTN1+(WTN0-WTN1)*r*r*exp(2-2*r);
  }

  //---------------------------------------------------------------------------
  // Hole relaxation time for EBM
  PetscScalar HoleEnergyRelaxTime(const PetscScalar &Tp,const PetscScalar &Tl)
  {
    return WTPL;
  }
  AutoDScalar HoleEnergyRelaxTime(const AutoDScalar &Tp,const AutoDScalar &Tl)
  {
    return WTPL;
  }

  // end of energy relax time

private:
  // [Schottky]
  PetscScalar ARICHN;
  PetscScalar ARICHP;
  PetscScalar VSURFN;   // Thermionic emission velocity of electron
  PetscScalar VSURFP;

  void   Schottky_Init()
  {
    ARICHN = 6.285700e+00*A/(K*cm)/(K*cm);
    ARICHP = 1.050000e+02*A/(K*cm)/(K*cm);
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("ARICHN",    PARA("ARICHN",    "Effective Richardson constant for electrons.", "A K^-2 cm^-2", A/(K*cm)/(K*cm), &ARICHN)) );
    parameter_map.insert(para_item("ARICHP",    PARA("ARICHP",    "Effective Richardson constant for holes.", "A K^-2 cm^-2", A/(K*cm)/(K*cm), &ARICHP)) );
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
  PetscScalar BB_Tunneling(const PetscScalar &Tl,const  PetscScalar &E)
  {
     return A_BTBT*E*E/sqrt(Eg(Tl))*exp(-B_BTBT*std::pow(Eg(Tl),PetscScalar(1.5))/(E+1*V/cm));
  }
  AutoDScalar BB_Tunneling(const AutoDScalar &Tl,const  AutoDScalar &E)
  {
     return A_BTBT*E*E/sqrt(Eg(Tl))*exp(-B_BTBT*adtl::pow(Eg(Tl),PetscScalar(1.5))/(E+1*V/cm));
  }


// constructor and destructor
public:
  GSS_AlGaAs_BandStructure(const PMIS_Environment &env):PMIS_BandStructure(env)
  {
    T300 = 300.0*K;
    Eg_Init();
    Lifetime_Init();
    Recomb_Init();
    RelaxTime_Init();
    Schottky_Init();
    BBTunneling_Init();
  }

  ~GSS_AlGaAs_BandStructure()
  {}
}
;


extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_BandStructure*  PMIS_AlGaAs_BandStructure_Default (const PMIS_Environment& env)
  {
    return new GSS_AlGaAs_BandStructure(env);
  }
}
