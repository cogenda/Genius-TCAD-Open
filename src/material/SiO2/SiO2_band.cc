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
// Material Type: SiO2


#include "PMI.h"

class GSS_SiO2_BandStructure : public PMII_BandStructure
{
private:
  PetscScalar BANDGAP;   // The bandgap for the material.

  PetscScalar HCI_ECN;   // critical electric field for electron scattering in the insulator
  PetscScalar HCI_ECP;   // critical electric field for hole scattering in the insulator
  PetscScalar HCI_BARLN; // represents barrier lowering effects to electron due to the image field
  PetscScalar HCI_TUNLN; // accounts for the possibility of tunneling to electron
  PetscScalar HCI_BARLP; // represents barrier lowering effects to hole due to the image field
  PetscScalar HCI_TUNLP; // accounts for the possibility of tunneling to hole
  PetscScalar HCI_MFP;   // mean-free-path for electrons in silicon dioxide

  PetscScalar FN_A;      // Coefficient of the pre-exponential term for the Fowler-Nordheim  tunneling model
  PetscScalar FN_B;      // Coefficient of the exponential term for the Fowler-Nordheim tunneling model.

  PetscScalar ELECMASS;   // The relative effective mass of electron
  PetscScalar CBET_alpha; // Coefficient of the Conduction band electron tunneling
  PetscScalar VBHT_alpha; // Coefficient of the Valence band hole tunneling
  PetscScalar VBET_alpha; // Coefficient of the Valence band electron tunneling


  PetscScalar GN_DG;   // Fit parameter for Density Gradient model
  PetscScalar GP_DG;   // Fit parameter for Density Gradient model

  PetscScalar TrapA;         // trapA density
  PetscScalar TrapAElecCS;   // trapA cross section of electron
  PetscScalar TrapAHoleCS;   // trapA cross section of hole

  PetscScalar TrapB;         // trapB density
  PetscScalar TrapBElecCS;   // trapB cross section of electron
  PetscScalar TrapBHoleCS;   // trapB cross section of hole
  PetscScalar TrapBHIonR;    // trapB H+ release rate

  PetscScalar InterfaceState;    // interface state density
  PetscScalar HIonInterfaceCS;   // H+ reaction cross section with Si-H at interface

  void   Band_Init()
  {
    BANDGAP  = 9.0*eV;

    HCI_ECN  = 1.650000e+05*V/cm;
    HCI_ECP  = 1.650000e+05*V/cm;
    HCI_BARLN = 2.590000E-04*std::pow(V*cm, 0.5);
    HCI_TUNLN = 3.000000E-05*std::pow(V*cm*cm, 1.0/3.0);
    HCI_BARLP = 2.590000E-04*std::pow(V*cm, 0.5);
    HCI_TUNLP = 3.000000E-05*std::pow(V*cm*cm, 1.0/3.0);
    HCI_MFP   = 3.200000E-07*cm;

    FN_A      = 6.320000E-07*A/(V*V);
    FN_B      = 2.210000E+08*V/cm;

    ELECMASS   =  0.4*me;
    CBET_alpha = 1.0;
    VBHT_alpha = 1.0;
    VBET_alpha = 1.0;

    GN_DG      = 5.0;
    GP_DG      = 20.0;

    TrapA       = 5e18*std::pow(cm, -3);
    TrapAElecCS = 1e-13*cm*cm;
    TrapAHoleCS = 1e-14*cm*cm;

    TrapB       = 1e16*std::pow(cm, -3);
    TrapBElecCS = 1e-12*cm*cm;
    TrapBHoleCS = 1e-14*cm*cm;
    TrapBHIonR  = 1e-5/s;

    InterfaceState  = 1e13*std::pow(cm, -2);
    HIonInterfaceCS = 1e-11*cm*cm;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("BANDGAP",   PARA("BANDGAP",  "The bandgap for the material", "eV", eV, &BANDGAP)) );
    parameter_map.insert(para_item("ELECMASS",  PARA("ELECMASS", "The relative effective mass of electron", "electron mass", me, &ELECMASS)) );

    parameter_map.insert(para_item("FN.A",PARA("FN.A", "Coefficient of the pre-exponential term for the Fowler-Nordheim  tunneling model", "A/V^2", A/(V*V), &FN_A)) );
    parameter_map.insert(para_item("FN.B",PARA("FN.B", "Coefficient of the exponential term for the Fowler-Nordheim tunneling model", "V/cm", V/cm, &FN_B)) );

    parameter_map.insert(para_item("CBET.alpha",PARA("CBET.alpha", "Coefficient of the Conduction band electron tunneling", "-", 1.0, &CBET_alpha)) );
    parameter_map.insert(para_item("VBHT.alpha",PARA("VBHT.alpha", "Coefficient of the Valence band hole tunneling", "-", 1.0, &VBHT_alpha)) );
    parameter_map.insert(para_item("VBET.alpha",PARA("VBET.alpha", "Coefficient of the Valence band electron tunneling", "-", 1.0, &VBET_alpha)) );
    parameter_map.insert(para_item("Gamman",    PARA("Gamman",    "Electron fit parameter for Density Gradient model", "-", 1.0, &GN_DG)) );
    parameter_map.insert(para_item("Gammap",    PARA("Gammap",    "Hole fit parameter for Density Gradient model", "-", 1.0, &GP_DG)) );

    parameter_map.insert(para_item("TrapADensity",PARA("TrapADensity","Trap A density", "cm^-3", 1.0/(cm*cm*cm), &TrapA)) );
    parameter_map.insert(para_item("TrapAElecCS",PARA("TrapAElecCS","Positive trap A capture cross section of electron", "cm^2", cm*cm, &TrapAElecCS)) );
    parameter_map.insert(para_item("TrapAHoleCS",PARA("TrapAHoleCS","Trap A capture cross section of hole", "cm^2", cm*cm, &TrapAHoleCS)) );

    parameter_map.insert(para_item("TrapBDensity",PARA("TrapBDensity","Trap B density", "cm^-3", 1.0/(cm*cm*cm), &TrapB)) );
    parameter_map.insert(para_item("TrapBElecCS",PARA("TrapBElecCS","Positive trap B capture cross section of electron", "cm^2", cm*cm, &TrapBElecCS)) );
    parameter_map.insert(para_item("TrapBHoleCS",PARA("TrapBHoleCS","Trap B capture cross section of hole", "cm^2", cm*cm, &TrapBHoleCS)) );
    parameter_map.insert(para_item("TrapBHIonRelease",PARA("TrapBHIonRelease","Trap B H+ release rate", "1/s", 1.0/s, &TrapBHIonR)) );

    parameter_map.insert(para_item("InterfaceState",PARA("InterfaceState","Interface state density", "cm^-2", 1.0/(cm*cm), &InterfaceState)) );
    parameter_map.insert(para_item("HIonInterfaceCS",PARA("HIonInterfaceCS","H+ reaction cross section with Si-H at interface", "cm^2", cm*cm, &HIonInterfaceCS)) );

#endif
  }


public:

  PetscScalar Eg            (const PetscScalar &Tl) const { return BANDGAP;  }
  PetscScalar EffecElecMass (const PetscScalar &Tl) const { return ELECMASS; }
  PetscScalar EffecHoleMass  (const PetscScalar &Tl) const { return ELECMASS; }

  PetscScalar TrapADensity() const { return TrapA; }
  PetscScalar TrapACaptureElecCS(const PetscScalar &Tl) const { return TrapAElecCS; }
  PetscScalar TrapACaptureHoleCS(const PetscScalar &Tl) const { return TrapAHoleCS; }

  PetscScalar TrapBDensity() const { return TrapB; }
  PetscScalar TrapBCaptureElecCS(const PetscScalar &Tl) const { return TrapBElecCS; }
  PetscScalar TrapBCaptureHoleCS(const PetscScalar &Tl) const { return TrapBHoleCS; }
  PetscScalar TrapBReleaseHIonRate(const PetscScalar &Tl) const { return TrapBHIonR; }

  PetscScalar InterfaceStateDensity() const { return InterfaceState; }
  PetscScalar HIonInterfaceTrapCS(const PetscScalar &Tl) const { return HIonInterfaceCS; }

  PetscScalar ARichardson() const { return 1.100000e+02*A/(K*cm)/(K*cm); }

  /**
   * @return electron fit parameter of Density Gradient solver
   */
  PetscScalar Gamman         () const {return GN_DG;}

  /**
   * @return hole fit parameter of Density Gradient solver
   */
  PetscScalar Gammap         () const {return GP_DG;}


  PetscScalar HCI_Probability_Insulator_n(const PetscScalar &t_ins, const PetscScalar &E_ins) const
  {
    if(E_ins > 0)
    {
      PetscScalar alpha = std::min(t_ins/HCI_MFP, sqrt(HCI_ECN/E_ins) );
      if( alpha > 30.0 ) return 0;
      return exp( -alpha  );
    }
    return exp( -t_ins/HCI_MFP  );
  }

  PetscScalar HCI_Probability_Insulator_p(const PetscScalar &t_ins, const PetscScalar &E_ins) const
  {
    if(E_ins < 0)
    {
      PetscScalar alpha = std::min(t_ins/HCI_MFP, sqrt(HCI_ECN/(-E_ins)) );
      if( alpha > 30.0 ) return 0;
      return exp( -alpha  );
    }
    return exp( -t_ins/HCI_MFP  );
  }

  PetscScalar HCI_Barrier_n(const PetscScalar &affinity_semi, const PetscScalar &, const PetscScalar &affinity_ins,
                            const PetscScalar &t_ins, const PetscScalar &E_ins) const
  {
    if(E_ins > 0)
      return (affinity_semi - affinity_ins) - HCI_BARLN*std::pow(E_ins, 0.5) - HCI_TUNLN*std::pow(E_ins, 2.0/3.0) ;
    return (affinity_semi - affinity_ins) - E_ins*t_ins;
  }

  PetscScalar HCI_Barrier_p(const PetscScalar &affinity_semi, const PetscScalar & Eg_semi, const PetscScalar &affinity_ins,
                            const PetscScalar &t_ins, const PetscScalar &E_ins) const
  {
    if(E_ins < 0)
      return (affinity_ins + BANDGAP - affinity_semi -  Eg_semi) - HCI_BARLP*std::pow(fabs(E_ins), 0.5) - HCI_TUNLP*std::pow(fabs(E_ins), 2.0/3.0) ;
    return (affinity_ins + BANDGAP - affinity_semi -  Eg_semi) + E_ins*t_ins;
  }


  PetscScalar J_FN_Tunneling(const PetscScalar &E_ins, const PetscScalar &alpha) const
  {
    PetscScalar E = fabs(E_ins) + 1*V/cm;
    if(alpha == 1.0)
    {
      if( FN_B/E > 70.0 ) return 0.0;
      return FN_A*E_ins*E_ins*exp( - FN_B/E );
    }
    else
    {
      if( FN_B/E*(1-std::pow(1.0-alpha, 1.5)) > 70.0 ) return 0.0;
      return FN_A*E_ins*E_ins/std::pow(1.0-sqrt(1.0-alpha), 2.0)*exp( - FN_B/E*(1-std::pow(1.0-alpha, 1.5)) );
    }
  }

public:

  PetscScalar J_CBET_Tunneling(const PetscScalar &m, const PetscScalar &Tl,
                               const PetscScalar &Efn1, const PetscScalar &Efn2,
                               const PetscScalar &Ec1,  const PetscScalar &Ec2,
                               const PetscScalar &B1,   const PetscScalar &B2,
                               const PetscScalar &t) const
  {
    std::vector<PetscScalar> JE;
    const PetscScalar dE = 0.02*kb*Tl/e;
    PetscScalar E=std::max(Ec1, Ec2);
    // direct tunneling
    for(; E<std::min(B1, B2); E+=dE)
    {
       PetscScalar phi1 = B1 - E;
       PetscScalar phi2 = B2 - E;
       PetscScalar TE = exp( - 4.0/3.0*sqrt(2*ELECMASS)/hbar*t*(phi1 + sqrt(phi1*phi2)  + phi2 )/(sqrt(phi1) + sqrt(phi2)) );
       PetscScalar f1E = 1.0/(1.0+exp((E-Efn1)/(kb*Tl)));
       PetscScalar f2E = 1.0/(1.0+exp((E-Efn2)/(kb*Tl)));
       PetscScalar g1E = (E >= Ec1 ? 1.0 : 0.0);
       PetscScalar g2E = (E >= Ec2 ? 1.0 : 0.0);
       JE.push_back(TE*(f1E-f2E)*g1E*g2E);
    }
    // FN tunneling
    for(; E<std::max(B1, B2); E+=dE)
    {
      PetscScalar phi = std::max(B1, B2) - E;
      PetscScalar tbar = phi/std::abs(B1-B2)*t;
      PetscScalar TE = exp( - 4.0/3.0*sqrt(2*ELECMASS)/hbar*tbar*(sqrt(phi)) );
      PetscScalar f1E = 1.0/(1.0+exp((E-Efn1)/(kb*Tl)));
      PetscScalar f2E = 1.0/(1.0+exp((E-Efn2)/(kb*Tl)));
      PetscScalar g1E = (E >= Ec1 ? 1.0 : 0.0);
      PetscScalar g2E = (E >= Ec2 ? 1.0 : 0.0);
      JE.push_back(TE*(f1E-f2E)*g1E*g2E);
    }

    JE.push_back(JE[JE.size()-1]); // prevent negative value of JE.size()-1

    PetscScalar J_CBET = 0.0;
    for(unsigned int n=0; n<JE.size()-1; n++)
    {
      J_CBET += 0.5*(JE[n]+JE[n+1])*dE;
    }
    return J_CBET*CBET_alpha*e*m*kb*Tl/(2*hbar*hbar*hbar);
  }


  AutoDScalar J_CBET_Tunneling(const PetscScalar &m, const AutoDScalar &Tl,
                               const AutoDScalar &Efn1, const AutoDScalar &Efn2,
                               const AutoDScalar &Ec1,  const AutoDScalar &Ec2,
                               const AutoDScalar &B1,   const AutoDScalar &B2,
                               const PetscScalar &t) const
  {
    std::vector<AutoDScalar> JE;
    const AutoDScalar dE = 0.02*kb*Tl/e;
    AutoDScalar E=adtl::fmax(Ec1, Ec2);

    // direct tunneling
    for(; E<adtl::fmin(B1, B2); E+=dE)
    {
      AutoDScalar phi1 = B1 - E;
      AutoDScalar phi2 = B2 - E;
      AutoDScalar TE = exp( - 4.0/3.0*sqrt(2*ELECMASS)/hbar*t*(phi1 + sqrt(phi1*phi2)  + phi2 )/(sqrt(phi1) + sqrt(phi2)) );
      AutoDScalar f1E = 1.0/(1.0+exp((E-Efn1)/(kb*Tl)));
      AutoDScalar f2E = 1.0/(1.0+exp((E-Efn2)/(kb*Tl)));
      PetscScalar g1E = (E >= Ec1 ? 1.0 : 0.0);
      PetscScalar g2E = (E >= Ec2 ? 1.0 : 0.0);

      JE.push_back(TE*(f1E-f2E)*g1E*g2E);
    }
    // FN tunneling
    for(; E<adtl::fmax(B1, B2); E+=dE)
    {
      AutoDScalar phi = std::max(B1, B2) - E;
      AutoDScalar tbar = phi/adtl::fabs(B1-B2)*t;
      AutoDScalar TE = exp( - 4.0/3.0*sqrt(2*ELECMASS)/hbar*tbar*(sqrt(phi)) );
      AutoDScalar f1E = 1.0/(1.0+exp((E-Efn1)/(kb*Tl)));
      AutoDScalar f2E = 1.0/(1.0+exp((E-Efn2)/(kb*Tl)));
      PetscScalar g1E = (E >= Ec1 ? 1.0 : 0.0);
      PetscScalar g2E = (E >= Ec2 ? 1.0 : 0.0);
      JE.push_back(TE*(f1E-f2E)*g1E*g2E);
    }

    JE.push_back(JE[JE.size()-1]); // prevent negative value of JE.size()-1

    AutoDScalar J_CBET = 0.0;
    for(unsigned int n=0; n<JE.size()-1; n++)
    {
      J_CBET += 0.5*(JE[n]+JE[n+1])*dE;
    }
    return J_CBET*CBET_alpha*e*m*kb*Tl/(2*hbar*hbar*hbar);
  }


  PetscScalar J_VBHT_Tunneling(const PetscScalar &m, const PetscScalar &Tl,
                               const PetscScalar &Efp1, const PetscScalar &Efp2,
                               const PetscScalar &Ev1,  const PetscScalar &Ev2,
                               const PetscScalar &B1,   const PetscScalar &B2,
                               const PetscScalar &t) const
  {
    std::vector<PetscScalar> JE;
    const PetscScalar dE = 0.02*kb*Tl/e;
    PetscScalar E=std::min(Ev1, Ev2);
    // direct tunneling
    for(; E>std::max(B1, B2); E-=dE)
    {
      PetscScalar phi1 = E - B1;
      PetscScalar phi2 = E - B2;
      PetscScalar TE = exp( - 4.0/3.0*sqrt(2*ELECMASS)/hbar*t*(phi1 + sqrt(phi1*phi2)  + phi2 )/(sqrt(phi1) + sqrt(phi2)) );
      PetscScalar f1E = 1.0/(1.0+exp((Efp1-E)/(kb*Tl)));
      PetscScalar f2E = 1.0/(1.0+exp((Efp2-E)/(kb*Tl)));
      PetscScalar g1E = (E <= Ev1 ? 1.0 : 0.0);
      PetscScalar g2E = (E <= Ev2 ? 1.0 : 0.0);
      JE.push_back(TE*(f1E-f2E)*g1E*g2E);
    }

    // FN tunneling
    for(; E>std::min(B1, B2); E-=dE)
    {
      PetscScalar phi = E - std::min(B1, B2);
      PetscScalar tbar = phi/std::abs(B1-B2)*t;
      PetscScalar TE = exp( - 4.0/3.0*sqrt(2*ELECMASS)/hbar*tbar*(sqrt(phi)) );
      PetscScalar f1E = 1.0/(1.0+exp((Efp1-E)/(kb*Tl)));
      PetscScalar f2E = 1.0/(1.0+exp((Efp2-E)/(kb*Tl)));
      PetscScalar g1E = (E <= Ev1 ? 1.0 : 0.0);
      PetscScalar g2E = (E <= Ev2 ? 1.0 : 0.0);
      JE.push_back(TE*(f1E-f2E)*g1E*g2E);
    }

    JE.push_back(JE[JE.size()-1]); // prevent negative value of JE.size()-1

    PetscScalar J_VBET = 0.0;
    for(unsigned int n=0; n<JE.size()-1; n++)
    {
      J_VBET += 0.5*(JE[n]+JE[n+1])*dE;
    }
    return J_VBET*VBHT_alpha*e*m*kb*Tl/(2*hbar*hbar*hbar);
  }

  AutoDScalar J_VBHT_Tunneling(const PetscScalar &m, const AutoDScalar &Tl,
                               const AutoDScalar &Efp1, const AutoDScalar &Efp2,
                               const AutoDScalar &Ev1,  const AutoDScalar &Ev2,
                               const AutoDScalar &B1,   const AutoDScalar &B2,
                               const PetscScalar &t) const
  {
    std::vector<AutoDScalar> JE;
    const AutoDScalar dE = 0.02*kb*Tl/e;
    AutoDScalar E=adtl::fmin(Ev1, Ev2);
    // direct tunneling
    for(; E>adtl::fmax(B1, B2); E-=dE)
    {
      AutoDScalar phi1 = E - B1;
      AutoDScalar phi2 = E - B2;
      AutoDScalar TE = exp( - 4.0/3.0*sqrt(2*ELECMASS)/hbar*t*(phi1 + sqrt(phi1*phi2)  + phi2 )/(sqrt(phi1) + sqrt(phi2)) );
      AutoDScalar f1E = 1.0/(1.0+exp((Efp1-E)/(kb*Tl)));
      AutoDScalar f2E = 1.0/(1.0+exp((Efp2-E)/(kb*Tl)));
      PetscScalar g1E = (E <= Ev1 ? 1.0 : 0.0);
      PetscScalar g2E = (E <= Ev2 ? 1.0 : 0.0);
      JE.push_back(TE*(f1E-f2E)*g1E*g2E);
    }

    // FN tunneling
    for(; E>adtl::fmin(B1, B2); E-=dE)
    {
      AutoDScalar phi = E - adtl::fmin(B1, B2);
      AutoDScalar tbar = phi/adtl::fabs(B1-B2)*t;
      AutoDScalar TE = exp( - 4.0/3.0*sqrt(2*ELECMASS)/hbar*tbar*(sqrt(phi)) );
      AutoDScalar f1E = 1.0/(1.0+exp((Efp1-E)/(kb*Tl)));
      AutoDScalar f2E = 1.0/(1.0+exp((Efp2-E)/(kb*Tl)));
      PetscScalar g1E = (E <= Ev1 ? 1.0 : 0.0);
      PetscScalar g2E = (E <= Ev2 ? 1.0 : 0.0);
      JE.push_back(TE*(f1E-f2E)*g1E*g2E);
    }

    JE.push_back(JE[JE.size()-1]); // prevent negative value of JE.size()-1

    AutoDScalar J_VBET = 0.0;
    for(unsigned int n=0; n<JE.size()-1; n++)
    {
      J_VBET += 0.5*(JE[n]+JE[n+1])*dE;
    }
    return J_VBET*VBHT_alpha*e*m*kb*Tl/(2*hbar*hbar*hbar);
  }




  PetscScalar J_VBET_Tunneling(const PetscScalar &m, const PetscScalar &Tl,
                               const PetscScalar &Efn1, const PetscScalar &Efn2,
                               const PetscScalar &Ec1,  const PetscScalar &Ec2,
                               const PetscScalar &Ev1,  const PetscScalar &Ev2,
                               const PetscScalar &B1,   const PetscScalar &B2,
                               const PetscScalar &t) const
  {
    std::vector<PetscScalar> JE;
    const PetscScalar dE = 0.02*kb*Tl/e;
    PetscScalar E=std::min(Ev1, Ev2);
    // direct tunneling
    for(; E<std::min(B1, B2); E+=dE)
    {
      PetscScalar phi1 = B1 - E;
      PetscScalar phi2 = B2 - E;
      PetscScalar TE = exp( - 4.0/3.0*sqrt(2*ELECMASS)/hbar*t*(phi1 + sqrt(phi1*phi2)  + phi2 )/(sqrt(phi1) + sqrt(phi2)) );
      PetscScalar f1E = 1.0/(1.0+exp((E-Efn1)/(kb*Tl)));
      PetscScalar f2E = 1.0/(1.0+exp((E-Efn2)/(kb*Tl)));
      PetscScalar g   = ((E >= Ec1 && E <= Ev2) ||  (E >= Ec2 && E <= Ev1)) ? 1.0 : 0.0;
      JE.push_back(TE*(f1E-f2E)*g);
    }

    JE.push_back(JE[JE.size()-1]); // prevent negative value of JE.size()-1

    PetscScalar J_VBET = 0.0;
    for(unsigned int n=0; n<JE.size()-1; n++)
    {
      J_VBET += 0.5*(JE[n]+JE[n+1])*dE;
    }
    return J_VBET*VBET_alpha*e*m*kb*Tl/(2*hbar*hbar*hbar);
  }



  AutoDScalar J_VBET_Tunneling(const PetscScalar &m, const AutoDScalar &Tl,
                               const AutoDScalar &Efn1, const AutoDScalar &Efn2,
                               const AutoDScalar &Ec1,  const AutoDScalar &Ec2,
                               const AutoDScalar &Ev1,  const AutoDScalar &Ev2,
                               const AutoDScalar &B1,   const AutoDScalar &B2,
                               const PetscScalar &t) const
  {
    std::vector<AutoDScalar> JE;
    const AutoDScalar dE = 0.02*kb*Tl/e;
    AutoDScalar E=adtl::fmin(Ev1, Ev2);
    // direct tunneling
    for(; E<adtl::fmin(B1, B2); E+=dE)
    {
      AutoDScalar phi1 = B1 - E;
      AutoDScalar phi2 = B2 - E;
      AutoDScalar TE = exp( - 4.0/3.0*sqrt(2*ELECMASS)/hbar*t*(phi1 + sqrt(phi1*phi2)  + phi2 )/(sqrt(phi1) + sqrt(phi2)) );
      AutoDScalar f1E = 1.0/(1.0+exp((E-Efn1)/(kb*Tl)));
      AutoDScalar f2E = 1.0/(1.0+exp((E-Efn2)/(kb*Tl)));
      PetscScalar g   = ((E >= Ec1 && E <= Ev2) ||  (E >= Ec2 && E <= Ev1)) ? 1.0 : 0.0;
      JE.push_back(TE*(f1E-f2E)*g);
    }

    JE.push_back(JE[JE.size()-1]); // prevent negative value of JE.size()-1

    AutoDScalar J_VBET = 0.0;
    for(unsigned int n=0; n<JE.size()-1; n++)
    {
      J_VBET += 0.5*(JE[n]+JE[n+1])*dE;
    }
    return J_VBET*VBET_alpha*e*m*kb*Tl/(2*hbar*hbar*hbar);
  }

public:
  GSS_SiO2_BandStructure(const PMII_Environment &env):PMII_BandStructure(env)
  {
    Band_Init();
  }
  ~GSS_SiO2_BandStructure()
  {
  }
}
;

extern "C"
{
  DLL_EXPORT_DECLARE  PMII_BandStructure* PMII_SiO2_BandStructure_Default (const PMII_Environment& env)
  {
    return new GSS_SiO2_BandStructure(env);
  }
}
