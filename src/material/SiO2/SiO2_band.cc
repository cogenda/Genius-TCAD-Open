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
  PetscScalar AFFINITY;  // The electron affinity for the material.
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

  void   Band_Init()
  {
    AFFINITY = 1.070000e+00*eV;
    BANDGAP  = 9.0*eV;

    HCI_ECN  = 1.650000e+05*V/cm;
    HCI_ECP  = 1.650000e+05*V/cm;
    HCI_BARLN = 2.590000E-04*pow(V*cm, 0.5);
    HCI_TUNLN = 3.000000E-05*pow(V*cm*cm, 1.0/3.0);
    HCI_BARLP = 2.590000E-04*pow(V*cm, 0.5);
    HCI_TUNLP = 3.000000E-05*pow(V*cm*cm, 1.0/3.0);
    HCI_MFP   = 3.200000E-07*cm;

    FN_A      = 6.320000E-07*A/(V*V);
    FN_B      = 2.210000E+08*V/cm;

    ELECMASS   =  0.4*me;
    CBET_alpha = 1.0;
    VBHT_alpha = 1.0;
    VBET_alpha = 1.0;
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("AFFINITY",  PARA("AFFINITY", "The electron affinity for the material", "eV", eV, &AFFINITY)) );
    parameter_map.insert(para_item("BANDGAP",   PARA("BANDGAP",  "The bandgap for the material", "eV", eV, &BANDGAP)) );
    parameter_map.insert(para_item("ELECMASS",  PARA("ELECMASS", "The relative effective mass of electron", "electron mass", me, &ELECMASS)) );
    parameter_map.insert(para_item("CBET.alpha",PARA("CBET.alpha", "Coefficient of the Conduction band electron tunneling", "-", 1.0, &CBET_alpha)) );
    parameter_map.insert(para_item("VBHT.alpha",PARA("VBHT.alpha", "Coefficient of the Valence band hole tunneling", "-", 1.0, &VBHT_alpha)) );
    parameter_map.insert(para_item("VBET.alpha",PARA("VBET.alpha", "Coefficient of the Valence band electron tunneling", "-", 1.0, &VBET_alpha)) );
#endif
  }


public:

  PetscScalar Eg            (const PetscScalar &Tl) const { return BANDGAP;  }

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

  PetscScalar HCI_Barrier_n(const PetscScalar &affinity_semi, const PetscScalar &,
                            const PetscScalar &t_ins, const PetscScalar &E_ins) const
  {
    if(E_ins > 0)
      return (affinity_semi - AFFINITY) - HCI_BARLN*pow(E_ins, 0.5) - HCI_TUNLN*pow(E_ins, 2.0/3.0) ;
    return (affinity_semi - AFFINITY) - E_ins*t_ins;
  }

  PetscScalar HCI_Barrier_p(const PetscScalar &affinity_semi, const PetscScalar & Eg_semi,
                            const PetscScalar &t_ins, const PetscScalar &E_ins) const
  {
    if(E_ins < 0)
      return (AFFINITY + BANDGAP - affinity_semi -  Eg_semi) - HCI_BARLP*pow(fabs(E_ins), 0.5) - HCI_TUNLP*pow(fabs(E_ins), 2.0/3.0) ;
    return (AFFINITY + BANDGAP - affinity_semi -  Eg_semi) + E_ins*t_ins;
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
