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
#include <algorithm>


#include "PMI.h"
#include "fvm_node_data.h"






class GSS_Si_BandStructure_Sdevice : public PMIS_BandStructure
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

    //parameter_map.insert(para_item("N0.BGN",   PARA("N0.BGN",   "The concentration parameter used in band-gap narrowing model", "cm^-3", std::pow(cm,-3), &N0_BGN)) );
    //parameter_map.insert(para_item("E0.BGN",   PARA("E0.BGN",   "The energy parameter used in band-gap narrowing model", "eV", eV, &E0_BGN)) );

    parameter_map.insert(para_item("BGN.MODEL",PARA("BGN.MODEL","Specific bandgap narrowing model, support Bennett(default), Slotboom, OldSlotboom, delAlamo", &BGN_Model)) );

#endif

  }
public:
  //---------------------------------------------------------------------------
  // procedure of Bandgap
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


  //particle energy to elec-hole pare generation rate
  PetscScalar ParticleQuantumEffect(const PetscScalar &Tl) { return 3.6*eV; }


  //end of Bandgap
  
  //[incomplete ionization]
public:
  //
  int IonType( const std::string & ion_string )
  {
    // convert ion_string to lower case
    std::string ion_string_lower_case(ion_string);
    std::transform(ion_string_lower_case.begin(), ion_string_lower_case.end(), ion_string_lower_case.begin(), ::tolower);

    // exactly match
    if( species_map.find(ion_string) != species_map.end() )
      return species_map.find(ion_string)->second.ion;

    // partly match species
    std::map<std::string, Species>::const_iterator it = species_map.begin();
    for(; it != species_map.end(); ++it)
      if(ion_string_lower_case.find(it->first)!=std::string::npos)
        return it->second.ion;

    return 0;
  }



  PetscScalar Na_II(const PetscScalar &p, const PetscScalar &Tl, bool fermi)
  {
    PetscScalar Ni = ReadDopingNa() + ReadDopingNd();
    PetscScalar Nv = this->Nv(Tl);

    PetscScalar gamma = 1.0;
    if( fermi ) gamma = gamma_f(p/Nv);

    PetscScalar Na_eff = 0.0;

    std::map<std::string, Species>::const_iterator it = species_map.begin();
    for( ; it != species_map.end(); ++it)
    {
      if(it->second.ion < 0)
      {
        const std::string & ion_string = it->first;
        if( HasVariable(ion_string) )
        {
          PetscScalar Na = ReadRealVariable(VariableIndex(ion_string));
          PetscScalar Na_crit = it->second.N_crit;
          PetscScalar E0      = it->second.E0;
          PetscScalar g       = it->second.GB;
          PetscScalar alpha   = it->second.alpha;
          if( Na > 0.0  )
          {
            if(Na < Na_crit)
            {
              PetscScalar dEa = E0 - alpha*std::pow(Ni, 1.0/3.0);
              PetscScalar p1 = gamma*Nv*exp(-dEa/(kb*Tl));
              Na_eff += Na/( 1 + g*p/p1 );
            }
            else
              Na_eff += Na;
          }
        }
      }
    }
    return Na_eff;
  }


  AutoDScalar Na_II(const AutoDScalar &p, const AutoDScalar &Tl, bool fermi)
  {
    PetscScalar Ni = ReadDopingNa() + ReadDopingNd();
    AutoDScalar Nv = this->Nv(Tl);

    AutoDScalar gamma = 1.0;
    if( fermi ) gamma = gamma_f(p/Nv);

    AutoDScalar Na_eff = 0.0;

    std::map<std::string, Species>::const_iterator it = species_map.begin();
    for( ; it != species_map.end(); ++it)
    {
      if(it->second.ion < 0)
      {
        const std::string & ion_string = it->first;
        if( HasVariable(ion_string) )
        {
          PetscScalar Na = ReadRealVariable(VariableIndex(ion_string));
          PetscScalar Na_crit = it->second.N_crit;
          PetscScalar E0      = it->second.E0;
          PetscScalar g       = it->second.GB;
          PetscScalar alpha   = it->second.alpha;
          if( Na > 0.0  )
          {
            if(Na < Na_crit)
            {
              PetscScalar dEa = E0 - alpha*std::pow(Ni, 1.0/3.0);
              AutoDScalar p1 = gamma*Nv*exp(-dEa/(kb*Tl));
              Na_eff += Na/( 1 + g*p/p1 );
            }
            else
              Na_eff += Na;
          }
        }
      }
    }

    return Na_eff;
  }

  PetscScalar Nd_II(const PetscScalar &n, const PetscScalar &Tl, bool fermi)
  {
    PetscScalar Ni = ReadDopingNa() + ReadDopingNd();
    PetscScalar Nc = this->Nc(Tl);

    PetscScalar gamma = 1.0;
    if( fermi ) gamma = gamma_f(n/Nc);

    PetscScalar Nd_eff = 0.0;

    std::map<std::string, Species>::const_iterator it = species_map.begin();
    for( ; it != species_map.end(); ++it)
    {
      if(it->second.ion > 0)
      {
        const std::string & ion_string = it->first;
        if( HasVariable(ion_string) )
        {
          PetscScalar Nd = ReadRealVariable(VariableIndex(ion_string));
          PetscScalar Nd_crit = it->second.N_crit;
          PetscScalar E0      = it->second.E0;
          PetscScalar g       = it->second.GB;
          PetscScalar alpha   = it->second.alpha;
          if( Nd > 0.0  )
          {
            if(Nd < Nd_crit)
            {
              PetscScalar dEd = E0 - alpha*std::pow(Ni, 1.0/3.0);
              PetscScalar n1 = gamma*Nc*exp(-dEd/(kb*Tl));
              Nd_eff += Nd/( 1 + g*n/n1 );
            }
            else
              Nd_eff += Nd;
          }
        }
      }
    }

    return Nd_eff;
  }

  AutoDScalar Nd_II(const AutoDScalar &n, const AutoDScalar &Tl, bool fermi)
  {
    PetscScalar Ni = ReadDopingNa() + ReadDopingNd();
    AutoDScalar Nc = this->Nc(Tl);

    AutoDScalar gamma = 1.0;
    if( fermi ) gamma = gamma_f(n/Nc);

    AutoDScalar Nd_eff = 0.0;

    std::map<std::string, Species>::const_iterator it = species_map.begin();
    for( ; it != species_map.end(); ++it)
    {
      if(it->second.ion > 0)
      {
        const std::string & ion_string = it->first;
        if( HasVariable(ion_string) )
        {
          PetscScalar Nd = ReadRealVariable(VariableIndex(ion_string));
          PetscScalar Nd_crit = it->second.N_crit;
          PetscScalar E0      = it->second.E0;
          PetscScalar g       = it->second.GB;
          PetscScalar alpha   = it->second.alpha;
          if( Nd > 0.0  )
          {
            if(Nd < Nd_crit)
            {
              PetscScalar dEd = E0 - alpha*std::pow(Ni, 1.0/3.0);
              AutoDScalar n1 = gamma*Nc*exp(-dEd/(kb*Tl));
              Nd_eff += Nd/( 1 + g*n/n1 );
            }
            else
              Nd_eff += Nd;
          }
        }
      }
    }

    return Nd_eff;
  }


private:

  struct Species
  {
    std::string name;
    int ion;             // ion type, -1 for p type and +1 for n type
    PetscScalar E0;      // The constant term used in the calculation of the band ionization energy
    PetscScalar GB;      // The band degeneracy factor
    PetscScalar alpha;   // The prefactor for the doping dependent term used in the calculation of the band ionization energy
    PetscScalar beta;    // The prefactor the temperature dependent term used in the calculation of the band ionization energy.
    PetscScalar gamma;   // The exponent of temperature used in the calculation of the band ionization energy.
    PetscScalar N_crit;  // The impurity concentration form which the doping transition from incomplete ionization to complete ionization
  };

  // predefined and user specified species
  std::map<std::string, Species> species_map;

  // Init value
  void IncompleteIonization_Init()
  {
    // p type
    Species boron    = {"boron",    -1, 0.045*eV, 4.0, 3.037e-8*eV*cm, 200.0, 0.950, 1e22*std::pow(cm, -3)  };
    Species aluminum = {"aluminum", -1, 0.067*eV, 4.0, 3.037e-8*eV*cm, 200.0, 0.950, 1e22*std::pow(cm, -3)  };
    Species gallium  = {"gallium",  -1, 0.072*eV, 4.0, 3.037e-8*eV*cm, 200.0, 0.950, 1e22*std::pow(cm, -3)  };
    Species indium   = {"indium",   -1, 0.160*eV, 4.0, 3.037e-8*eV*cm, 200.0, 0.950, 1e22*std::pow(cm, -3)  };

    // n type
    Species nitrogen   = {"nitrogen",   1, 0.045*eV, 2.0, 3.100e-8*eV*cm, 200.0, 1.000, 1e22*std::pow(cm, -3)  };
    Species phosphorus = {"phosphorus", 1, 0.045*eV, 2.0, 3.100e-8*eV*cm, 200.0, 1.000, 1e22*std::pow(cm, -3)  };
    Species arsenic    = {"arsenic",    1, 0.054*eV, 2.0, 3.100e-8*eV*cm, 200.0, 1.000, 1e22*std::pow(cm, -3)  };
    Species antimony   = {"antimony",   1, 0.039*eV, 2.0, 3.100e-8*eV*cm, 200.0, 1.000, 1e22*std::pow(cm, -3)  };

    species_map["boron"     ] = boron;
    species_map["aluminum"  ] = aluminum;
    species_map["gallium"   ] = gallium;
    species_map["indium"    ] = indium;
    species_map["nitrogen"  ] = nitrogen;
    species_map["phosphorus"] = phosphorus;
    species_map["arsenic"   ] = arsenic;
    species_map["antimony"  ] = antimony;
    
    //alias
    species_map["BoronActive"     ] = boron;
    species_map["AluminumActive"  ] = aluminum;
    species_map["GalliumActive"   ] = gallium;
    species_map["IndiumActive"    ] = indium;
    species_map["NitrogenActive"  ] = nitrogen;
    species_map["PhosphorusActive"] = phosphorus;
    species_map["ArsenicActive"   ] = arsenic;
    species_map["AntimonyActive"  ] = antimony;    
  }


  void IncompleteIonization_Setup(std::vector<Parser::Parameter> & pmi_parameters)
  {
    // check if any user defind species
    std::string species_name;
    for(std::vector<Parser::Parameter>::iterator it = pmi_parameters.begin(); it != pmi_parameters.end();)
    {
      if( it->type() == Parser::STRING && it->name() == "species" )
      {
        species_name = it->get_string();
        it = pmi_parameters.erase(it);
      }
      else
        ++it;
    }

    if( !species_name.empty() )
    {
      int ion=0 ;
      PetscScalar E0=0.0;
      PetscScalar GB=0;
      PetscScalar alpha=0.0;
      PetscScalar beta=0.0;
      PetscScalar gamma=1.0;
      PetscScalar N_crit=1e22;

      for(std::vector<Parser::Parameter>::iterator it = pmi_parameters.begin(); it != pmi_parameters.end();)
      {
        if( it->type() == Parser::INTEGER && it->name() == "ion" )
        { ion = it->get_int(); it = pmi_parameters.erase(it); }
        else if( it->type() == Parser::REAL && it->name() == "eb0" )
        { E0 = it->get_real(); it = pmi_parameters.erase(it); }
        else if( it->type() == Parser::REAL && it->name() == "gb" )
        { GB = it->get_real(); it = pmi_parameters.erase(it); }
        else if( it->type() == Parser::REAL && it->name() == "alpha" )
        { alpha = it->get_real(); it = pmi_parameters.erase(it); }
        else if( it->type() == Parser::REAL && it->name() == "beta" )
        { beta = it->get_real(); it = pmi_parameters.erase(it); }
        else if( it->type() == Parser::REAL && it->name() == "gamma" )
        { gamma = it->get_real(); it = pmi_parameters.erase(it); }
        else if( it->type() == Parser::REAL && it->name() == "ncrit" )
        { N_crit = it->get_real(); it = pmi_parameters.erase(it); }
        else
          ++it;
      }

      Species species = { species_name, ion, E0*eV, GB, alpha*eV*cm, beta, gamma, N_crit*std::pow(cm, -3) };
      species_map[species_name] = species;
    }


  }


private:
  //[Lifetime]
  //Lifetimes
  PetscScalar TAUNMIN;         // The Shockley-Read-Hall electron lifetime.
  PetscScalar TAUPMIN;         // The Shockley-Read-Hall hole lifetime.
  PetscScalar TAUNMAX;         // The Shockley-Read-Hall electron lifetime.
  PetscScalar TAUPMAX;         // The Shockley-Read-Hall hole lifetime.
  PetscScalar STAUN;         // The electron surface recombination velocity.
  PetscScalar STAUP;         // The hole surface recombination velocity.
  //Concentration-Dependent Lifetimes
  PetscScalar NSRHN;         // The Shockley-Read-Hall concentration parameter for electrons.
  PetscScalar NSRHP;         // The Shockley-Read-Hall concentration parameter for holes.
  // Lattice Temperature-Dependent Lifetimes
  PetscScalar TAUNTALH;       // The exponent of lattice temperature dependent electron lifetime.
  PetscScalar TAUPTALH;       // The exponent of lattice temperature dependent hole lifetime.

  //Init value
  void Lifetime_Init()
  {
    //Lifetimes
    TAUNMIN = 0.0000e+00*s;         // The Shockley-Read-Hall electron lifetime.
    TAUPMIN = 0.0000e+00*s;         // The Shockley-Read-Hall hole lifetime.
    TAUNMAX = 1.0000e-05*s;         // The Shockley-Read-Hall electron lifetime.
    TAUPMAX = 3.0000e-06*s;         // The Shockley-Read-Hall hole lifetime.

    STAUN   = 0.000000e+00*cm/s;         // The electron surface recombination velocity.
    STAUP   = 0.000000e+00*cm/s;         // The hole surface recombination velocity.
    //Concentration-Dependent Lifetimes
    NSRHN  = 1.0000e+16*std::pow(cm,-3);         // The Shockley-Read-Hall concentration parameter for electrons.
    NSRHP  = 1.0000e+16*std::pow(cm,-3);         // The Shockley-Read-Hall concentration parameter for holes.
    // Lattice Temperature-Dependent Lifetimes
    TAUNTALH = -1.5000e+00;       // The exponent of lattice temperature dependent electron lifetime.
    TAUPTALH = -1.5000e+00;       // The exponent of lattice temperature dependent hole lifetime.


#ifdef __CALIBRATE__
    parameter_map.insert(para_item("TAUNMIN",    PARA("TAUNMIN",    "The Shockley-Read-Hall electron lifetime", "s", s, &TAUNMIN)) );
    parameter_map.insert(para_item("TAUPMIN",    PARA("TAUPMIN",    "The Shockley-Read-Hall hole lifetime", "s", s, &TAUPMIN)) );
    parameter_map.insert(para_item("TAUNMAX",    PARA("TAUNMAX",    "The Shockley-Read-Hall electron lifetime", "s", s, &TAUNMAX)) );
    parameter_map.insert(para_item("TAUPMAX",    PARA("TAUPMAX",    "The Shockley-Read-Hall hole lifetime", "s", s, &TAUPMAX)) );

    parameter_map.insert(para_item("STAUN",    PARA("STAUN",    "The electron surface recombination velocity", "cm/s", cm/s, &STAUN)) );
    parameter_map.insert(para_item("STAUP",    PARA("STAUP",    "The hole surface recombination velocity", "cm/s", cm/s, &STAUP)) );

    parameter_map.insert(para_item("NSRHN", PARA("NSRHN", "The Shockley-Read-Hall concentration parameter for electrons", "cm^-3", std::pow(cm,-3), &NSRHN)) );
    parameter_map.insert(para_item("NSRHP", PARA("NSRHP", "The Shockley-Read-Hall concentration parameter for holes", "cm^-3", std::pow(cm,-3), &NSRHP)) );

    parameter_map.insert(para_item("TAUNTALH",    PARA("TAUNTALH",    "The exponent of lattice temperature dependent electron lifetime", "-", 1.0, &TAUNTALH)) );
    parameter_map.insert(para_item("TAUPTALH",    PARA("TAUPTALH",    "The exponent of lattice temperature dependent hole lifetime", "-", 1.0, &TAUPTALH)) );
#endif

  }




public:
  //---------------------------------------------------------------------------
  // electron lift time for SHR Recombination
  PetscScalar TAUN (const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar taun = TAUNMIN + ( TAUNMAX - TAUNMIN ) / ( 1 + ( (Na+Nd)/NSRHN ));
    return taun * std::pow(Tl/T300, TAUNTALH);
  }

  AutoDScalar TAUN (const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar taun = TAUNMIN + ( TAUNMAX - TAUNMIN ) / ( 1 + ( (Na+Nd)/NSRHN ));
    return taun * adtl::pow(Tl/T300, TAUNTALH);
  }

  //---------------------------------------------------------------------------
  // hole lift time for SHR Recombination
  PetscScalar TAUP (const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar taup = TAUPMIN + ( TAUPMAX - TAUPMIN ) / ( 1 + ( (Na+Nd)/NSRHP ));
    return taup * std::pow(Tl/T300, TAUPTALH);
  }
  AutoDScalar TAUP (const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar taup = TAUPMIN + ( TAUPMAX - TAUPMIN ) / ( 1 + ( (Na+Nd)/NSRHP ));
    return taup * adtl::pow(Tl/T300, TAUPTALH);
  }
  // End of Lifetime

  //[the fit parameter for density-gradient solver]
private:
  PetscScalar GN_DG;   // Fit parameter for Density Gradient model
  PetscScalar GP_DG;   // Fit parameter for Density Gradient model

  void DG_Init()
  {
    GN_DG = 3.6;
    GP_DG = 5.6;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("Gamman",    PARA("Gamman",    "Electron fit parameter for Density Gradient model", "-", 1.0, &GN_DG)) );
    parameter_map.insert(para_item("Gammap",    PARA("Gammap",    "Hole fit parameter for Density Gradient model", "-", 1.0, &GP_DG)) );
#endif

  }

public:

  PetscScalar Gamman         () {return GN_DG;}
  PetscScalar Gammap         () {return GP_DG;}


private:
  //[Recombination]
  // SRH, Auger, and Direct Recombination
  PetscScalar ETRAP;         // The trap level (Et - Ei) used in determining the Shockley-Read-Hall recombination rate.

  PetscScalar AUGAN;          // The Auger coefficient for electrons.
  PetscScalar AUGAP;          // The Auger coefficient for holes.
  PetscScalar AUGBN;          // The Auger coefficient for electrons.
  PetscScalar AUGBP;          // The Auger coefficient for holes.
  PetscScalar AUGCN;          // The Auger coefficient for electrons.
  PetscScalar AUGCP;          // The Auger coefficient for holes.
  PetscScalar AUGHN;          // The Auger coefficient for electrons.
  PetscScalar AUGHP;          // The Auger coefficient for holes.
  PetscScalar AUGN0;          // The Auger coefficient for electrons.
  PetscScalar AUGP0;          // The Auger coefficient for holes.

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

    AUGAN   = 6.7000e-32*std::pow(cm,6)/s;
    AUGAP   = 7.2000e-32*std::pow(cm,6)/s;
    AUGBN   = 2.4500e-31*std::pow(cm,6)/s;
    AUGBP   = 4.5000e-33*std::pow(cm,6)/s;
    AUGCN   = -2.2000e-32*std::pow(cm,6)/s;
    AUGCP   = 2.6300e-32*std::pow(cm,6)/s;
    AUGHN   = 3.46667;
    AUGHP   = 8.25688;
    AUGN0   =1.0000e+18*std::pow(cm,-3);
    AUGP0   =1.0000e+18*std::pow(cm,-3);


    C_DIRECT = 0.000000e+00*std::pow(cm,3)/s;
    M_RTUN   = 2.500000e-01;
    S_RTUN   = 2.500000e+00;
    B_RTUN   = 4.000000e+14*std::pow(cm,S_RTUN-3)*std::pow(V,-S_RTUN)/s;
    E_RTUN   = 1.900000e+07*V/cm;
#ifdef __CALIBRATE__
    parameter_map.insert(para_item("ETRAP",    PARA("ETRAP",    "The trap level (Et - Ei) used in determining the Shockley-Read-Hall recombination rate", "eV", eV, &ETRAP)) );

    parameter_map.insert(para_item("AUGAN",    PARA("AUGAN",     "The Auger coefficient A for electrons", "cm^6/s", std::pow(cm,6)/s, &AUGAN)) );
    parameter_map.insert(para_item("AUGAP",    PARA("AUGAP",     "The Auger coefficient A for holes", "cm^6/s", std::pow(cm,6)/s, &AUGAP)) );
    parameter_map.insert(para_item("AUGBN",    PARA("AUGBN",     "The Auger coefficient B for electrons", "cm^6/s", std::pow(cm,6)/s, &AUGBN)) );
    parameter_map.insert(para_item("AUGBP",    PARA("AUGBP",     "The Auger coefficient B for holes", "cm^6/s", std::pow(cm,6)/s, &AUGBP)) );
    parameter_map.insert(para_item("AUGCN",    PARA("AUGCN",     "The Auger coefficient C for electrons", "cm^6/s", std::pow(cm,6)/s, &AUGCN)) );
    parameter_map.insert(para_item("AUGCP",    PARA("AUGCP",     "The Auger coefficient C for holes", "cm^6/s", std::pow(cm,6)/s, &AUGCP)) );
    parameter_map.insert(para_item("AUGHN",    PARA("AUGHN",     "The Auger coefficient H for electrons", "", 1.0, &AUGHN)) );
    parameter_map.insert(para_item("AUGHP",    PARA("AUGHP",     "The Auger coefficient H for holes", "", 1.0, &AUGHP)) );
    parameter_map.insert(para_item("AUGN0",    PARA("AUGN0",     "The Auger coefficient N0 for electrons", "cm^-3", std::pow(cm,-3), &AUGN0)) );
    parameter_map.insert(para_item("AUGp0",    PARA("AUGP0",     "The Auger coefficient N0 for holes", "cm^-3", std::pow(cm,-3), &AUGP0)) );

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

  PetscScalar AUGERN(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    return (AUGAN + AUGBN * (Tl/T300) + AUGCN * std::pow(Tl/T300, 2)) * (1 + AUGHN * exp(-n/AUGN0));
  }

  AutoDScalar AUGERN(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    return (AUGAN + AUGBN * (Tl/T300) + AUGCN * adtl::pow(Tl/T300, 2)) * (1 + AUGHN * exp(-n/AUGN0));
  }

  PetscScalar AUGERP(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    return (AUGAP + AUGBP * (Tl/T300) + AUGCP * std::pow(Tl/T300, 2)) * (1 + AUGHP * exp(-p/AUGP0));
  }

  AutoDScalar AUGERP(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    return (AUGAP + AUGBP * (Tl/T300) + AUGCP * adtl::pow(Tl/T300, 2)) * (1 + AUGHP * exp(-p/AUGP0));
  }

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
    PetscScalar ni = nie(p, n, Tl);
    PetscScalar dn = p*n-ni*ni;
    return (AUGERN(p,n,Tl)*n + AUGERP(p,n,Tl)*p)*dn;
  }
  AutoDScalar R_Auger     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni = nie(p, n, Tl);
    AutoDScalar dn = p*n-ni*ni;
    return (AUGERN(p,n,Tl)*n + AUGERP(p,n,Tl)*p)*dn;
  }

  //---------------------------------------------------------------------------
  // Electron Auger Recombination
  PetscScalar R_Auger_N     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni = nie(p, n, Tl);
    PetscScalar dn = p*n-ni*ni;
    return (AUGERN(p,n,Tl)*n)*dn;
  }
  AutoDScalar R_Auger_N     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni = nie(p, n, Tl);
    AutoDScalar dn = p*n-ni*ni;
    return (AUGERN(p,n,Tl)*n)*dn;
  }
  //---------------------------------------------------------------------------
  // Hole Auger Recombination
  PetscScalar R_Auger_P     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni = nie(p, n, Tl);
    PetscScalar dn = p*n-ni*ni;
    return (AUGERP(p,n,Tl)*p)*dn;
  }
  AutoDScalar R_Auger_P     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni = nie(p, n, Tl);
    AutoDScalar dn = p*n-ni*ni;
    return (AUGERP(p,n,Tl)*p)*dn;
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
    PetscScalar Raug = (AUGERN(p,n,Tl)*n + AUGERP(p,n,Tl)*p)*dn;
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
    AutoDScalar Raug = (AUGERN(p,n,Tl)*n + AUGERP(p,n,Tl)*p)*dn;
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
  GSS_Si_BandStructure_Sdevice(const PMIS_Environment &env):PMIS_BandStructure(env)
  {
    T300 = 300.0*K;
    PMI_Info = "This is the Sdevice compatable model for band structure parameters of Silicon";
    Eg_Init();
    IncompleteIonization_Init();
    Lifetime_Init();
    Recomb_Init();
    RelaxTime_Init();
    Schottky_Init();
    HCI_Init();
    DG_Init();
    BBTunneling_Init();
  }

  ~GSS_Si_BandStructure_Sdevice()
  {}


  // set parameters for each band model
  int calibrate(std::vector<Parser::Parameter> & pmi_parameters)
  {

    IncompleteIonization_Setup(pmi_parameters);

    return PMI_Server::calibrate(pmi_parameters);
  }


  void post_calibrate_process()
  {
    Bandgap_Setup();
  }

private:
  
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

  /*-----------------------------------------------------------------------
   *   GAMMA calculates f1/2(eta)/exp(eta) according to the approximate
   *   formula in casey's book,dummy arguement x=f1/2(eta).
   */
  PetscScalar gamma_f(PetscScalar x)
  {
    const PetscScalar a=3.53553e-1,b=4.95009e-3,c=1.48386e-4;
    const PetscScalar d=4.42563e-6,pi1=1.772453851e0,pi2=9.869604401e0;
    const PetscScalar VerySmallNumericValue = 1.0e-30;
    const PetscScalar MaximumExponent = 76.0;

    PetscScalar temx;
    if(x>1.0e1)
    {
      temx=sqrt(std::pow(7.5e-1*pi1*x,PetscScalar(4.e0/3.e0))-pi2/6.e0);
      if(x > MaximumExponent)
        return VerySmallNumericValue;
      else
        return x/exp(temx);
    }
    else if(x>0.0)
    {
      temx=x*(a+x*(-b+x*(c-x*d)));
      return 1.0/exp(temx);
    }
    else
      return 1.0;
  }


  AutoDScalar gamma_f(const AutoDScalar &x)
  {
    const PetscScalar a=3.53553e-1,b=4.95009e-3,c=1.48386e-4;
    const PetscScalar d=4.42563e-6,pi1=1.772453851e0,pi2=9.869604401e0;
    const PetscScalar VerySmallNumericValue = 1.0e-30;
    const PetscScalar MaximumExponent = 76.0;

    AutoDScalar temx;
    if(x>1.0e1)
    {
      temx=sqrt(adtl::pow(7.5e-1*pi1*x,PetscScalar(4.e0/3.e0))-pi2/6.e0);
      if(x > MaximumExponent)
        return VerySmallNumericValue;
      else
        return x/exp(temx);
    }
    else if(x>0.0)
    {
      temx=x*(a+x*(-b+x*(c-x*d)));
      return 1.0/exp(temx);
    }
    else
      return 1.0;
  }  
};


extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_BandStructure*  PMIS_Si_BandStructure_Sdevice (const PMIS_Environment& env)
  {
    return new GSS_Si_BandStructure_Sdevice(env);
  }
}
