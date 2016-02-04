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
// Material Type: 4H-SiC
#include <algorithm>

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

    ELECMASS  = 0.39*me;
    HOLEMASS  = 0.82*me;
    N0_BGN    = 1.000000e+17*std::pow(cm,-3);
    V0_BGN    = 9.000000E-03*eV;
    CON_BGN   = 5.000000e-01;

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("TREF",   PARA("TREF",   "The reference temperature for bandgap model", "K", K, &TREF)) );
    parameter_map.insert(para_item("EGREF",  PARA("EGREF",  "The energy bandgap of the material at TREF", "eV", eV, &EGREF)) );
    parameter_map.insert(para_item("EGALPH", PARA("EGALPH", "The value of alpha used in calculating the temperature depended energy bandgap", "eV/K", eV/K, &EGALPH)) );
    parameter_map.insert(para_item("EGBETA", PARA("EGBETA", "The value of beta used in calculating the temperature depended energy bandgap", "K", K, &EGBETA)) );

    parameter_map.insert(para_item("ELECMASS", PARA("ELECMASS", "The relative effective mass of electron", "electron mass", me, &ELECMASS)) );
    parameter_map.insert(para_item("HOLEMASS", PARA("HOLEMASS", "The relative effective mass of hole", "electron mass", me, &HOLEMASS)) );

    parameter_map.insert(para_item("N0.BGN",   PARA("N0.BGN",   "The concentration parameter used in Slotboom's band-gap narrowing model", "cm^-3", std::pow(cm,-3), &N0_BGN)) );
    parameter_map.insert(para_item("V0.BGN",   PARA("V0.BGN",   "The voltage parameter used in Slotboom's band-gap narrowing model", "V", V, &V0_BGN)) );
    parameter_map.insert(para_item("CON.BGN",  PARA("CON.BGN",  "The const parameter used in Slotboom's band-gap narrowing model", "-", 1.0, &CON_BGN)) );
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
    return 2*4*std::pow(2*pi*EffecElecMass(Tl)*kb*Tl/(h*h),1.5);
  }
  AutoDScalar Nc (const AutoDScalar &Tl)
  {
    return 2*4*adtl::pow(2*pi*EffecElecMass(Tl)*kb*Tl/(h*h),1.5);
  }

  PetscScalar Nv (const PetscScalar &Tl)
  {
    return 2*std::pow(2*pi*EffecHoleMass(Tl)*kb*Tl/(h*h),1.5);
  }
  AutoDScalar Nv (const AutoDScalar &Tl)
  {
    return 2*adtl::pow(2*pi*EffecHoleMass(Tl)*kb*Tl/(h*h),1.5);
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
public:  
  //[incomplete ionization]
  
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
          PetscScalar Eh      = it->second.Eh;
          PetscScalar Ec      = it->second.Ec;
          PetscScalar g       = it->second.GB;
          PetscScalar alpha   = it->second.alpha;
          if( Na > 0.0  )
          {
            if(Na < Na_crit)
            {
              PetscScalar dEh = Eh - alpha*std::pow(Ni, 1.0/3.0);
              PetscScalar dEc = Ec - alpha*std::pow(Ni, 1.0/3.0);
              PetscScalar ph = gamma*Nv*exp(-dEh/(kb*Tl));
              PetscScalar pc = gamma*Nv*exp(-dEc/(kb*Tl));
              Na_eff += 0.5*Na/( 1 + g*p/ph ) + 0.5*Na/( 1 + g*p/pc );
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
          PetscScalar Eh      = it->second.Eh;
          PetscScalar Ec      = it->second.Ec;
          PetscScalar g       = it->second.GB;
          PetscScalar alpha   = it->second.alpha;
          if( Na > 0.0  )
          {
            if(Na < Na_crit)
            {
              PetscScalar dEh = Eh - alpha*std::pow(Ni, 1.0/3.0);
              PetscScalar dEc = Ec - alpha*std::pow(Ni, 1.0/3.0);
              AutoDScalar ph = gamma*Nv*exp(-dEh/(kb*Tl));
              AutoDScalar pc = gamma*Nv*exp(-dEc/(kb*Tl));
              Na_eff += 0.5*Na/( 1 + g*p/ph ) + 0.5*Na/( 1 + g*p/pc );
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
          PetscScalar Eh      = it->second.Eh;
          PetscScalar Ec      = it->second.Ec;
          PetscScalar g       = it->second.GB;
          PetscScalar alpha   = it->second.alpha;
          if( Nd > 0.0  )
          {
            if(Nd < Nd_crit)
            {
              PetscScalar dEh = Eh - alpha*std::pow(Ni, 1.0/3.0);
              PetscScalar dEc = Ec - alpha*std::pow(Ni, 1.0/3.0);
              PetscScalar nh = gamma*Nc*exp(-dEh/(kb*Tl));
              PetscScalar nc = gamma*Nc*exp(-dEc/(kb*Tl));
              Nd_eff += 0.5*Nd/( 1 + g*n/nh ) + 0.5*Nd/( 1 + g*n/nc );
            }
            else
              Nd_eff += Nd;
          }
        }//if( HasVariable(ion_string) )
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
          PetscScalar Eh      = it->second.Eh;
          PetscScalar Ec      = it->second.Ec;
          PetscScalar g       = it->second.GB;
          PetscScalar alpha   = it->second.alpha;
          if( Nd > 0.0  )
          {
            if(Nd < Nd_crit)
            {
              PetscScalar dEh = Eh - alpha*std::pow(Ni, 1.0/3.0);
              PetscScalar dEc = Ec - alpha*std::pow(Ni, 1.0/3.0);
              AutoDScalar nh = gamma*Nc*exp(-dEh/(kb*Tl));
              AutoDScalar nc = gamma*Nc*exp(-dEc/(kb*Tl));
              Nd_eff += 0.5*Nd/( 1 + g*n/nh ) + 0.5*Nd/( 1 + g*n/nc );
            }
            else
              Nd_eff += Nd;
          }
        }//if( HasVariable(ion_string) )
      }
    }

    return Nd_eff;
  }


private:

  struct Species
  {
    std::string name;
    int ion;             // ion type, -1 for p type and +1 for n type
    PetscScalar Eh;      // The constant term used in the calculation of the band ionization energy at hexagonal site
    PetscScalar Ec;      // The constant term used in the calculation of the band ionization energy at cubic site
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
    Species boron    = {"boron",    -1, 0.300*eV, 0.300*eV, 4.0, 0.0*eV*cm,      200.0, 0.950, 1e22*std::pow(cm, -3)  };
    Species aluminum = {"aluminum", -1, 0.230*eV, 0.230*eV, 4.0, 0.0*eV*cm,      200.0, 0.950, 1e22*std::pow(cm, -3)  };
    
    // n type
    Species nitrogen   = {"nitrogen",   1, 0.050*eV, 0.100*eV, 2.0, 0.0*eV*cm, 200.0, 1.000, 1e22*std::pow(cm, -3)  };
    Species phosphorus = {"phosphorus", 1, 0.045*eV, 0.100*eV, 2.0, 0.0*eV*cm, 200.0, 1.000, 1e22*std::pow(cm, -3)  };

    species_map["boron"     ] = boron;
    species_map["aluminum"  ] = aluminum;
    species_map["nitrogen"  ] = nitrogen;
    species_map["phosphorus"] = phosphorus;
    
    // alias
    species_map["BoronActive"     ] = boron;
    species_map["AluminumActive"  ] = aluminum;
    species_map["NitrogenActive"  ] = nitrogen;
    species_map["PhosphorusActive"] = phosphorus;
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
      PetscScalar Eh=0.0;
      PetscScalar Ec=0.0;
      PetscScalar GB=0;
      PetscScalar alpha=0.0;
      PetscScalar beta=0.0;
      PetscScalar gamma=1.0;
      PetscScalar N_crit=1e22;

      for(std::vector<Parser::Parameter>::iterator it = pmi_parameters.begin(); it != pmi_parameters.end();)
      {
        if( it->type() == Parser::INTEGER && it->name() == "ion" )
        { ion = it->get_int(); it = pmi_parameters.erase(it); }
        else if( it->type() == Parser::REAL && it->name() == "ebh" )
        { Eh = it->get_real(); it = pmi_parameters.erase(it); }
        else if( it->type() == Parser::REAL && it->name() == "ebc" )
        { Ec = it->get_real(); it = pmi_parameters.erase(it); }
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

      Species species = { species_name, ion, Eh*eV, Ec*eV, GB, alpha*eV*cm, beta, gamma, N_crit*std::pow(cm, -3) };
      species_map[species_name] = species;
    }
  }
  
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
    IncompleteIonization_Init();
    Lifetime_Init();
    Recomb_Init();
    RelaxTime_Init();
    Schottky_Init();
    BBTunneling_Init();
  }

  ~GSS_SiC4H_BandStructure()
  {}
  
  // set parameters for each band model
  int calibrate(std::vector<Parser::Parameter> & pmi_parameters)
  {
    IncompleteIonization_Setup(pmi_parameters);

    return PMI_Server::calibrate(pmi_parameters);
  }  
  
private:
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
  
  
  
}
;


extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_BandStructure*  PMIS_SiC4H_BandStructure_Default (const PMIS_Environment& env)
  {
    return new GSS_SiC4H_BandStructure(env);
  }
}
