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
#include "mathfunc.h"

class GSS_Si_BandStructure_Schenk : public PMIS_BandStructure
{
private:
  PetscScalar T300;
  //[Bandgap]
  // Bandgap and Effective Density of States
  PetscScalar EG0;       // The energy bandgap of the material at 0 K.
  PetscScalar EG300;     // The energy bandgap of the material at 300 K.
  PetscScalar EGALPH;    // The value of alpha used in calculating the temperature depended energy bandgap.
  PetscScalar EGBETA;    // The value of beta  used in calculating the temperature depended energy bandgap.
  PetscScalar ELECMASS;  // The relative effective mass of electron
  PetscScalar HOLEMASS;  // The relative effective mass of hole
  PetscScalar NC300;     // The effective density of states in the conduction band at 300K.
  PetscScalar NV300;     // The effective density of states in the valence band at 300K.
  PetscScalar NC_F;      // The parameter for temperature depended effective density of states in the conduction band.
  PetscScalar NV_F;      // The parameter for temperature depended effective density of states in the valence band.

  // Schenk Model of Bandgap Narrowing

  // Table 1 of Schenk's paper.
  PetscScalar me_eff;  // effective electron mass me/m0
  PetscScalar mh_eff;  // effective hole mass mh/m0
  PetscScalar ge;      // degeneracy factor of electron
  PetscScalar gh;      // degeneracy factor of hole
  PetscScalar alphae;
  PetscScalar alphah;
  PetscScalar Ryex;
  PetscScalar aex;
  PetscScalar epss;

  // Table 2
  PetscScalar be;
  PetscScalar bh;
  PetscScalar ce;
  PetscScalar ch;
  PetscScalar de;
  PetscScalar dh;
  PetscScalar pe;
  PetscScalar ph;

  // Table 3
  PetscScalar he;
  PetscScalar hh;
  PetscScalar je;
  PetscScalar jh;
  PetscScalar ke;
  PetscScalar kh;
  PetscScalar qe;
  PetscScalar qh;

  PetscScalar n_scale;

  // Init value
  void Eg_Init()
  {
    // Use parameters from Green (JAP 67, p.2945, 1990) for
    // silicon bandgap and densities of states
    EG0       = 1.16964*eV;
    EG300     = 1.1241*eV;
    EGALPH    = 2.73E-4*eV/K;
    EGBETA    = 0.0*K;

    ELECMASS  = 1.0903*me;
    HOLEMASS  = 1.1525*me;
    NC300     = 2.86E19*std::pow(cm,-3);
    NV300     = 3.10E19*std::pow(cm,-3);
    NC_F      = 1.58;
    NV_F      = 1.85;

    // Table 1 of Schenk's paper.
    me_eff = 0.321;  // effective electron mass me/m0
    mh_eff = 0.346;
    alphae = mh_eff / (me_eff+mh_eff);
    alphah = me_eff / (me_eff+mh_eff);
    ge = 12.0;
    gh = 4.0;
    Ryex = 16.55e-3 * eV; // Excitonic Rydberg
    aex = 37.19e-8 * cm;  // Excitonic Bohr radius
    epss = 11.7;

    // Table 2
    be = 8.0;
    bh = 1.0;
    ce = 1.3346;
    ch = 1.2365;
    de = 0.893;
    dh = 1.153;
    pe = 7.0 / 30.0;
    ph = 7.0 / 30.0;

    // Table 3
    he = 3.91;
    hh = 4.20;
    je = 2.8585;
    jh = 2.9307;
    ke = 0.012;
    kh = 0.19;
    qe = 3.0 / 4.0;
    qh = 1.0 / 4.0;

    n_scale = std::pow(aex, -3.0);

#ifdef __CALIBRATE__
    parameter_map.insert(para_item("EG0",    PARA("EG0",    "The energy bandgap of the material at 0 K", "eV", eV, &EG0)) );
    parameter_map.insert(para_item("EG300",  PARA("EG300",  "The energy bandgap of the material at 300 K", "eV", eV, &EG300)) );
    parameter_map.insert(para_item("EGALPH", PARA("EGALPH", "The value of alpha used in calculating the temperature depended energy bandgap", "eV/K", eV/K, &EGALPH)) );
    parameter_map.insert(para_item("EGBETA", PARA("EGBETA", "The value of beta used in calculating the temperature depended energy bandgap", "K", K, &EGBETA)) );

    parameter_map.insert(para_item("ELECMASS", PARA("ELECMASS", "The relative effective mass of electron", "electron mass", me, &ELECMASS)) );
    parameter_map.insert(para_item("HOLEMASS", PARA("HOLEMASS", "The relative effective mass of hole", "electron mass", me, &HOLEMASS)) );
    parameter_map.insert(para_item("NC300",    PARA("NC300",    "The effective density of states in the conduction band at 300K", "cm^-3", std::pow(cm,-3), &NC300)) );
    parameter_map.insert(para_item("NV300",    PARA("NV300",    "The effective density of states in the valence band at 300K", "cm^-3", std::pow(cm,-3), &NV300)) );
    parameter_map.insert(para_item("NC.F",     PARA("NC.F",     "The parameter for temperature depended effective density of states in the conduction band", "-", 1.0, &NC_F)) );
    parameter_map.insert(para_item("NV.F",     PARA("NV.F",     "The parameter for temperature depended effective density of states in the valence band", "-", 1.0, &NV_F)) );

    parameter_map.insert(para_item("BGN.ME", PARA("BGN.ME", "Density of state effective mass for electrons in bandgap narrowing model.", "-", 1.0, &me_eff)));
    parameter_map.insert(para_item("BGN.MH", PARA("BGN.MH", "Density of state effective mass for holes in bandgap narrowing model.", "-", 1.0, &mh_eff)));
    parameter_map.insert(para_item("BGN.GE", PARA("BGN.GE", "Degeneracy factor for electrons in bandgap narrowing model.", "-", 1.0, &ge )));
    parameter_map.insert(para_item("BGN.GH", PARA("BGN.GH", "Degeneracy factor for holes in bandgap narrowing model.", "-", 1.0, &gh )));
    parameter_map.insert(para_item("BGN.RYEX", PARA("BGN.RYEX", "Exitonic Rydberg", "eV", eV, &Ryex )));
    parameter_map.insert(para_item("BGN.AEX", PARA("BGN.AEX", "Exitonic Bohr radius", "cm", cm, &aex)));

    parameter_map.insert(para_item("BGN.BE", PARA("BGN.BE", "First fitting parameter for electrons in carrier-carrier bandgap narrowing model", "-", 1.0, &be)));
    parameter_map.insert(para_item("BGN.BH", PARA("BGN.BH", "First fitting parameter for holes in carrier-carrier bandgap narrowing model", "-", 1.0, &bh)));
    parameter_map.insert(para_item("BGN.CE", PARA("BGN.CE", "Second fitting parameter for electrons in carrier-carrier bandgap narrowing model", "-", 1.0, &ce)));
    parameter_map.insert(para_item("BGN.CH", PARA("BGN.CH", "Second fitting parameter for holes in carrier-carrier bandgap narrowing model", "-", 1.0, &ch)));
    parameter_map.insert(para_item("BGN.DE", PARA("BGN.DE", "Third fitting parameter for electrons in carrier-carrier bandgap narrowing model", "-", 1.0, &de)));
    parameter_map.insert(para_item("BGN.DH", PARA("BGN.DH", "Third fitting parameter for holes in carrier-carrier bandgap narrowing model", "-", 1.0, &dh)));
    parameter_map.insert(para_item("BGN.PE", PARA("BGN.PE", "Forth fitting parameter for electrons in carrier-carrier bandgap narrowing model", "-", 1.0, &pe)));
    parameter_map.insert(para_item("BGN.PH", PARA("BGN.PH", "Forth fitting parameter for holes in carrier-carrier bandgap narrowing model", "-", 1.0, &ph)));

    parameter_map.insert(para_item("BGN.HE", PARA("BGN.HE", "First fitting parameter for electrons in ionic bandgap narrowing model", "-", 1.0, &he)));
    parameter_map.insert(para_item("BGN.HH", PARA("BGN.HH", "First fitting parameter for holes in ionic bandgap narrowing model", "-", 1.0, &hh)));
    parameter_map.insert(para_item("BGN.JE", PARA("BGN.JE", "Second fitting parameter for electrons in ionic bandgap narrowing model", "-", 1.0, &je)));
    parameter_map.insert(para_item("BGN.JH", PARA("BGN.JH", "Second fitting parameter for holes in ionic bandgap narrowing model", "-", 1.0, &jh)));
    parameter_map.insert(para_item("BGN.KE", PARA("BGN.KE", "Third fitting parameter for electrons in ionic bandgap narrowing model", "-", 1.0, &ke)));
    parameter_map.insert(para_item("BGN.KH", PARA("BGN.KH", "Third fitting parameter for holes in ionic bandgap narrowing model", "-", 1.0, &kh)));
    parameter_map.insert(para_item("BGN.QE", PARA("BGN.QE", "Forth fitting parameter for electrons in ionic bandgap narrowing model", "-", 1.0, &qe)));
    parameter_map.insert(para_item("BGN.QH", PARA("BGN.QH", "Forth fitting parameter for holes in ionic bandgap narrowing model", "-", 1.0, &qh)));
#endif

  }

  void Eg_PostCalibrateProcess()
  {
    alphae = mh_eff / (me_eff+mh_eff);
    alphah = me_eff / (me_eff+mh_eff);
  }

public:


  //---------------------------------------------------------------------------
  // procedure of Bandgap
  PetscScalar Eg (const PetscScalar &Tl)
  {
    //return EG0 - EGALPH*Tl*Tl / (EGBETA + Tl);
    return EG300+EGALPH*(T300*T300/(T300+EGBETA) - Tl*Tl/(Tl+EGBETA));
  }
  AutoDScalar Eg (const AutoDScalar &Tl)
  {
    //return EG0 - EGALPH*Tl*Tl / (EGBETA + Tl);
    return EG300+EGALPH*(T300*T300/(T300+EGBETA) - Tl*Tl/(Tl+EGBETA));
  }


  //---------------------------------------------------------------------------
  // procedure of Bandgap Narrowing by Schenk model
  PetscScalar EgNarrow(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    return EgNarrowToEc(p, n, Tl) + EgNarrowToEv(p, n, Tl);
  }

  PetscScalar EgNarrowToEc   (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    const PetscScalar pi = 3.1415927;

    PetscScalar ndS = ReadDopingNd() / n_scale;
    PetscScalar naS = ReadDopingNa() / n_scale;
    PetscScalar niS = ndS + naS;

    PetscScalar neS = std::abs(n) / n_scale;
    PetscScalar nhS = std::abs(p) / n_scale;

    PetscScalar F = kb * Tl / Ryex;    // [adim]

    PetscScalar nsigma = neS + nhS;
    PetscScalar np = alphae * neS + alphah * nhS;

    PetscScalar nsigma2 = ndS + naS;                // For SCR:  nsigma = niS
    PetscScalar np2 = alphae * ndS + alphah * naS;  // For SCR: np = alphae * NDOP + alphah * PDOP

    PetscScalar Ui = nsigma * nsigma / pow(F,3.0);

    // BGN due to exchange-correlation (carrier-carrier)
    PetscScalar De_xc = -( pow(4.0 * pi, 3.0) * nsigma*nsigma * ( pow(48.0 * neS / pi / ge, 1.0/3.0) + ce * log(1.0 + de * pow(np, pe))  )
                           +(8.0 * pi * alphae/ge)* neS * F*F
                           +sqrt(8.0 *pi*nsigma) * pow(F, 5.0 / 2.0)
                         ) / (pow(4.0*pi, 3.0) * nsigma*nsigma + pow(F, 3.0) + be * sqrt(nsigma) * F*F + 40.0 * pow(nsigma, 3.0/2.0) * F);

    PetscScalar De_i = -( niS * (1.0 + Ui))
                       /( sqrt(F * nsigma2 / (2.0 * pi)) * (1.0 + he * log(1.0 + sqrt(nsigma2) / F) )
                          +je * Ui * pow(np2, 3.0/4.0) * (1.0 + ke * pow(np2, qe))
                        );
#if defined(HAVE_FENV_H) && defined(DEBUG)
    genius_assert( !fetestexcept(FE_INVALID) );
#endif
    return  - (De_xc + De_i)*Ryex;
  }

  PetscScalar EgNarrowToEv   (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    const PetscScalar pi = 3.1415927;

    PetscScalar ndS = ReadDopingNd() / n_scale;
    PetscScalar naS = ReadDopingNa() / n_scale;
    PetscScalar niS = ndS + naS;

    PetscScalar neS = std::abs(n) / n_scale;
    PetscScalar nhS = std::abs(p) / n_scale;

    PetscScalar F = kb * Tl / Ryex;    // [adim]

    PetscScalar nsigma = neS + nhS;
    PetscScalar np = alphae * neS + alphah * nhS;

    PetscScalar nsigma2 = ndS + naS;                // For SCR:  nsigma = niS
    PetscScalar np2 = alphae * ndS + alphah * naS;  // For SCR: np = alphae * NDOP + alphah * PDOP

    PetscScalar Ui = nsigma * nsigma / pow(F,3.0);

    PetscScalar Dh_xc = -( pow(4.0 * pi, 3.0) * nsigma*nsigma * ( pow(48.0 * nhS / pi / gh, 1.0/3.0) + ch * log(1.0 + dh * pow(np, ph))  )
                           +(8.0 * pi * alphah/gh)* nhS * F*F
                           +sqrt(8.0 *pi*nsigma) * pow(F, 5.0 / 2.0)
                         ) / (pow(4.0*pi, 3.0) * nsigma*nsigma + pow(F,3.0) + bh * sqrt(nsigma) * F*F + 40.0 * pow(nsigma, 3.0/2.0) * F);
    PetscScalar Dh_i = -( niS * (1.0 + Ui))
                       /( sqrt(F * nsigma2 / (2.0 * pi)) * (1.0 + hh * log(1.0 + sqrt(nsigma2) / F) )
                          +jh * Ui * pow(np2, 3.0/4.0) * (1.0 + kh * pow(np2, qh))
                        );

    return  - (Dh_xc +  Dh_i) * Ryex;
  }


  AutoDScalar EgNarrow(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    return EgNarrowToEc(p, n, Tl) + EgNarrowToEv(p, n, Tl);
  }

  AutoDScalar EgNarrowToEc   (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    const PetscScalar pi = 3.1415927;

    PetscScalar ndS = ReadDopingNd() / n_scale;
    PetscScalar naS = ReadDopingNa() / n_scale;
    PetscScalar niS = ndS + naS;

    AutoDScalar neS = adtl::fabs(n) / n_scale;
    AutoDScalar nhS = adtl::fabs(p) / n_scale;

    AutoDScalar F = kb * Tl / Ryex;    // [adim]

    AutoDScalar nsigma = neS + nhS;
    AutoDScalar np = alphae * neS + alphah * nhS;

    PetscScalar nsigma2 = ndS + naS;                // For SCR:  nsigma = niS
    PetscScalar np2 = alphae * ndS + alphah * naS;  // For SCR: np = alphae * NDOP + alphah * PDOP

    AutoDScalar Ui = nsigma * nsigma / pow(F,3.0);

    // BGN due to exchange-correlation (carrier-carrier)
    AutoDScalar De_xc = -( pow(4.0 * pi, 3.0) * nsigma*nsigma * ( pow(48.0 * neS / pi / ge, 1.0/3.0) + ce * log(1.0 + de * pow(np, pe))  )
                           +(8.0 * pi * alphae/ge)* neS * F*F
                           +sqrt(8.0 *pi*nsigma) * pow(F, 5.0 / 2.0)
                         ) / (pow(4.0*pi, 3.0) * nsigma*nsigma + pow(F, 3.0) + be * sqrt(nsigma) * F*F + 40.0 * pow(nsigma, 3.0/2.0) * F);

    AutoDScalar De_i = -( niS * (1.0 + Ui))
                       /( sqrt(F * nsigma2 / (2.0 * pi)) * (1.0 + he * log(1.0 + sqrt(nsigma2) / F) )
                          +je * Ui * pow(np2, 3.0/4.0) * (1.0 + ke * pow(np2, qe))
                        );
    return  - (De_xc + De_i) * Ryex;
  }

  AutoDScalar EgNarrowToEv   (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    const PetscScalar pi = 3.1415927;

    PetscScalar ndS = ReadDopingNd() / n_scale;
    PetscScalar naS = ReadDopingNa() / n_scale;
    PetscScalar niS = ndS + naS;

    AutoDScalar neS = adtl::fabs(n) / n_scale;
    AutoDScalar nhS = adtl::fabs(p) / n_scale;

    AutoDScalar F = kb * Tl / Ryex;    // [adim]

    AutoDScalar nsigma = neS + nhS;
    AutoDScalar np = alphae * neS + alphah * nhS;

    PetscScalar nsigma2 = ndS + naS;                // For SCR:  nsigma = niS
    PetscScalar np2 = alphae * ndS + alphah * naS;  // For SCR: np = alphae * NDOP + alphah * PDOP

    AutoDScalar Ui = nsigma * nsigma / pow(F,3.0);

    AutoDScalar Dh_xc = -( pow(4.0 * pi, 3.0) * nsigma*nsigma * ( pow(48.0 * nhS / pi / gh, 1.0/3.0) + ch * log(1.0 + dh * pow(np, ph))  )
                           +(8.0 * pi * alphah/gh)* nhS * F*F
                           +sqrt(8.0 *pi*nsigma) * pow(F, 5.0 / 2.0)
                         ) / (pow(4.0*pi, 3.0) * nsigma*nsigma + pow(F,3.0) + bh * sqrt(nsigma) * F*F + 40.0 * pow(nsigma, 3.0/2.0) * F);
    AutoDScalar Dh_i = -( niS * (1.0 + Ui))
                       /( sqrt(F * nsigma2 / (2.0 * pi)) * (1.0 + hh * log(1.0 + sqrt(nsigma2) / F) )
                          +jh * Ui * pow(np2, 3.0/4.0) * (1.0 + kh * pow(np2, qh))
                        );

    return  - (Dh_xc +  Dh_i) * Ryex;
  }


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
    return sqrt(Nc*Nv)*exp(-bandgap/(2*kb*Tl))*exp(EgNarrow(p, n, Tl));
  }

  //end of Bandgap

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


  //[incomplete ionization]
  PetscScalar Na_II(const PetscScalar &p, const PetscScalar &Tl, bool fermi)
  {
    PetscScalar Ni = ReadDopingNa() + ReadDopingNd();
    PetscScalar Nv = this->Nv(Tl);

    PetscScalar gamma = 1.0;
    if( fermi ) gamma = gamma_f(p/Nv);

    PetscScalar Na_eff = 0.0;

    for(unsigned int i=0; i<p_ions.size(); ++i)
    {
      PetscScalar Na = ReadRealVariable(p_ions[i].first);
      PetscScalar Na_crit = p_ions[i].second.N_crit;
      PetscScalar E0      = p_ions[i].second.E0;
      PetscScalar g       = p_ions[i].second.GB;
      PetscScalar alpha   = p_ions[i].second.alpha;
      PetscScalar gamma   = p_ions[i].second.gamma;
      if( Na > 0.0  )
      {
        if(Na < Na_crit)
        {
          PetscScalar dEa = E0 - alpha*pow(Ni, 1.0/3.0);
          PetscScalar p1 = gamma*Nv*exp(-dEa/(kb*Tl));
          Na_eff += Na/( 1 + g*p/p1 );
        }
        else
          Na_eff += Na;
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

    for(unsigned int i=0; i<p_ions.size(); ++i)
    {
      PetscScalar Na = ReadRealVariable(p_ions[i].first);
      PetscScalar Na_crit = p_ions[i].second.N_crit;
      PetscScalar E0      = p_ions[i].second.E0;
      PetscScalar g       = p_ions[i].second.GB;
      PetscScalar alpha   = p_ions[i].second.alpha;
      PetscScalar gamma   = p_ions[i].second.gamma;
      if( Na > 0.0  )
      {
        if(Na < Na_crit)
        {
          PetscScalar dEa = E0 - alpha*pow(Ni, 1.0/3.0);
          AutoDScalar p1 = gamma*Nv*exp(-dEa/(kb*Tl));
          Na_eff += Na/( 1 + g*p/p1 );
        }
        else
          Na_eff += Na;
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

    for(unsigned int i=0; i<n_ions.size(); ++i)
    {
      PetscScalar Nd = ReadRealVariable(n_ions[i].first);
      PetscScalar Nd_crit = n_ions[i].second.N_crit;
      PetscScalar E0      = n_ions[i].second.E0;
      PetscScalar g       = n_ions[i].second.GB;
      PetscScalar alpha   = n_ions[i].second.alpha;
      PetscScalar gamma   = n_ions[i].second.gamma;
      if( Nd > 0.0  )
      {
        if(Nd < Nd_crit)
        {
          PetscScalar dEd = E0 - alpha*pow(Ni, 1.0/3.0);
          PetscScalar n1 = gamma*Nc*exp(-dEd/(kb*Tl));
          Nd_eff += Nd/( 1 + g*n/n1 );
        }
        else
          Nd_eff += Nd;
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

    for(unsigned int i=0; i<n_ions.size(); ++i)
    {
      PetscScalar Nd = ReadRealVariable(n_ions[i].first);
      PetscScalar Nd_crit = n_ions[i].second.N_crit;
      PetscScalar E0      = n_ions[i].second.E0;
      PetscScalar g       = n_ions[i].second.GB;
      PetscScalar alpha   = n_ions[i].second.alpha;
      PetscScalar gamma   = n_ions[i].second.gamma;
      if( Nd > 0.0  )
      {
        if(Nd < Nd_crit)
        {
          PetscScalar dEd = E0 - alpha*pow(Ni, 1.0/3.0);
          AutoDScalar n1 = gamma*Nc*exp(-dEd/(kb*Tl));
          Nd_eff += Nd/( 1 + g*n/n1 );
        }
        else
          Nd_eff += Nd;
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
    Species boron    = {"boron",    -1, 0.045*eV, 4.0, 3.037e-8, 200.0, 0.950, 1e22*std::pow(cm, -3)  };
    Species aluminum = {"aluminum", -1, 0.067*eV, 4.0, 3.037e-8, 200.0, 0.950, 1e22*std::pow(cm, -3)  };
    Species gallium  = {"gallium",  -1, 0.072*eV, 4.0, 3.037e-8, 200.0, 0.950, 1e22*std::pow(cm, -3)  };
    Species indium   = {"indium",   -1, 0.160*eV, 4.0, 3.037e-8, 200.0, 0.950, 1e22*std::pow(cm, -3)  };
    Species pdopant  = {"na",       -1, 0.045*eV, 4.0, 3.037e-8, 200.0, 0.950, 1e22*std::pow(cm, -3)  };

    // n type
    Species nitrogen   = {"nitrogen",   1, 0.045*eV, 2.0, 3.100e-8, 200.0, 1.000, 1e22*std::pow(cm, -3)  };
    Species phosphorus = {"phosphorus", 1, 0.045*eV, 2.0, 3.100e-8, 200.0, 1.000, 1e22*std::pow(cm, -3)  };
    Species arsenic    = {"arsenic",    1, 0.054*eV, 2.0, 3.100e-8, 200.0, 1.000, 1e22*std::pow(cm, -3)  };
    Species antimony   = {"antimony",   1, 0.039*eV, 2.0, 3.100e-8, 200.0, 1.000, 1e22*std::pow(cm, -3)  };
    Species ndopant    = {"nd",         1, 0.045*eV, 2.0, 3.100e-8, 200.0, 1.000, 1e22*std::pow(cm, -3)  };

    species_map["boron"     ] = boron;
    species_map["aluminum"  ] = aluminum;
    species_map["gallium"   ] = gallium;
    species_map["indium"    ] = indium;
    species_map["na"        ] = pdopant;
    species_map["nitrogen"  ] = nitrogen;
    species_map["phosphorus"] = phosphorus;
    species_map["arsenic"   ] = arsenic;
    species_map["antimony"  ] = antimony;
    species_map["nd"        ] = ndopant;

    // fill n_ions and p_ions
    std::map<std::string, Species>::const_iterator it = species_map.begin();
    for( ; it != species_map.end(); ++it)
    {
      const std::string & ion_string = it->second.name;
      if( HasVariable(ion_string) )
      {
        if(it->second.ion > 0)
          n_ions.push_back( std::make_pair( VariableIndex(ion_string),  (it->second)) );
        if(it->second.ion < 0)
          p_ions.push_back( std::make_pair( VariableIndex(ion_string),  (it->second)) );
      }
    }
  }

  std::vector< std::pair<unsigned int, Species> > n_ions;
  std::vector< std::pair<unsigned int, Species> > p_ions;

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

      Species species = { species_name, ion, E0*eV, GB, alpha, beta, gamma, N_crit*std::pow(cm, -3) };
      species_map[species_name] = species;
    }

    // update n_ions and p_ions
    n_ions.clear();
    p_ions.clear();
    for( std::map<std::string, Species>::const_iterator it = species_map.begin(); it != species_map.end(); ++it)
    {
      const std::string & ion_string = it->second.name;
      if( HasVariable(ion_string) )
      {
        if(it->second.ion > 0)
          n_ions.push_back( std::make_pair( VariableIndex(ion_string),  (it->second)) );
        if(it->second.ion < 0)
          p_ions.push_back( std::make_pair( VariableIndex(ion_string),  (it->second)) );
      }
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

    HCI_Fiegna_A = 4.87E+02*m/s/pow(eV, 2.5);
    HCI_Fiegna_X = 1.30E+08*pow(V/(cm*eV*eV), 1.5);

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
    return HCI_Fiegna_A/(3*HCI_Fiegna_X)*pow(Eeff, 1.5)/sqrt(phin)*exp(-HCI_Fiegna_X*pow(phin, 3.0)/pow(Eeff, 1.5));
  }


  PetscScalar HCI_Integral_Fiegna_p(const PetscScalar &phip, const PetscScalar &Eeff)
  {
    if( HCI_Fiegna_X > 30*Eeff  ) return 0;
    return HCI_Fiegna_A/(3*HCI_Fiegna_X)*pow(Eeff, 1.5)/sqrt(phip)*exp(-HCI_Fiegna_X*pow(phip, 3.0)/pow(Eeff, 1.5));
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
  GSS_Si_BandStructure_Schenk(const PMIS_Environment &env):PMIS_BandStructure(env)
  {
    T300 = 300.0*K;
    PMI_Info = "This is the Schenk model for band structure parameters of Silicon";
    Eg_Init();
    Lifetime_Init();
    Recomb_Init();
    RelaxTime_Init();
    Schottky_Init();
    BBTunneling_Init();
  }

  ~GSS_Si_BandStructure_Schenk()
  {}


public:
  int calibrate(std::vector<Parser::Parameter> & pmi_parameters)
  {
    IncompleteIonization_Setup(pmi_parameters);
    return  PMIS_BandStructure::calibrate(pmi_parameters);
  }

  void post_calibrate_process()
  {
    Eg_PostCalibrateProcess();
  }

};


extern "C"
{
  DLL_EXPORT_DECLARE  PMIS_BandStructure*  PMIS_Si_BandStructure_Schenk (const PMIS_Environment& env)
  {
    return new GSS_Si_BandStructure_Schenk(env);
  }
}
