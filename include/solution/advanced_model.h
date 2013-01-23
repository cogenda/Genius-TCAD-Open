/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/

//  $Id: advanced_model.h,v 1.7 2008/07/09 05:58:16 gdiso Exp $

#ifndef __advanced_model_h__
#define __advanced_model_h__


#include <string>
#include <vector>

#include "enum_advanced_model.h"


/**
 * region advanced model specification
 * mainly for semiconductor region
 */
struct AdvancedModel
{

  /**
   * constructor, set default values
   */
  AdvancedModel()
  {
    Carrier = ModelSpecify::PN_Both;

    ESurface          = true;
    TransverseMobilityInBoundaryLayer = true;
    HighFieldMobility = true;
    HighFieldMobilitySelfConsistently = true;
    Mob_Force = ModelSpecify::ESimple;
    QuasiFermiCarrierTruc = 1e-2;

    HotCarrierInjection = false;
    DIRTunneling = false;
    FNTunneling = false;
    TunnelingSelfConsistently = false;
    BandBandTunneling = false;

    ImpactIonization = false;
    II_Force = ModelSpecify::GradQf;

    Fermi=false;
    IncompleteIonization=false;

    Trap=false;

    QNEnabled = false;
    QPEnabled = false;
    QNFactor = 1.0;
    QPFactor = 1.0;
    QMinConcentration = 1.0;

    EB_Level = ModelSpecify::NONE;
  }



  /**
   * specify which carrier should be considered during the simulation
   */
  ModelSpecify::Carrier Carrier;


  //-----------------------------------------------------------
  // parameters for high field mobility
  //-----------------------------------------------------------

  /**
   * specify if effective surface elecical field should be used
   */
  bool    ESurface;

  /**
   * specify if transverse mobility only be considered near the insulator interface
   */
  bool    TransverseMobilityInBoundaryLayer;

  /**
   * specify if high field mobility should be used
   */
  bool    HighFieldMobility;

  /**
   * specify if we should consider the High Field Mobility self-consistently
   * if it is false, the previous E field will be used to evaluate the mobility
   */
  bool    HighFieldMobilitySelfConsistently;

  /**
   * how to evaluate the driving force in the High Field Mobility
   */
  ModelSpecify::MobilityForce Mob_Force;


  /**
   * density truncation for Quasi-Fermi level evaluation
   */
  double  QuasiFermiCarrierTruc;

  //-----------------------------------------------------------
  // parameters for hot carrier injection
  //-----------------------------------------------------------

  /**
   * specify if hot carrier injection should be supported
   */
  bool HotCarrierInjection;


  //-----------------------------------------------------------
  // parameters for tunneling
  //-----------------------------------------------------------

  /**
   * specify if direct tunneling should be supported
   */
  bool    DIRTunneling;

  /**
   * specify if Fowler-Nordheim tunneling should be supported
   */;
  bool    FNTunneling;


  /**
   * tunneling(including dir and fn) self consistently
   */
  bool    TunnelingSelfConsistently;


  //-----------------------------------------------------------
  // parameters for band band tunneling
  //-----------------------------------------------------------

  /**
   * specify if band band tunneling should be supported
   */
  bool    BandBandTunneling;


  //-----------------------------------------------------------
  // parameters for impact ionization
  //-----------------------------------------------------------

  /**
   * specify if impact ionization should be supported
   */
  bool    ImpactIonization;

  /**
   * how to evaluate the driving force in the impact ionization
   */
  ModelSpecify::IIForce II_Force;


  //-----------------------------------------------------------
  // parameters for Fermi statistics and Incomplete Ionization
  //-----------------------------------------------------------

  /**
   * specify if Fermi statistics should be supported
   */
  bool    Fermi;

  /**
   * specify if Incomplete Ionization should be supported
   */
  bool    IncompleteIonization;


  //------------------------------------------------------
  // parameters for Charge Trapping
  //------------------------------------------------------
  /**
   * specify if charge trapping should be supported
   */
  bool    Trap;

  //------------------------------------------------------
  // parameters for DG-DDM simulation
  //------------------------------------------------------

  /**
   * flag for electron quantum potential
   */
  bool    QNEnabled;

  /**
   * flag for hole quantum potential
   */
  bool    QPEnabled;

  /**
   * damping factor for electron quantum potential
   */
  double  QNFactor;

  /**
   * damping factor for hole quantum potential
   */
  double  QPFactor;

  /**
   * min concentration for DG solver
   * this is a factor to nie, default is 1.0
   */
  double  QMinConcentration;


  //------------------------------------------------------
  // parameters for EBM simulation
  //------------------------------------------------------

  /**
   *  enum value for specify levels for energy balance solver
   */
  ModelSpecify::EBMLevel       EB_Level;

  /**
   * @return true if EBM solver considers electron temperature
   */
  bool enable_Tn() const
    { return (EB_Level==ModelSpecify::Tn || EB_Level==ModelSpecify::TnTp ||EB_Level==ModelSpecify::TnTl || EB_Level==ModelSpecify::ALL); }

  /**
   * @return true if EBM solver considers hole temperature
   */
  bool enable_Tp() const
    { return (EB_Level==ModelSpecify::Tp || EB_Level==ModelSpecify::TnTp || EB_Level==ModelSpecify::TpTl || EB_Level==ModelSpecify::ALL); }

  /**
   * @return true if EBM solver considers lattice temperature
   */
  bool enable_Tl() const
    { return (EB_Level==ModelSpecify::Tl || EB_Level==ModelSpecify::TnTl || EB_Level==ModelSpecify::TpTl || EB_Level==ModelSpecify::ALL); }


  /**
   * force to consider temperature
   */
  void force_temperature_usage()
  {
    if (EB_Level == ModelSpecify::NONE) EB_Level=ModelSpecify::Tl;
    if (EB_Level == ModelSpecify::Tn)   EB_Level=ModelSpecify::TnTl;
    if (EB_Level == ModelSpecify::Tp)   EB_Level=ModelSpecify::TpTl;
    if (EB_Level == ModelSpecify::TnTp) EB_Level=ModelSpecify::ALL;
  }

  /**
   *  @return the level of SG discrete scheme of electron current equation
   */
  unsigned int Jn_level() const
  {
    switch(EB_Level)
    {
    case ModelSpecify::NONE :
    case ModelSpecify::Tp   : return 1;
    case ModelSpecify::Tl   :
    case ModelSpecify::TpTl : return 2;
    case ModelSpecify::Tn   :
    case ModelSpecify::TnTl :
    case ModelSpecify::TnTp :
    case ModelSpecify::ALL  : return 3;
    }
    return 0; //prevent compiler warning
  }

  /**
   *  @return the level of SG discrete scheme of hole current equation
   */
  unsigned int Jp_level() const
  {
    switch(EB_Level)
    {
    case ModelSpecify::NONE :
    case ModelSpecify::Tn   : return 1;
    case ModelSpecify::Tl   :
    case ModelSpecify::TnTl : return 2;
    case ModelSpecify::Tp   :
    case ModelSpecify::TpTl :
    case ModelSpecify::TnTp :
    case ModelSpecify::ALL  : return 3;
    }
    return 0; //prevent compiler warning
  }

  /**
   *  @return the level of joule heating of electron
   */
  unsigned int Hn_level() const
  {
    switch(EB_Level)
    {
    case ModelSpecify::NONE : return 0;
    case ModelSpecify::Tp   : return 0;
    case ModelSpecify::Tl   : return 1;
    case ModelSpecify::TpTl : return 1;
    case ModelSpecify::Tn   : return 2;
    case ModelSpecify::TnTl : return 2;
    case ModelSpecify::TnTp : return 2;
    case ModelSpecify::ALL  : return 2;
    }
    return 0; //prevent compiler warning
  }

  /**
   *  @return the level of joule heating of hole
   */
  unsigned int Hp_level() const
  {
    switch(EB_Level)
    {
    case ModelSpecify::NONE : return 0;
    case ModelSpecify::Tp   : return 2;
    case ModelSpecify::Tl   : return 1;
    case ModelSpecify::TpTl : return 2;
    case ModelSpecify::Tn   : return 0;
    case ModelSpecify::TnTl : return 1;
    case ModelSpecify::TnTp : return 2;
    case ModelSpecify::ALL  : return 2;
    }
    return 0; //prevent compiler warning
  }

};



#endif //#define __advanced_model_h__
