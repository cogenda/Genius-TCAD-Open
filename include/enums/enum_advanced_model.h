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

//  $Id: enum_advanced_model.h,v 1.3 2008/07/09 05:58:15 gdiso Exp $



#ifndef __enum_adv_model_h__
#define __enum_adv_model_h__



// ------------------------------------------------------------
// enum Model Specify

namespace ModelSpecify
{

  /**
   * define carriers should be considered
   */
  enum Carrier
  {
    N_Type=0,
    P_Type,
    PN_Both,
    PN_None
  };

  /**
   * for EB solver, specify different levels
   */
  enum EBMLevel
  {
    NONE=0,
    Tl,
    Tn,
    Tp,
    TnTp,
    TnTl,
    TpTl,
    ALL
  };

  /**
   * define how to evaluate high field mobility force
   */
  enum MobilityForce
  {
    EJ   = 0,
    ESimple,
    EQF,
    CarrierTemp
  };


  /**
   * define how to evaluate Impact Ionization force
   */
  enum IIForce
  {
    IIForce_EdotJ   = 0,
    GradQf,
    EVector,
    ESide,
    SoftII,
    HardII,
    TempII
  };


}




#endif // #define __enum_adv_model_h__
