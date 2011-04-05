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



#ifndef __enum_bc_h__
#define __enum_bc_h__


/**
 * the Boundary Condition Type genius support
 */
enum BCType
{

  /**
   * Neumann Boundary is the most general boundary,
   * which has ZERO flux through this boundary
   */
  NeumannBoundary             = 0x0001,

  /**
   * The interface of Semiconductor region to Vacuum region.
   * In most of the case, it is considered as Neumann Boundary
   */
  IF_Semiconductor_Vacuum     = 0x0002,

  /**
   * The interface of Insulator region to Vacuum region.
   * In most of the case, it is considered as Neumann Boundary
   */
  IF_Insulator_Vacuum         = 0x0003,

  /**
   * The interface of Electrode region to Vacuum region.
   * In most of the case, it is considered as Neumann Boundary
   */
  IF_Electrode_Vacuum         = 0x0004,


  /**
   * The interface of resistive metal region to Vacuum region.
   * In most of the case, it is considered as Neumann Boundary
   */
  IF_Metal_Vacuum        = 0x0005,

  /**
   * The interface of PML region to Vacuum/Insulator region.
   */
  IF_PML_Scatter              = 0x0006,

  /**
   * The interface of PML region to another PML region.
   */
  IF_PML_PML                  = 0x0007,

  /**
   * The interface of Electrode region to Insulator region.
   * we assume potential and temperature continuous on this boundary
   */
  IF_Electrode_Insulator      = 0x0010,


  /**
   * The interface of Semiconductor region to Insulator region.
   * this is important in MOS simulation.
   */
  IF_Insulator_Semiconductor  = 0x0012,

  /**
   * The interface of Insulator region to Insulator region.
   * Since some MOS has a nitride-oxide gate. we have to deal with this
   * boundary condition
   */
  IF_Insulator_Insulator      = 0x0013,

  /**
   * The interface of Electrode region to Electrode region.
   * nearly meaningless.
   */
  IF_Electrode_Electrode      = 0x0014,


  /**
   * The interface of Electrode region to resistive metal region.
   * nearly meaningless.
   */
  IF_Electrode_Metal          = 0x0015,


  /**
   * The interface of resistive metal region to Insulator region.
   * we assume div(J)=0 on the boundary.
   * NOTE, it can be a resistat gate boundary
   */
  IF_Insulator_Metal          = 0x0016,

  /**
   * The interface of resistive metal region to resistive metal region.
   * i.e. n-poly and p-poly contact
   */
  IF_Metal_Metal              = 0x0017,

  /**
   * The interface of Semiconductor region to semiconductor region with same material.
   * nearly meaningless.
   */
  HomoInterface               = 0x0018,

  /**
   * The interface of Semiconductor region to semiconductor region with different material.
   * which forms heterojunction.
   */
  HeteroInterface             = 0x0019,

  /**
   * The interface of Semiconductor region to Electrode region.
   * by default, it forms OhmicContact.
   */
  IF_Electrode_Semiconductor  = 0x001a,

  /**
   * The interface of Semiconductor region to resistive metal region.
   * by default, it forms OhmicContact.
   */
  IF_Metal_Semiconductor  = 0x001b,

  /**
   * Electrode which froms  Ohmic Contact.
   */
  OhmicContact        = 0x0101,

  /**
   * Electrode which froms resistive metal Ohmic Contact.
   */
  IF_Metal_Ohmic = 0x0102,

  /**
   * Electrode which froms  Schottky Contact.
   */
  SchottkyContact     = 0x0103,

  /**
   * Electrode which froms resistive metal Schottky Contact.
   */
  IF_Metal_Schottky  = 0x0104,

  /**
   * Simple Boundary condition for MOS gate
   */
  SimpleGateContact   = 0x0105,

  /**
   * Boundary condition for MOS gate
   */
  GateContact         = 0x0106,

  /**
   * Boundary condition for float metal with charge,
   * useful for EEPROM simulation
   */
  ChargedContact      = 0x0107,

  /**
   * Boundary condition for float metal with charge,
   * useful for EEPROM simulation
   */
  ChargeIntegral      = 0x0108,


  /**
   * boundary condition for electrode pad. pad sould be on the surface of resistance region
   */
  SolderPad           = 0x0109,

  /**
   * Boundary condition for inter connect of electrodes
   * useful for IC cell i.e. Inverter or SRAM simulation
   */
  InterConnect        = 0x0110,


  /**
   * Absorbing boundary for electromagnetic simulation
   */
  AbsorbingBoundary   = 0x1001,

  /**
   * boundary for adding light wave.
   */
  SourceBoundary      = 0x1002,


  INVALID_BC_TYPE     = 0xffff           // should always be last

};


#endif // #define __enum_bc_h__




