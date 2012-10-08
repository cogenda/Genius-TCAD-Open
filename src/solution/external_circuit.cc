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



#include "parser_card.h"
#include "physical_unit.h"

#include "external_circuit.h"
#include "external_circuit_rcl.h"
#include "external_circuit_pi.h"
#include "external_circuit_rct.h"

using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::C;
using PhysicalUnit::s;
using PhysicalUnit::um;

ExternalCircuit * ExternalCircuit::build_default()
{ return new ExternalCircuitRCL(); }


ExternalCircuit * ExternalCircuit::build(const Parser::Card &c)
{
  std::string type = c.get_string("ext.circuit", "rcl");

  if( type == "rcl" )
  {
    Real res = std::max(1e-6, c.get_real ( "res", 0.0 )) *V/A;
    Real cap = c.get_real ( "cap", 0.0 ) *C/V;
    Real ind = c.get_real ( "ind", 0.0 ) *V*s/A;
    return new ExternalCircuitRCL(res, cap, ind);
  }

  if( type == "pi" )
  {
    Real r    = std::max(1e-3, c.get_real ( "res.load", 0.0 )) *V/A;
    Real R    = std::max(1e-3, c.get_real ( "res", 0.0 )) *V/A;
    Real cap1 = c.get_real ( "cap1", 0.0 ) *C/V;
    Real cap2 = c.get_real ( "cap2", 0.0 ) *C/V;
    return new ExternalCircuitPI(r, R, cap1, cap2);
  }

  if( type == "tline" )
  {
    Real r      = std::max(1e-3, c.get_real ( "res.load", 0.0 )) *V/A;
    Real Rl     = c.get_real ( "res.um", 0.0 ) *V/A/um;
    Real capl   = c.get_real ( "cap.um", 0.0 ) *C/V/um;
    Real length = c.get_real ( "tl.length", 0.0 ) *um;
    int  div    = c.get_int  ( "tl.div", 100);
    return new ExternalCircuitRCTLine(r, Rl, capl, length, div);
  }

  return 0;
}


