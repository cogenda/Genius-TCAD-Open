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




#include "external_circuit_rcl.h"
#include "physical_unit.h"

std::string ExternalCircuitRCL::format() const
{
  std::stringstream ss;
  ss  <<"string<ext.circuit>=rcl"<<" "
      <<"real<res>="<<_res/(PhysicalUnit::V/PhysicalUnit::A)<<" "
      <<"real<cap>="<<_cap/(PhysicalUnit::C/PhysicalUnit::V)<<" "
      <<"real<ind>="<<_ind/(PhysicalUnit::V*PhysicalUnit::s/PhysicalUnit::A)<<" "
      <<"real<potential>="<<this->potential()/PhysicalUnit::V<<" "
      <<"real<potential.old>="<<this->potential_old()/PhysicalUnit::V;
  return ss.str();
}


Real ExternalCircuitRCL::inter_connect_resistance() const
{ return std::max(1e-3*PhysicalUnit::Ohm, _res); }


Real ExternalCircuitRCL::mna_function(Real dt)
{
  if( this->is_voltage_driven() )
    return  (_potential-_Vapp) + (_ind/dt+_res)*_cap/dt*_potential - (_ind/dt+_res)*_cap/dt*_potential_old - _ind/dt*(_current_old+_cap_current_old);
  else if( this->is_current_driven() )
    return _cap_current_old - _Iapp;
  else
    return 0.0;
}


Real ExternalCircuitRCL:: mna_jacobian(Real dt)
{
  if( this->is_voltage_driven() )
    return 1+(_ind/dt+_res)*_cap/dt;
  else if( this->is_current_driven() )
    return 0.0;
  else
    return 0.0;
}

