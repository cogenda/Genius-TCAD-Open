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




#include "external_circuit_pi.h"
#include "physical_unit.h"

std::string ExternalCircuitPI::format() const
{
  std::stringstream ss;
  ss  <<"string<ext.circuit>=pi"<<" "
      <<"real<res.load>="<<_r_app/(PhysicalUnit::V/PhysicalUnit::A)<<" "
      <<"real<res>="<<_res/(PhysicalUnit::V/PhysicalUnit::A)<<" "
      <<"real<cap1>="<<_cap1/(PhysicalUnit::C/PhysicalUnit::V)<<" "
      <<"real<cap2>="<<_cap2/(PhysicalUnit::C/PhysicalUnit::V)<<" "
      <<"real<potential>="<<this->potential()/PhysicalUnit::V<<" "
      <<"real<potential.old>="<<this->potential_old()/PhysicalUnit::V;
  return ss.str();
}


Real ExternalCircuitPI::inter_connect_resistance() const
{ return _r_app+_res; }



Real ExternalCircuitPI::mna_function(Real dt)
{
  if( this->is_voltage_driven() )
  {
    //  (V0-V1)/R + C2*dV0/dt - I         = 0
    //  (V1-Va)/r + C1*dV1/dt + (V1-V0)/R = 0
    //| A B | V0 = a + I
    //| C D | V1 = b
    Real A = 1/_res + _cap2/dt;
    Real B = -1/_res;
    Real C = -1/_res;
    Real D = 1/_res + 1/_r_app + _cap1/dt;
    Real a =  _cap2*_V0_last/dt;
    Real b = _Vapp/_r_app + _cap1*_V1_last/dt;

    _V1 = (b-C*_V0)/D;
    Real I = (A-B/D*C)*_V0 + B/D*b - a;
    return I;
  }

  if( this->is_current_driven() )
  {
    //  (V0-V1)/R + C2*dV0/dt - I  = 0
    //  -Ia + C1*dV1/dt + (V1-V0)/R = 0
    //| A B | V0 = a + I
    //| C D | V1 = b
    Real A = 1/_res + _cap2/dt;
    Real B = -1/_res;
    Real C = -1/_res;
    Real D = 1/_res + _cap1/dt;
    Real a =  _cap2*_V0_last/dt;
    Real b = _Iapp + _cap1*_V1_last/dt;

    _V1 = (b-C*_V0)/D;
    Real I = (A-B/D*C)*_V0 + B/D*b - a;
    return I;
  }

  return 0.0;
}


Real ExternalCircuitPI::mna_jacobian(Real dt)
{
  if( this->is_voltage_driven() )
  {
    //  (V0-V1)/R + C2*dV0/dt - I         = 0
    //  (V1-Va)/r + C1*dV1/dt + (V1-V0)/R = 0
    //| A B | V0 = a + I
    //| C D | V1 = b
    Real A = 1/_res + _cap2/dt;
    Real B = -1/_res;
    Real C = -1/_res;
    Real D = 1/_res + 1/_r_app + _cap1/dt;
    //Real a =  _cap2*_V0_last/dt;
    //Real b = _Vapp/_r_app + _cap1*_V1_last/dt;

    //_V1 = (b-C*_V0)/D;
    //Real I = (A-B/D*C)*_V0 + B/D*b - a;
    return (A-B/D*C); //dI/dV0
  }

  if( this->is_current_driven() )
  {
    //  (V0-V1)/R + C2*dV0/dt - I  = 0
    //  -Ia + C1*dV1/dt + (V1-V0)/R = 0
    //| A B | V0 = a + I
    //| C D | V1 = b
    Real A = 1/_res + _cap2/dt;
    Real B = -1/_res;
    Real C = -1/_res;
    Real D = 1/_res + _cap1/dt;
    //Real a =  _cap2*_V0_last/dt;
    //Real b = _Iapp + _cap1*_V1_last/dt;

    //_V1 = (b-C*_V0)/D;
    //Real I = (A-B/D*C)*_V0 + B/D*b - a;
    return (A-B/D*C); //dI/dV0
  }

  return 0.0;
}


std::complex <Real> ExternalCircuitPI::mna_ac_jacobian(Real omega)
{
  //  (V0-V1)/R + C2*j*omega*V0 - I     = 0
  //  (V1-Va)/r + C1*j*omega*V1 + (V1-V0)/R = 0
  //| A B | V0 = a + I
  //| C D | V1 = b
  std::complex <Real> j(0.0, 1.0);
  std::complex <Real> A = 1/_res + j*omega*_cap2;
  std::complex <Real> B = -1/_res;
  std::complex <Real> C = -1/_res;
  std::complex <Real> D = 1/_res + 1/_r_app + j*omega*_cap1;
  //std::complex <Real> a =  0.0;
  //std::complex <Real> b = _Vac/_r_app;

  //Real I = (A-B/D*C)*_V0 + B/D*b - a;
  return (A-B/D*C); //dI/dV0
}


void ExternalCircuitPI::tran_op_init()
{
  ExternalCircuit::tran_op_init();

  if( this->is_voltage_driven() )
  {
    //  (V0-V1)/R + C2*dV0/dt - I         = 0
    //  (V1-Va)/r + C1*dV1/dt + (V1-V0)/R = 0
    //| A B | V0 = a + I
    //| C D | V1 = b
    Real A = 1/_res;
    Real B = -1/_res;
    Real C = -1/_res;
    Real D = 1/_res + 1/_r_app;
    Real a =  0;
    Real b = _Vapp/_r_app;

    _V1 = _V1_last = (b-C*_V0)/D;

  }

  if( this->is_current_driven() )
  {
    //  (V0-V1)/R + C2*dV0/dt - I  = 0
    //  -Ia + C1*dV1/dt + (V1-V0)/R = 0
    //| A B | V0 = a + I
    //| C D | V1 = b
    Real A = 1/_res;
    Real B = -1/_res;
    Real C = -1/_res;
    Real D = 1/_res;
    Real a =  0;
    Real b = _Iapp;

    _V1 = _V1_last = (b-C*_V0)/D;
  }
}



