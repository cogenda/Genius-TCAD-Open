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



#ifndef __external_circuit_pi_h__
#define __external_circuit_pi_h__

#include "external_circuit.h"


  // For voltage driven
  //
  //         _____   V1   _____     V0
  //    ----|_____|------|_____|--------> to electrode (Ve, I)
  //    |      r     |       R      |
  //   Vapp         === C1         === C2
  //    |____________|______________|
  //                GND
  //
  //  (V0-V1)/R + C2*dV0/dt - I         = 0
  //  (V1-Va)/r + C1*dV1/dt + (V1-V0)/R = 0
  


  // And for current driven
  //                       _____   Ve
  //    ------------------|_____|--------> to electrode (Ve, I)
  //    |            |       R      |
  //   Iapp         === C1         === C2
  //    |____________|______________|
  //                GND
  //
  //  (V0-V1)/R + C2*dV0/dt - I  = 0
  //  -Ia + C1*dV1/dt + (V1-V0)/R = 0



class ExternalCircuitPI : public ExternalCircuit
{
public:
  ExternalCircuitPI(Real r, Real R, Real C1=0.0, Real C2=0.0)
  :_r_app(r), _res(R),_cap1(C1),_cap2(C2),_V0(_potential),_V0_last(_potential_old)
  {
    _V1      = 0.0;
    _V1_last = 0.0;
  }

  virtual ~ExternalCircuitPI() {}


  /**
   * type of external circuit
   */
  virtual std::string type() const { return "PI"; }

  /**
   * format string for export
   */
  virtual std::string format() const;

  /**
   * @return the lumped inductance
   */
  Real r() const    {return _r_app;}

  /**
   * @return writable reference to lumped inductance
   */
  Real & r()        {return _r_app;}

  /**
   * @return the lumped resistor
   */
  Real R() const    {return _res;}

  /**
   * @return writable reference to lumped resistor
   */
  Real & R()        {return _res;}

  /**
   * @return the lumped capacitance
   */
  Real C1() const    {return _cap1;}

  /**
   * @return writable reference to lumped capacitance
   */
  Real & C1()        {return _cap1;}

  /**
   * @return the lumped capacitance
   */
  Real C2() const    {return _cap2;}

  /**
   * @return writable reference to lumped capacitance
   */
  Real & C2()        {return _cap2;}

  /**
   * value of serial resistance for inter connect
   */
  virtual Real inter_connect_resistance() const;

  /**
   * value of serial resistance to electrode
   */
  virtual Real serial_resistance() const { return _r_app+_res; }


  /**
   * value of serial resistance to electrode
   */
  virtual void set_serial_resistance(Real res) { _r_app = res-_res; }


  /**
   * use this value as scaling to electrode
   */
  virtual Real electrode_scaling(Real dt) const { return std::min(_r_app,_res); }


  /**
   * use this value as scaling to Modified Nodal Analysis
   */
  virtual Real mna_scaling(Real dt) const
  {
    return 1.0;
  }

  /**
   * calculate function of Modified Nodal Analysis
   */
  virtual Real mna_function(Real dt);

  /**
   * calculate jacobian of Modified Nodal Analysis
   */
  virtual Real mna_jacobian(Real dt);


  /**
   * use this value as scaling to Modified Nodal Analysis in AC model
   */
  virtual std::complex <Real> mna_ac_scaling(Real omega) { return std::complex <Real>(1.0, 0.0); }


  /**
   * calculate AC jacobian of Modified Nodal Analysis
   */
  virtual std::complex <Real> mna_ac_jacobian(Real omega);

  /**
   * when a solution is achieved, update the potential/current
   */
  virtual void update()
  {
    ExternalCircuit::update();
    _V1_last = _V1;
  }

  /**
   * roll back to previous solution
   */
  virtual void rollback()
  {
    ExternalCircuit::rollback();
    _V1 = _V1_last;
  }

  /**
   * init op state before transient simulation
   */
  virtual void tran_op_init();

private:

  Real _r_app;

  Real _res;

  Real _cap1;

  Real _cap2;

  // template variable
  Real &_V0, &_V0_last;

  Real _V1, _V1_last;


};

#endif
