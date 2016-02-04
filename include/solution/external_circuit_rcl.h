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



#ifndef __external_circuit_rcl_h__
#define __external_circuit_rcl_h__

#include "external_circuit.h"


  // For voltage driven
  //
  //          _____                Ve
  //    -----|_____|----/\/\/\/\-------> to electrode (Ve, I)
  //    |       R          L       |
  //   Vapp                     C ===
  //    |__________________________|
  //           GND
  //
  // And for current driven
  //                               Ve
  //    --------------------------------> to electrode (Ve, I)
  //    |                          |
  //   Iapp                     C ===
  //    |__________________________|
  //           GND


class ExternalCircuitRCL : public ExternalCircuit
{
public:
  ExternalCircuitRCL(Real R=0.0, Real C=0.0, Real L=0.0)
  :_res(R), _cap(C), _ind(L),_cap_current(0.0),_cap_current_old(0.0)
  {}

  virtual ~ExternalCircuitRCL() {}


  /**
   * type of external circuit
   */
  virtual std::string type() const { return "RCL"; }

  /**
   * format string for export
   */
  virtual std::string format() const;


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
  Real C() const    {return _cap;}

  /**
   * @return writable reference to lumped capacitance
   */
  Real & C()        {return _cap;}

  /**
   * @return the lumped inductance
   */
  Real L() const    {return _ind;}

  /**
   * @return writable reference to lumped inductance
   */
  Real & L()        {return _ind;}


  /**
   * value of serial resistance to electrode
   */
  virtual Real serial_resistance() const { return _res; }


  /**
   * value of serial resistance for inter connect
   */
  virtual Real inter_connect_resistance() const;


  /**
   * value of serial resistance to electrode
   */
  virtual void set_serial_resistance(Real res) { _res = res; }


  /**
   * use this value as scaling to electrode
   */
  virtual Real electrode_scaling(Real dt) const { return 1.0/(1.0+_res); }


  /**
   * use this value as scaling to Modified Nodal Analysis
   */
  virtual Real mna_scaling(Real dt) const;

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
  virtual std::complex <Real> mna_ac_scaling(Real omega);


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
    _cap_current_old = _cap_current;
  }

  /**
   * roll back to previous solution
   */
  virtual void rollback()
  {
    ExternalCircuit::rollback();
    _cap_current = _cap_current_old;
  }

  /**
   * init op state before transient simulation
   */
  virtual void tran_op_init()
  {
    ExternalCircuit::tran_op_init();
    _cap_current = _cap_current_old = 0.0;
  }

private:

  Real _res;

  Real _cap;

  Real _ind;

  /**
   * the current which flow through lumped capacitance to ground
   */
  Real      _cap_current;

  /**
   * the current which flow through lumped capacitance to ground for last step
   */
  Real      _cap_current_old;

};

#endif
