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



#ifndef __external_circuit_h__
#define __external_circuit_h__

#include <string>
#include <complex>

#include "genius_common.h"


namespace Parser{ class Card; }


/**
 * this is a class for external circuit for electrode
 */
class ExternalCircuit
{
public:

  enum DRIVEN{VDRIVEN, IDRIVEN, FLOAT};

  /**
   * constructor
   */
  ExternalCircuit()
  : _Vapp(0), _Iapp(0), _drv(VDRIVEN),
    _potential(0), _current(0),
    _current_displacement(0), _current_conductance(0),
    _current_electron(0), _current_hole(0),
    _Vac(0.0026)
  {}

  /**
   * destructor
   */
  virtual ~ExternalCircuit() {}

  /**
   * type of external circuit
   */
  virtual std::string type() const = 0;

  /**
   * format string for export
   */
  virtual std::string format() const = 0;

  /**
   * create default ExternalCircuit (RCL model)
   */
  static ExternalCircuit * build_default();

  /**
   * create ExternalCircuit by input
   */
  static ExternalCircuit * build(const Parser::Card &c);

  // Vapp or Iapp settings 
public:
  /**
   * @return the application voltage.
   */
  Real Vapp() const
  { return _Vapp;}

  /**
   * @return writable reference to application voltage.
   */
  Real & Vapp()
  { return _Vapp;}

  /**
   * @return the application current.
   */
  Real Iapp() const
  { return _Iapp;}

  /**
   * @return writable reference to application current.
   */
  Real & Iapp()
  { return _Iapp;}

  /**
   * force to ground this electrode
   */
  void ground()
  {
    _Vapp = 0.0;
    _Iapp = 0.0;
  }

  /**
   * indicate this electrode is stimulated by voltage source
   */
  void set_voltage_driven()
  { _drv = VDRIVEN; }

  /**
   * indicate this electrode is stimulated by current source
   */
  void set_current_driven()
  { _drv = IDRIVEN; }

  /**
   * indicate this electrode is float
   */
  void set_float()
  { _drv = FLOAT; }

  /**
   * @return true iff voltage driven
   */
  bool is_voltage_driven() const
  { return _drv == VDRIVEN; }

  /**
   * @return true iff current driven
   */
  bool is_current_driven() const
  { return _drv == IDRIVEN; }

  /**
   * @return true iff electrode is float
   */
  bool is_float() const
  { return _drv == FLOAT; }

  /**
   * @return the driven state of electrode
   */
  DRIVEN driven_state() const
  { return _drv; }

  /**
   * value of serial resistance for inter connect
   */
  virtual Real inter_connect_resistance() const=0;

  /**
   * value of serial resistance to electrode
   */
  virtual Real serial_resistance() const = 0;

  /**
   * value of serial resistance to electrode
   */
  virtual void set_serial_resistance(Real res) = 0;

protected:

  /**
   * application voltage
   */
  Real   _Vapp;

  /**
   * application current
   */
  Real   _Iapp;

  /**
   * indicator this electrode is driven by voltage source or
   * current source
   */
  DRIVEN        _drv;


  // ckt load
public:

  /**
   * @return the electrode potential.
   */
  Real potential() const
  { return _potential;}

  /**
   * @return writable reference to electrode potential.
   */
  Real & potential()
  { return _potential;}

  /**
   * @return electrode potential of last step.
   */
  Real  potential_old() const
  { return _potential_old;}

  /**
   * @return writable reference to electrode potential of last step.
   */
  Real & potential_old()
  { return _potential_old;}

  /**
   * @return the electrode current
   */
  Real current() const
  { return _current;}

  /**
   * @return the writable reference to electrode current
   */
  Real & current()
  { return _current;}

  /**
   * @return electrode current of last step.
   * @note this variable is read only for other application
   */
  Real  current_old() const
  { return _current_old;}

  /**
   * use this value as scaling to electrode
   */
  virtual Real electrode_scaling(Real dt) const { return 1.0; }

  /**
   * use this value as scaling to Modified Nodal Analysis
   */
  virtual Real mna_scaling(Real dt) const { return 1.0; }

  /**
   * calculate function of Modified Nodal Analysis
   */
  virtual Real mna_function(Real dt) = 0;

  /**
   * calculate jacobian of Modified Nodal Analysis
   */
  virtual Real mna_jacobian(Real dt) = 0;

  /**
   * use this value as scaling to Modified Nodal Analysis in AC model
   */
  virtual std::complex <Real> mna_ac_scaling(Real omega) = 0;

  /**
   * calculate AC jacobian of Modified Nodal Analysis
   */
  virtual std::complex <Real> mna_ac_jacobian(Real omega) = 0;

  /**
   * when a solution is achieved, update the potential/current
   */
  virtual void update()
  {
    _potential_old = _potential;
    _current_old = _current;
  }

  /**
   * roll back to previous solution
   */
  virtual void rollback()
  {
    _potential = _potential_old;
    _current   = _current_old;
  }

  /**
   * init op state before transient simulation
   */
  virtual void tran_op_init()
  {
    _potential_old = _potential;
    _current_old = _current;
  }


protected:
  /**
   * the potential of this electrode
   */
  Real      _potential;

  /**
   * the potential for last step
   */
  Real      _potential_old;

  /**
   * the current flow out of this electrode
   */
  Real      _current;

  /**
   * the current for last step
   */
  Real      _current_old;


  // current statistic
public:

  /**
   * @return the displacement current.
   */
  Real current_displacement() const
  { return _current_displacement;}

  /**
   * @return writable reference to displacement current.
   */
  Real & current_displacement()
  { return _current_displacement;}


  /**
   * @return the conductance current.
   */
  Real current_conductance() const
  { return _current_conductance;}

  /**
   * @return writable reference to conductance current.
   */
  Real & current_conductance()
  { return _current_conductance;}

  /**
   * @return the electron current.
   */
  Real current_electron() const
  { return _current_electron;}

  /**
   * @return writable reference to electron current.
   */
  Real & current_electron()
  { return _current_electron;}

  /**
   * @return the hole current.
   */
  Real current_hole() const
  { return _current_hole;}

  /**
   * @return writable reference to hole current.
   */
  Real & current_hole()
  { return _current_hole;}

protected:

  /**
   * the displacement current flow out of this electrode
   */
  Real      _current_displacement;

  /**
   * the conductance current flow out of this electrode
   */
  Real      _current_conductance;

  /**
   * the electron current flow out of this electrode
   */
  Real      _current_electron;

  /**
   * the hole current flow out of this electrode
   */
  Real      _current_hole;

  // ac settings 
public:
  /**
   * @return the stimulate voltage for AC scan.
   */
  Real Vac() const
  { return _Vac; }

  /**
   * @return the writable reference to stimulate voltage for AC scan.
   */
  Real & Vac()
  { return _Vac; }

  /**
   * @return the electrode potential of AC scan.
   */
  std::complex<Real> potential_ac() const
  { return _potential_ac;}

  /**
   * @return writable reference to electrode potential of AC scan.
   */
  std::complex<Real> & potential_ac()
  { return _potential_ac;}


  /**
   * @return the electrode current of AC scan.
   */
  std::complex<Real> current_ac() const
  { return _current_ac;}

  /**
   * @return writable reference to electrode current of AC scan.
   */
  std::complex<Real> & current_ac()
  { return _current_ac;}


protected:
  /**
   * the application voltage for AC sweep.
   */
  Real      _Vac;

  /**
   * the (complex) electrode potential for AC sweep.
   */
  std::complex<Real>      _potential_ac;

  /**
   * the (complex) electrode current for AC sweep.
   */
  std::complex<Real>      _current_ac;


};


#endif //define __external_circuit_h__
