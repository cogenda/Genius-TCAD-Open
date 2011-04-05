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

//  $Id: external_circuit.h,v 1.10 2008/07/09 05:58:16 gdiso Exp $


#ifndef __external_circuit_h__
#define __external_circuit_h__

#include <complex>


#include "solver_specify.h"

/**
 * this is a class for external circuit for electrode bcs
 */
class ExternalCircuit
{
public:

  enum DRIVEN{VDRIVEN, IDRIVEN, FLOAT};

  /**
   * constructor
   */
  ExternalCircuit()
  : _res(0), _cap(0), _ind(0), _Vapp(0), _Iapp(0),
    _potential(0), _current(0), _current_displacement(0),
    _current_conductance(0), _cap_current(0),
    _Vac(0.0026),
    _drv(VDRIVEN)
  {}

  /**
   * destructor
   */
  virtual ~ExternalCircuit(){}


  /**
   * @return the lumped resistor
   */
  virtual PetscScalar R() const
  {return _res;}

  /**
   * @return writable reference to lumped resistor
   */
  virtual PetscScalar & R()
  {return _res;}

  /**
   * @return the lumped capacitance
   */
  virtual PetscScalar C() const
  {return _cap;}

  /**
   * @return writable reference to lumped capacitance
   */
  virtual PetscScalar & C()
  {return _cap;}

  /**
   * @return the lumped inductance
   */
  virtual PetscScalar L() const
  {return _ind;}

  /**
   * @return writable reference to lumped inductance
   */
  virtual PetscScalar & L()
  {return _ind;}




  /**
   * @return the application voltage.
   */
  virtual PetscScalar Vapp() const
  { return _Vapp;}

  /**
   * @return writable reference to application voltage.
   */
  virtual PetscScalar & Vapp()
  { return _Vapp;}

  /**
   * @return the application current.
   */
  virtual PetscScalar Iapp() const
  { return _Iapp;}

  /**
   * @return writable reference to application current.
   */
  virtual PetscScalar & Iapp()
  { return _Iapp;}




  /**
   * @return the electrode potential.
   */
  virtual PetscScalar potential() const
  { return _potential;}

  /**
   * @return writable reference to electrode potential.
   */
  virtual PetscScalar & potential()
  { return _potential;}

  /**
   * @return electrode potential of last step.
   * @note this variable is read only for other application
   */
  virtual PetscScalar  potential_old() const
  { return _potential_old;}

  /**
   * @return electrode potential of current iterating.
   * @note write only
   */
  virtual PetscScalar  & potential_itering()
  { return _potential_itering;}


  /**
   * @return the electrode current.
   */
  virtual PetscScalar current() const
  { return _current;}

  /**
   * @return writable reference to electrode current.
   */
  virtual PetscScalar & current()
  { return _current;}


  /**
   * @return the displacement current.
   */
  virtual PetscScalar current_displacement() const
  { return _current_displacement;}

  /**
   * @return writable reference to displacement current.
   */
  virtual PetscScalar & current_displacement()
  { return _current_displacement;}

  /**
   * @return the conductance current.
   */
  virtual PetscScalar current_conductance() const
  { return _current_conductance;}

  /**
   * @return writable reference to conductance current.
   */
  virtual PetscScalar & current_conductance()
  { return _current_conductance;}


  /**
   * @return electrode current of last step.
   * @note this variable is read only for other application
   */
  virtual PetscScalar  current_old() const
  { return _current_old;}

  /**
   * @return electrode current of current iterating.
   * @note write only
   */
  virtual PetscScalar  & current_itering()
  { return _current_itering;}

  /**
   * @return cap current
   * @note this variable is read only for other application
   */
  virtual PetscScalar  cap_current() const
  { return _cap_current;}



  /**
   * @return the stimulate voltage for AC scan.
   */
  virtual PetscScalar Vac() const
  { return _Vac; }

  /**
   * @return the writable reference to stimulate voltage for AC scan.
   */
  virtual PetscScalar & Vac()
  { return _Vac; }

  /**
   * @return the electrode potential of AC scan.
   */
  virtual std::complex<PetscScalar> potential_ac() const
  { return _potential_ac;}

  /**
   * @return writable reference to electrode potential of AC scan.
   */
  virtual std::complex<PetscScalar> & potential_ac()
  { return _potential_ac;}


  /**
   * @return the electrode current of AC scan.
   */
  virtual std::complex<PetscScalar> current_ac() const
  { return _current_ac;}

  /**
   * @return writable reference to electrode current of AC scan.
   */
  virtual std::complex<PetscScalar> & current_ac()
  { return _current_ac;}


  /**
   * when a solution is achieved, update the potential/current
   */
  void update()
  {
    _potential_old = _potential;
    _potential = _potential_itering;

    _current_old = _current;
    _current = _current_itering;

    _cap_current = _cap*(_potential-_potential_old)/SolverSpecify::dt;
  }



  /**
   * force to ground this electrode
   */
  virtual void ground()
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


private:

  /**
   * external lumped resistor
   */
  PetscScalar   _res;

  /**
   * external lumped capacitance
   */
  PetscScalar   _cap;

  /**
   * external lumped inductance
   */
  PetscScalar   _ind;

  /**
   * application voltage
   */
  PetscScalar   _Vapp;

  /**
   * application current
   */
  PetscScalar   _Iapp;


  /**
   * the potential of this electrode
   */
  PetscScalar      _potential;

  /**
   * the potential for this iterative cycle
   */
  PetscScalar      _potential_itering;

  /**
   * the potential for last step
   */
  PetscScalar      _potential_old;

  /**
   * the current flow out of this electrode
   */
  PetscScalar      _current;

  /**
   * the displacement current flow out of this electrode
   */
  PetscScalar      _current_displacement;

  /**
   * the conductance current flow out of this electrode
   */
  PetscScalar      _current_conductance;

  /**
   * the current for this iterative cycle
   */
  PetscScalar      _current_itering;

  /**
   * the current for last step
   */
  PetscScalar      _current_old;

  /**
   * the current which flow through lumped capacitance to ground
   */
  PetscScalar      _cap_current;

  /**
   * the application voltage for AC sweep.
   */
  PetscScalar      _Vac;

  /**
   * the (complex) electrode potential for AC sweep.
   */
  std::complex<PetscScalar>      _potential_ac;

  /**
   * the (complex) electrode current for AC sweep.
   */
  std::complex<PetscScalar>      _current_ac;

  /**
   * indicator this electrode is driven by voltage source or
   * current source
   */
  DRIVEN        _drv;

};


#endif //define __external_circuit_h__
