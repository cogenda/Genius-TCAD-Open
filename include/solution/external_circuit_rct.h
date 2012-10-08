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



#ifndef __external_circuit_rct_h__
#define __external_circuit_rct_h__

#include <vector>

#include "external_circuit.h"


  // For voltage driven
  //
  //         _____        _____     V0
  //    ----|_____|------|_____|--------> to electrode (Ve, I)
  //    |      r     |       R
  //   Vapp         === C
  //    |____________|______________
  //                GND
  //
  //  (V0-V1)/R  - I         = 0
  //  (Vi-Vi_1)/R + (Vi-Vi+1)/R + C*dVi/dt = 0
  //  (Vn-Vapp)/r + (Vn-Vn_1)/R + C*dVn/dt =0


  // And for current driven
  //                       _____   Ve
  //    ------------------|_____|--------> to electrode (Ve, I)
  //    |            |       R
  //   Iapp         === C
  //    |____________|______________
  //                GND
  //
  //  (V0-V1)/R  - I         = 0
  //  (Vi-Vi_1)/R + (Vi-Vi+1)/R + C*dVi/dt = 0
  //  -Ia + (Vn-Vn_1)/R + C*dVn/dt =0



class ExternalCircuitRCTLine : public ExternalCircuit
{
public:
  /**
   * @param r      serial resistance
   * @param Rl     resistance per um
   * @param Cl     capacitance per um
   * @param length the length of the transmition line
   */
  ExternalCircuitRCTLine(Real r, Real Rl, Real Cl, Real length, int div);

  virtual ~ExternalCircuitRCTLine() {}

  /**
   * type of external circuit
   */
  virtual std::string type() const { return "RCT"; }

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
   * value of serial resistance for inter connect
   */
  virtual Real inter_connect_resistance() const;

  /**
   * value of serial resistance to electrode
   */
  virtual Real serial_resistance() const { return _r_app+_r_per_um*_length; }


  /**
   * value of serial resistance to electrode
   */
  virtual void set_serial_resistance(Real res) { _r_app = res-_r_per_um*_length; }


  /**
   * use this value as scaling to electrode
   */
  virtual Real electrode_scaling(Real dt) const { return std::min(_r_app,_r_per_um*_length); }


  /**
   * use this value as scaling to Modified Nodal Analysis
   */
  virtual Real mna_scaling(Real dt) const  { return 1.0;  }

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
    _v_last = _v;
  }

  /**
   * roll back to previous solution
   */
  virtual void rollback()
  {
    ExternalCircuit::rollback();
    _v = _v_last;
  }

  /**
   * init op state before transient simulation
   */
  virtual void tran_op_init();

private:

  Real _r_app;

  Real _r_per_um;

  Real _c_per_um;

  Real _length;

  /**
   * N segments of TL
   */
  int N;

  std::vector<Real> _v;

  std::vector<Real> _v_last;


  void solveMatrix (int n, std::vector<Real> a, std::vector<Real> b, std::vector<Real> c, std::vector<Real> v, std::vector<Real> &x)
  {
   /**
    * n - number of equations
    * a - sub-diagonal (means it is the diagonal below the main diagonal) -- indexed from 0..n-2
    * b - the main diagonal
    * c - sup-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
    * v - right part
    * x - the answer
    */
    for (int i = 1; i < n; i++)
    {
      Real m = a[i-1]/b[i-1];
      b[i] = b[i] - m * c[i - 1];
      v[i] = v[i] - m*v[i-1];
    }

    x[n-1] = v[n-1]/b[n-1];

    for (int i = n - 2; i >= 0; --i)
      x[i] = (v[i] - c[i] * x[i+1]) / b[i];
  }


  void solveMatrix2 (int n, const std::vector<Real> &a, const std::vector<Real> &b, const std::vector<Real> &c, const std::vector<Real> &v, std::vector<Real> &x)
  {
   /**
    * n - number of equations
    * a - sub-diagonal (means it is the diagonal below the main diagonal) -- indexed from 0..n-2
    * b - the main diagonal
    * c - sup-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
    * v - right part
    * x - the answer
    */
    std::vector<Real> d(n);
    std::vector<Real> y(n);
    for (int i = 1; i < n; i++)
    {
      Real m = a[i-1]/b[i-1];
      d[i] = b[i] - m * c[i - 1];
      y[i] = v[i] - m*v[i-1];
    }

    x[n-1] = y[n-1]/d[n-1];

    for (int i = n - 2; i >= 0; --i)
      x[i] = (y[i] - c[i] * x[i+1]) / d[i];
  }

};

#endif
