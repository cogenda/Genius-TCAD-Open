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




#include "external_circuit_rct.h"
#include "physical_unit.h"

ExternalCircuitRCTLine::ExternalCircuitRCTLine(Real r, Real Rl, Real Cl, Real length, int div)
  :_r_app(r),_r_per_um(Rl),_c_per_um(Cl),_length(length), N(div)
{
  _v.resize(N, 0.0);
  _v_last.resize(N, 0.0);
}


std::string ExternalCircuitRCTLine::format() const
{
  std::stringstream ss;
  ss  <<"string<ext.circuit>=tline"<<" "
      <<"real<res.load>="<<_r_app/(PhysicalUnit::V/PhysicalUnit::A)<<" "
      <<"real<res.um>="<<_r_per_um/(PhysicalUnit::V/PhysicalUnit::A/PhysicalUnit::um)<<" "
      <<"real<cap.um>="<<_c_per_um/(PhysicalUnit::C/PhysicalUnit::V/PhysicalUnit::um)<<" "
      <<"real<tl.length>="<<_length/(PhysicalUnit::um)<<" "
      <<"int<tl.div>="<<N<<" "
      <<"real<potential>="<<this->potential()/PhysicalUnit::V<<" "
      <<"real<potential.old>="<<this->potential_old()/PhysicalUnit::V;
  return ss.str();
}


Real ExternalCircuitRCTLine::inter_connect_resistance() const
{ return _r_app+_r_per_um*_length; }



Real ExternalCircuitRCTLine::mna_function(Real dt)
{
  Real res = _r_per_um*_length/N;
  Real cap = _c_per_um*_length/N;

  Real _V0 = _potential;

  std::vector<Real> A(N-1, 0.0); // low diag
  std::vector<Real> B(N,   0.0); // diag
  std::vector<Real> C(N-1, 0.0); // up diag
  std::vector<Real> r(N, 0.0);   // rhs
  for(int i=0; i<N; i++)
  {
    if(i==0)
    {
      B[i] = 1.0/res + 1.0/res + cap/dt;
      C[i] = -1.0/res;
      r[i] = _V0/res+cap*_v_last[i]/dt;
    }
    else if(i==N-1)
    {
      B[i] = 1.0/res + 1.0/_r_app + cap/dt;
      A[i-1] = -1.0/res;
      r[i] = _Vapp/_r_app+cap*_v_last[i]/dt;
    }
    else
    {
      B[i]   = 1.0/res + 1.0/res + cap/dt;
      C[i]   = -1.0/res;
      A[i-1] = -1.0/res;
      r[i] = cap*_v_last[i]/dt;
    }
  }

  solveMatrix(N, A, B, C, r, _v);

  return (_V0-_v[0])/res;
}


Real ExternalCircuitRCTLine::mna_jacobian(Real dt)
{
  Real res = _r_per_um*_length/N;
  Real cap = _c_per_um*_length/N;

  Real A = 1.0/res;
  std::vector<Real> B(N, 0.0);
  B[0] = -1.0/res;
  std::vector<Real> C(N, 0.0);
  C[0] = -1.0/res;

  std::vector<Real> Da(N-1, -1.0/res);
  std::vector<Real> Db(N, 2.0/res+cap/dt);
  std::vector<Real> Dc(N-1,-1.0/res);
  Db[N-1] = 1.0/res + 1.0/_r_app + cap/dt;

  //Real I = (A-B*D^-1*C)*_V0 + B/D*b - a
  // we need dI/d_V0 = (A-B*D^-1*C)

  // assume D^-1*C=M we solve DM=C
  std::vector<Real> M(N);
  solveMatrix(N, Da, Db, Dc, C, M);

  //Q=B*M
  Real Q = 0;
  for(int i=0; i<N; i++)
    Q += B[i]*M[i];

  return A-Q;
}


std::complex <Real> ExternalCircuitRCTLine::mna_ac_jacobian(Real omega)
{

  return 0.0;
}


void ExternalCircuitRCTLine::tran_op_init()
{
  ExternalCircuit::tran_op_init();

  Real res = _r_per_um*_length/N;

  Real _V0 = _potential;

  std::vector<Real> A(N-1, 0.0); // low diag
  std::vector<Real> B(N,   0.0); // diag
  std::vector<Real> C(N-1, 0.0); // up diag
  std::vector<Real> r(N, 0.0);   // rhs
  for(int i=0; i<N; i++)
  {
    if(i==0)
    {
      B[i] = 1.0/res + 1.0/res;
      C[i] = -1.0/res;
      r[i] = _V0/res;
    }
    else if(i==N-1)
    {
      B[i] = 1.0/res + 1.0/_r_app;
      A[i-1] = -1.0/res;
      r[i] = _Vapp/_r_app;
    }
    else
    {
      B[i]   = 1.0/res + 1.0/res;
      C[i]   = -1.0/res;
      A[i-1] = -1.0/res;
      r[i] = 0.0;
    }
  }

  solveMatrix(N, A, B, C, r, _v);

  _v_last = _v;
}



