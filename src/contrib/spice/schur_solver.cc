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
#include <cassert>


#include "../klu/klu.h"
#include "schur_solver.h"





void SchurSolver::build(unsigned int N, const std::vector<unsigned int> & schur_block, const std::vector<unsigned int> & non_schur_block)
{
  unsigned int perm_order = 0;
  col_perm.resize(N);
  row_perm.resize(N);
  for(std::vector<unsigned int>::const_iterator it = schur_block.begin() ; it!=schur_block.end(); ++it)
    row_perm[*it] = col_perm[*it] = perm_order++;

  for(std::vector<unsigned int>::const_iterator it = non_schur_block.begin() ; it!=non_schur_block.end(); ++it)
    row_perm[*it] = col_perm[*it] = perm_order++;

  A11.m = A11.n = schur_block.size();
  A22.m = A22.n = non_schur_block.size();

  A12.m = A21.n = non_schur_block.size();
  A12.n = A21.m = schur_block.size();

  b1.n = schur_block.size();
  b1.v.resize( schur_block.size() );

  b2.n = non_schur_block.size();
  b2.v.resize( non_schur_block.size() );

  x1.n = schur_block.size();
  x1.v.resize( schur_block.size() );

  x2.n = non_schur_block.size();
  x2.v.resize( non_schur_block.size() );

  M_schur.m = M_schur.n = schur_block.size();
  v_schur.n = schur_block.size();
  v_schur.v.resize( schur_block.size() );
}


SchurSolver::~SchurSolver()
{}



void SchurSolver::MatZero()
{
  A11.zero();
  A12.zero();
  A21.zero();
  A22.zero();
  M_schur.zero();
}


void SchurSolver::RHSZero()
{
  b1.zero();
  b2.zero();
  v_schur.zero();

  x1.zero();
  x2.zero();
}


bool SchurSolver::MatSetValue(unsigned int i, unsigned int j, double v, bool add)
{
  unsigned int row = row_perm[i];
  unsigned int col = col_perm[j];

  if(row < A11.m && col < A11.n)  //A11
    A11.set_value(row, col, v, add);
  if(row < A11.m && col >= A11.n) //A12
    A12.set_value(row, col-A11.n, v, add);
  if(row >= A11.m && col < A11.n) //A21
    A21.set_value(row-A11.m, col, v, add);
  if(row >= A11.m && col >= A11.n) //A22
    A22.set_value(row-A11.m, col-A11.n, v, add);

  return true;
}


bool SchurSolver::RHSSetValue(unsigned int i, double v, bool add)
{
  unsigned int row = row_perm[i];
  if( row < A11.m ) //b1
  {
    b1.set_value(row, v, add);
  }
  else //b2
  {
    b2.set_value(row-A11.m, v, add);
  }
  return true;
}


bool SchurSolver::SchurSolve()
{
  /*
  std::cout<<"Matrix"<<std::endl;
  std::cout<<A11;
  std::cout<<A12;
  std::cout<<A21;
  std::cout<<A22;
  
  std::cout<<"rhs"<<std::endl;
  std::cout<<b1;
  std::cout<<b2;
  */

  // A11 - A12*A22^-1*A21
  // image A22^-1*A21=K, solve A22*K=A21

  std::vector<int> A22p;
  std::vector<int> A22i;
  std::vector<double> A22x;
  A22.triple_to_csc(A22p, A22i, A22x);

  std::vector<std::vector<double> > cvs;
  A21.triple_to_cvs(cvs);

  // image A22^-1*b2 = M, solve A22*M=b2;
  std::vector<double> b = b2.v;

  // use KLU
  klu_symbolic *Symbolic;
  klu_numeric *Numeric;
  klu_common Common;
  klu_defaults (&Common);

  Symbolic = klu_analyze(A22p.size()-1, &A22p[0], &A22i[0], &Common) ;
  assert(Common.status==KLU_OK);

  Numeric = klu_factor(&A22p[0], &A22i[0], &A22x[0], Symbolic, &Common) ;
  if(Common.status!=KLU_OK)
  {
    std::cout<<"A22 error with singular col " << Common.singular_col << std::endl;
  }

  for(unsigned int n=0; n<cvs.size(); ++n)
    klu_solve (Symbolic, Numeric, cvs[n].size(), 1, &(cvs[n][0]), &Common) ;

  klu_solve (Symbolic, Numeric, b.size(), 1, &(b[0]), &Common) ;// now b=M

  klu_free_symbolic (&Symbolic, &Common) ;
  klu_free_numeric (&Numeric, &Common) ;

  std::vector<std::vector<double> > rvs;
  A12.triple_to_rvs(rvs);

  // now K=cvs and b=M

  // A11-A12*K
  M_schur.triple_entries = A11.triple_entries;

  assert(rvs.size() == cvs.size());
  for(unsigned int i=0; i<rvs.size(); ++i)
  {
    const std::vector<double> & row = rvs[i];
    for(unsigned int j=0; j<cvs.size(); ++j)
    {
      const std::vector<double> & col = cvs[j];
      for(unsigned int k=0; k<row.size(); ++k)
        M_schur.set_value(i, j, -row[k]*col[k], true);
    }
  }

  // b1 - A12*A22^-1*b2
  v_schur.v = b1.v;
  for(unsigned int i=0; i<rvs.size(); ++i)
  {
    const std::vector<double> & row = rvs[i];
    assert(b.size() == row.size());
    double r = 0.0;
    for(unsigned int j=0; j<b.size(); ++j)
      r += row[j]*b[j];
    v_schur.set_value(i, -r, true);
  }

  /*
  std::cout<<"M schur"<<std::endl;
  std::cout<<M_schur;
  std::cout<<"v schur"<<std::endl;
  std::cout<<v_schur;
  */

  return true;
}



void SchurSolver::SchurMatrix(std::vector<double > & mat_entries) const
{
  for(unsigned int i=0; i<M_schur.m; i++)
    for(unsigned int j=0; j<M_schur.n; j++)
      mat_entries.push_back( M_schur.get_value(i, j) );
}



void SchurSolver::SchueSolveX1(std::vector<double> & x)
{
  x = v_schur.v;

  std::vector<int> Ap;
  std::vector<int> Ai;
  std::vector<double> Ax;
  M_schur.triple_to_csc(Ap, Ai, Ax);

  // use KLU
  klu_symbolic *Symbolic;
  klu_numeric *Numeric;
  klu_common Common;
  klu_defaults (&Common);
  Symbolic = klu_analyze(Ap.size()-1, &Ap[0], &Ai[0], &Common) ;
  Numeric = klu_factor(&Ap[0], &Ai[0], &Ax[0], Symbolic, &Common) ;

  klu_solve (Symbolic, Numeric, x.size(), 1, &(x[0]), &Common) ;

  klu_free_symbolic (&Symbolic, &Common) ;
  klu_free_numeric (&Numeric, &Common) ;
}




void SchurSolver::SchueSolveX( const std::vector<double> & x)
{
  //we can solve x2 by  A22*x2 = b2 - A21*x1
  x1.v = x;

  std::vector<int> A22p;
  std::vector<int> A22i;
  std::vector<double> A22x;
  A22.triple_to_csc(A22p, A22i, A22x);

  std::vector<std::vector<double> > rvs;
  A21.triple_to_rvs(rvs);

  std::vector<double> b = b2.v;
  for(unsigned int i=0; i<rvs.size(); ++i)
  {
    const std::vector<double> & row = rvs[i];
    assert(row.size() == x.size());
    for(unsigned int j=0; j<x.size(); ++j)
    {
      b[i] -= row[j]*x[j];
    }
  }
  
  // use KLU
  klu_symbolic *Symbolic;
  klu_numeric *Numeric;
  klu_common Common;
  klu_defaults (&Common);
  Symbolic = klu_analyze(A22p.size()-1, &A22p[0], &A22i[0], &Common) ;
  Numeric = klu_factor(&A22p[0], &A22i[0], &A22x[0], Symbolic, &Common) ;

  klu_solve (Symbolic, Numeric, b.size(), 1, &(b[0]), &Common) ;

  klu_free_symbolic (&Symbolic, &Common) ;
  klu_free_numeric (&Numeric, &Common) ;

  x2.v = b;
}


void SchurSolver::XGetValue( std::vector<double> & x) const
{
  x.resize(x1.n + x2.n);
  for(unsigned int i=0; i<x.size(); i++)
  {
    unsigned int col = col_perm[i];
    if( col < A11.m ) //x1
    {
      x[i] = x1.v[col];
    }
    else //x2
    {
      x[i] = x2.v[col-A11.m];
    }
  }
}



//---------------------------------------------------------------

void SchurSolver::Mat::set_value(unsigned int i, unsigned int j, double value, bool add)
{
  if(add)
    triple_entries[std::make_pair(i,j)] += value;
  else
    triple_entries[std::make_pair(i,j)] =  value;
}

void SchurSolver::Mat::zero()
{
  std::map<std::pair<unsigned int, unsigned int>, double >::iterator it = triple_entries.begin();
  std::map<std::pair<unsigned int, unsigned int>, double >::iterator it_end = triple_entries.end();
  for(; it != it_end; ++it)
    it->second = 0.0;
}


void SchurSolver::Mat::triple_to_csc(std::vector<int> &Ap, std::vector<int> &Ai, std::vector<double> &Ax) const
{
  Ap.resize(n+1);
  Ai.resize(triple_entries.size());
  Ax.resize(triple_entries.size());

  //first scan
  std::vector<int> col_entry_size(n);
  std::map<std::pair<unsigned int, unsigned int>, double >::const_iterator it = triple_entries.begin();
  std::map<std::pair<unsigned int, unsigned int>, double >::const_iterator it_end = triple_entries.end();
  for(; it != it_end; ++it)
    col_entry_size[it->first.second]++;

  for(unsigned int i=0; i<col_entry_size.size(); ++i)
    Ap[i+1] = Ap[i]+col_entry_size[i];

  //second scan
  std::vector<int> col_entry_count(n);
  for(it = triple_entries.begin(); it != it_end; ++it)
  {
    unsigned int col = it->first.second;
    int col_offsite = Ap[col] + col_entry_count[col];
    Ai[col_offsite] = it->first.first;
    Ax[col_offsite] = it->second;
    col_entry_count[col]++;
  }
}


void SchurSolver::Mat::triple_to_rvs(std::vector<std::vector<double> > &rvs) const
{
  rvs.resize(n);
  for(unsigned int i=0; i<n; i++)
    rvs[i].resize(m, 0.0);

  std::map<std::pair<unsigned int, unsigned int>, double >::const_iterator it = triple_entries.begin();
  std::map<std::pair<unsigned int, unsigned int>, double >::const_iterator it_end = triple_entries.end();
  for(; it != it_end; ++it)
  {
    unsigned int row = it->first.first;
    unsigned int col = it->first.second;
    rvs[row][col] = it->second;
  }
}

void SchurSolver::Mat::triple_to_cvs(std::vector<std::vector<double> > &cvs) const
{
  cvs.resize(m);
  for(unsigned int i=0; i<m; i++)
    cvs[i].resize(n, 0.0);

  std::map<std::pair<unsigned int, unsigned int>, double >::const_iterator it = triple_entries.begin();
  std::map<std::pair<unsigned int, unsigned int>, double >::const_iterator it_end = triple_entries.end();
  for(; it != it_end; ++it)
  {
    unsigned int row = it->first.first;
    unsigned int col = it->first.second;
    cvs[col][row] = it->second;
  }
}


double SchurSolver::Mat::get_value(unsigned int i, unsigned int j) const
{
  std::map<std::pair<unsigned int, unsigned int>, double >::const_iterator it = triple_entries.find( std::make_pair(i, j));
  if( it == triple_entries.end() ) return 0.0;
  return it->second;
}


void SchurSolver::Mat::print(std::ostream& os) const
{
  std::vector<std::vector<double> > rvs;
  triple_to_rvs(rvs);

  for(unsigned int i=0; i< rvs.size(); ++i)
  {  
    for(unsigned int j=0; j<rvs[i].size(); ++j)
      os << rvs[i][j] << " ";
    os << std::endl;
  }
  os << std::endl;
}


//---------------------------------------------------------------

void SchurSolver::Vec::set_value(unsigned int i, double value, bool add)
{
  if(add)
    v[i] += value;
  else
    v[i] = value;
}

double SchurSolver::Vec::get_value(unsigned int i) const
{
  return v[i];
}

void SchurSolver::Vec::zero()
{
  std::fill(v.begin(), v.end(), 0.0);
}

void SchurSolver::Vec::print(std::ostream& os) const
{
  for(unsigned int i=0; i<v.size(); ++i)
    os << v[i] << " ";
  os << std::endl;
  os << std::endl;
}

