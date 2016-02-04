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
 
#ifndef __schur_solver_h__
#define __schur_solver_h__

#include <iostream>
#include <vector>
#include <map>

/**
 * Schur complement operator
 *
 * | A11 A12 | x1     b1
 * |         |     =  
 * | A21 A22 | x2     b2
 *
 * SchurM   = A11 - A12*A22^-1*A21
 * SchurRHS = b1 - A12*A22^-1*b2
 * We can solve x1 as
 * SchurM*x1 = SchurRHS
 */
class SchurSolver
{
public:

  /**
   * empty constructor
   */
  SchurSolver() {}


  /**
   * constructor, set matrix size and schur block A11 and non schur block A22
   */
  SchurSolver(unsigned int N, const std::vector<unsigned int> &schur_block, const std::vector<unsigned int> &non_schur_block)
  { build(N, schur_block, non_schur_block); }


  ~SchurSolver();


  /**
   * set matrix size and schur block A11 and non schur block A22
   */
  void build(unsigned int N, const std::vector<unsigned int> &schur_block, const std::vector<unsigned int> &non_schur_block);

  /**
   * zero matrix entries
   */
  void MatZero();

  /**
   * zero rhs entries
   */
  void RHSZero();

  /**
   * set value to mat entry
   */
  bool MatSetValue(unsigned int i, unsigned int j, double v, bool add=false);

  /**
   * set value to vector entry
   */
  bool RHSSetValue(unsigned int i, double v, bool add=false);

  /**
   * do schur solve!
   */
  bool SchurSolve();


  /**
   * get SchurMatrix A11 - A12*A22^-1*A21 as dense matrix
   */
  void SchurMatrix(std::vector<double > & mat_entries) const;

  /**
   * get SchurRHS as b1 - A12*A22^-1*b2
   */
  void SchueRHS( std::vector<double> &vec_entries) const
  { vec_entries = v_schur.v; }


  /**
   * solve x1! 
   */
  void SchueSolveX1(std::vector<double> & x1);


  /**
   * by given x1, then we can solve x2 by  A22*x2 = b2 - A21*x1
   */
  void SchueSolveX( const std::vector<double> & x1);

  /**
   * get X solved by SchueSolveX
   */
  void XGetValue( std::vector<double> & x) const;


private:

  // mat
  struct Mat
  {
    /// row and col size
    unsigned int m,n;
    ///triple mat entries
    std::map<std::pair<unsigned int, unsigned int>, double > triple_entries;
    /// set value
    void set_value(unsigned int i, unsigned int j, double v, bool add);
    /// zero entries
    void zero();
    /// convert triple to csc
    void triple_to_csc(std::vector<int> &Ap, std::vector<int> &Ai, std::vector<double> &Ax) const;
    /// convert to row vectors
    void triple_to_rvs(std::vector<std::vector<double> > &rvs) const;
    /// convert to col vectors
    void triple_to_cvs(std::vector<std::vector<double> > &cvs) const;
    /// get value
    double get_value(unsigned int i, unsigned int j) const;
    ///Formatted print to \p std::cout.
    void print(std::ostream& os) const;
    ///
    friend std::ostream& operator << (std::ostream& os, const Mat& m)
    {
      m.print(os);
      return os;
    }
  };


  struct Vec
  {
    /// vec size
    unsigned int n;
    /// vec entries
    std::vector<double> v;
    /// set value
    void set_value(unsigned int i, double v, bool add);
    /// get value
    double get_value(unsigned int i) const;
    /// zero
    void zero();
    ///Formatted print to \p std::cout.
    void print(std::ostream& os) const;
    ///
    friend std::ostream& operator << (std::ostream& os, const Vec& v)
    {
      v.print(os);
      return os;
    }
  };


  Mat A11;
  Mat A12;
  Mat A21;
  Mat A22;


  Vec x1;
  Vec x2;

  Vec b1;
  Vec b2;


  Mat M_schur;
  Vec v_schur;

  /**
   * Permutation vector, row_perm[old_index] = new_index
   * build when schur block given.
   * this makes sure schur block is at the begin of matrix 
   */
  std::vector<unsigned int> row_perm;


  /**
   * Permutation vector, col_perm[old_index] = new_index
   * build when schur block given.
   * this makes sure schur block is at the begin of matrix
   */
  std::vector<unsigned int> col_perm;

};


#endif
