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

#ifndef __boundary_condition_resistance_ohmic_h__
#define __boundary_condition_resistance_ohmic_h__


#include "boundary_condition.h"


/**
 * The Resistance Ohmic Boundary Condition
 */
class IF_Metal_OhmicBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  IF_Metal_OhmicBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~IF_Metal_OhmicBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
  { return IF_Metal_Ohmic; }

  /**
   * @return boundary condition type in string
   */
  virtual std::string bc_type_name() const
  { return "IF_Metal_Ohmic"; }

  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
  { return INTERFACE; }


  /**
   * @return true iff this boundary has a current flow
   */
  virtual bool has_current_flow() const
  { return true; }


  /**
   * @return average psi of this boundary.
   */
  virtual PetscScalar psi() const
  {return _average_psi;}

  /**
   * @return writable reference to average psi of this boundary.
   */
  virtual PetscScalar & psi()
  {return _average_psi;}


  /**
   * @return the current flow of this boundary.
   */
  virtual PetscScalar current() const
  {return _current_flow;}

  /**
   * @return writable reference to current flow of this boundary
   */
  virtual PetscScalar & current()
  { return _current_flow;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * current flow in/out
   */
  PetscScalar   _current_flow;

  /**
   * average psi
   */
  PetscScalar   _average_psi;

  /**
   * @return true when any of the electron/hole recomb velocity is infinity
   */
  bool _infinity_recombination() const
  {
    return scalar("elec.recomb.velocity") == std::numeric_limits<PetscScalar>::infinity() ||
           scalar("hole.recomb.velocity") == std::numeric_limits<PetscScalar>::infinity();
  }

  /**
   * Ohmic BC for DDM1, assuming infinity recombination rate
   */
  void _DDM1_Function_Infinite_Recombination(PetscScalar * , Vec , InsertMode &);

  /**
   * Ohmic BC for DDM1, assuming infinity recombination rate
   */
  void _DDM1_Jacobian_Infinite_Recombination(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * Ohmic BC for DDM1, assuming limited recombination rate
   */
  void _DDM1_Function_Limited_Recombination(PetscScalar * , Vec , InsertMode &);

  /**
   * Ohmic BC for DDM1, assuming limited recombination rate
   */
  void _DDM1_Jacobian_Limited_Recombination(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * Ohmic BC for DDM1, assuming infinity recombination rate
   */
  void _DDM1R_Function_Infinite_Recombination(PetscScalar * , Vec , InsertMode &);

  /**
   * Ohmic BC for DDM1, assuming infinity recombination rate
   */
  void _DDM1R_Jacobian_Infinite_Recombination(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * Ohmic BC for DDM1, assuming limited recombination rate
   */
  void _DDM1R_Function_Limited_Recombination(PetscScalar * , Vec , InsertMode &);

  /**
   * Ohmic BC for DDM1, assuming limited recombination rate
   */
  void _DDM1R_Jacobian_Limited_Recombination(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);


  /**
   * Ohmic BC for DDM2, assuming infinity recombination rate
   */
  void _DDM2_Function_Infinite_Recombination(PetscScalar * , Vec , InsertMode &);

  /**
   * Ohmic BC for DDM2, assuming infinity recombination rate
   */
  void _DDM2_Jacobian_Infinite_Recombination(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * Ohmic BC for DDM2, assuming limited recombination rate
   */
  void _DDM2_Function_Limited_Recombination(PetscScalar * , Vec , InsertMode &);

  /**
   * Ohmic BC for DDM2, assuming limited recombination rate
   */
  void _DDM2_Jacobian_Limited_Recombination(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * Ohmic BC for EBM3, assuming infinity recombination rate
   */
  void _EBM3_Function_Infinite_Recombination(PetscScalar * , Vec , InsertMode &);

  /**
   * Ohmic BC for EBM3, assuming infinity recombination rate
   */
  void _EBM3_Jacobian_Infinite_Recombination(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * Ohmic BC for EBM3, assuming limited recombination rate
   */
  void _EBM3_Function_Limited_Recombination(PetscScalar * , Vec , InsertMode &);

  /**
   * Ohmic BC for EBM3, assuming limited recombination rate
   */
  void _EBM3_Jacobian_Limited_Recombination(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);


  // function buffer
  std::vector<PetscScalar> _current_buffer;

  // Jacobian buffer
  std::vector< PetscInt > _buffer_rows;
  std::vector< std::vector<PetscInt> >    _buffer_cols;
  std::vector< std::vector<PetscScalar> > _buffer_jacobian_entries;

public:

#ifdef TCAD_SOLVERS

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data and scaling constant into petsc vector
   */
  virtual void Poissin_Fill_Value(Vec x, Vec L);

  /**
   * preprocess Jacobian function for poisson solver
   */
  virtual void Poissin_Function_Preprocess(PetscScalar *, Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * preprocess Jacobian Matrix for poisson solver
   */
  virtual void Poissin_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * update solution value of poisson's equation
   */
  virtual void Poissin_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data and scaling constant into petsc vector
   */
  virtual void DDM1_Fill_Value(Vec x, Vec L);

  /**
   * preprocess function for level 1 DDM solver
   */
  virtual void DDM1_Function_Preprocess(PetscScalar *, Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * preprocess Jacobian Matrix of level 1 DDM equation.
   */
  virtual void DDM1_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *);




  //////////////////////////////////////////////////////////////////////////////////
  //-------------Function and Jacobian evaluate for Density Gradient--------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of Density Gradient equation.
   */
  virtual void DG_Fill_Value(Vec x, Vec L);

  /**
   * preprocess function for Density Gradient solver
   */
  virtual void DG_Function_Preprocess(PetscScalar *, Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for Density Gradient solver
   */
  virtual void DG_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * preprocess Jacobian Matrix of Density Gradient equation.
   */
  virtual void DG_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for Density Gradient solver
   */
  virtual void DG_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * update solution data of Density Gradient solver.
   */
  virtual void DG_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * preprocess function for level 2 DDM solver
   */
  virtual void DDM2_Function_Preprocess(PetscScalar * ,Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * preprocess Jacobian Matrix of level 2 DDM equation.
   */
  virtual void DDM2_Jacobian_Preprocess(PetscScalar *,SparseMatrix<PetscScalar> *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * update solution data of DDML2 solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * preprocess function for level 3 EBM solver
   */
  virtual void EBM3_Function_Preprocess(PetscScalar *,Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * preprocess Jacobian Matrix of level 3 EBM equation.
   */
  virtual void EBM3_Jacobian_Preprocess(PetscScalar * ,SparseMatrix<PetscScalar> *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);



  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const PetscScalar omega, InsertMode & add_value_flag );

  /**
   * update solution value for ddm ac solver
   */
  virtual void DDMAC_Update_Solution(const PetscScalar * lxx , const Mat, const PetscScalar omega);


#ifdef COGENDA_COMMERCIAL_PRODUCT
  //////////////////////////////////////////////////////////////////////////////////
  //----------------- functions for Gummel DDML1 solver --------------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for preprocess build RHS and matrix for gummel equation.
   */
  virtual void DDM1_Gummel_Carrier_Preprocess(const std::string & carrier, Mat A, Vec r,
                                              std::vector<PetscInt> &src,  std::vector<PetscInt> &dst, std::vector<PetscInt> &clear);

  /**
   * function for build RHS and matrix for gummel equation.
   */
  virtual void DDM1_Gummel_Carrier(const std::string & carrier, PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);


  /**
   * function for reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Half_Implicit_Current_Reserve(Mat A, InsertMode &add_value_flag);

  /**
   * function for preprocess build RHS and matrix for half implicit current continuity equation.
   */
  virtual void DDM1_Half_Implicit_Current_Preprocess(Vec f, Mat A, std::vector<PetscInt> &src,  std::vector<PetscInt> &dst, std::vector<PetscInt> &clear);

  /**
   * function for build RHS and matrix for half implicit current continuity equation.
   */
  virtual void DDM1_Half_Implicit_Current(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);

  /**
   * function for update solution value for half implicit current continuity equation.
   */
  virtual void DDM1_Half_Implicit_Current_Update_Solution(PetscScalar *);

  /**
   * function for preprocess build RHS and matrix for half implicit poisson correction equation.
   */
  virtual void DDM1_Half_Implicit_Poisson_Correction_Preprocess(Vec f, std::vector<PetscInt> &src,  std::vector<PetscInt> &dst, std::vector<PetscInt> &clear);

  /**
   * function for build RHS and matrix for half implicit poisson correction equation.
   */
  virtual void DDM1_Half_Implicit_Poisson_Correction(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);

#endif

#endif


};





#endif

