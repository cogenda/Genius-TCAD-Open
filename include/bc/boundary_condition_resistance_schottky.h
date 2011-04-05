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

#ifndef __boundary_condition_resistance_schottky_h__
#define __boundary_condition_resistance_schottky_h__


#include "boundary_condition.h"


/**
 * The Resistance Schottky Boundary Condition
 */
class IF_Metal_SchottkyBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  IF_Metal_SchottkyBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~IF_Metal_SchottkyBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
  { return IF_Metal_Schottky; }

  /**
   * @return boundary condition type in string
   */
  virtual std::string bc_type_name() const
  { return "MetalSchottky"; }

  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
  { return INTERFACE; }


  /**
   * @return true when it has external circuit
   */
  virtual bool is_electrode()  const
  {return ext_circuit()!=NULL;}

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

  PetscScalar _current_flow;

  /**
   * average psi
   */
  PetscScalar   _average_psi;

public:

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
  virtual void Poissin_Function_Preprocess(Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * preprocess Jacobian Matrix for poisson solver
   */
  virtual void Poissin_Jacobian_Preprocess(Mat *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &);

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
  virtual void DDM1_Function_Preprocess(Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * preprocess Jacobian Matrix of level 1 DDM equation.
   */
  virtual void DDM1_Jacobian_Preprocess(Mat *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * preprocess function for level 2 DDM solver
   */
  virtual void DDM2_Function_Preprocess(Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * preprocess Jacobian Matrix of level 2 DDM equation.
   */
  virtual void DDM2_Jacobian_Preprocess(Mat *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * preprocess function for level 3 EBM solver
   */
  virtual void EBM3_Function_Preprocess(Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * preprocess Jacobian Matrix of level 3 EBM equation.
   */
  virtual void EBM3_Jacobian_Preprocess(Mat *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);



  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag ){}

  /**
   * update solution value for ddm ac solver
   */
  virtual void DDMAC_Update_Solution(const PetscScalar * lxx , const Mat, const double omega){}

  //////////////////////////////////////////////////////////////////////////////////
  //----------------- functions for Fast Hydrodynamic solver  --------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for evaluating boundary for HDM method
   */
  virtual void HDM_Boundary( const PetscScalar * , Vec /*x*/, InsertMode &  ){}

  //////////////////////////////////////////////////////////////////////////////////
  //-----------------  functions for Linear Poissin solver   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for reserve none zero pattern in petsc matrix.
   */
  virtual void LinearPoissin_Reserve(Mat A, InsertMode &){}

  /**
   * function for build matrix of linear poisson's equation.
   */
  virtual void LinearPoissin_Matrix(Mat A, InsertMode &){}

  /**
   * function for build RHS vector of linear poisson's equation.
   */
  virtual void LinearPoissin_RHS(Vec b, InsertMode &){}

};





#endif

