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

#ifndef __boundary_condition_solderpad_h__
#define __boundary_condition_solderpad_h__

#include "boundary_condition.h"

/**
 * The pad Boundary Condition
 */
class SolderPadBC : public BoundaryCondition
{
public:

  /**
   * constructor, set default value
   */
  SolderPadBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~SolderPadBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return SolderPad; }

  /**
   * @return boundary condition type in string
   */
  virtual std::string bc_type_name() const
  { return "SolderPad"; }

  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return BOUNDARY; }

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
   * @return the heat transfer rate of this boundary
   */
  virtual PetscScalar Heat_Transfer() const
  {return _Heat_Transfer;}

  /**
   * @return writable reference to heat transfer rate of this boundary
   */
  virtual PetscScalar & Heat_Transfer()
  {return _Heat_Transfer;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * heat transfer rate of this boundary
   */
  PetscScalar   _Heat_Transfer;

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
   * fill solution data into petsc vector of level 1 DDM equation.
   */
  virtual void DDM1_Fill_Value(Vec x, Vec L);

  /**
   * preprocess function for level 1 DDM solver
   */
  virtual void DDM1_Function_Preprocess(Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * preprocess Jacobian Matrix of level 1 DDM equation.
   */
  virtual void DDM1_Jacobian_Preprocess(Mat *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag);

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *x);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed DDML1-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 1 DDM equation.
   */
  virtual void MixA_DDM1_Fill_Value(Vec x, Vec L);

  /**
   * preprocess function for Advanced Mixed DDML1 solver
   */
  virtual void MixA_DDM1_Function_Preprocess(Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for Advanced Mixed DDML1 solver
   */
  virtual void MixA_DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * preprocess Jacobian Matrix for Advanced Mixed type level 1 DDM equation.
   */
  virtual void MixA_DDM1_Jacobian_Preprocess(Mat *, std::vector<PetscInt> &,  std::vector<PetscInt> &, std::vector<PetscInt> &);

  /**
   * build function and its jacobian for Advanced Mixed DDML1 solver
   */
  virtual void MixA_DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 2 DDM equation.
   */
  virtual void DDM2_Fill_Value(Vec x, Vec L);

  /**
   * preprocess function for level 2 DDM solver
   */
  virtual void DDM2_Function_Preprocess(Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * preprocess Jacobian Matrix of level 2 DDM equation.
   */
  virtual void DDM2_Jacobian_Preprocess(Mat *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag);

  /**
   * update solution data of level 2 DDM solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *x);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed DDML2-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 2 DDM equation.
   */
  virtual void MixA_DDM2_Fill_Value(Vec x, Vec L);

  /**
   * preprocess function for Advanced Mixed DDML2 solver
   */
  virtual void MixA_DDM2_Function_Preprocess(Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for Advanced Mixed DDML2 solver
   */
  virtual void MixA_DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * preprocess Jacobian Matrix for Advanced Mixed type level 2 DDM equation.
   */
  virtual void MixA_DDM2_Jacobian_Preprocess(Mat *, std::vector<PetscInt> &,  std::vector<PetscInt> &, std::vector<PetscInt> &);

  /**
   * build function and its jacobian for Advanced Mixed DDML2 solver
   */
  virtual void MixA_DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 3 EBM equation.
   */
  virtual void EBM3_Fill_Value(Vec x, Vec L);

  /**
   * preprocess function for level 3 EBM solver
   */
  virtual void EBM3_Function_Preprocess(Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag);

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
  virtual void EBM3_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag);

  /**
   * update solution data of EBM3 solver.
   */
  virtual void EBM3_Update_Solution(PetscScalar *x);



  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed EBM3 -------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of EBM3 equation.
   */
  virtual void MixA_EBM3_Fill_Value(Vec x, Vec L);

  /**
   * preprocess function for Advanced Mixed EBM3 solver
   */
  virtual void MixA_EBM3_Function_Preprocess(Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for Advanced Mixed EBM3 solver
   */
  virtual void MixA_EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_EBM3_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * preprocess Jacobian Matrix for Advanced Mixed type EBM3 equation.
   */
  virtual void MixA_EBM3_Jacobian_Preprocess(Mat *, std::vector<PetscInt> &,  std::vector<PetscInt> &, std::vector<PetscInt> &);

  /**
   * build function and its jacobian for Advanced Mixed EBM3 solver
   */
  virtual void MixA_EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);



  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );

  /**
   * update solution value for ddm ac solver
   */
  virtual void DDMAC_Update_Solution(const PetscScalar * lxx , const Mat, const double omega);

};






#endif

