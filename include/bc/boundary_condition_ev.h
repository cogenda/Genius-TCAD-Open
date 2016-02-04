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

#ifndef __boundary_condition_ev_h__
#define __boundary_condition_ev_h__


#include "boundary_condition.h"


/**
 * The electrode-to-vacuum interface
 */
class ElectrodeVacuumInterfaceBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  ElectrodeVacuumInterfaceBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~ElectrodeVacuumInterfaceBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
  { return IF_Electrode_Vacuum; }

  /**
   * @return boundary condition type in string
   */
  virtual std::string bc_type_name() const
  { return "IF_Electrode_Vacuum"; }

  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTERFACE; }

  /**
   * @return false
   */
  virtual bool is_electrode()  const
    {return true;}

  /**
   * @return true iff this boundary has a current flow
   */
  virtual bool has_current_flow() const
  { return true; }

  /**
   * @return the current flow of this boundary.
   */
  virtual PetscScalar current() const
  {return 0.0;}

  /**
   * find nearest points in gate region
   * find elements near the semiconductor-insulator interface
   */
  virtual void prepare_for_use() {}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;


public:

#ifdef TCAD_SOLVERS

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &)
  { /* no thing to do for Neumann boundary */ }

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &)
  { /* no thing to do for Neumann boundary */ }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &)
  { /* no thing to do for Neumann boundary */ }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &)
  { /* no thing to do for Neumann boundary */ }



  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);


  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);



  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &);


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

#endif

#ifdef IDC_SOLVERS
#ifdef COGENDA_COMMERCIAL_PRODUCT
  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for RIC Solver-------------------//
  //////////////////////////////////////////////////////////////////////////////////


  /**
   * fill solution data into petsc vector of RIC equation.
   */
  virtual void RIC_Fill_Value(Vec x, Vec L);

  /**
   * preprocess function for RIC solver
   */
  virtual void RIC_Function_Preprocess(PetscScalar *, Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for RIC solver
   */
  virtual void RIC_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * preprocess Jacobian Matrix of RIC equation.
   */
  virtual void RIC_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for RIC solver
   */
  virtual void RIC_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for DICTAT Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////


  /**
   * fill solution data into petsc vector of DICTAT equation.
   */
  virtual void DICTAT_Fill_Value(Vec x, Vec L);

  /**
   * preprocess function for DICTAT solver
   */
  virtual void DICTAT_Function_Preprocess(PetscScalar *, Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DICTAT solver
   */
  virtual void DICTAT_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * preprocess Jacobian Matrix of DICTAT equation.
   */
  virtual void DICTAT_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DICTAT solver
   */
  virtual void DICTAT_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &);

#endif
#endif
};





#endif

