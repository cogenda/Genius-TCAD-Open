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

#ifndef __boundary_condition_charge_integral_h__
#define __boundary_condition_charge_integral_h__


#include "boundary_condition.h"


/**
 * The charge integral Boundary Condition
 */
class ChargeIntegralBC : public BoundaryCondition
{
public:

  /**
   * constructor, set default value
   */
  ChargeIntegralBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~ChargeIntegralBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
  { return ChargeIntegral; }


  /**
   * @return boundary condition type in string
   */
  virtual std::string bc_type_name() const
  { return "ChargeIntegral"; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTER_CONNECT; }

  /**
   * @return the psi of float metal
   */
  virtual PetscScalar psi() const
  { return _psi; }

  /**
   * @return writable reference to psi of float metal
   */
  virtual PetscScalar & psi()
  { return _psi; }


  /**
   * indicate that this bc is an electrode
   * @return false
   */
  virtual bool is_electrode() const
    {return false;}

  /**
   * @return true iff this boundary has a current flow
   */
  virtual bool has_current_flow() const
  { return false; }

  /**
   * @return true if this bc is the hub of inter-connect layer
   */
  virtual bool is_inter_connect_hub() const
    { return true; }

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * fermi level
   */
  PetscScalar      _psi;

private:

  /**
   * fill initial value to Charge Integral node
   */
  void charge_integral_fill_value(Vec x, Vec L);

  /**
   * set governing equation for Charge Integral node
   */
  void charge_integral_function(PetscScalar *x , Vec f, InsertMode &add_value_flag);

  /**
   * reserve matrix entries
   */
  void charge_integral_reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * set jacobian matrix entries
   */
  void charge_integral_jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag);

  /**
   * update solution data
   */
  void charge_integral_update_solution(PetscScalar *x);


public:
  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 1 DDM equation.
   */
  virtual void Poissin_Fill_Value(Vec x, Vec L)
  { charge_integral_fill_value(x, L); }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void Poissin_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { charge_integral_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { charge_integral_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { charge_integral_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of DDML1 solver.
   */
  virtual void Poissin_Update_Solution(PetscScalar *x)
  { charge_integral_update_solution(x); }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 1 DDM equation.
   */
  virtual void DDM1_Fill_Value(Vec x, Vec L)
  { charge_integral_fill_value(x, L); }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { charge_integral_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { charge_integral_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { charge_integral_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *x)
  { charge_integral_update_solution(x); }

  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 2 DDM equation.
   */
  virtual void DDM2_Fill_Value(Vec x, Vec L)
  { charge_integral_fill_value(x, L); }

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { charge_integral_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { charge_integral_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { charge_integral_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of level 2 DDM solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *x)
  { charge_integral_update_solution(x); }

  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 3 EBM equation.
   */
  virtual void EBM3_Fill_Value(Vec x, Vec L)
  { charge_integral_fill_value(x, L); }

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { charge_integral_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { charge_integral_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { charge_integral_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of EBM3 solver.
   */
  virtual void EBM3_Update_Solution(PetscScalar *x)
  { charge_integral_update_solution(x); }


#ifdef COGENDA_COMMERCIAL_PRODUCT
  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for RIC   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of RIC equation.
   */
  virtual void RIC_Fill_Value(Vec x, Vec L)
  { charge_integral_fill_value(x, L); }

  /**
   * build function and its jacobian for RIC solver, nothing to do
   */
  virtual void RIC_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { charge_integral_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void RIC_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { charge_integral_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for RIC solver, nothing to do
   */
  virtual void RIC_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { charge_integral_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of RIC solver.
   */
  virtual void RIC_Update_Solution(PetscScalar *x)
  { charge_integral_update_solution(x); }

#endif

};





#endif

