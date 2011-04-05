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

#ifndef __boundary_condition_electrode_interconnect_h__
#define __boundary_condition_electrode_interconnect_h__


#include "boundary_condition.h"


/**
 * The inter connect Boundary Condition
 */
class ElectrodeInterConnectBC : public BoundaryCondition
{
public:

  /**
   * constructor, set default value
   */
  ElectrodeInterConnectBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~ElectrodeInterConnectBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return InterConnect; }

  /**
   * @return boundary condition type in string
   */
  virtual std::string bc_type_name() const
  { return "InterConnect"; }

  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTER_CONNECT; }

  /**
   * indicate that this bc is an electrode
   * @return true
   */
  virtual bool is_electrode() const
  {return ext_circuit()!=NULL;}

  /**
   * @return true iff this boundary has a current flow
   */
  virtual bool has_current_flow() const
  { return true; }

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
   * fill initial value to inter-connect node
   */
  void inter_connect_fill_value(Vec x, Vec L);

  /**
   * set governing equation for inter-connect node. use nodal analysis method
   */
  void inter_connect_function(PetscScalar *x , Vec f, InsertMode &add_value_flag);

  /**
   * reserve matrix entries
   */
  void inter_connect_reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * set jacobian matrix entries
   */
  void inter_connect_jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag);

  /**
   * update solution data
   */
  void inter_connect_update_solution(PetscScalar *x);


public:
  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &)
  { /* no thing to do for poisson solver */ }

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &)
  { /* no thing to do for poisson solver */ }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 1 DDM equation.
   */
  virtual void DDM1_Fill_Value(Vec x, Vec L)
  { inter_connect_fill_value(x, L); }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { inter_connect_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { inter_connect_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { inter_connect_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *x)
  { inter_connect_update_solution(x); }

  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 2 DDM equation.
   */
  virtual void DDM2_Fill_Value(Vec x, Vec L)
  { inter_connect_fill_value(x, L); }

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { inter_connect_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { inter_connect_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { inter_connect_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of level 2 DDM solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *x)
  { inter_connect_update_solution(x); }

  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * fill solution data into petsc vector of level 3 EBM equation.
   */
  virtual void EBM3_Fill_Value(Vec x, Vec L)
  { inter_connect_fill_value(x, L); }

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { inter_connect_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { inter_connect_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { inter_connect_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of EBM3 solver.
   */
  virtual void EBM3_Update_Solution(PetscScalar *x)
  { inter_connect_update_solution(x); }

};





#endif

