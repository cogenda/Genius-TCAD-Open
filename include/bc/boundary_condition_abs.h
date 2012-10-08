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

#ifndef __boundary_condition_abs_h__
#define __boundary_condition_abs_h__


#include "boundary_condition.h"



/**
 * The Neumann Boundary Condition
 */
class AbsorbingBC : public BoundaryCondition
{
public:

  /**
   * constructor, set default value
   */
  AbsorbingBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~AbsorbingBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return AbsorbingBoundary; }

  /**
   * @return boundary condition type in string
   */
  virtual std::string bc_type_name() const
  { return "AbsorbingBoundary"; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return BOUNDARY; }


  /**
   * @return true iff this boundary has a current flow
   */
  virtual bool has_current_flow() const
  { return false; }

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &)
  { /* no thing to do for AbsorbingBC boundary */ }

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &)
  { /* no thing to do for AbsorbingBC boundary */ }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &)
  { /* no thing to do for AbsorbingBC boundary */ }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &)
  { /* no thing to do for AbsorbingBC boundary */ }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &)
  { /* no thing to do for AbsorbingBC boundary */ }

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &)
  { /* no thing to do for AbsorbingBC boundary */ }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &)
  { /* no thing to do for AbsorbingBC boundary */ }

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &)
  { /* no thing to do for AbsorbingBC boundary */ }


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat , Vec , const Mat , const double , InsertMode & )
  { /* no thing to do for AbsorbingBC boundary */ }


};






#endif

