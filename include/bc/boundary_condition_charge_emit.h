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

#ifndef __boundary_condition_charge_emit_h__
#define __boundary_condition_charge_emit_h__


#include "boundary_condition.h"

/**
 * The Charge Emit Boundary Condition
 */
class ChargeEmitBC : public BoundaryCondition
{
public:

  /**
   * constructor, set default value
   */
  ChargeEmitBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~ChargeEmitBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
  { return ChargeEmit; }

  /**
   * @return boundary condition type in string
   */
  virtual std::string bc_type_name() const
  { return "ChargeEmit"; }

  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return BOUNDARY; }

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
  { return true; }


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

    // current buffer
    std::vector<PetscScalar> _current_buffer;
    std::vector<PetscScalar> _electron_current_buffer;
    std::vector<PetscScalar> _hole_current_buffer;

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
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &)
  { /* nothing to do for Neumann boundary */ }


  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &)
  { /* nothing to do for Neumann boundary */ }



  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &)
  { /* nothing to do for Neumann boundary */ }


  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &)
  { /* nothing to do for Neumann boundary */ }


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const PetscScalar omega, InsertMode & add_value_flag )
  { /* nothing to do for Neumann boundary */ }

#endif

};





#endif

