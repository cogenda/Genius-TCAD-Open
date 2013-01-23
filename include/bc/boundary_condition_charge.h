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

#ifndef __boundary_condition_charge_h__
#define __boundary_condition_charge_h__

#include "boundary_condition.h"


/**
 * The float metal with free charge
 */
class ChargedContactBC : public BoundaryCondition
{
public:
  /**
   * constructor
   */
  ChargedContactBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~ChargedContactBC(){}

  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return ChargedContact; }

  /**
   * @return boundary condition type in string
   */
  virtual std::string bc_type_name() const
  { return "ChargedContact"; }

  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTERFACE; }

  /**
   * @return the psi of float metal
   */
  virtual PetscScalar psi() const
    {return _psi;}

  /**
   * @return writable reference to psi of float metal
   */
  virtual PetscScalar & psi()
  {return _psi;}

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
   * @return true iff this boundary has a current flow
   */
  virtual bool has_current_flow() const
  { return true; }

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * float metal potential
   */
  PetscScalar      _psi;


  /**
   * current flow in/out ( HotCarrierInjection, FNTunneling, etc )
   */
  PetscScalar   _current_flow;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * preprocess Jacobian function for poisson solver
   */
  virtual void Poissin_Function_Preprocess(PetscScalar *, Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * preprocess Jacobian Matrix for poisson solver
   */
  virtual void Poissin_Jacobian_Preprocess(PetscScalar *, Mat *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

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
   * preprocess function for level 1 DDM solver
   */
  virtual void DDM1_Function_Preprocess(PetscScalar *, Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat * , InsertMode &);

  /**
   * preprocess Jacobian Matrix of level 1 DDM equation.
   */
  virtual void DDM1_Jacobian_Preprocess(PetscScalar *, Mat *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *);

  /**
   * function for pre process of level 1 DDM equation.
   */
  virtual void DDM1_Pre_Process();


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
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat * , InsertMode &);

  /**
   * preprocess Jacobian Matrix of level 2 DDM equation.
   */
  virtual void DDM2_Jacobian_Preprocess(PetscScalar *,Mat *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of DDML2 solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *);

  /**
   * function for pre process of level 2 DDM equation.
   */
  virtual void DDM2_Pre_Process();

  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for EBM   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void EBM3_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for EBM3 solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for EBM3 solver
   */
  virtual void EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of level 3 EBM solver.
   */
  virtual void EBM3_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );


#ifdef COGENDA_COMMERCIAL_PRODUCT
  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for RIC Solver-------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * preprocess function for RIC solver
   */
  virtual void RIC_Function_Preprocess(PetscScalar *, Vec, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for RIC solver
   */
  virtual void RIC_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void RIC_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * preprocess Jacobian Matrix of RIC equation.
   */
  virtual void RIC_Jacobian_Preprocess(PetscScalar *, Mat *, std::vector<PetscInt>&, std::vector<PetscInt>&, std::vector<PetscInt>&);

  /**
   * build function and its jacobian for RIC solver
   */
  virtual void RIC_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of RIC solver.
   */
  virtual void RIC_Update_Solution(PetscScalar *);

  /**
   * function for pre process of RIC equation.
   */
  virtual void RIC_Pre_Process();

#endif

}
;





#endif

