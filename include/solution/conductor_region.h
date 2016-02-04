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

//  $Id: conductor_region.h,v 1.15 2008/07/09 05:58:16 gdiso Exp $

#ifndef __conductor_region_h__
#define __conductor_region_h__

#include <vector>
#include <string>
#include "fvm_node_info.h"
#include "material.h"
#include "simulation_region.h"

class Elem;

/**
 * the data and support function for conductor material
 */
class ElectrodeSimulationRegion :  public SimulationRegion
{
public:

  ElectrodeSimulationRegion(const std::string &name, const std::string &material, const double T, const unsigned int dim, const double z);

  virtual ~ElectrodeSimulationRegion()
  { delete mt; }

  /**
   * @return the region property
   */
  virtual SimulationRegionType type() const
  { return ElectrodeRegion; }

  /**
   * @return the region property in string
   */
  virtual std::string type_name() const
  { return "ElectrodeRegion";}


  /**
   * insert local mesh element into the region, only copy the pointer
   * and create cell data
   */
  virtual void insert_cell (const Elem * e);

  /**
   * @note only node belongs to current processor and ghost node
   * own FVM_NodeData
   */
  virtual void insert_fvm_node(FVM_Node * fn);

  /**
   * init node data for this region
   */
  virtual void init(PetscScalar T_external);

  /**
   * re-init region data after import solution from data file
   */
  virtual void reinit_after_import();

  /**
   * @return the pointer to material data
   */
  Material::MaterialConductor * material() const
  {return mt;}

  /**
   * @return the base class of material database
   */
  virtual Material::MaterialBase * get_material_base() const
  { return (Material::MaterialBase *)mt; }

  /**
   * @return the optical refraction index of the region
   */
  virtual Complex get_optical_refraction(double lamda) const
  { 
    std::complex<PetscScalar> r = material()->optical->RefractionIndex(lamda, T_external()); 
    return Complex(r.real(), r.imag());
  }

  /**
   * @return relative permittivity of material
   */
  virtual double get_eps() const
  { return mt->basic->Permittivity(); }

  /**
   * @return maretial density [g cm^-3]
   */
  virtual double get_density() const
  { return mt->basic->Density(T_external()); }

  /**
   * @return affinity of material
   */
  virtual double get_affinity() const
  { return mt->basic->Affinity(T_external()); }

   /**
   * @return material conductance [A/V/cm]
    */
  virtual double get_conductance() const
  { return mt->basic->Conductance(); }

  /**
   * virtual function for set different model, calibrate parameters to PMI
   */
  virtual void set_pmi(const std::string &type, const std::string &model_name, std::vector<Parser::Parameter> & pmi_parameters);

  /**
   * get atom fraction of region material
   */
  virtual void atom_fraction(std::vector<Atom> &atoms, std::vector<double> & fraction) const
  {
    atoms.clear();
    fraction.clear();
    mt->basic->G4Material(atoms, fraction);
  }

  /**
   * set the variables for this region
   */
  virtual void set_region_variables();

private:
  /**
   * the pointer to material database
   */
  Material::MaterialConductor *mt;


public:

#ifdef TCAD_SOLVERS

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of poisson's equation.
   */
  virtual void Poissin_Fill_Value(Vec x, Vec L);

  /**
   * build function and its jacobian for poisson solver
   */
  virtual void Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for poisson solver
   */
  virtual void Poissin_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for poisson solver (experimental)
   */
  virtual void Poissin_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for poisson solver (experimental)
   */
  virtual void Poissin_Jacobian_Hanging_Node(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * update solution data of FVM_NodeData by petsc vector of poisson's equation.
   */
  virtual void Poissin_Update_Solution(PetscScalar *lxx);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of L1 DDM.
   */
  virtual void DDM1_Fill_Value(Vec x, Vec L);

  /**
   * build function and its jacobian for L1 DDM
   */
  virtual void DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for L1 DDM
   */
  virtual void DDM1_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L1 DDM, do nothing here
   */
  virtual void DDM1_Time_Dependent_Function(PetscScalar * , Vec , InsertMode &) {}

  /**
   * build time derivative term and its jacobian for L1 DDM, do nothing here
   */
  virtual void DDM1_Time_Dependent_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &) {}

  /**
   * update solution data of FVM_NodeData by petsc vector of L1 DDM.
   */
  virtual void DDM1_Update_Solution(PetscScalar *lxx);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for L1 HALL DDM------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of L1 HALL DDM.
   */
  virtual void HALL_Fill_Value(Vec x, Vec L);

  /**
   * build function and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Function(const VectorValue<PetscScalar> & B, PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Jacobian(const VectorValue<PetscScalar> & B, PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Time_Dependent_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * update solution data of FVM_NodeData by petsc vector of L1 HALL DDM
   */
  virtual void HALL_Update_Solution(PetscScalar *lxx);


  //////////////////////////////////////////////////////////////////////////////////
  //-------------Function and Jacobian evaluate for Density Gradient--------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for get nodal variable number
   */
  virtual unsigned int dg_n_variables() const;

  /**
   * function for get offset of nodal variable
   */
  virtual unsigned int dg_variable_offset(SolutionVariable var) const;

  /**
   * function for fill vector of Density Gradient equation.
   */
  virtual void DG_Fill_Value(Vec x, Vec L);

  /**
   * function for evaluating Density Gradient equation.
   */
  virtual void DG_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * function for evaluating Jacobian of Density Gradient equation.
   */
  virtual void DG_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * function for evaluating time derivative term of Density Gradient equation.
   */
  virtual void DG_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * function for evaluating Jacobian of time derivative term of Density Gradient equation.
   */
  virtual void DG_Time_Dependent_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * function for update solution value of Density Gradient equation.
   */
  virtual void DG_Update_Solution(PetscScalar *lxx);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of L2 DDM.
   */
  virtual void DDM2_Fill_Value(Vec x, Vec L);

  /**
   * build function and its jacobian for L2 DDM
   */
  virtual void DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for L2 DDM
   */
  virtual void DDM2_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L2 DDM
   */
  virtual void DDM2_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L2 DDM
   */
  virtual void DDM2_Time_Dependent_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * update solution data of FVM_NodeData by petsc vector of L2 DDM.
   */
  virtual void DDM2_Update_Solution(PetscScalar *lxx);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * get nodal variable number, which depedent on EBM Level
   */
  virtual unsigned int ebm_n_variables() const;

  /**
   * get offset of nodal variable, which depedent on EBM Level
   */
  virtual unsigned int ebm_variable_offset(SolutionVariable var) const;

  /**
   * filling solution data from FVM_NodeData into petsc vector of L3 EBM.
   */
  virtual void EBM3_Fill_Value(Vec x, Vec L);

  /**
   * build function and its jacobian for L3 EBM
   */
  virtual void EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for L3 EBM
   */
  virtual void EBM3_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L3 EBM
   */
  virtual void EBM3_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L3 EBM
   */
  virtual void EBM3_Time_Dependent_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag);

  /**
   * update solution data of FVM_NodeData by petsc vector of L3 EBM.
   */
  virtual void EBM3_Update_Solution(PetscScalar *lxx);


  //////////////////////////////////////////////////////////////////////////////////
  //-----------------   Vec and Matrix evaluate for EBM AC   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of DDMAC.
   */
  virtual void DDMAC_Fill_Value(Vec x, Vec L) const;

  /**
   * fill matrix of DDMAC equation
   */
  virtual void DDMAC_Fill_Matrix_Vector(Mat A,  Vec b, const Mat J, const PetscScalar omega, InsertMode &add_value_flag) const;

  /**
   * filling AC transformation matrix (as preconditioner) entry by Jacobian matrix
   */
  virtual void DDMAC_Fill_Transformation_Matrix(Mat T, const Mat J, const PetscScalar omega, InsertMode &add_value_flag) const;

  /**
   * fill matrix of DDMAC equation for fvm_node
   */
  virtual void DDMAC_Fill_Nodal_Matrix_Vector(const FVM_Node *fvm_node, Mat A, Vec b, const Mat J, const PetscScalar omega, InsertMode & add_value_flag,
                                              const SimulationRegion * adjacent_region=NULL, const FVM_Node * adjacent_fvm_node=NULL) const;

  /**
   * fill matrix of DDMAC equation for variable of fvm_node
   */
  virtual void DDMAC_Fill_Nodal_Matrix_Vector(const FVM_Node *fvm_node, const SolutionVariable var,
                                              Mat A, Vec b, const Mat J, const PetscScalar omega, InsertMode & add_value_flag,
                                              const SimulationRegion * adjacent_region=NULL, const FVM_Node * adjacent_fvm_node=NULL) const;

  /**
   * filling AC matrix entry by force variable of FVM_Node1 equals to FVM_Node2
   */
  virtual void DDMAC_Force_equal(const FVM_Node *fvm_node, Mat A, InsertMode & add_value_flag,
                                 const SimulationRegion * adjacent_region=NULL,
                                 const FVM_Node * adjacent_fvm_node=NULL) const;

  /**
   * filling AC matrix entry by force given variable of FVM_Node1 equals to FVM_Node2
   */
  virtual void DDMAC_Force_equal(const FVM_Node *fvm_node, const SolutionVariable var,
                                 Mat A, InsertMode & add_value_flag,
                                 const SimulationRegion * adjacent_region=NULL,
                                 const FVM_Node * adjacent_fvm_node=NULL) const;

  /**
   * update solution value of DDMAC equation
   */
  virtual void DDMAC_Update_Solution(PetscScalar *lxx);


#ifdef COGENDA_COMMERCIAL_PRODUCT
  //////////////////////////////////////////////////////////////////////////////////
  //----------------- functions for Gummel DDML1 solver --------------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for build RHS and matrix for half implicit current continuity equation.
   */
  virtual void DDM1_Half_Implicit_Current(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);

  /**
   * function for build RHS and matrix for half implicit current poisson correction equation.
   */
  virtual void DDM1_Half_Implicit_Poisson_Correction(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);

#endif



#endif


#ifdef IDC_SOLVERS
#ifdef COGENDA_COMMERCIAL_PRODUCT
  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for RIC   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of RIC equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   */
  virtual void RIC_Fill_Value(Vec , Vec );

  /**
   * function for evaluating RIC equation.
   */
  virtual void RIC_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * function for evaluating Jacobian of RIC equation.
   */
  virtual void RIC_Jacobian(PetscScalar *, SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * function for evaluating time derivative term of RIC equation.
   */
  virtual void RIC_Time_Dependent_Function(PetscScalar *, Vec , InsertMode &) {}

  /**
   * function for evaluating Jacobian of time derivative term of RIC equation.
   */
  virtual void RIC_Time_Dependent_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &) {}

  /**
   * function for update solution value of RIC equation.
   */
  virtual void RIC_Update_Solution(PetscScalar *);

  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for DICTAT---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of DICTAT equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   */
  virtual void DICTAT_Fill_Value(Vec , Vec );

  /**
   * function for evaluating DICTAT equation.
   */
  virtual void DICTAT_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * function for evaluating Jacobian of DICTAT equation.
   */
  virtual void DICTAT_Jacobian(PetscScalar *, SparseMatrix<PetscScalar> *, InsertMode &);

  /**
   * function for evaluating time derivative term of DICTAT equation.
   */
  virtual void DICTAT_Time_Dependent_Function(PetscScalar *, Vec , InsertMode &){}

  /**
   * function for evaluating Jacobian of time derivative term of DICTAT equation.
   */
  virtual void DICTAT_Time_Dependent_Jacobian(PetscScalar * , SparseMatrix<PetscScalar> *, InsertMode &){}

  /**
   * function for update solution value of DICTAT equation.
   */
  virtual void DICTAT_Update_Solution(PetscScalar *);

#endif
#endif

};


#endif //#ifndef __conductor_region_h__
