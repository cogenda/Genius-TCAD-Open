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


#ifndef __vacuum_region_h__
#define __vacuum_region_h__

#include <vector>
#include <string>
#include "fvm_node_info.h"
#include "material.h"
#include "simulation_region.h"

class Elem;

/**
 * the data and support function for vacuum
 */
class VacuumSimulationRegion :  public SimulationRegion
{
public:

  /**
   * constructor
   */
  VacuumSimulationRegion(const std::string &name, const std::string &material, const double T, const double z);

  /**
   * destructor
   */
  virtual ~VacuumSimulationRegion()
  { delete mt; }

  /**
   * @return the region type
   */
  virtual SimulationRegionType type() const
  { return VacuumRegion; }

  /**
   * @return the region property in string
   */
  virtual std::string type_name() const
  { return "VacuumRegion";}

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
  Material::MaterialVacuum * material() const
  {return mt;}


  /**
   * @return the base class of material database
   */
  virtual Material::MaterialBase * get_material_base() const
  { return (Material::MaterialBase *)mt; }


  /**
   * @return the optical refraction index of the region
   */
  virtual Complex get_optical_refraction(double lamda)
  { return material()->optical->RefractionIndex(lamda, T_external()); }

  /**
   * @return relative permittivity of material
   */
  virtual double get_eps() const
  { return mt->basic->Permittivity(); }

  /**
   * @return maretial density [g cm^-3]
   */
  virtual double get_density(PetscScalar T) const
  { return mt->basic->Density(T); }

  /**
   * get atom fraction of region material
   */
  virtual void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction) const
  {
    atoms.clear();
    fraction.clear();
    mt->basic->atom_fraction(atoms, fraction);
  }

  /**
   * set the variables for this region
   */
  virtual void set_region_variables();

private:
  /**
   * the pointer to material database
   */
  Material::MaterialVacuum *mt;


public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of poisson's equation.
   */
  virtual void Poissin_Fill_Value(Vec , Vec ) {}

  /**
   * build function and its jacobian for poisson solver
   */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &) {}

  /**
   * build function and its jacobian for poisson solver
   */
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * process function and its jacobian at hanging node for poisson solver (experimental)
   */
  virtual void Poissin_Function_Hanging_Node(PetscScalar * , Vec , InsertMode &) {}

  /**
   * process function and its jacobian at hanging node for poisson solver (experimental)
   */
  virtual void Poissin_Jacobian_Hanging_Node(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * update solution data of FVM_NodeData by petsc vector of poisson's equation.
   */
  virtual void Poissin_Update_Solution(PetscScalar *) {}

  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of L1 DDM.
   */
  virtual void DDM1_Fill_Value(Vec , Vec ) {}

  /**
   * build function and its jacobian for L1 DDM
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &) {}

  /**
   * build function and its jacobian for L1 DDM
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * build time derivative term and its jacobian for L1 DDM, do nothing here
   */
  virtual void DDM1_Time_Dependent_Function(PetscScalar * , Vec , InsertMode &) {}

  /**
   * build time derivative term and its jacobian for L1 DDM, do nothing here
   */
  virtual void DDM1_Time_Dependent_Jacobian(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * process function and its jacobian at hanging node for L1 DDM (experimental)
   */
  virtual void DDM1_Function_Hanging_Node(PetscScalar * , Vec , InsertMode &) {}

  /**
   * process function and its jacobian at hanging node for L1 DDM (experimental)
   */
  virtual void DDM1_Jacobian_Hanging_Node(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * update solution data of FVM_NodeData by petsc vector of L1 DDM.
   */
  virtual void DDM1_Update_Solution(PetscScalar *) {}


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for L1 HALL DDM------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of L1 HALL DDM.
   */
  virtual void HALL_Fill_Value(Vec , Vec ) {}

  /**
   * build function and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Function(const VectorValue<double> & , PetscScalar * , Vec , InsertMode &) {}

  /**
   * build function and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Jacobian(const VectorValue<double> & , PetscScalar * , Mat *, InsertMode &) {}

  /**
   * process function and its jacobian at hanging node for L1 HALL DDM (experimental)
   */
  virtual void HALL_Function_Hanging_Node(PetscScalar * , Vec , InsertMode &) {}

  /**
   * process function and its jacobian at hanging node for L1 HALL DDM (experimental)
   */
  virtual void HALL_Jacobian_Hanging_Node(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * build time derivative term and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Time_Dependent_Function(PetscScalar * , Vec , InsertMode &) {}

  /**
   * build time derivative term and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Time_Dependent_Jacobian(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * update solution data of FVM_NodeData by petsc vector of L1 HALL DDM
   */
  virtual void HALL_Update_Solution(PetscScalar *) {}



  //////////////////////////////////////////////////////////////////////////////////
  //-------------Function and Jacobian evaluate for Density Gradient--------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for fill vector of Density Gradient equation.
   */
  virtual void DG_Fill_Value(Vec x, Vec L){}

  /**
   * function for evaluating Density Gradient equation.
   */
  virtual void DG_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag){}

  /**
   * function for evaluating Jacobian of Density Gradient equation.
   */
  virtual void DG_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag){}

  /**
   * function for evaluating time derivative term of Density Gradient equation.
   */
  virtual void DG_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag){}

  /**
   * function for evaluating Jacobian of time derivative term of Density Gradient equation.
   */
  virtual void DG_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag){}

  /**
   * function for update solution value of Density Gradient equation.
   */
  virtual void DG_Update_Solution(PetscScalar *lxx){}


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for get nodal variable number
   */
  virtual unsigned int dg_n_variables() const { return 0; }

  /**
   * function for get offset of nodal variable
   */
  virtual unsigned int dg_variable_offset(SolutionVariable var) const {return 0;}


  /**
   * filling solution data from FVM_NodeData into petsc vector of L2 DDM.
   */
  virtual void DDM2_Fill_Value(Vec , Vec ) {}

  /**
   * build function and its jacobian for L2 DDM
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &) {}

  /**
   * build function and its jacobian for L2 DDM
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * build time derivative term and its jacobian for L2 DDM
   */
  virtual void DDM2_Time_Dependent_Function(PetscScalar * , Vec , InsertMode &) {}

  /**
   * build time derivative term and its jacobian for L2 DDM
   */
  virtual void DDM2_Time_Dependent_Jacobian(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * process function and its jacobian at hanging node for L2 DDM (experimental)
   */
  virtual void DDM2_Function_Hanging_Node(PetscScalar * , Vec , InsertMode &) {}

  /**
   * process function and its jacobian at hanging node for L2 DDM (experimental)
   */
  virtual void DDM2_Jacobian_Hanging_Node(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * update solution data of FVM_NodeData by petsc vector of L2 DDM.
   */
  virtual void DDM2_Update_Solution(PetscScalar *) {}


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * get nodal variable number, which depedent on EBM Level
   */
  virtual unsigned int ebm_n_variables() const { return 0; }

  /**
   * get offset of nodal variable, which depedent on EBM Level
   */
  virtual unsigned int ebm_variable_offset(SolutionVariable ) const { return 0; }

  /**
   * filling solution data from FVM_NodeData into petsc vector of L3 EBM.
   */
  virtual void EBM3_Fill_Value(Vec , Vec ) {}

  /**
   * build function and its jacobian for L3 EBM
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &) {}

  /**
   * build function and its jacobian for L3 EBM
   */
  virtual void EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * process function and its jacobian at hanging node for L3 EBM (experimental)
   */
  virtual void EBM3_Function_Hanging_Node(PetscScalar * , Vec , InsertMode &) {}

  /**
   * process function and its jacobian at hanging node for L3 EBM (experimental)
   */
  virtual void EBM3_Jacobian_Hanging_Node(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * build time derivative term and its jacobian for L3 EBM
   */
  virtual void EBM3_Time_Dependent_Function(PetscScalar * , Vec , InsertMode &) {}

  /**
   * build time derivative term and its jacobian for L3 EBM
   */
  virtual void EBM3_Time_Dependent_Jacobian(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * update solution data of FVM_NodeData by petsc vector of L3 EBM.
   */
  virtual void EBM3_Update_Solution(PetscScalar *) {}


  //////////////////////////////////////////////////////////////////////////////////
  //-----------------   Vec and Matrix evaluate for EBM AC   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of DDMAC.
   */
  virtual void DDMAC_Fill_Value(Vec , Vec ) const {}

  /**
   * fill matrix of DDMAC equation
   */
  virtual void DDMAC_Fill_Matrix_Vector(Mat ,  Vec , const Mat , const double , InsertMode &) const {}

  /**
   * filling AC transformation matrix (as preconditioner) entry by Jacobian matrix
   */
  virtual void DDMAC_Fill_Transformation_Matrix(Mat , const Mat , const double , InsertMode &) const {}

  /**
   * fill matrix of DDMAC equation for fvm_node
   */
  virtual void DDMAC_Fill_Nodal_Matrix_Vector(const FVM_Node *, Mat , Vec , const Mat , const double , InsertMode & ,
                                              const SimulationRegion * =NULL, const FVM_Node * =NULL) const {}

  /**
   * fill matrix of DDMAC equation for variable of fvm_node
   */
  virtual void DDMAC_Fill_Nodal_Matrix_Vector(const FVM_Node *, const SolutionVariable ,
                                              Mat , Vec , const Mat , const double , InsertMode & ,
                                              const SimulationRegion * =NULL, const FVM_Node * =NULL) const {}
  /**
   * filling AC matrix entry by force variable of FVM_Node1 equals to FVM_Node2
   */
  virtual void DDMAC_Force_equal(const FVM_Node *fvm_node, Mat A, InsertMode & add_value_flag,
                                 const SimulationRegion * adjacent_region=NULL,
                                 const FVM_Node * adjacent_fvm_node=NULL) const {}

  /**
   * filling AC matrix entry by force given variable of FVM_Node1 equals to FVM_Node2
   */
  virtual void DDMAC_Force_equal(const FVM_Node *fvm_node, const SolutionVariable var,
                                 Mat A, InsertMode & add_value_flag,
                                 const SimulationRegion * adjacent_region=NULL,
                                 const FVM_Node * adjacent_fvm_node=NULL) const {}

  /**
   * update solution value of DDMAC equation
   */
  virtual void DDMAC_Update_Solution(PetscScalar *) {}


  //////////////////////////////////////////////////////////////////////////////////
  //-----------------  functions for Linear Poissin solver   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for build matrix of linear poisson's equation.
   */
  virtual void LinearPoissin_Matrix(Mat, InsertMode &) {}

  /**
   * function for build RHS vector of linear poisson's equation.
   */
  virtual void LinearPoissin_RHS(Vec, InsertMode & ) {}


  /**
   * function for update solution value of linear poisson's equation.
   */
  virtual void LinearPoissin_Update_Solution(const PetscScalar * ) {}

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
  virtual void RIC_Jacobian(PetscScalar *, Mat *, InsertMode &);

  /**
   * function for evaluating time derivative term of RIC equation.
   */
  virtual void RIC_Time_Dependent_Function(PetscScalar *, Vec , InsertMode &) {}

  /**
   * function for evaluating Jacobian of time derivative term of RIC equation.
   */
  virtual void RIC_Time_Dependent_Jacobian(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * function for update solution value of RIC equation.
   */
  virtual void RIC_Update_Solution(PetscScalar *);

#endif

};




#endif //#ifndef __vacuum_region_h__
