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
#include "fvm_node_data_conductor.h"
#include "fvm_cell_data_conductor.h"
#include "material.h"
#include "simulation_region.h"

class Elem;

/**
 * the data and support function for conductor material
 */
class ConductorSimulationRegion :  public SimulationRegion
{
public:

  ConductorSimulationRegion(const std::string &name, const std::string &material, SimulationSystem & system);

  virtual ~ConductorSimulationRegion()
  { delete mt; }

  virtual SimulationRegionType type() const
  { return ConductorRegion; }

  /**
   * insert local mesh element into the region, only copy the pointer
   * and create cell data
   */
  virtual void insert_cell (const Elem * e);

  /**
   * @note only node belongs to current processor and ghost node
   * own FVM_NodeData
   */
  virtual void insert_fvm_node(FVM_Node * fn)
  {

    // node (or ghost node) belongs to this processor
    // we should build FVM_NodeData structure for it

    if  ( fn->root_node()->on_local() )
      fn->hold_node_data( new FVM_Conductor_NodeData );

    _region_node[fn->root_node()->id()] = fn;
  }

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
  virtual Material::MaterialBase * get_material_base()
  { return (Material::MaterialBase *)mt; }

  /**
   * @return the optical refraction index of the region
   */
  virtual Complex get_optical_refraction(double lamda)
  { return material()->optical->RefractionIndex(lamda); }

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

private:
  /**
   * the pointer to material database
   */
  Material::MaterialConductor *mt;


public:

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
  virtual void Poissin_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for poisson solver (experimental)
   */
  virtual void Poissin_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for poisson solver (experimental)
   */
  virtual void Poissin_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

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
  virtual void DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

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
  virtual void DDM1_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for L1 DDM (experimental)
   */
  virtual void DDM1_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

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
  virtual void HALL_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for L1 HALL DDM (experimental)
   */
  virtual void HALL_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for L1 HALL DDM (experimental)
   */
  virtual void HALL_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * update solution data of FVM_NodeData by petsc vector of L1 HALL DDM
   */
  virtual void HALL_Update_Solution(PetscScalar *lxx);



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
  virtual void DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L2 DDM
   */
  virtual void DDM2_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L2 DDM
   */
  virtual void DDM2_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for L2 DDM (experimental)
   */
  virtual void DDM2_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for L2 DDM (experimental)
   */
  virtual void DDM2_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

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
  virtual void EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for L3 EBM (experimental)
   */
  virtual void EBM3_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * process function and its jacobian at hanging node for L3 EBM (experimental)
   */
  virtual void EBM3_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L3 EBM
   */
  virtual void EBM3_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L3 EBM
   */
  virtual void EBM3_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

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
  virtual void DDMAC_Fill_Matrix_Vector(Mat A,  Vec b, const Mat J, const double omega, InsertMode &add_value_flag) const;

  /**
   * fill matrix of DDMAC equation for fvm_node
   */
  virtual void DDMAC_Fill_Nodal_Matrix_Vector(const FVM_Node *fvm_node, Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag,
                                              const SimulationRegion * adjacent_region=NULL, const FVM_Node * adjacent_fvm_node=NULL) const;

  /**
   * fill matrix of DDMAC equation for variable of fvm_node
   */
  virtual void DDMAC_Fill_Nodal_Matrix_Vector(const FVM_Node *fvm_node, const SolutionVariable var,
                                              Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag,
                                              const SimulationRegion * adjacent_region=NULL, const FVM_Node * adjacent_fvm_node=NULL) const;

  /**
   * update solution value of DDMAC equation
   */
  virtual void DDMAC_Update_Solution(PetscScalar *lxx);

};


#endif //#ifndef __conductor_region_h__
