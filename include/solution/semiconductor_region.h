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

//  $Id: semiconductor_region.h,v 1.20 2008/07/09 05:58:16 gdiso Exp $

#ifndef __semiconductor_region_h__
#define __semiconductor_region_h__

#include <vector>
#include <string>


#include "genius_common.h"
#include "enum_advanced_model.h"
#include "fvm_node_info.h"
#include "material.h"
#include "simulation_region.h"

class Elem;
class GateContactBC;
class SimpleGateContactBC;
class ResistanceInsulatorBC;
class InsulatorSemiconductorInterfaceBC;



/**
 * the data and support function for semiconductor material
 */
class SemiconductorSimulationRegion :  public SimulationRegion
{
public:

  /**
   * constructor
   */
  SemiconductorSimulationRegion(const std::string &name, const std::string &material, const PetscScalar T);

  /**
   * destructor
   */
  virtual ~SemiconductorSimulationRegion()
  {
    delete mt;
  }

  /**
   * @return the region type
   */
  virtual SimulationRegionType type() const
    { return SemiconductorRegion; }

  /**
   * @return the region property in string
   */
  virtual std::string type_name() const
    { return "SemiconductorRegion";}

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
   * init data of semiconductor region,
   * call this function after doping profile done!
   */
  virtual void init(PetscScalar T_external);

  /**
   * re-init region data after import solution from data file
   */
  virtual void reinit_after_import();

  /**
   * clear stored data
   */
  virtual void clear();

  /**
   * @return the pointer to material data
   */
  Material::MaterialSemiconductor * material() const
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
   * @return the energy bandgap for optical simulation
   */
  virtual double get_optical_Eg(double T)
  { return material()->band->Eg(T); }

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
   * @return affinity of material
   */
  virtual double get_affinity(PetscScalar T) const
    { return mt->basic->Affinity(T); }

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

#if defined(HAVE_UNORDERED_MAP)
  typedef std::unordered_multimap<unsigned int, std::pair<unsigned int, SimulationRegion *> > _multimap_elem_on_interface_type;
#elif defined(HAVE_TR1_UNORDERED_MAP) || defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER)
  typedef std::tr1::unordered_multimap<unsigned int, std::pair<unsigned int, SimulationRegion *> > _multimap_elem_on_interface_type;
#else
  typedef std::multimap<unsigned int, std::pair<unsigned int, SimulationRegion *> >  _multimap_elem_on_interface_type;
#endif
  /**
   * record all the elem-side pair which on insulator interface, we need to do special process of electrical field
   * with these elements
   */
  _multimap_elem_on_interface_type  _elem_on_insulator_interface;

  /**
   * function to fill _elem_on_insulator_interface
   */
  void find_elem_on_insulator_interface();


  /**
   * record elem in the mos channel. they have a location close enough to the insulator region, and has a gate region at opposite side of insulator region
   */
  std::set<const Elem *> _elem_in_mos_channel;

  // these bc will fill _elem_in_mos_channel
  friend class GateContactBC;
  friend class SimpleGateContactBC;
  friend class ResistanceInsulatorBC;
  friend class InsulatorSemiconductorInterfaceBC;

  /**
   * normal of nearest insulator interface to elem, the distance of the elem to insulator interface is less than 0.1um
   */
  std::map<const Elem *, Point> _nearest_interface_normal;

  /**
   * function to fill _nearest_interface_normal
   */
  void find_nearest_interface_normal();

  /**
   * function to set node current to zero
   */
  void zero_node_current();

public:

  /**
   * @return true if elem on insulator interface
   */
  bool is_elem_on_insulator_interface(const Elem *elem) const;

  /**
   * get side and pointer to SimulationRegion when the neighbor of elem side is in insulator region
   */
  void elem_on_insulator_interface(const Elem *elem, std::vector<unsigned int> & sides,
                                   std::vector<SimulationRegion *> &regions) const;

  /**
   * @return true if elem in mos channel
   */
  bool is_elem_in_mos_channel(const Elem * elem) const
  { return _elem_in_mos_channel.find(elem) != _elem_in_mos_channel.end(); }

private:

  /**
   * the pointer to material database
   */
  Material::MaterialSemiconductor *mt;


private:

  bool highfield_mobility() const;


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
   * build time derivative term and its jacobian for L1 DDM
   */
  virtual void DDM1_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build time derivative term and its jacobian for L1 DDM
   */
  virtual void DDM1_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

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
  virtual void HALL_Function(const VectorValue<double> & B, PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for L1 HALL DDM
   */
  virtual void HALL_Jacobian(const VectorValue<double> & B, PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

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
   * filling AC transformation matrix (as preconditioner) entry by Jacobian matrix
   */
  virtual void DDMAC_Fill_Transformation_Matrix(Mat T, const Mat J, const double omega, InsertMode &add_value_flag) const;

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



  //////////////////////////////////////////////////////////////////////////////////
  //----------------- functions for Fast Hydrodynamic solver  --------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for fill vector of Hydrodynamic equation.
   */
  virtual void HDM_Fill_Value(Vec x, Vec vol);

  /**
   * function for evaluating flux of Hydrodynamic equation.
   */
  virtual void HDM_Flux(const PetscScalar * x, Vec flux, Vec t);

  /**
   * function for evaluating flux of Hydrodynamic equation.
   */
  virtual void HDM_Source(const PetscScalar * lx, const PetscScalar *lt, Vec x);

  /**
   * function for update solution value of Hydrodynamic equation.
   */
  virtual void HDM_Update_Solution(const PetscScalar * x);


  //////////////////////////////////////////////////////////////////////////////////
  //-----------------  functions for Linear Poissin solver   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for build matrix of linear poisson's equation.
   */
  virtual void LinearPoissin_Matrix(Mat A, InsertMode &add_value_flag);


  /**
   * function for build RHS vector of linear poisson's equation.
   */
  virtual void LinearPoissin_RHS(Vec b, InsertMode &add_value_flag);


  /**
   * function for update solution value of linear poisson's equation.
   */
  virtual void LinearPoissin_Update_Solution(const PetscScalar * x);
};



#endif //#ifndef __semiconductor_region_h__
