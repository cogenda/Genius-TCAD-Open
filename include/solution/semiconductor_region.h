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
  SemiconductorSimulationRegion(const std::string &name, const std::string &material, const double T, const double z);

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
   * @return the quality if each fvm cell
   */
  virtual Real fvm_cell_quality() const;

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
   * virtual function for set different model, calibrate parameters to PMI
   */
  virtual void set_pmi(const std::string &type, const std::string &model_name, std::vector<Parser::Parameter> & pmi_parameters);

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

#if defined(HAVE_UNORDERED_SET)
  typedef std::unordered_set<const Elem *> _elem_in_mos_channel_type;
#elif defined(HAVE_TR1_UNORDERED_SET) || defined(HAVE_TR1_UNORDERED_SET_WITH_STD_HEADER)
  typedef std::tr1::unordered_set<const Elem *> _elem_in_mos_channel_type;
#else
  typedef std::set<const Elem *>  _elem_in_mos_channel_type;
#endif

  /**
   * record elem in the mos channel. they have a location close enough to the insulator region, and has a gate region at opposite side of insulator region
   */
  _elem_in_mos_channel_type _elem_in_mos_channel;

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
   * record element has side/edge/node on boundary
   */
  std::set<const Elem *> _elem_touch_boundary;

  /**
   * function to fill _elem_touch_boundary
   */
  void find_elem_touch_boundary();

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

  /**
   * @return true if elem touch the boundary (has side, edge or node on boundary)
   */
  bool is_elem_touch_boundary(const Elem * elem) const
  { return _elem_touch_boundary.find(elem) != _elem_touch_boundary.end(); }


private:

  /**
   * the pointer to material database
   */
  Material::MaterialSemiconductor *mt;


private:

  /**
   * @return true when highfield mobility enabled
   */
  bool highfield_mobility() const;

  /**
   * @return truncated partial area associated with edge ne of elem
   */
  Real truncated_partial_area(const Elem * elem, unsigned int ne) const;

  /**
   * elem has its circumcircle center outside the region
   */
  //void bad_elem_info() const;

private:

  void DDM1_Gummel_Carrier_Electron(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);

  void DDM1_Gummel_Carrier_Hole(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);

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
   * function for evaluating pseudo time step of level 1 DDM equation.
   */
  virtual void DDM1_Pseudo_Time_Step_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * function for evaluating Jacobian of pseudo time step of level 1 DDM equation.
   */
  virtual void DDM1_Pseudo_Time_Step_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * function for convergence test of pseudo time step of level 1 DDM equation.
   */
  virtual int DDM1_Pseudo_Time_Step_Convergence_Test(PetscScalar * x);

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
  virtual void DG_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * function for evaluating time derivative term of Density Gradient equation.
   */
  virtual void DG_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * function for evaluating Jacobian of time derivative term of Density Gradient equation.
   */
  virtual void DG_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

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
   * function for evaluating pseudo time step for L2 DDM
   */
  virtual void DDM2_Pseudo_Time_Step_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * function for evaluating Jacobian of pseudo time step for L2 DDM
   */
  virtual void DDM2_Pseudo_Time_Step_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * function for convergence test of pseudo time step for L2 DDM
   */
  virtual int DDM2_Pseudo_Time_Step_Convergence_Test(PetscScalar * x);

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

#ifdef COGENDA_COMMERCIAL_PRODUCT
  //////////////////////////////////////////////////////////////////////////////////
  //----------------- functions for Gummel DDML1 solver --------------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for build RHS and matrix for gummel equation.
   */
  virtual void DDM1_Gummel_Carrier(const std::string & carrier, PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);

  /**
   * function for build RHS and matrix for gummel carrier equation.
   */
  virtual void DDM1_Implicit_Gummel_Carrier_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag);

  /**
   * function for evaluating Jacobian for gummel carrier equation.
   */
  virtual void DDM1_Implicit_Gummel_Carrier_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag);

  /**
   * function for build RHS and matrix for half implicit current continuity equation.
   */
  virtual void DDM1_Half_Implicit_Current(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);

  /**
   * function for build RHS and matrix for half implicit poisson correction equation.
   */
  virtual void DDM1_Half_Implicit_Poisson_Correction(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);

  /**
   * function for build RHS and matrix for half implicit poisson correction with Polsky's method.
   */
  virtual void DDM1_Half_Implicit_Poisson_Correction_Polsky(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag);

#endif

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


  //////////////////////////////////////////////////////////////////////////////////
  //-----------------  functions for mobility evaluation     ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * function for evaluating mobility of edge.
   *
   * the moblity is ordered as region's edge, and weighted by partial area of the edge
   *
   */
  virtual void Mob_Evaluation( std::vector< std::pair<unsigned int, unsigned int> > &edge,
                               std::vector< std::pair<double, double> > & mob,
                               std::vector< double > & weight) const;

#ifdef COGENDA_COMMERCIAL_PRODUCT
  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for RIC   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * filling solution data from FVM_NodeData into petsc vector of RIC equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   */
  virtual void RIC_Fill_Value(Vec , Vec ) {}

  /**
   * function for evaluating RIC equation.
   */
  virtual void RIC_Function(PetscScalar * , Vec , InsertMode &) {}

  /**
   * function for evaluating Jacobian of RIC equation.
   */
  virtual void RIC_Jacobian(PetscScalar *, Mat *, InsertMode &) {}

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
  virtual void RIC_Update_Solution(PetscScalar *) {}

#endif
};



#endif //#ifndef __semiconductor_region_h__
