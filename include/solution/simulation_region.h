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

//  $Id: simulation_region.h,v 1.45 2008/07/09 05:58:16 gdiso Exp $

#ifndef __simulation_region_h__
#define __simulation_region_h__

#include <vector>
#include <string>


#include "fvm_node_info.h"
#include "fvm_node_data.h"
#include "fvm_cell_data.h"
#include "petscvec.h"
#include "petscmat.h"
#include "advanced_model.h"
#include "enum_region.h"

namespace Parser {
  class Parameter;
}
class SolverBase;
class SimulationSystem;
class BoundaryCondition;
namespace Material {
  class MaterialBase;
}

/**
 * this is the most important data structure for numerical simulation.
 * it contains FVM Cell (it can also be used as first order FEM Cell)
 * and Node in this region, pointer to material database,
 * boundary and interface information, as well as parallel mechanism
 */
class SimulationRegion
{
  friend class BoundaryCondition;

public:

  /**
   * constructor
   */
  SimulationRegion(const std::string &name, const std::string &material, SimulationSystem & system);

  /**
   * destructor
   */
  virtual ~SimulationRegion();

  /**
   * @return the region property
   */
  virtual SimulationRegionType type() const=0;

  /**
   * virtual function for region data init,
   * take T_external as parameter
   */
  virtual void init(PetscScalar T_external)=0;

  /**
   * virtual function for re-init region data after import solution from data file
   */
  virtual void reinit_after_import()=0;

  /**
   * get the external temperature, as region's default temperature
   */
  PetscScalar T_external() const;

  /**
   * set subdomain id of the region
   */
  void set_subdomain_id(unsigned int sub_id)
  { _subdomain_id = sub_id; }

  /**
   * set subdomain id to region pointer map
   */
  void set_subdomain_id_to_region_map(const std::map<unsigned int,  SimulationRegion *> &_map)
  { _subdomain_id_to_region_map = _map; }

  /**
   * insert local mesh element into the region, only copy the pointer
   * and create cell data
   */
  virtual void insert_cell (const Elem * e)=0;

  /**
   * insert FVM node and corresponding FVM Node data into region
   * pure virtual function
   */
  virtual void insert_fvm_node(FVM_Node *)=0;

  /**
   * @return the total cell number in this region
   * @note includes on processor cell and ghost cell
   */
  unsigned int n_cell() const
  { return _region_cell.size(); }

  /**
   * @return the on processor cell number in this region
   */
  unsigned int n_on_processor_cell() const;

  /**
   * @return nth region elem
   */
  const Elem * get_region_elem(unsigned int n) const
  { return _region_cell[n]; }

  /**
   * @return nth region elem data
   */
  const FVM_CellData * get_region_elem_data(unsigned int n) const
  { return _region_cell_data[n]; }

  /**
   * @return nth region elem data
   */
  FVM_CellData * get_region_elem_Data(unsigned int n)
  { return _region_cell_data[n]; }

  /**
   * @return the FVM Node number in this region
   * @note all the node are stored. no matter which processor_id the node is.
   */
  unsigned int n_node() const
    { return _region_node.size(); }

  /**
   * @return the on processor FVM Node number in this region
   */
  unsigned int n_on_processor_node() const;

  /**
   * @return the fvm_node pointer by Node *
   * if no find, NULL is returned
   */
  FVM_Node * region_fvm_node(const Node* node) const
  {
    const_node_iterator it = _region_node.find( node->id() );
    if( it!=_region_node.end() )
      return (*it).second;
    return NULL;
  }

  /**
   * @return the fvm_node pointer by Node id
   * if no find, NULL is returned
   */
  FVM_Node * region_fvm_node(unsigned int id) const
  {
    const_node_iterator it = _region_node.find( id );
    if( it!=_region_node.end() )
      return (*it).second;
    return NULL;
  }

  /**
   * @return the node data pointer by Node *
   * if no find, NULL is returned
   */
  FVM_NodeData * region_node_data(const Node* node) const
  {
    const_node_iterator it = _region_node.find( node->id() );
    if( it!=_region_node.end() )
      return (*it).second->node_data();
    return NULL;
  }

  /**
   * @return the node data pointer by Node *
   * if no find, NULL is returned
   */
  FVM_NodeData * region_node_data(unsigned int id) const
  {
    const_node_iterator it = _region_node.find( id );
    if( it!=_region_node.end() )
      return (*it).second->node_data();
    return NULL;
  }

  /**
   * @return region's name
   */
  const std::string & name() const
  { return _region_name; }

  /**
   * @return region's material
   */
  const std::string & material() const
    { return _region_material; }

  /**
   * @return const reference of element in this region
   */
  const std::vector<const Elem *> & cells() const
    { return _region_cell; }

  /**
   * clear stored data
   */
  void clear();


  typedef std::vector<const Elem*>::iterator             element_iterator;

  typedef std::vector<const Elem*>::const_iterator       const_element_iterator;

  /**
   * region elements default begin() accessor
   */
  element_iterator elements_begin ()
  {
    return _region_cell.begin();
  }

  /**
   * region elements default end() accessor
   */
  element_iterator elements_end   ()
  {
    return _region_cell.end();
  }

  /**
   * region elements default const begin() accessor
   */
  const_element_iterator elements_begin () const
  {
    return _region_cell.begin();
  }

  /**
   * region elements default const end() accessor
   */
  const_element_iterator elements_end   () const
  {
    return _region_cell.end();
  }

  /**
   * typedef node_iterator
   */
  typedef std::map< unsigned int, FVM_Node * >::iterator             node_iterator;

  /**
   * typedef const_node_iterator
   */
  typedef std::map< unsigned int, FVM_Node * >::const_iterator       const_node_iterator;

  /**
   * node default begin() accessor
   */
  node_iterator nodes_begin        ()
  {
    return _region_node.begin();
  }

  /**
   * node default end() accessor
   */
  node_iterator nodes_end          ()
  {
    return _region_node.end();
  }

  /**
   * node default const begin() accessor
   */
  const_node_iterator nodes_begin        () const
  {
    return _region_node.begin();
  }

  /**
   * node default const end() accessor
   */
  const_node_iterator nodes_end          () const
  {
    return _region_node.end();
  }

  /**
   * for some post process
   */
  void prepare_for_use();

  /**
   * setting some physical model to region.
   */
  AdvancedModel & advanced_model()
  { return  _advanced_model; }

  /**
   * get reference to advanced physical model
   */
  const AdvancedModel & advanced_model() const
  { return  _advanced_model;}

  /**
   * get pointer to advanced physical model
   */
  const AdvancedModel * get_advanced_model() const
  { return  & _advanced_model;}

  /**
   * @return the base class of material database
   */
  virtual Material::MaterialBase * get_material_base()=0;

  /**
   * @return the optical refraction index of the region
   */
  virtual Complex get_optical_refraction(double )
  { return Complex(0.0, 0.0); }

  /**
   * @return the energy bandgap for optical simulation
   * energy requited for generate (n-p) electron pair from vacuum
   * 0.510998902MeV*2
   */
  virtual double get_optical_Eg(double)
  { return 0.510998902*1e6*2; }

  /**
   * @return relative permittivity of material
   */
  virtual double get_eps() const
  { return 1.0; }

  /**
   * @return maretial density [g cm^-3]
   */
  virtual double get_density(PetscScalar ) const
  { return 0.0; }

  /**
   * get atom fraction of region material
   */
  virtual void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction)const=0;

  /**
   * virtual function for set different model, calibrate parameters to PMI
   */
  virtual void set_pmi(const std::string &type, const std::string &model_name,
                       const std::vector<Parser::Parameter> & pmi_parameters);

  /**
   * get an information string of the PMI models
   */
  virtual std::string get_pmi_info(const std::string& type, const int verbosity);

  /**
   * @return the bound box of this region, which is described by two point
   */
  std::pair<Point, Point> region_bound_box() const
  {
    Point p1(1e30,1e30,1e30), p2(-1e30,-1e30,-1e30);

    for(const_node_iterator it=_region_node.begin(); it!=_region_node.end(); ++it )
    {
      (*it).second->root_node()->assign_min_to(p1);
      (*it).second->root_node()->assign_max_to(p2);
    }
    return std::pair<Point, Point>(p1,p2);
  }


  /**
   * @return true if 2D hanging node exist
   */
  bool is_2d_hanging_node() const;


  /**
   * @return true if 3D hanging node exist
   */
  bool is_3d_hanging_node() const;


  /**
   * store hanging node on elem side for later use
   */
  void add_hanging_node_on_side(const Node * node, const Elem * elem, unsigned int s)
  {
    node_iterator it = _region_node.find(node->id());
    genius_assert( it!=_region_node.end() );

    const FVM_Node * fvm_node = (*it).second;
    _hanging_node_on_elem_side[fvm_node] = std::pair<const Elem *, unsigned int>(elem, s);
  }

  /**
   * store hanging node on elem edge for later use
   */
  void add_hanging_node_on_edge(const Node * node, const Elem * elem, unsigned int e)
  {
    node_iterator it = _region_node.find(node->id());
    genius_assert( it!=_region_node.end() );

    const FVM_Node * fvm_node = (*it).second;
    _hanging_node_on_elem_edge[fvm_node] = std::pair<const Elem *, unsigned int>(elem, e);

  }

  /**
   * @return the number of hanging nodes on element side
   */
  unsigned int n_hanging_node_on_elem_side() const
  { return _hanging_node_on_elem_side.size(); }


  /**
   * @return the number of hanging nodes on element edge
   */
  unsigned int n_hanging_node_on_elem_edge() const
  { return _hanging_node_on_elem_edge.size(); }


  typedef std::map<const FVM_Node *, std::pair<const Elem *, unsigned int> >::const_iterator       hanging_node_on_elem_side_iterator;
  typedef std::map<const FVM_Node *, std::pair<const Elem *, unsigned int> >::const_iterator       hanging_node_on_elem_edge_iterator;

  /**
   * const hanging_node_on_elem_side begin() accessor
   */
  hanging_node_on_elem_side_iterator hanging_node_on_elem_side_begin        () const
  {
    return _hanging_node_on_elem_side.begin();
  }

  /**
   * const hanging_node_on_elem_side end() accessor
   */
  hanging_node_on_elem_side_iterator hanging_node_on_elem_side_end          () const
  {
    return _hanging_node_on_elem_side.end();
  }

  /**
   * const hanging_node_on_elem_edge begin() accessor
   */
  hanging_node_on_elem_edge_iterator hanging_node_on_elem_edge_begin        () const
  {
    return _hanging_node_on_elem_edge.begin();
  }

  /**
   * const hanging_node_on_elem_edge end() accessor
   */
  hanging_node_on_elem_edge_iterator hanging_node_on_elem_edge_end          () const
  {
    return _hanging_node_on_elem_edge.end();
  }

  /**
   * save a pointer of solver in region level
   */
  void set_current_solver(SolverBase * sv)
  { solver = sv; }

protected:

  /**
   * the region's name, defined by user
   */
  std::string                    _region_name;


  /**
   * the region's material, defined by user
   */
  std::string                    _region_material;


  /**
   * every region hold this data, then one can find neighbor region
   */
  std::map<unsigned int,  SimulationRegion *>  _subdomain_id_to_region_map;

  /**
   * record the elements which belong to this region,
   * only local element (on processor and ghost element) are recorded.
   */
  std::vector<const Elem *>      _region_cell;

  /**
   * cell data
   */
  std::vector<FVM_CellData *>    _region_cell_data;

  /**
   *  the node belongs to this region. stored as \< node_id, FVM_Node *\>
   *  all the nodes (on and off processor) are stored.
   *  I am afraid that map may have efficent problem, use vector instead?
   */
  std::map< unsigned int, FVM_Node * > _region_node;


  /**
   * sub domain index
   */
  unsigned int             _subdomain_id;

  /**
   * region advanced models
   */
  AdvancedModel            _advanced_model;

  /**
   * stores the hanging node which lies on the side center of of an element.
   * for 2D case, the side is an edge. hanging node lies on edge center
   * for 3D case, only the Quad4 side has center hanging node
   */
  std::map<const FVM_Node *, std::pair<const Elem *, unsigned int> >  _hanging_node_on_elem_side;


  /**
   * stores the hanging node which lies on the edge center of of an element.
   * this is for 3D case only.
   * maybe several elem shares the same hanging node. we record only one
   */
  std::map<const FVM_Node *, std::pair<const Elem *, unsigned int> >  _hanging_node_on_elem_edge;

  /**
   * the reference to corresponding SimulationSystem
   */
  SimulationSystem    & _system;


  /**
   * a pointer to solver, each solver can set this pointer point to itself
   * then region functions can access some solver level cariable
   */
  SolverBase * solver;

public:


  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for fill vector of poisson's equation.
   *
   * filling solution data from FVM_NodeData into petsc vector of poisson's equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   *
   * @param x                global solution vector
   * @param L                the left scaling vector, usually contains the cell volumn
   * @note fill items of global solution vector with belongs to local processor
   * each derived region should override it
   */
  virtual void Poissin_Fill_Value(Vec x, Vec L)=0;

  /**
   * @brief virtual function for evaluating poisson's equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of poisson's equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void Poissin_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating hanging node for poisson's equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void Poissin_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian for hanging node of poisson's equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void Poissin_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for update solution value of poisson's equation.
   *
   * update solution data of FVM_NodeData by petsc vector of poisson's equation.
   *
   * @param x                global solution vector
   *
   * @note each derived region should override it
   */
  virtual void Poissin_Update_Solution(PetscScalar *lxx)=0;




  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////



  /**
   * @brief virtual function for fill vector of level 1 DDM equation.
   *
   * filling solution data from FVM_NodeData into petsc vector of level 1 DDM equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   *
   * @param x                global solution vector
   * @param L                the left scaling vector, usually contains the cell volumn
   * @note fill items of global solution vector with belongs to local processor
   * each derived region should override it
   */
  virtual void DDM1_Fill_Value(Vec x, Vec L)=0;

  /**
   * @brief virtual function for evaluating level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating time derivative term of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of time derivative term of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating hanging node for level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of hanging node for level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;


  /**
   * @brief virtual function for update solution value of level 1 DDM equation.
   *
   * update solution data of FVM_NodeData by petsc vector of level 1 DDM equation.
   *
   * @param x                global solution vector
   *
   * @note each derived region should override it
   */
  virtual void DDM1_Update_Solution(PetscScalar *lxx)=0;




  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for L1 Hall DDM------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for fill vector of level 1 DDM equation with hall correction.
   *
   * filling solution data from FVM_NodeData into petsc vector of level 1 DDM equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   *
   * @param x                global solution vector
   * @param L                the left scaling vector, usually contains the cell volumn
   * @note fill items of global solution vector with belongs to local processor
   * each derived region should override it
   */
  virtual void HALL_Fill_Value(Vec x, Vec L)=0;

  /**
   * @brief virtual function for evaluating level 1 DDM equation with hall correction.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void HALL_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of level 1 DDM equation with hall correction.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void HALL_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating time derivative term of level 1 DDM equation with hall correction.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void HALL_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of time derivative term of level 1 DDM equation with hall correction.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void HALL_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating hanging node for level 1 DDM equation with hall correction.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void HALL_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of hanging node for level 1 DDM equation with hall correction.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void HALL_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;


  /**
   * @brief virtual function for update solution value of level 1 DDM equation with hall correction.
   *
   * update solution data of FVM_NodeData by petsc vector of level 1 DDM equation with hall correction.
   *
   * @param x                global solution vector
   *
   * @note each derived region should override it
   */
  virtual void HALL_Update_Solution(PetscScalar *lxx)=0;




  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////



  /**
   * @brief virtual function for fill vector of level 2 DDM equation.
   *
   * filling solution data from FVM_NodeData into petsc vector of level 2 DDM equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   *
   * @param x                global solution vector
   * @param L                the left scaling vector, usually contains the cell volumn
   * @note fill items of global solution vector with belongs to local processor
   * each derived region should override it
   */
  virtual void DDM2_Fill_Value(Vec x, Vec L)=0;

  /**
   * @brief virtual function for evaluating level 2 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of level 2 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating time derivative term of level 2 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM2_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of time derivative term of level 2 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM2_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating hanging node for level 2 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM2_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of hanging node for level 2 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM2_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for update solution value of level 2 DDM equation.
   *
   * update solution data of FVM_NodeData by petsc vector of level 2 DDM equation.
   *
   * @param x                global solution vector
   *
   * @note each derived region should override it
   */
  virtual void DDM2_Update_Solution(PetscScalar *lxx)=0;



  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////



  /**
   * @brief virtual function for get nodal variable number
   * @note depedent on EBM Level
   * each derived region should override it
   */
  virtual unsigned int ebm_n_variables() const=0;

  /**
   * @brief virtual function for get offset of nodal variable
   * @return local offset of variable var to global offset
   * @note depedent on EBM Level
   * each derived region should override it
   */
  virtual unsigned int ebm_variable_offset(SolutionVariable var) const=0;

  /**
   * @brief virtual function for fill vector of level 3 EBM equation.
   *
   * filling solution data from FVM_NodeData into petsc vector of level 3 EBM equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   *
   * @param x                global solution vector
   * @param L                the left scaling vector, usually contains the cell volumn
   * @note fill items of global solution vector with belongs to local processor
   * each derived region should override it
   */
  virtual void EBM3_Fill_Value(Vec x, Vec L)=0;

  /**
   * @brief virtual function for evaluating level 3 EBM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of level 3 EBM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;


  /**
   * @brief virtual function for evaluating hanging node for level 3 EBM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void EBM3_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of hanging node for level 3 EBM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void EBM3_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;


  /**
   * @brief virtual function for evaluating time derivative term of level 3 EBM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void EBM3_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of time derivative term of level 3 EBM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void EBM3_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for update solution value of level 3 EBM equation.
   *
   * update solution data of FVM_NodeData by petsc vector of level 3 EBM equation.
   *
   * @param x                global solution vector
   *
   * @note each derived region should override it
   */
  virtual void EBM3_Update_Solution(PetscScalar *lxx)=0;


  //////////////////////////////////////////////////////////////////////////////////
  //-----------------   Vec and Matrix evaluate for DDMAC   ----------------------//
  //////////////////////////////////////////////////////////////////////////////////


  /**
   * @brief virtual function for fill vector of DDMAC equation.
   *
   * filling solution data from FVM_NodeData into petsc vector of level 3 EBM equation.
   *
   * @param x                global solution vector
   * @param L                the left scaling vector, usually contains the cell volumn
   * @note fill items of global solution vector with belongs to local processor
   * each derived region should override it
   */
  virtual void DDMAC_Fill_Value(Vec x, Vec L) const=0;

  /**
   * @brief virtual function for fill matrix of DDMAC equation.
   *
   * filling AC matrix entry by Jacobian matrix of level 3 EBM equation.
   *
   * @param A                AC matrix
   * @param b                AC rhs vector
   * @param J                Jacobian matrix
   * @param omega            AC frequency
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note fill entry of AC matrix A with Jacobian matrix J at frequency omega
   * each derived region should override it
   */
  virtual void DDMAC_Fill_Matrix_Vector(Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag) const=0;

  /**
   * @brief virtual function for fill matrix of DDMAC equation.
   *
   * filling AC matrix entry by Jacobian matrix of level 3 EBM equation for a specified fvm_node.
   * @param fvm_node           process this node
   * @param A                  AC matrix
   * @param J                  Jacobian matrix
   * @param omega              AC frequency
   * @param add_value_flag     flag for last operator is ADD_VALUES
   * @param adjacent_region    the AC matrix entry add to adjacent_fvm_node in this region
   * @param adjacent_fvm_node  the AC matrix entry add to this node in adjacent_region region
   *
   * @note fill entry of AC matrix A with Jacobian matrix J at frequency omega
   * each derived region should override it
   */
  virtual void DDMAC_Fill_Nodal_Matrix_Vector(const FVM_Node *fvm_node, Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag,
                                              const SimulationRegion * adjacent_region=NULL,
                                              const FVM_Node * adjacent_fvm_node=NULL) const=0;

  /**
   * @brief virtual function for fill matrix of DDMAC equation.
   *
   * filling AC matrix entry by Jacobian matrix of level 3 EBM equation for specified variable of fvm_node.
   * @param fvm_node           process this node
   * @param var                process the var of this node
   * @param A                  AC matrix
   * @param J                  Jacobian matrix
   * @param omega              AC frequency
   * @param add_value_flag     flag for last operator is ADD_VALUES
   * @param adjacent_region    the AC matrix entry add to adjacent_fvm_node in this region
   * @param adjacent_fvm_node  the AC matrix entry add to this node in adjacent_region region
   *
   * @note fill entry of AC matrix A with Jacobian matrix J at frequency omega
   * each derived region should override it
   */
  virtual void DDMAC_Fill_Nodal_Matrix_Vector(const FVM_Node *fvm_node, const SolutionVariable var,
                                              Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag,
                          const SimulationRegion * adjacent_region=NULL,
                                              const FVM_Node * adjacent_fvm_node=NULL) const=0;

  /**
   * @brief virtual function for update solution value of DDMAC equation.
   *
   *
   * @param x                global solution vector
   *
   * @note each derived region should override it
   */
  virtual void DDMAC_Update_Solution(PetscScalar *lxx)=0;

};



#endif //#ifndef __simulation_region_h__
