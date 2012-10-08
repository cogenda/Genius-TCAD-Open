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

#include "variable_define.h"
#include "fvm_node_info.h"
#include "fvm_node_data.h"
#include "fvm_cell_data.h"
#include "petscvec.h"
#include "petscmat.h"
#include "advanced_model.h"
#include "enum_region.h"

#if defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
#elif defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER) || defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
#endif

#if defined(HAVE_TR1_UNORDERED_SET)
# include <tr1/unordered_set>
#elif defined(HAVE_TR1_UNORDERED_SET_WITH_STD_HEADER) || defined(HAVE_UNORDERED_SET)
# include <unordered_set>
#else
# include <set>
#endif

namespace Parser {
  class Parameter;
}
class SolverBase;
class MeshBase;
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
  SimulationRegion(const std::string &name, const std::string &material, const double T, const double z);

  /**
   * destructor
   */
  virtual ~SimulationRegion();

  /**
   * @return the region property
   */
  virtual SimulationRegionType type() const=0;

  /**
   * @return the region property in string
   */
  virtual std::string type_name() const=0;

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
  double T_external() const;

  /**
   * @return the dimension in Z direction.
   * @note only used for 2D mesh which lies on xy plane
   */
  double z_width() const;

  /**
   * set subdomain id of the region
   */
  void set_subdomain_id(unsigned int sub_id)
  { _subdomain_id = sub_id; }

  /**
   * get subdomain id of the region
   */
  unsigned int subdomain_id() const
  { return _subdomain_id; }

  /**
   * add boundary condition belongs to this region
   */
  void add_boundary(BoundaryCondition * bc);

  /**
   * set subdomain id to region pointer map
   */
  static void set_subdomain_id_to_region_map(const std::map<unsigned int,  SimulationRegion *> &_map)
  { _subdomain_id_to_region_map = _map; }

  /**
   * reseve memory for region data block
   */
  void reserve_data_block(unsigned int n_cell_data, unsigned int n_node_data);

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
   * @return the total edge number in this region
   */
  unsigned int n_edge() const
  { return _region_edges.size(); }

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
  FVM_CellData * get_region_elem_data(unsigned int n)
  { return _region_cell_data[n]; }

  /**
   * @return the FVM Node number in this region
   * @note all the node are count. no matter which processor_id the node is.
   */
  unsigned int n_node() const
  { return _n_region_node; }

  /**
   * @return the on processor FVM Node number in this region
   */
  unsigned int n_on_processor_node() const;

  /**
   * get all the region node ids by order,
   * must executed in parallel.
   */
  void region_node(std::vector<unsigned int> & nodes) const;

  /**
   * @return the fvm_node pointer by Node *
   * if no find, NULL is returned
   */
  FVM_Node * region_fvm_node(const Node* node) const;

  /**
   * @return the fvm_node pointer by Node id
   * if no find, NULL is returned
   */
  FVM_Node * region_fvm_node(unsigned int id) const;

  /**
   * @return the node data pointer by Node *
   * if no find, NULL is returned
   */
  FVM_NodeData * region_node_data(const Node* node) const;

  /**
   * @return the node data pointer by Node *
   * if no find, NULL is returned
   */
  FVM_NodeData * region_node_data(unsigned int id) const;

  /**
   * @return region's name
   */
  const std::string & name() const
  { return _region_name; }


  /**
   * @return region's label(name)
   */
  const std::string & label() const
  { return _region_name; }


  /**
   * @return region's material
   */
  const std::string & material() const
    { return _region_material; }

  /**
   * @return the boundingbox of the region
   */
  const std::pair<Point, Point> & boundingbox() const
  { return _region_bounding_box; }


  /**
   * @return the boundingbox of the region
   */
  std::pair<Point, Real> boundingsphere() const
  {
    const double  diag = (_region_bounding_box.second - _region_bounding_box.first).size();
    const Point   cent = (_region_bounding_box.second + _region_bounding_box.first)*0.5;
    return std::pair<Point, Real>(cent, diag);
  }


  /**
   * @return const reference of element in this region
   */
  const std::vector<const Elem *> & cells() const
    { return _region_cell; }


  const std::map<short int, BoundaryCondition *> & region_boundaries() const
  { return _region_boundaries; }

  /**
   * clear stored data
   */
  virtual void clear();

  /**
   * @return the quality if each fvm cell in (0-poor, 1-fine]
   */
  virtual Real fvm_cell_quality() const;


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
   * typedef local_node_iterator
   */
  typedef std::vector<FVM_Node * >::iterator             local_node_iterator;

  /**
   * typedef const_local_node_iterator
   */
  typedef std::vector<FVM_Node * >::const_iterator       const_local_node_iterator;

  /**
   * local node default begin() accessor
   */
  local_node_iterator on_local_nodes_begin        ()
  {
    return _region_local_node.begin();
  }

  /**
   * local node default end() accessor
   */
  local_node_iterator on_local_nodes_end          ()
  {
    return _region_local_node.end();
  }

  /**
   * local node default const begin() accessor
   */
  const_local_node_iterator on_local_nodes_begin        () const
  {
    return _region_local_node.begin();
  }

  /**
   * local node default const end() accessor
   */
  const_local_node_iterator on_local_nodes_end          () const
  {
    return _region_local_node.end();
  }


  /**
   * typedef processor_node_iterator
   */
  typedef std::vector<FVM_Node * >::iterator             processor_node_iterator;

  /**
   * typedef const_processor_node_iterator
   */
  typedef std::vector<FVM_Node * >::const_iterator       const_processor_node_iterator;

  /**
   * processor node default begin() accessor
   */
  processor_node_iterator on_processor_nodes_begin        ()
  {
    return _region_processor_node.begin();
  }

  /**
   * processor node default end() accessor
   */
  processor_node_iterator on_processor_nodes_end          ()
  {
    return _region_processor_node.end();
  }

  /**
   * processor node default const begin() accessor
   */
  const_processor_node_iterator on_processor_nodes_begin        () const
  {
    return _region_processor_node.begin();
  }

  /**
   * processor node default const end() accessor
   */
  const_processor_node_iterator on_processor_nodes_end          () const
  {
    return _region_processor_node.end();
  }



  /**
   * typedef edge_iterator
   */
  typedef std::vector< std::pair<FVM_Node *, FVM_Node *> >::iterator             edge_iterator;


  /**
   * typedef const_edge_iterator
   */
  typedef std::vector< std::pair<FVM_Node *, FVM_Node *> >::const_iterator       const_edge_iterator;


  /**
   * processor edge default begin() accessor
   */
  edge_iterator edges_begin        ()
  {
    return _region_edges.begin();
  }

  /**
   * processor edge default end() accessor
   */
  edge_iterator edges_end          ()
  {
    return _region_edges.end();
  }

  /**
   * processor edge default const begin() accessor
   */
  const_edge_iterator edges_begin        () const
  {
    return _region_edges.begin();
  }

  /**
   * processor edge default const end() accessor
   */
  const_edge_iterator edges_end        () const
  {
    return _region_edges.end();
  }

  /**
   * @return the corresponding location of an element's edge in _region_edges
   * by given an element pointer, and the local index of the edge
   */
  unsigned int elem_edge_index(const Elem* elem, unsigned int e) const
  { return _region_elem_edge_in_edges_index.find(elem)->second[e]; }

  /**
   * (re)build _region_local_node and _region_processor_node for fast iteration
   */
  void rebuild_region_fvm_node_list();

  /**
   * for some pre process
   */
  virtual void prepare_for_use();

  /**
   * for some pre process, which should be executed in parallel
   * call it after prepare_for_use
   */
  virtual void prepare_for_use_parallel();

  /**
   * delete fvm_node NOT on this processor, dangerous
   */
  void remove_remote_object();

  /**
   * @return true if we are neighbor
   */
  bool is_neighbor(const SimulationRegion *r) const;

  /**
   * @return region neighbors
   */
  const std::vector<SimulationRegion *> & region_neighbors() const
  { return _region_neighbors; }


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
  virtual Material::MaterialBase * get_material_base() const=0;

  /**
   * set the variables for this region
   */
  virtual void set_region_variables()=0;

  /**
   * add a solution variable by full define, also allocate memory when the variable_valid flag is true
   * @return variable_index
   */
  unsigned int add_variable(const SimulationVariable &v);

  /**
   * add a predefined solution variable by it's name define, also allocate memory
   * @return true for success
   */
  bool add_variable(const std::string &v, DataLocation);

  /**
   * @return true when the variable exist
   */
  bool has_variable(const std::string &v, DataLocation) const;

  /**
   * get a SimulationVariable by its name and location
   * should make sure the variable exist
   */
  const SimulationVariable & get_variable(const std::string &v, DataLocation) const;

  /**
   * get a SimulationVariable by its name and location, @return true when this variable exist
   */
  bool get_variable(const std::string &v, DataLocation, SimulationVariable & ) const;

  /**
   * get all the user_defined SimulationVariable by data location and data type
   */
  void get_user_defined_variable(DataLocation, DataType, std::vector<SimulationVariable> & ) const;

  /**
   * region level data access functions, gather the variable in to the vector
   * must executed in parallel. the data is ordered by node id or cell id,
   * thus has the same order as _region_node or _region_cell
   * @return true for success
   */
  template <typename T>
  bool get_variable_data(const std::string &v, DataLocation, std::vector<T> &) const;


  /**
   * set all the variables v to value d
   * must executed in parallel.
   * @return true for success
   */
  template <typename T>
  bool set_variable_data(const std::string &v, DataLocation, const T d) ;

  /**
   * @return the region node based variables
   */
  const std::map<std::string, SimulationVariable> & region_point_variables() const
  { return _region_point_variables; }

  /**
   * @return the region cell based variables
   */
  const std::map<std::string, SimulationVariable> & region_cell_variables() const
  { return _region_cell_variables; }

  /**
   * sync point variable with ghost node
   * @return true for success
   */
  template <typename T>
  bool sync_point_variable(const std::string &v);

  /**
   * @return the optical refraction index of the region
   */
  virtual Complex get_optical_refraction(double lamda)
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
   * set material conductance [A/V/cm]
   */
  virtual void set_conductance(double )  {}

  /**
   * @return material conductance [A/V/cm]
   */
  virtual double get_conductance() const
  { return 0.0; }

  /**
   * @return material density [g cm^-3]
   */
  virtual double get_density(PetscScalar ) const
  { return 0.0; }

  /**
   * @return affinity of material
   */
  virtual double get_affinity(PetscScalar T) const
  { return 0.0; }

  /**
   * get atom fraction of region material
   */
  virtual void atom_fraction(std::vector<std::string> &atoms, std::vector<double> & fraction)const=0;

  /**
   * virtual function for set different model, calibrate parameters to PMI
   */
  virtual void set_pmi(const std::string &type, const std::string &model_name, std::vector<Parser::Parameter> & pmi_parameters);

  /**
   * get an information string of the PMI models
   */
  virtual std::string get_pmi_info(const std::string& type, const int verbosity);

  /**
   * @return true if 2D hanging node exist
   */
  bool has_2d_hanging_node() const
  { return _hanging_node_on_elem_side_flag; }


  /**
   * @return true if 3D hanging node exist
   */
  bool has_3d_hanging_node() const
  { return _hanging_node_on_elem_edge_flag; }

  /**
   * store hanging node on elem side for later use
   */
  void add_hanging_node_on_side(const Node * node, const Elem * elem, unsigned int s);

  /**
   * store hanging node on elem edge for later use
   */
  void add_hanging_node_on_edge(const Node * node, const Elem * elem, unsigned int e);

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
   * approx memory usage
   */
  virtual size_t memory_size() const;

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
   * region's default temperature
   */
  double                        _T_external;


  /**
   * when 2D mesh is used, we may need the dimension in z direction
   */
  double                        _z_width;

  /**
   * set subdomain_id_to_region_map as static member
   */
  static std::map<unsigned int,  SimulationRegion *>  _subdomain_id_to_region_map;

  /**
   * neighbor regions
   */
  std::vector<SimulationRegion *> _region_neighbors;

  /**
   * the region boundaries
   */
  std::map<short int, BoundaryCondition *> _region_boundaries;

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
   * data block for cell based value
   */
  DataStorage _cell_data_storage;


  /**
   * all the nodes belongs to this region
   */
  unsigned int _n_region_node;

  /**
   *  the node belongs to this region. stored as \< node_id, FVM_Node *\>
   */
  std::map< unsigned int, FVM_Node * > _region_node;


  /**
   * on local nodes (on processor nodes + ghost nodes) belong to this region.
   */
  std::vector<FVM_Node *>    _region_local_node;

  /**
   * on processor nodes belong to this region.
   */
  std::vector<FVM_Node *>    _region_processor_node;

  /**
   * ghost nodes belong to this region.
   */
  std::vector<FVM_Node *>    _region_ghost_node;

  /**
   * on processor nodes which is the ghost nodes in other regions, build when required.
   */
  std::vector<FVM_Node *>    _region_image_node;

  /**
   * data block for node based value
   */
  DataStorage _node_data_storage;

  /**
   * the edges belongs to this regon, for fast FVM integral
   * the two fvm_node of this edge is ordered as id(1) \< id(2)
   */
  std::vector< std::pair<FVM_Node *, FVM_Node *> > _region_edges;

  /**
   * the corresponding location of an element's edge in _region_edges
   * by given an element pointer, and the local index of the edge
   * use unordered_map when possible
   */
#if defined(HAVE_UNORDERED_MAP)
    std::unordered_map<const Elem *, std::vector<unsigned int> > _region_elem_edge_in_edges_index;
#elif defined(HAVE_TR1_UNORDERED_MAP) || defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER)
    std::tr1::unordered_map<const Elem *, std::vector<unsigned int> > _region_elem_edge_in_edges_index;
#else
    std::map<const Elem *, std::vector<unsigned int> > _region_elem_edge_in_edges_index;
#endif


  /**
   * the boundingbox of the region
   */
  std::pair<Point, Point>    _region_bounding_box;

  /**
   * sub domain index
   */
  unsigned int             _subdomain_id;

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


  bool _hanging_node_on_elem_side_flag;

  bool _hanging_node_on_elem_edge_flag;


  /**
   * region advanced models
   */
  AdvancedModel            _advanced_model;


  /**
   * the point based solution variables defined for region, it should be filled by derived class
   */
  std::map<std::string, SimulationVariable>  _region_point_variables;

  /**
   * the cell based solution variables defined for region, it should be filled by derived class
   */
  std::map<std::string, SimulationVariable>  _region_cell_variables;


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
   * @brief virtual function for evaluating pseudo time step of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region can override it
   */
  virtual void DDM1_Pseudo_Time_Step_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag) {}

  /**
   * @brief virtual function for evaluating Jacobian of pseudo time step of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region can override it
   */
  virtual void DDM1_Pseudo_Time_Step_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag) {}

  /**
   * @brief virtual function for convergence test of pseudo time step of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @return                 the number of nodes not convergenced yet
   *
   * @note each derived region can override it
   */
  virtual int DDM1_Pseudo_Time_Step_Convergence_Test(PetscScalar * x) { return 0; }

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
  //---------------Function and Jacobian evaluate for new L1 DDM------------------//
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
  virtual void DDM1R_Fill_Value(Vec x, Vec L)
  { this->DDM1_Fill_Value(x, L); }

  /**
   * @brief virtual function for evaluating level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1R_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
  { this->DDM1_Function(x, f, add_value_flag); }

  /**
   * @brief virtual function for evaluating Jacobian of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1R_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
  { this->DDM1_Jacobian(x, jac, add_value_flag); }

  /**
   * @brief virtual function for evaluating time derivative term of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1R_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
  { this->DDM1_Time_Dependent_Function(x, f, add_value_flag); }

  /**
   * @brief virtual function for evaluating Jacobian of time derivative term of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1R_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
  { this->DDM1_Time_Dependent_Jacobian(x, jac, add_value_flag); }

  /**
   * @brief virtual function for evaluating pseudo time step of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region can override it
   */
  virtual void DDM1R_Pseudo_Time_Step_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
  { this->DDM1_Pseudo_Time_Step_Function(x, f, add_value_flag);}

  /**
   * @brief virtual function for evaluating Jacobian of pseudo time step of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region can override it
   */
  virtual void DDM1R_Pseudo_Time_Step_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
  { this->DDM1_Pseudo_Time_Step_Jacobian(x, jac, add_value_flag);}

  /**
   * @brief virtual function for convergence test of pseudo time step of level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @return                 the number of nodes not convergenced yet
   *
   * @note each derived region can override it
   */
  virtual int DDM1R_Pseudo_Time_Step_Convergence_Test(PetscScalar * x)
  { return DDM1_Pseudo_Time_Step_Convergence_Test(x); }

  /**
   * @brief virtual function for evaluating hanging node for level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1R_Function_Hanging_Node(PetscScalar * x, Vec f, InsertMode &add_value_flag)
  {this->DDM1_Function_Hanging_Node(x, f, add_value_flag);}

  /**
   * @brief virtual function for evaluating Jacobian of hanging node for level 1 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1R_Jacobian_Hanging_Node(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
  { this->DDM1_Jacobian_Hanging_Node(x, jac, add_value_flag);}

  /**
   * @brief virtual function for update solution value of level 1 DDM equation.
   *
   * update solution data of FVM_NodeData by petsc vector of level 1 DDM equation.
   *
   * @param x                global solution vector
   *
   * @note each derived region should override it
   */
  virtual void DDM1R_Update_Solution(PetscScalar *lxx)
  { this->DDM1_Update_Solution(lxx); }

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
   * @param B                magnetic field
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void HALL_Function(const VectorValue<double> & B, PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of level 1 DDM equation with hall correction.
   *
   * @param B                magnetic field
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void HALL_Jacobian(const VectorValue<double> & B, PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

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
   * @brief virtual function for evaluating pseudo time step of level 2 DDM equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region can override it
   */
  virtual void DDM2_Pseudo_Time_Step_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag) {}

  /**
   * @brief virtual function for evaluating Jacobian of pseudo time step of level 2 DDM equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region can override it
   */
  virtual void DDM2_Pseudo_Time_Step_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag) {}

  /**
   * @brief virtual function for convergence test of pseudo time step of level 2 DDM equation.
   *
   * @param x                local unknown vector
   * @return                 the number of nodes not convergenced yet
   *
   * @note each derived region can override it
   */
  virtual int DDM2_Pseudo_Time_Step_Convergence_Test(PetscScalar * x) { return 0; }

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
   * @brief virtual function for fill transport matrix of DDMAC equation.
   *
   * filling AC transformation matrix (as preconditioner) entry by Jacobian matrix
   *
   * @param T                AC transformation matrix
   * @param J                Jacobian matrix
   * @param omega            AC frequency
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note fill entry of AC transpose matrix T with Jacobian matrix J at frequency omega
   * each derived region should override it
   */
  virtual void DDMAC_Fill_Transformation_Matrix(Mat T, const Mat J, const double omega, InsertMode & add_value_flag) const=0;

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
   * @brief virtual function for fill matrix of DDMAC equation.
   *
   * filling AC matrix entry by force variable of FVM_Node1 equals to FVM_Node2
   * @param fvm_node           process this node
   * @param A                  AC matrix
   * @param add_value_flag     flag for last operator is ADD_VALUES
   * @param adjacent_region    the AC matrix entry add to adjacent_fvm_node in this region
   * @param adjacent_fvm_node  the AC matrix entry add to this node in adjacent_region region
   * each derived region should override it
   */
  virtual void DDMAC_Force_equal(const FVM_Node *fvm_node, Mat A, InsertMode & add_value_flag,
                                 const SimulationRegion * adjacent_region=NULL,
                                 const FVM_Node * adjacent_fvm_node=NULL) const=0;

  /**
   * @brief virtual function for fill matrix of DDMAC equation.
   *
   * filling AC matrix entry by force given variable of FVM_Node1 equals to FVM_Node2
   * @param fvm_node           process this node
   * @param var                process the var of this node
   * @param A                  AC matrix
   * @param add_value_flag     flag for last operator is ADD_VALUES
   * @param adjacent_region    the AC matrix entry add to adjacent_fvm_node in this region
   * @param adjacent_fvm_node  the AC matrix entry add to this node in adjacent_region region
   * each derived region should override it
   */
  virtual void DDMAC_Force_equal(const FVM_Node *fvm_node, const SolutionVariable var,
                                 Mat A, InsertMode & add_value_flag,
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



#ifdef COGENDA_COMMERCIAL_PRODUCT
  //////////////////////////////////////////////////////////////////////////////////
  //----------------- functions for Gummel DDML1 solver --------------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for build RHS and matrix for gummel carrier equation.
   *
   * @param carrier          carrier type
   * @param x                local unknown vector
   * @param A                petsc matrix as dF/dx
   * @param r                petsc vector F(x)
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note only semiconductor region should override it
   */
  virtual void DDM1_Gummel_Carrier(const std::string & carrier, PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag) {}


  /**
   * @brief virtual function for build RHS and matrix for gummel carrier equation.
   *
   * @param x                local unknown vector
   * @param f                petsc global function vector
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note only semiconductor region should override it
   */
  virtual void DDM1_Implicit_Gummel_Carrier_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag) {}

  /**
   * @brief virtual function for evaluating Jacobian for gummel carrier equation.
   *
   * @param x                local unknown vector
   * @param jac              petsc global jacobian matrix
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1_Implicit_Gummel_Carrier_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag) {}


  /**
   * @brief virtual function for build RHS and matrix for half implicit current continuity equation.
   *
   * @param x                local unknown vector
   * @param A                petsc matrix as dF/dx
   * @param r                petsc vector F(x)
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1_Half_Implicit_Current(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag) {}


  /**
   * @brief virtual function for build RHS and matrix for half implicit poisson correction equation.
   *
   * @param x                local unknown vector
   * @param A                petsc matrix as dF/dx
   * @param r                petsc vector F(x)
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note each derived region should override it
   */
  virtual void DDM1_Half_Implicit_Poisson_Correction(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag) {}


  /**
   * @brief virtual function for build RHS and matrix for half implicit poisson correction with Polsky's method.
   *
   * @param x                local unknown vector
   * @param A                petsc matrix as dF/dx
   * @param r                petsc vector F(x)
   * @param add_value_flag   flag for last operator is ADD_VALUES
   *
   * @note semiconductor region should override it
   */
  virtual void DDM1_Half_Implicit_Poisson_Correction_Polsky(PetscScalar * x, Mat A, Vec r, InsertMode &add_value_flag) {}

#endif

  //////////////////////////////////////////////////////////////////////////////////
  //----------------- functions for Fast Hydrodynamic solver  --------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for fill vector of Hydrodynamic equation.
   *
   * @param x                global solution vector
   * @param vol              the 1.0/volume of each CV
   * @note fill items of global solution vector with belongs to local processor
   * only semiconductor region need to override it
   */
  virtual void HDM_Fill_Value(Vec /*x*/, Vec /*vol*/) {}

  /**
   * @brief virtual function for evaluating flux of Hydrodynamic equation.
   *
   * @param lx               local solution array
   * @param flux             flux vector
   * @param t                local time step vector
   *
   * @note only semiconductor region need to override it
   */
  virtual void HDM_Flux(const PetscScalar * /*lx*/, Vec /*flux*/, Vec /*t*/) {}


  /**
   * @brief virtual function for evaluating flux of Hydrodynamic equation.
   *
   * @param lx               local solution array
   * @param lt               local time step vector
   * @param x                solution vector
   *
   * @note only semiconductor region need to override it
   */
  virtual void HDM_Source(const PetscScalar * /*lx*/, const PetscScalar * /*lt*/, Vec /*x*/) {}


  /**
   * @brief virtual function for update solution value of Hydrodynamic equation.
   *
   *
   * @param x                local solution vector
   *
   * @note only semiconductor region need to override it
   */
  virtual void HDM_Update_Solution(const PetscScalar * /*x*/) {}


  //////////////////////////////////////////////////////////////////////////////////
  //-----------------  functions for Linear Poissin solver   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for build matrix of linear poisson's equation.
   *
   * @param A                matrix
   *
   * @note each derived region should override it
   */
  virtual void LinearPoissin_Matrix(Mat A, InsertMode &add_value_flag) = 0;


  /**
   * @brief virtual function for build RHS vector of linear poisson's equation.
   *
   * @param b                RHS vector
   *
   * @note each derived region should override it
   */
  virtual void LinearPoissin_RHS(Vec b, InsertMode &add_value_flag) = 0;


  /**
   * @brief virtual function for update solution value of linear poisson's equation.
   *
   * @param x                local solution vector
   *
   * @note each derived region should override it
   */
  virtual void LinearPoissin_Update_Solution(const PetscScalar * x) = 0;


  //////////////////////////////////////////////////////////////////////////////////
  //-----------------  functions for mobility evaluation     ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for evaluating mobility of edge.
   *
   * the moblity is ordered as region's edge, and weighted by partial area of the edge
   *
   * @note only semiconductor region should override it
   */
  virtual void Mob_Evaluation( std::vector< std::pair<unsigned int, unsigned int> > &edge,
                               std::vector< std::pair<double, double> > & mob,
                               std::vector< double > & weight) const {}


};



#endif //#ifndef __simulation_region_h__
