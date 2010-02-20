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

//  $Id: boundary_condition.h,v 1.39 2008/07/09 09:10:08 gdiso Exp $

#ifndef __boundary_condition_h__
#define __boundary_condition_h__


#include <vector>
#include <string>

#include "genius_env.h"
#include "log.h"
#include "node.h"
#include "elem.h"
#include "fvm_node_info.h"
#include "external_circuit.h"
#include "material_define.h"
#include "physical_unit.h"
#include "petscvec.h"
#include "petscmat.h"
#include "enum_region.h"

//predefine
class FVM_Node;

/**
 *  Boundary type is the geometry type of boundary,
 *  which can be a "boundary" or "interface" of two region with
 *  different subdomain id.
 *  Please don't be confused with boundary condition type,
 *  which only has mathematic/physical property.
 *  The INTER_CONNECT is a special bc type for electrode inter-connect
 */
enum BoundaryType {BOUNDARY, INTERFACE, MIXED_BOUNDARY_INTERFACE, INTER_CONNECT};


/**
 *  Boundary Condition Type
 */
enum BCType
{

  /**
   * Neumann Boundary is the most general boundary,
   * which has ZERO flux through this boundary
   */
  NeumannBoundary             = 0x0001,

  /**
   * The interface of Semiconductor region to Vacuum region.
   * In most of the case, it is considered as Neumann Boundary
   */
  IF_Semiconductor_Vacuum     = 0x0002,

  /**
   * The interface of Insulator region to Vacuum region.
   * In most of the case, it is considered as Neumann Boundary
   */
  IF_Insulator_Vacuum         = 0x0003,

  /**
   * The interface of Electrode region to Vacuum region.
   * In most of the case, it is considered as Neumann Boundary
   */
  IF_Electrode_Vacuum         = 0x0004,

  /**
   * The interface of PML region to Vacuum/Insulator region.
    */
  IF_PML_Scatter              = 0x0005,

  /**
   * The interface of PML region to another PML region.
    */
  IF_PML_PML                  = 0x0006,

  /**
   * The interface of Electrode region to Insulator region.
   * we assume potential and temperature continuous on this boundary
   */
  IF_Electrode_Insulator      = 0x0011,

  /**
   * The interface of Semiconductor region to Insulator region.
   * this is important in MOS simulation.
   */
  IF_Insulator_Semiconductor  = 0x0012,

  /**
   * The interface of Insulator region to Insulator region.
   * Since some MOS has a nitride-oxide gate. we have to deal with this
   * boundary condition
   */
  IF_Insulator_Insulator      = 0x0013,

  /**
   * The interface of Electrode region to Electrode region.
   * nearly meaningless.
   */
  IF_Electrode_Electrode      = 0x0014,

  /**
   * The interface of Semiconductor region to semiconductor region with same material.
   * nearly meaningless.
   */
  HomoInterface               = 0x0015,

  /**
   * The interface of Semiconductor region to semiconductor region with different material.
   * which forms heterojunction.
   */
  HeteroInterface             = 0x0016,

  /**
   * The interface of Semiconductor region to Electrode region.
   * by default, it forms OhmicContact.
   */
  IF_Electrode_Semiconductor  = 0x0017,

  /**
   * Electrode which froms  Ohmic Contact.
   */
  OhmicContact        = 0x0101,

  /**
   * Electrode which froms  Schottky Contact.
   */
  SchottkyContact     = 0x0102,

  /**
   * Simple Boundary condition for MOS gate
   */
  SimpleGateContact   = 0x0103,

  /**
   * Boundary condition for MOS gate
   */
  GateContact         = 0x0104,

  /**
   * Boundary condition for float metal with charge,
   * useful for EEPROM simulation
   */
  ChargedContact      = 0x0105,

  /**
   * Boundary condition for inter connect of electrodes
   * useful for IC cell i.e. Inverter or SRAM simulation
   */
  InterConnect        = 0x0110,

  /**
   * Absorbing boundary for electromagnetic simulation
   */
  AbsorbingBoundary   = 0x1001,

  /**
   * boundary for adding light wave.
   */
  SourceBoundary      = 0x1002,


  INVALID_BC_TYPE     = 0xffff           // should always be last

};



/**
 *  @return enum BCType by string
 */
extern BCType BC_string_to_enum(const std::string & str);

/**
 *  @return string BCType by enum
 */
extern std::string BC_enum_to_string(BCType);


/**
 * @return the BCType of an interface boundary by two subdomain material type
 */
extern BCType determine_bc_by_subdomain(const std::string & mat1, const std::string mat2);



// predefine
class SimulationSystem;
class SimulationRegion;

/**
 * The base class of Boundary Condition
 */
class BoundaryCondition
{

public:

  /**
   * constructor
   */
  BoundaryCondition(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~BoundaryCondition()
  {
    _bd_nodes.clear();
    delete _ext_circuit;
  }

  /**
   * @return the const reference to label
   */
  const std::string & label() const
    { return _boundary_name; }

  /**
   * @return writable reference to label
   */
  std::string & label()
  { return _boundary_name; }


  /**
   * @return the electrode region label. only the "main" bc can own this value.
   * i.e. an electrode Anode has the ohmic contact with region Silicon,
   * then the ohmic bc can be load by Anode instead of its formal name Anode_to_Silicon
   */
  const std::string & electrode_label() const
    { return _electrode_name; }

  /**
   * @return writable reference to electrode region label
   */
  std::string & electrode_label()
  { return _electrode_name; }


  /**
   * add node* with this boundary type into vector _bd_nodes
   */
  void add_node(const Node * node)
  { _bd_nodes.push_back(node); }

  /**
   * set corresponding pointer to (region and FVM_Node) of a boundary Node
   */
  void insert(const Node * node, SimulationRegion *, FVM_Node *);

  /**
   * @return const reference to boundary nodes vector
   */
  const std::vector<const Node *> & nodes() const
    { return _bd_nodes; }

  /**
   * @return boundary nodes number
   */
  unsigned int  n_nodes() const
    { return _bd_nodes.size(); }



  /**
   * add elem with this boundary type into vector _bd_elems
   */
  void add_elem(const Elem * elem, unsigned int side)
  { _bd_elems.push_back(std::make_pair(elem, side)); }

  /**
   * @return boundary elements number
   */
  unsigned int  n_elems() const
  { return _bd_elems.size(); }

  /**
   * get nth element and its boundary face index
   */
  std::pair<const Elem *, unsigned int> get_elem(unsigned int n) const
  { return  _bd_elems[n]; }




  typedef std::vector<const Node *>::const_iterator const_node_iterator;


  /**
   * node default const begin() accessor
   */
  const_node_iterator nodes_begin        () const
  {
    return _bd_nodes.begin();
  }

  /**
   * node default const end() accessor
   */
  const_node_iterator nodes_end          () const
  {
    return _bd_nodes.end();
  }

  /**
   * @return the number of FVM nodes with Node n as its root_node
   */
  unsigned int n_region_node_with_root_node(const Node * n) const
    { return (*_bd_fvm_nodes.find(n)).second.size(); }

  /**
   * @return true if the node on external boundary
   */
  bool node_on_boundary(const Node * n) const
    { return n_region_node_with_root_node(n) == 1; }

  /**
   * @return true if the node on internal interface
   */
  bool node_on_interface(const Node * n) const
    { return n_region_node_with_root_node(n) > 1; }

  typedef std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> >::iterator region_node_iterator;

  /**
   * begin() accessor of all the (region and corresponding FVM_Node) of a Node
   */
  region_node_iterator region_node_begin( const Node * n )
  { return  (*_bd_fvm_nodes.find(n)).second.begin(); }

  /**
   * end() accessor of all the (region and corresponding FVM_Node) of a Node
   */
  region_node_iterator region_node_end( const Node * n )
  { return  (*_bd_fvm_nodes.find(n)).second.end(); }


  typedef std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> >::reverse_iterator  region_node_reverse_iterator;

  /**
   * rbegin() accessor of all the (region and corresponding FVM_Node) of a Node
   */
  region_node_reverse_iterator region_node_rbegin( const Node * n )
  { return  (*_bd_fvm_nodes.find(n)).second.rbegin(); }

  /**
   * rend() accessor of all the (region and corresponding FVM_Node) of a Node
   */
  region_node_reverse_iterator region_node_rend( const Node * n )
  { return  (*_bd_fvm_nodes.find(n)).second.rend(); }

  typedef std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> >::const_iterator const_region_node_iterator;

  /**
   * const begin() accessor of all the (region and corresponding FVM_Node) of a Node
   */
  const_region_node_iterator region_node_begin( const Node * n ) const
    { return  (*_bd_fvm_nodes.find(n)).second.begin(); }

  /**
   * const end() accessor of all the (region and corresponding FVM_Node) of a Node
   */
  const_region_node_iterator region_node_end( const Node * n ) const
    { return  (*_bd_fvm_nodes.find(n)).second.end(); }


  typedef std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> >::const_reverse_iterator  const_region_node_reverse_iterator;

  /**
   * rbegin() accessor of all the (region and corresponding FVM_Node) of a Node
   */
  const_region_node_reverse_iterator const_region_node_rbegin( const Node * n ) const
    { return  (*_bd_fvm_nodes.find(n)).second.rbegin(); }

  /**
   * rend() accessor of all the (region and corresponding FVM_Node) of a Node
   */
  const_region_node_reverse_iterator const_region_node_rend( const Node * n ) const
    { return  (*_bd_fvm_nodes.find(n)).second.rend(); }


  /**
   * find the FVM_Node by Node and its region pointer.
   */
  FVM_Node * get_region_fvm_node(const Node * n, SimulationRegion * region)
  {
    region_node_iterator reg_it     = region_node_begin(n);
    region_node_iterator reg_it_end = region_node_end(n);
    for(; reg_it!=reg_it_end; ++reg_it )
      if( (*reg_it).second.first == region )
        return (*reg_it).second.second;

    return NULL;
  }

  /**
   * @return true if boundary node associated with a region with specified SimulationRegionType
   */
  bool has_associated_region(const Node * n, SimulationRegionType rt) const
  {
    typedef std::map<const Node *, std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> > >::const_iterator It;
    It it = _bd_fvm_nodes.find(n);
    return it->second.find(rt)!=it->second.end();
  }

  /**
   * find the FVM_Node by Node and its region type. the region type should be unique in this multimap!
   */
  FVM_Node * get_region_fvm_node(const Node * n, SimulationRegionType type)
  {
    genius_assert( (*_bd_fvm_nodes.find(n)).second.count(type) == 1 );
    return (*(*_bd_fvm_nodes.find(n)).second.find(type)).second.second;
  }

  /**
   * find the SimulationRegion by Node and its region type. the region type should be unique in this multimap!
   */
  SimulationRegion * get_fvm_node_region(const Node * n, SimulationRegionType type)
  {
    genius_assert( (*_bd_fvm_nodes.find(n)).second.count(type) == 1 );
    return (*(*_bd_fvm_nodes.find(n)).second.find(type)).second.first;
  }

  /**
   * @return the node neighbor number
   * @Note node n must be on processor node
   */
  unsigned int n_node_neighbors(const Node * n) const
  {
    assert (n->processor_id()==Genius::processor_id());

    std::set<const Node *> node_set;

    const_region_node_iterator  reg_it     = region_node_begin(n);
    const_region_node_iterator  reg_it_end = region_node_end(n);
    for(; reg_it!=reg_it_end; ++reg_it )
    {
      const FVM_Node * fvm_node = (*reg_it).second.second;
      FVM_Node::fvm_neighbor_node_iterator it = fvm_node->neighbor_node_begin();
      FVM_Node::fvm_neighbor_node_iterator it_end = fvm_node->neighbor_node_end();
      for(; it!=it_end; ++it)
        node_set.insert((*it).first);
    }
    return  node_set.size();
  }

  /**
   * @return boundary type, one of BOUNDARY, INTERFACE, MIXED_BOUNDARY_INTERFACE and INTER_CONNECT
   */
  virtual BoundaryType boundary_type() const=0;


  /**
   * set boundary type
   */
  virtual void set_boundary_type(BoundaryType )
  { genius_error(); }


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const=0;

  /**
   * @return reference to system
   */
  SimulationSystem    & system()
  { return _system; }


  /**
   * @return the temperature of external entironment.
   */
  virtual PetscScalar T_external() const
    {return _T_Ext;}

  /**
   * @return writable reference to temperature of external entironment
   * however, we should never reach here
   */
  virtual PetscScalar & T_external()
  {return _T_Ext;}

  /**
   * @return a flag to show whether a boundary is full reflection
   * default is false
   */
  virtual bool reflection() const
  { return false; }

  /**
   * @return writable reference to a reflection flag
   * however, we should never reach here
   */
  virtual bool & reflection()
  { genius_error(); return _bool_dummy_;}

  /**
   * @return the heat transfer rate of this boundary
   * however, we should never reach here
   */
  virtual PetscScalar Heat_Transfer() const
    {genius_error(); return _dummy_;}

  /**
   * @return writable reference to heat transfer rate of this boundary
   * however, we should never reach here
   */
  virtual PetscScalar & Heat_Transfer()
  {genius_error(); return _dummy_;}

  /**
   * @return the work function of electrode material
   * however, we should never reach here
   */
  virtual PetscScalar Work_Function() const
    {genius_error(); return _dummy_;}

  /**
   * @return writable reference to work function of electrode material
   * however, we should never reach here
   */
  virtual PetscScalar & Work_Function()
  {genius_error(); return _dummy_;}


  /**
   * @return the thichness of gate material
   * however, we should never reach here
   */
  virtual PetscScalar Thickness() const
    {genius_error(); return _dummy_;}

  /**
   * @return writable reference to thichness of gate material
   * however, we should never reach here
   */
  virtual PetscScalar & Thickness()
  {genius_error(); return _dummy_;}

  /**
   * @return the electric constant of gate material
   * however, we should never reach here
   */
  virtual PetscScalar eps() const
    {genius_error(); return _dummy_;}

  /**
   * @return writable reference to electric constant of gate material
   * however, we should never reach here
   */
  virtual PetscScalar & eps()
  {genius_error(); return _dummy_;}

  /**
   * @return the free charge density.
   * @note it has differente meaning in different BCs
   * however, we should never reach here
   */
  virtual PetscScalar Qf() const
    {genius_error(); return _dummy_;}

  /**
   * @return writable reference to free charge density
   * @note it has differente meaning in different BCs
   * however, we should never reach here
   */
  virtual PetscScalar & Qf()
  {genius_error(); return _dummy_;}


  /**
   * @return the psi of float metal.
   * however, we should never reach here
   */
  virtual PetscScalar psi() const
    {genius_error(); return _dummy_;}

  /**
   * @return writable reference to psi of float metal
   * however, we should never reach here
   */
  virtual PetscScalar & psi()
  {genius_error(); return _dummy_;}

  /**
   * @return true iff this boundary is an electrode
   */
  virtual bool is_electrode() const=0;

  /**
   * @return the width in z direction
   * for 2D mesh, z_width is the device dimension in Z direction;
   * for 3D mesh, z_width is always 1.0
   */
  virtual PetscScalar z_width() const
  { return _z_width;}

  /**
   * @return the writable reference to width in z direction
   * for 2D mesh, z_width is the device dimension in Z direction;
   * for 3D mesh, z_width is always 1.0
   */
  virtual PetscScalar & z_width()
  { return _z_width;}


  /**
   * @return true if this electrode belongs to inter-connect layer
   */
  bool is_inter_connect_electrode() const
  {
    genius_assert(this->is_electrode());
    return !_inter_connect.empty();
  }

  /**
   * @return const reference of electrodes belongs to this inter_connect
   */
  const std::vector<BoundaryCondition * > & inter_connect() const
  { return _inter_connect; }

  /**
   * @return reference of electrodes belongs to this inter_connect
   */
  std::vector<BoundaryCondition * > & inter_connect()
  { return _inter_connect; }

  /**
   * set inter_connect electrodes
   */
  void set_inter_connect(const std::set<BoundaryCondition * > & bcs)
  {
     std::set<BoundaryCondition * >::const_iterator it;
     for(it=bcs.begin(); it!=bcs.end(); ++it)
       _inter_connect.push_back(*it);
  }

  /**
   * @return true if this bc is the hub of inter-connect layer
   */
  virtual bool is_inter_connect_hub() const
  { return false; }

  /**
   * @return pointer to inter_connect_hub
   */
  BoundaryCondition * inter_connect_hub()
  {return _inter_connect_hub;}

  /**
   * set inter_connect_hub
   */
  void set_inter_connect_hub(BoundaryCondition * hub)
  { _inter_connect_hub = hub; }


  /**
   * @return the offset of nodal solution data in global petsc vector
   */
  unsigned int global_offset () const
    { return _global_offset; }

  /**
   * function for set global offset
   */
  void set_global_offset (unsigned int pos )
  { _global_offset = pos; }


  /**
   * @return the offset of nodal solution data in local vector
   */
  unsigned int local_offset () const
    { return _local_offset; }

  /**
   * function for set local offset
   */
  void set_local_offset (unsigned int pos )
  { _local_offset = pos; }

  /**
   * @return the offset of nodal solution data in array derived from VecGetArray
   */
  unsigned int array_offset () const
  { return _array_offset; }

  /**
   * function for set array offset
   */
  void set_array_offset (unsigned int pos )
  { _array_offset = pos; }


  /**
   * let this bc hold a pointer of External Circuit
   */
  void build_ext_circuit(ExternalCircuit * ckt)
  { _ext_circuit = ckt;};


  /**
   * @return the const pointer of External Circuit
   */
  const ExternalCircuit * ext_circuit() const
  {
    return  _ext_circuit;
  }

  /**
   * @return the writable pointer of External Circuit
   */
  ExternalCircuit * ext_circuit()
  {
    return  _ext_circuit;
  }

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const
  { return ""; }

private:

  /**
   * the nodes this boundary/interface has
   * @Note: please make sure the nodes are sorted by their id
   */
  std::vector<const Node *> _bd_nodes;

  /**
   * the element and corresponding side belong to this boundary
   */
  std::vector<std::pair<const Elem *, unsigned int> > _bd_elems;

  /**
   * the glpbal node to region node map. the regions are sorted by SimulationRegionType
   */
  std::map<const Node *, std::multimap<SimulationRegionType, std::pair<SimulationRegion *, FVM_Node *> > > _bd_fvm_nodes;

  /**
   * the reference to corresponding SimulationSystem
   * since bc contains physical equations, it is important
   * for bc having the ability to access region information
   */
  SimulationSystem    & _system;

  /**
   * the boundary name given by user
   */
  std::string _boundary_name;

  /**
   * the electrode region name, which can be used to specify the
   * electrode boundary
   */
  std::string _electrode_name;

  /**
   * pointer to External Circuit, only electrode owns this data
   */
  ExternalCircuit * _ext_circuit;

  /**
   * the width in z direction
   * for 2D mesh, z_width is the device dimension in Z direction;
   * for 3D mesh, z_width is always 1.0
   * @Note system level also has z.width variable.
   * However, it will be override by boundary level z.width
   */
  PetscScalar _z_width;


  /**
   * temperature of external entironment
   */
  PetscScalar   _T_Ext;

  /**
   * An inter-connect layer of IC, it can connect several electrodes
   * every electrodes belongs to this inter-connect layer owns the same _inter_connect structure
   */
  std::vector<BoundaryCondition * > _inter_connect;

  /**
   * pointer to _inter_connect_hub
   * every electrodes belongs to this inter-connect layer owns this pointer
   */
  BoundaryCondition * _inter_connect_hub;

  /**
   * the offset of bc equation in global petsc vector
   * this value should be set on each processor
   * different solver may set this variable with different value
   */
  unsigned int _global_offset;

  /**
   * the offset of bc equation in local petsc vector
   * (scattered from global petsc vector with VecScatterBegin/VecScatterEnd)
   * this value should be set on each processor
   * different solver may set this variable with different value
   */
  unsigned int _local_offset;

  /**
   * some times, we need to access bc variable in a \p global vector
   * with array generated by VecGetArray/VecRestoreArray
   */
  unsigned int _array_offset;

private:

  /**
   * dummy to prevent compile problem
   */
  static PetscScalar _dummy_;

  /**
   * the same purpose with dummy
   */
   static bool _bool_dummy_;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for fill vector of poisson's equation.
   *
   * filling solution data into petsc vector of poisson's equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   *
   * @param Vec              global solution vector
   * @param Vec              the left scaling vector
   * @note each derived boundary condition can override it
   */
  virtual void Poissin_Fill_Value(Vec , Vec ) {}


  /**
   * @brief virtual function for reserve none zero pattern in petsc matrix.
   *
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition can override it
   */
  virtual void Poissin_Jacobian_Reserve(Mat *, InsertMode &) {}


  /**
   * @brief virtual function for evaluating poisson's equation.
   *
   * @param PetscScalar*    local unknown vector
   * @param Vec             petsc global function vector
   * @param InsertMode&     flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of poisson's equation.
   *
   * @param PetscScalar*     local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void Poissin_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for update solution value of poisson's equation.
   *
   *
   * @param PetscScalar*     global solution vector
   *
   * @note do nothing, each derived boundary condition can override it
   */
  virtual void Poissin_Update_Solution(PetscScalar *) {}



  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////


  /**
   * @brief virtual function for fill vector of level 1 DDM equation.
   *
   * filling solution data into petsc vector of level 1 DDM equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   *
   * @param Vec              global solution vector
   * @param Vec              the left scaling vector
   * @note each derived boundary condition can override it
   */
  virtual void DDM1_Fill_Value(Vec , Vec ) {}

  /**
   * @brief virtual function for reserve none zero pattern in petsc matrix.
   *
   * @param Mat*              petsc global jacobian matrix
   * @param InsertMode&       flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition can override it
   */
  virtual void DDM1_Jacobian_Reserve(Mat *, InsertMode &) {}


  /**
   * @brief virtual function for evaluating level 1 DDM equation.
   *
   * @param PetscScalar*    local unknown vector
   * @param Vec             petsc global function vector
   * @param InsertMode&     flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of level 1 DDM equation.
   *
   * @param PetscScalar*     local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating trace parameter of level 1 DDM equation.
   *
   * @param Vec              local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param Vec              vector for dI/dx
   * @param Vec              vector for dF/dV
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void DDM1_Electrode_Trace(Vec, Mat *, Vec , Vec) {}

  /**
   * @brief virtual function for update solution value of level 1 DDM equation.
   *
   * @param PetscScalar*     global solution vector
   *
   * @note each derived boundary condition can override it
   */
  virtual void DDM1_Update_Solution(PetscScalar *) {}



  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed DDML1------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for evaluating Mixed type level 1 DDM equation.
   *
   * @param PetscScalar*    local unknown vector
   * @param Vec             petsc global function vector
   * @param InsertMode&     flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_DDM1_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { this->DDM1_Function(x , f, add_value_flag ); }

  /**
   * @brief virtual function for reserve none zero pattern in petsc matrix.
   *
   * @param Mat *              petsc global jacobian matrix
   * @param InsertMode&        flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_DDM1_Jacobian_Reserve(Mat * jac, InsertMode &add_value_flag)
  { this->DDM1_Jacobian_Reserve(jac, add_value_flag ); }

  /**
   * @brief virtual function for evaluating Mixed type Jacobian of level 1 DDM equation.
   *
   * @param PetscScalar*     local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_DDM1_Jacobian(PetscScalar * x, Mat * jac, InsertMode &add_value_flag)
  { this->DDM1_Jacobian(x, jac, add_value_flag); }


  /**
   * @brief virtual function for evaluating electrode Load of Mixed simulation.
   *
   * @param Vec              petsc global solution vector
   * @param Mat*             petsc global jacobian matrix
   * @param double&          electrode current
   * @param PetscScalar&     pdI_pdV of electrode
   * @param Vec&             pdI_pdx of electrode
   * @param Vec&             pdf(x)_pdV of electrode
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_DDM1_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &) {}


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed DDML1-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for evaluating Advanced Mixed type level 1 DDM equation.
   *
   * @param PetscScalar*    local unknown vector
   * @param Vec             petsc global function vector
   * @param InsertMode&     flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void MixA_DDM1_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { this->DDM1_Function(x , f, add_value_flag ); }

  /**
   * @brief virtual function for reserve none zero pattern in petsc matrix.
   *
   * @param Mat *              petsc global jacobian matrix
   * @param InsertMode&        flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void MixA_DDM1_Jacobian_Reserve(Mat * jac, InsertMode &add_value_flag)
  { this->DDM1_Jacobian_Reserve(jac, add_value_flag ); }

  /**
   * @brief virtual function for evaluating Advanced Mixed type Jacobian of level 1 DDM equation.
   *
   * @param PetscScalar*     local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void MixA_DDM1_Jacobian(PetscScalar * x, Mat * jac, InsertMode &add_value_flag)
  { this->DDM1_Jacobian(x, jac, add_value_flag); }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////


  /**
   * @brief virtual function for fill vector of level 2 DDM equation.
   *
   * filling solution data into petsc vector of level 2 DDM equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   *
   * @param Vec              global solution vector
   * @param Vec              the left scaling vector
   * @note each derived boundary condition can override it
   */
  virtual void DDM2_Fill_Value(Vec , Vec ) {}

  /**
   * @brief virtual function for reserve none zero pattern in petsc matrix.
   *
   * @param Mat *              petsc global jacobian matrix
   * @param InsertMode&        flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition can override it
   */
  virtual void DDM2_Jacobian_Reserve(Mat *, InsertMode &) {}


  /**
   * @brief virtual function for evaluating level 2 DDM equation.
   *
   * @param PetscScalar*    local unknown vector
   * @param Vec             petsc global function vector
   * @param InsertMode&     flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating Jacobian of level 2 DDM equation.
   *
   * @param PetscScalar*     local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)=0;

  /**
   * @brief virtual function for evaluating trace parameter of level 2 DDM equation.
   *
   * @param Vec              local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param Vec              vector for dI/dx
   * @param Vec              vector for dF/dV
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void DDM2_Electrode_Trace(Vec, Mat *, Vec , Vec) {}

  /**
   * @brief virtual function for update solution value of level 2 DDM equation.
   *
   * @param PetscScalar*     global solution vector
   *
   * @note each derived boundary condition can override it
   */
  virtual void DDM2_Update_Solution(PetscScalar *) {}

  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed DDML2------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for evaluating Mixed type level 2 DDM equation.
   *
   * @param PetscScalar*    local unknown vector
   * @param Vec             petsc global function vector
   * @param InsertMode&     flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_DDM2_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { this->DDM2_Function(x , f, add_value_flag ); }

  /**
   * @brief virtual function for reserve none zero pattern in petsc matrix.
   *
   * @param Mat *              petsc global jacobian matrix
   * @param InsertMode&        flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_DDM2_Jacobian_Reserve(Mat * jac, InsertMode &add_value_flag)
  { this->DDM2_Jacobian_Reserve(jac, add_value_flag ); }

  /**
   * @brief virtual function for evaluating Mixed type Jacobian of level 2 DDM equation.
   *
   * @param PetscScalar*     local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_DDM2_Jacobian(PetscScalar * x, Mat * jac, InsertMode &add_value_flag)
  { this->DDM2_Jacobian(x, jac, add_value_flag); }


  /**
   * @brief virtual function for evaluating electrode Load of Mixed simulation.
   *
   * @param Vec              petsc global solution vector
   * @param Mat*             petsc global jacobian matrix
   * @param double&          electrode current
   * @param PetscScalar&     pdI_pdV of electrode
   * @param Vec&             pdI_pdx of electrode
   * @param Vec&             pdf(x)_pdV of electrode
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_DDM2_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &) {}


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed DDML2-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for evaluating Advanced Mixed type level 2 DDM equation.
   *
   * @param PetscScalar*    local unknown vector
   * @param Vec             petsc global function vector
   * @param InsertMode&     flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void MixA_DDM2_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { this->DDM2_Function(x , f, add_value_flag ); }

  /**
   * @brief virtual function for reserve none zero pattern in petsc matrix.
   *
   * @param Mat *              petsc global jacobian matrix
   * @param InsertMode&        flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void MixA_DDM2_Jacobian_Reserve(Mat * jac, InsertMode &add_value_flag)
  { this->DDM2_Jacobian_Reserve(jac, add_value_flag ); }

  /**
   * @brief virtual function for evaluating Advanced Mixed type Jacobian of level 2 DDM equation.
   *
   * @param PetscScalar*     local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void MixA_DDM2_Jacobian(PetscScalar * x, Mat * jac, InsertMode &add_value_flag)
  { this->DDM2_Jacobian(x, jac, add_value_flag); }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Hybrid EBM-----------------//
  //////////////////////////////////////////////////////////////////////////////////


  /**
   * @brief virtual function for fill vector of level 3 EBM equation.
   *
   * filling solution data into petsc vector of level 3 EBM equation.
   * can be used as initial data of nonlinear equation or diverged recovery.
   *
   * @param Vec              global solution vector
   * @param Vec              the left scaling vector
   * @note each derived boundary condition can override it
   */
  virtual void EBM3_Fill_Value(Vec , Vec ) {}

  /**
   * @brief virtual function for reserve none zero pattern in petsc matrix.
   *
   * @param Mat *              petsc global jacobian matrix
   * @param InsertMode&        flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition can override it
   */
  virtual void EBM3_Jacobian_Reserve(Mat *, InsertMode &) {}


  /**
   * @brief virtual function for evaluating level 3 EBM equation.
   *
   * @param PetscScalar*    local unknown vector
   * @param Vec             petsc global function vector
   * @param InsertMode&     flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void EBM3_Function(PetscScalar *, Vec , InsertMode &) {}

  /**
   * @brief virtual function for evaluating Jacobian of level 3 EBM equation.
   *
   * @param PetscScalar*     local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &) {}

  /**
   * @brief virtual function for evaluating trace parameter of level 3 EBM equation.
   *
   * @param Vec              local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param Vec              vector for dI/dx
   * @param Vec              vector for dF/dV
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note each derived boundary condition should override it
   */
  virtual void EBM3_Electrode_Trace(Vec, Mat *, Vec , Vec) {}

  /**
   * @brief virtual function for update solution value of level 3 EBM equation.
   *
   * @param PetscScalar*     global solution vector
   *
   * @note each derived boundary condition can override it
   */
  virtual void EBM3_Update_Solution(PetscScalar *) {}


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed EBM3 ------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for evaluating Mixed type level 3 EBM equation.
   *
   * @param PetscScalar*    local unknown vector
   * @param Vec             petsc global function vector
   * @param InsertMode&     flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_EBM3_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { this->EBM3_Function(x , f, add_value_flag ); }

  /**
   * @brief virtual function for reserve none zero pattern in petsc matrix.
   *
   * @param Mat *              petsc global jacobian matrix
   * @param InsertMode&        flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_EBM3_Jacobian_Reserve(Mat * jac, InsertMode &add_value_flag)
  { this->EBM3_Jacobian_Reserve(jac, add_value_flag ); }

  /**
   * @brief virtual function for evaluating Mixed type Jacobian of level 2 DDM equation.
   *
   * @param PetscScalar*     local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_EBM3_Jacobian(PetscScalar * x, Mat * jac, InsertMode &add_value_flag)
  { this->EBM3_Jacobian(x, jac, add_value_flag); }


  /**
   * @brief virtual function for evaluating electrode Load of Mixed simulation.
   *
   * @param Vec              petsc global solution vector
   * @param Mat*             petsc global jacobian matrix
   * @param double&          electrode current
   * @param PetscScalar&     pdI_pdV of electrode
   * @param Vec&             pdI_pdx of electrode
   * @param Vec&             pdf(x)_pdV of electrode
   *
   * @note only electrode boundary need to override it
   */
  virtual void Mix_EBM3_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &) {}


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed EBM3 -------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for evaluating Advanced Mixed type level 3 EBM equation.
   *
   * @param PetscScalar*    local unknown vector
   * @param Vec             petsc global function vector
   * @param InsertMode&     flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void MixA_EBM3_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { this->EBM3_Function(x , f, add_value_flag ); }

  /**
   * @brief virtual function for reserve none zero pattern in petsc matrix.
   *
   * @param Mat *              petsc global jacobian matrix
   * @param InsertMode&        flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void MixA_EBM3_Jacobian_Reserve(Mat * jac, InsertMode &add_value_flag)
  { this->EBM3_Jacobian_Reserve(jac, add_value_flag ); }

  /**
   * @brief virtual function for evaluating Advanced Mixed type Jacobian of level 3 EBM equation.
   *
   * @param PetscScalar*     local unknown vector
   * @param Mat*             petsc global jacobian matrix
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void MixA_EBM3_Jacobian(PetscScalar * x, Mat * jac, InsertMode &add_value_flag)
  { this->EBM3_Jacobian(x, jac, add_value_flag); }


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * @brief virtual function for evaluating matrix and rhs vector for ddm ac solver
   *
   * @param Mat              petsc AC matrix
   * @param Vec              rhs vector
   * @param Mat              petsc global jacobian matrix
   * @param double           AC frequency
   * @param InsertMode&      flag for last operator is ADD_VALUES
   *
   * @note only electrode boundary need to override it
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat,  Vec, const Mat, const double, InsertMode &  ) {}

  /**
   * @brief virtual function for update solution value for ddm ac solver
   *
   * @param PetscScalar*     local solution vector
   * @param Mat              petsc global jacobian matrix
   *
   * @note each derived boundary condition can override it
   */
  virtual void DDMAC_Update_Solution(const PetscScalar * , const Mat, const double) {}}
;




/**
 * The inter connect Boundary Condition
 */
class InterConnectBC : public BoundaryCondition
{
public:

  /**
   * constructor, set default value
   */
  InterConnectBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~InterConnectBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return InterConnect; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTER_CONNECT; }

  /**
   * indicate that this bc is an electrode
   * @return false
   */
  virtual bool is_electrode() const
    {return true;}

  /**
   * @return true if this bc is the hub of inter-connect layer
   */
  virtual bool is_inter_connect_hub() const
  { return true; }

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * fill initial value to inter-connect node
   */
  void inter_connect_fill_value(Vec x, Vec L);

  /**
   * set governing equation for inter-connect node. use nodal analysis method
   */
  void inter_connect_function(PetscScalar *x , Vec f, InsertMode &add_value_flag);

  /**
   * reserve matrix entries
   */
  void inter_connect_reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * set jacobian matrix entries
   */
  void inter_connect_jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag);

  /**
   * update solution data
   */
  void inter_connect_update_solution(PetscScalar *x);


public:
  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &)
  { /* no thing to do for poisson solver */ }

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &)
  { /* no thing to do for poisson solver */ }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of level 1 DDM equation.
   */
  virtual void DDM1_Fill_Value(Vec x, Vec L)
  { inter_connect_fill_value(x, L); }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { inter_connect_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { inter_connect_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { inter_connect_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *x)
  { inter_connect_update_solution(x); }

  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of level 2 DDM equation.
   */
  virtual void DDM2_Fill_Value(Vec x, Vec L)
  { inter_connect_fill_value(x, L); }

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { inter_connect_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { inter_connect_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { inter_connect_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of level 2 DDM solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *x)
  { inter_connect_update_solution(x); }

  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of level 3 EBM equation.
   */
  virtual void EBM3_Fill_Value(Vec x, Vec L)
  { inter_connect_fill_value(x, L); }

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
  { inter_connect_function(x, f, add_value_flag); }

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
  { inter_connect_reserve(jac, add_value_flag); }

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
  { inter_connect_jacobian(x, jac, add_value_flag); }

  /**
   * update solution data of EBM3 solver.
   */
  virtual void EBM3_Update_Solution(PetscScalar *x)
  { inter_connect_update_solution(x); }

};



/**
 * The Neumann Boundary Condition
 */
class NeumannBC : public BoundaryCondition
{
public:

  /**
   * constructor, set default value
   */
  NeumannBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~NeumannBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return NeumannBoundary; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return BOUNDARY; }

  /**
   * @return the heat transfer rate of this boundary
   */
  virtual PetscScalar Heat_Transfer() const
    {return _Heat_Transfer;}

  /**
   * @return writable reference to heat transfer rate of this boundary
   */
  virtual PetscScalar & Heat_Transfer()
  {return _Heat_Transfer;}

  /**
   * @return the reflection flag of this boundary
   */
  virtual bool reflection() const
  {return _reflection;}

  /**
   * @return writable reference to reflection flag of this boundary
   */
  virtual bool & reflection()
  {return _reflection;}

  /**
   * indicate that this bc is an electrode
   * @return false
   */
  virtual bool is_electrode() const
    {return false;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * heat transfer rate of this boundary
   */
  PetscScalar _Heat_Transfer;

  /**
   * full reflection flag......just used in ray tracing
   */
  bool _reflection;


public:

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
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &)
  { /* no thing to do for Neumann boundary */ }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &)
  { /* no thing to do for Neumann boundary */ }

  /**
   * build function and its jacobian for level 1 DDM solver, nothing to do
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &)
  { /* no thing to do for Neumann boundary */ }


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);


  /**
   * build function and its jacobian for level 2 DDM solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);



  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &);


  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );

};



/**
 * The Ohmic Boundary Condition the metal contact semiconductor to
 * form an Ohmic type electrode
 */
class OhmicContactBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  OhmicContactBC(SimulationSystem  & system, const std::string & label="", bool _interface=false);

  /**
   * destructor
   */
  virtual ~OhmicContactBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return OhmicContact; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return _bd_type; }

  /**
   * set boundary type
   */
  virtual void set_boundary_type(BoundaryType type)
  { _bd_type = type; }

  /**
   * @return the heat transfer rate of this boundary
   */
  virtual PetscScalar Heat_Transfer() const
    {return _Heat_Transfer;}

  /**
   * @return writable reference to heat transfer rate of this boundary
   */
  virtual PetscScalar & Heat_Transfer()
  {return _Heat_Transfer;}

  /**
   * @return the reflection flag of this boundary
   */
  virtual bool reflection() const
  {return true;}

  /**
   * @return writable reference to reflection flag of this boundary
   */
  virtual bool & reflection()
  {return _reflection;}

  /**
   * @return true
   */
  virtual bool is_electrode()  const
    {return true;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:
  /**
   * heat transfer rate of this boundary
   */
  PetscScalar   _Heat_Transfer;

  /**
   * the boundary type, can be interface or boundary.
   */
  BoundaryType  _bd_type;

  /**
   * full reflection flag......just used in ray tracing
   */
  bool _reflection;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

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

  /*
   * fill solution data into petsc vector of level 1 DDM equation.
   */
  virtual void DDM1_Fill_Value(Vec x, Vec L);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * build trace parameter for DDML1 solver
   */
  virtual void DDM1_Electrode_Trace(Vec, Mat *, Vec , Vec);

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed DDML1------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Mixed DDML1 solver
   */
  virtual void Mix_DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Mix_DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Mixed DDML1 solver
   */
  virtual void Mix_DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * evaluating electrode Load of Mixed simulation
   */
  virtual void Mix_DDM1_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed DDML1-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Advanced Mixed DDML1 solver
   */
  virtual void MixA_DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Advanced Mixed DDML1 solver
   */
  virtual void MixA_DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of level 2 DDM equation.
   */
  virtual void DDM2_Fill_Value(Vec x, Vec L);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * build trace parameter for DDML2 solver
   */
  virtual void DDM2_Electrode_Trace(Vec, Mat *, Vec , Vec);

  /**
   * update solution data of DDML2 solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed DDML2------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Mixed DDML2 solver
   */
  virtual void Mix_DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Mix_DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Mixed DDML2 solver
   */
  virtual void Mix_DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * evaluating electrode Load of Mixed simulation
   */
  virtual void Mix_DDM2_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed DDML2-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Advanced Mixed DDML2 solver
   */
  virtual void MixA_DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Advanced Mixed DDML2 solver
   */
  virtual void MixA_DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of level 3 EBM equation.
   */
  virtual void EBM3_Fill_Value(Vec x, Vec L);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * build trace parameter for level 3 EBM solver
   */
  virtual void EBM3_Electrode_Trace(Vec, Mat *, Vec , Vec);

  /**
   * update solution data of level 3 EBM solver.
   */
  virtual void EBM3_Update_Solution(PetscScalar *);



  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed EBM3 ------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Mixed EBM3 solver
   */
  virtual void Mix_EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Mix_EBM3_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Mixed DDML2 solver
   */
  virtual void Mix_EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * evaluating electrode Load of Mixed simulation
   */
  virtual void Mix_EBM3_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed EBM3-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Advanced Mixed EBM3 solver
   */
  virtual void MixA_EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_EBM3_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Advanced Mixed EBM3 solver
   */
  virtual void MixA_EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );

  /**
   * update solution value for ddm ac solver
   */
  virtual void DDMAC_Update_Solution(const PetscScalar * lxx , const Mat, const double omega);

};



/**
 * The Schottky Boundary Condition, the metal contact semiconductor to
 * form a Schottky type electrode
 */
class SchottkyContactBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  SchottkyContactBC(SimulationSystem  & system, const std::string & label="", bool _interface=false);

  /**
   * destructor
   */
  virtual ~SchottkyContactBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return SchottkyContact; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return _bd_type; }

  /**
   * set boundary type
   */
  virtual void set_boundary_type(BoundaryType type)
  { _bd_type = type; }

  /**
   * @return the heat transfer rate of this boundary
   */
  virtual PetscScalar Heat_Transfer() const
    {return _Heat_Transfer;}

  /**
   * @return writable reference to heat transfer rate of this boundary
   */
  virtual PetscScalar & Heat_Transfer()
  {return _Heat_Transfer;}


  /**
   * @return the work function of electrode material
   */
  virtual PetscScalar Work_Function() const
    {return _WorkFunction;}

  /**
   * @return writable reference to work function of electrode material
   */
  virtual PetscScalar & Work_Function()
  {return _WorkFunction;}

  /**
   * @return the reflection flag of this boundary
   */
  virtual bool reflection() const
  {return true;}

  /**
   * @return writable reference to reflection flag of this boundary
   */
  virtual bool & reflection()
  {return _reflection;}

  /**
   * @return true
   */
  virtual bool is_electrode()  const
    {return true;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:
  /**
   * heat transfer rate of this boundary
   */
  PetscScalar   _Heat_Transfer;

  /**
   * workfunction of electrode material
   */
  PetscScalar   _WorkFunction;

  /**
   * the boundary type, can be interface or boundary.
   */
  BoundaryType  _bd_type;

  /**
   * full reflection flag......just used in ray tracing
   */
  bool _reflection;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

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

  /*
   * fill solution data into petsc vector of level 1 DDM equation.
   */
  virtual void DDM1_Fill_Value(Vec x, Vec L);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * build trace parameter for DDML1 solver
   */
  virtual void DDM1_Electrode_Trace(Vec, Mat *, Vec , Vec);

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *);



  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed DDML1------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Mixed DDML1 solver
   */
  virtual void Mix_DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Mix_DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Mixed DDML1 solver
   */
  virtual void Mix_DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * evaluating electrode Load of Mixed simulation
   */
  virtual void Mix_DDM1_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed DDML1-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Advanced Mixed DDML1 solver
   */
  virtual void MixA_DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Advanced Mixed DDML1 solver
   */
  virtual void MixA_DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of level 2 DDM equation.
   */
  virtual void DDM2_Fill_Value(Vec , Vec );

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * build trace parameter for DDML2 solver
   */
  virtual void DDM2_Electrode_Trace(Vec, Mat *, Vec , Vec);

  /**
   * update solution data of DDML2 solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed DDML2------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Mixed DDML2 solver
   */
  virtual void Mix_DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Mix_DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Mixed DDML2 solver
   */
  virtual void Mix_DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * evaluating electrode Load of Mixed simulation
   */
  virtual void Mix_DDM2_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed DDML2-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Advanced Mixed DDML2 solver
   */
  virtual void MixA_DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Advanced Mixed DDML2 solver
   */
  virtual void MixA_DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of level 3 EBM equation.
   */
  virtual void EBM3_Fill_Value(Vec x, Vec L);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * build trace parameter for level 3 EBM solver
   */
  virtual void EBM3_Electrode_Trace(Vec, Mat *, Vec , Vec);

  /**
   * update solution data of level 3 EBM solver.
   */
  virtual void EBM3_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed EBM3 ------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Mixed EBM3 solver
   */
  virtual void Mix_EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Mix_EBM3_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Mixed DDML2 solver
   */
  virtual void Mix_EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * evaluating electrode Load of Mixed simulation
   */
  virtual void Mix_EBM3_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed EBM3 -------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Advanced Mixed EBM3 solver
   */
  virtual void MixA_EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_EBM3_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Advanced Mixed EBM3 solver
   */
  virtual void MixA_EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );

  /**
   * update solution value for ddm ac solver
   */
  virtual void DDMAC_Update_Solution(const PetscScalar * lxx , const Mat, const double omega);

};



/**
 * The Gate Boundary Condition, it is the metal(or polySi) contact with insulator
 * to form the MOS gate
 */
class GateContactBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  GateContactBC(SimulationSystem  & system, const std::string & label="", bool _interface=false);

  /**
   * destructor
   */
  virtual ~GateContactBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return GateContact; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return _bd_type; }

  /**
   * set boundary type
   */
  virtual void set_boundary_type(BoundaryType type)
  { _bd_type = type; }

  /**
   * @return the heat transfer rate of this boundary
   */
  virtual PetscScalar Heat_Transfer() const
    {return _Heat_Transfer;}

  /**
   * @return writable reference to heat transfer rate of this boundary
   */
  virtual PetscScalar & Heat_Transfer()
  {return _Heat_Transfer;}

  /**
   * @return the work function of electrode material
   */
  virtual PetscScalar Work_Function() const
    {return _WorkFunction;}

  /**
   * @return writable reference to work function of electrode material
   */
  virtual PetscScalar & Work_Function()
  {return _WorkFunction;}

  /**
   * @return true
   */
  virtual bool is_electrode()  const
    {return true;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * heat transfer rate of this boundary
   */
  PetscScalar   _Heat_Transfer;

  /**
   * workfunction of electrode material
   */
  PetscScalar   _WorkFunction;

  /**
   * the boundary type, can be interface or boundary.
   */
  BoundaryType  _bd_type;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat * , InsertMode &);

  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

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

  /*
   * fill solution data into petsc vector of level 1 DDM equation.
   */
  virtual void DDM1_Fill_Value(Vec , Vec );

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat * , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed DDML1------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Mixed DDML1 solver
   */
  virtual void Mix_DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Mix_DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Mixed DDML1 solver
   */
  virtual void Mix_DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * evaluating electrode Load of Mixed simulation
   */
  virtual void Mix_DDM1_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed DDML1-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Advanced Mixed DDML1 solver
   */
  virtual void MixA_DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Advanced Mixed DDML1 solver
   */
  virtual void MixA_DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of level 2 DDM equation.
   */
  virtual void DDM2_Fill_Value(Vec , Vec );

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat * , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of DDML2 solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *);



  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed DDML2------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Mixed DDML2 solver
   */
  virtual void Mix_DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Mix_DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Mixed DDML2 solver
   */
  virtual void Mix_DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * evaluating electrode Load of Mixed simulation
   */
  virtual void Mix_DDM2_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed DDML2-------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Advanced Mixed DDML2 solver
   */
  virtual void MixA_DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Advanced Mixed DDML2 solver
   */
  virtual void MixA_DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L3 EBM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of level 3 EBM equation.
   */
  virtual void EBM3_Fill_Value(Vec x, Vec L);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of level 3 EBM solver.
   */
  virtual void EBM3_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Function and Jacobian evaluate for Mixed EBM3 ------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Mixed EBM3 solver
   */
  virtual void Mix_EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Mix_EBM3_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Mixed DDML2 solver
   */
  virtual void Mix_EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * evaluating electrode Load of Mixed simulation
   */
  virtual void Mix_EBM3_Electrode_Load(const Vec , const Mat *, double &, PetscScalar &, Vec &, Vec &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------Function and Jacobian evaluate for Advanced Mixed EBM3 -------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * build function and its jacobian for Advanced Mixed EBM3 solver
   */
  virtual void MixA_EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void MixA_EBM3_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for Advanced Mixed EBM3 solver
   */
  virtual void MixA_EBM3_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );

  /**
   * update solution value for ddm ac solver
   */
  virtual void DDMAC_Update_Solution(const PetscScalar * lxx , const Mat, const double omega);

};





/**
 * The Simple Gate Boundary Condition, no gate need to be constructed,
 * instead, the thickness and eps parameter should be given for the SiO2 layer
 */
class SimpleGateContactBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  SimpleGateContactBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~SimpleGateContactBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return SimpleGateContact; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return BOUNDARY; }

  /**
   * @return the heat transfer rate of this boundary
   */
  virtual PetscScalar Heat_Transfer() const
    {return _Heat_Transfer;}

  /**
   * @return writable reference to heat transfer rate of this boundary
   */
  virtual PetscScalar & Heat_Transfer()
  {return _Heat_Transfer;}


  /**
   * @return the work function of electrode material
   */
  virtual PetscScalar Work_Function() const
    {return _WorkFunction;}

  /**
   * @return writable reference to work function of electrode material
   */
  virtual PetscScalar & Work_Function()
  {return _WorkFunction;}


  /**
   * @return the thichness of gate material
   */
  virtual PetscScalar Thickness() const
    {return _Thickness;}

  /**
   * @return writable reference to thichness of gate material
   */
  virtual PetscScalar & Thickness()
  {return _Thickness;}

  /**
   * @return the electric constant of gate material
   */
  virtual PetscScalar eps() const
    { return _eps;}

  /**
   * @return writable reference to electric constant of gate material
   */
  virtual PetscScalar & eps()
  { return _eps;}

  /**
   * @return the free charge density at Si/SiO2 interface
   */
  virtual PetscScalar Qf() const
    {return _Qf;}

  /**
   * @return writable reference to free charge density
   */
  virtual PetscScalar & Qf()
  {return _Qf;}

  /**
   * @return the reflection flag of this boundary
   */
  virtual bool reflection() const
  {return true;}

  /**
   * @return writable reference to reflection flag of this boundary
   */
  virtual bool & reflection()
  {return _reflection;}

  /**
   * @return true
   */
  virtual bool is_electrode()  const
    {return true;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * heat transfer rate of this boundary
   */
  PetscScalar   _Heat_Transfer;


  /**
   * workfunction of electrode material
   */
  PetscScalar   _WorkFunction;

  /**
   * the thickness of SiO2
   */
  PetscScalar   _Thickness;

  /**
   * the electric constant of gate material (mostly SiO2)
   */
  PetscScalar   _eps;

  /**
   * the free charge density at Si/SiO2 interface
   */
  PetscScalar   _Qf;

  /**
   * full reflection flag......just used in ray tracing
   */
  bool _reflection;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////


  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

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

  /*
   * fill solution data into petsc vector of level 1 DDM equation.
   */
  virtual void DDM1_Fill_Value(Vec , Vec );

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of DDML2 equation.
   */
  virtual void DDM2_Fill_Value(Vec , Vec );

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of DDML2 solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for EBM   ---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /*
   * fill solution data into petsc vector of level 3 EBM equation.
   */
  virtual void EBM3_Fill_Value(Vec x, Vec L);

  /**
   * build function and its jacobian for level 3 EBM solver
   */
  virtual void EBM3_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void EBM3_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for level 3 EBM solver
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

  /**
   * update solution value for ddm ac solver
   */
  virtual void DDMAC_Update_Solution(const PetscScalar * lxx , const Mat, const double omega );

};




/**
 * The insulator-to-semiconductor interface
 */
class InsulatorSemiconductorInterfaceBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  InsulatorSemiconductorInterfaceBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~InsulatorSemiconductorInterfaceBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return IF_Insulator_Semiconductor; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTERFACE; }

  /**
   * @return the free charge density at Si/SiO2 interface
   */
  virtual PetscScalar Qf() const
    {return _Qf;}

  /**
   * @return writable reference to free charge density
   */
  virtual PetscScalar & Qf()
  {return _Qf;}

  /**
   * @return false
   */
  virtual bool is_electrode()  const
    {return false;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * InsulatorInterface fixed charge density
   */
  PetscScalar      _Qf;


public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat *, InsertMode &);

  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);

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

  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );

};




/**
 * The interface between same semiconductor material
 */
class HomoInterfaceBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  HomoInterfaceBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~HomoInterfaceBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return HomoInterface; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTERFACE; }

  /**
   * @return false
   */
  virtual bool is_electrode()  const
    {return false;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat * , InsertMode &);

  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);


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


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );

};



/**
 * The interface between different semiconductor material
 */
class HeteroInterfaceBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  HeteroInterfaceBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~HeteroInterfaceBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return HeteroInterface; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTERFACE; }

  /**
    * @return the fixed charge density on the interface
    */
  virtual PetscScalar Qf() const
    {return _Qf;}

  /**
   * @return writable reference to fixed charge density on the interface
   */
  virtual PetscScalar & Qf()
  {return _Qf;}

  /**
   * @return false
   */
  virtual bool is_electrode()  const
    {return false;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * fixed charge density on the interface
   */
  PetscScalar      _Qf;


public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat * , InsertMode &);

  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);


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


  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );

};



/**
 * The electrode-to-insulator interface
 */
class ElectrodeInsulatorInterfaceBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  ElectrodeInsulatorInterfaceBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~ElectrodeInsulatorInterfaceBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return IF_Electrode_Insulator; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTERFACE; }

  /**
   * @return false
   */
  virtual bool is_electrode()  const
    {return false;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);


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



  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );

};




/**
 * The insulator-to-insulator interface, i.e. nitride SiO2-Si3N4-SiO2 gate
 */
class InsulatorInsulatorInterfaceBC : public BoundaryCondition
{
public:

  /**
   * constructor
   */
  InsulatorInsulatorInterfaceBC(SimulationSystem  & system, const std::string & label="");

  /**
   * destructor
   */
  virtual ~InsulatorInsulatorInterfaceBC(){}


  /**
   * @return boundary condition type
   */
  virtual BCType bc_type() const
    { return IF_Insulator_Insulator; }


  /**
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTERFACE; }

  /**
   * @return false
   */
  virtual bool is_electrode()  const
    {return false;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat *, InsertMode &);

  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for poisson solver, nothing to do
   */
  virtual void Poissin_Jacobian(PetscScalar * , Mat *, InsertMode &);



  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L1 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);


  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat *, InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);


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



  //////////////////////////////////////////////////////////////////////////////////
  //--------------Matrix and RHS Vector evaluate for DDM AC Solver----------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   *  evaluating matrix and rhs vector for ddm ac solver
   */
  virtual void DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag );

};



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
   * @return boundary type
   */
  virtual BoundaryType boundary_type() const
    { return INTERFACE; }

  /**
   * @return the free charge on the surface of float metal
   */
  virtual PetscScalar Qf() const
    {return _Qf;}

  /**
   * @return writable reference to free charge
   */
  virtual PetscScalar & Qf()
  {return _Qf;}

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
   * @return false
   */
  virtual bool is_electrode()  const
    {return false;}

  /**
   * @return the string which indicates the boundary condition
   */
  virtual std::string boundary_condition_in_string() const;

private:

  /**
   * float metal free charge
   */
  PetscScalar      _Qf;

  /**
   * float metal potential
   */
  PetscScalar      _psi;

public:

  //////////////////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for Poisson's Equation---------------------//
  //////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void Poissin_Jacobian_Reserve(Mat *, InsertMode &);

  /**
  * build function and its jacobian for poisson solver, nothing to do
  */
  virtual void Poissin_Function(PetscScalar * , Vec , InsertMode &);

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
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM1_Jacobian_Reserve(Mat * , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML1 solver
   */
  virtual void DDM1_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of DDML1 solver.
   */
  virtual void DDM1_Update_Solution(PetscScalar *);

  //////////////////////////////////////////////////////////////////////////////////
  //----------------Function and Jacobian evaluate for L2 DDM---------------------//
  //////////////////////////////////////////////////////////////////////////////////

  /**
   * reserve none zero pattern in petsc matrix.
   */
  virtual void DDM2_Jacobian_Reserve(Mat * , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Function(PetscScalar * , Vec , InsertMode &);

  /**
   * build function and its jacobian for DDML2 solver
   */
  virtual void DDM2_Jacobian(PetscScalar * , Mat *, InsertMode &);

  /**
   * update solution data of DDML2 solver.
   */
  virtual void DDM2_Update_Solution(PetscScalar *);

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

}
;


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

