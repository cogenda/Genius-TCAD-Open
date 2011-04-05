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


#ifndef __spice_ckt_h__
#define __spice_ckt_h__

#include <string>
#include <vector>
#include <map>

#include "config.h"

//for windows dll support
#ifdef CYGWIN
  class HINSTANCE__; // Forward or never
  typedef HINSTANCE__* HINSTANCE;
#endif

class BoundaryCondition;

/**
 * mixed type simulation with ngspice
 * the ngspice is compiled into a dll file
 */
class SPICE_CKT
{
public:
  /**
   * empty constructor
   */
  SPICE_CKT():dll_file(0)
  {
    _p_rhs = 0;
    _p_rhs_old = 0;
    _ckt_mode = 0;
    _spice_index_offset_to_global_index = invalid_uint;
    _spice_index_offset_to_local_index  = invalid_uint;
    _spice_index_offset_to_array_index  = invalid_uint;

    _init_dcop = 0;
    _init_dctrcurv = 0;
    _init_dctran = 0;
    _rotate_state_vectors = 0;
    _circuit_load = 0;

    matrix_reorder = false;
  }

  /**
   * constructor, take spice circuit file as input
   */
  SPICE_CKT(const std::string & ckt_file);

  /**
   * destructor, free the dll file pointer
   */
  ~SPICE_CKT();

  /**
   * sync circuit information between all the processors
   */
  void sync();

  /**
   * @return total node number in circuit
   * including voltage node and current branch node.
   * node 0 (ground) is also considered.
   */
  unsigned int n_ckt_nodes() const
  { return _n_nodes; }

  /**
   * @return the node's name used in spice circuit.
   */
  const std::string & ckt_node_name(unsigned int n) const
    { return _node_info_array[n].name; }

  /**
   * @return the electrode to spice node map
   */
  const std::map<std::string, unsigned int> & get_electrode_info() const
  { return _electrode_to_spice_node_map; }

  /**
   * link genius electrode to spice node
   */
  void link_electrode(unsigned int n, const BoundaryCondition * bc)
  {  _bc_to_spice_node_map[bc] = n; }

  /**
   * @return the spice node index by bc
   */
  unsigned int get_spice_node_by_bc(const BoundaryCondition * bc) const
  {
    if(_bc_to_spice_node_map.find(bc)!=_bc_to_spice_node_map.end())
      return _bc_to_spice_node_map.find(bc)->second;
    return invalid_uint;
  }

  /**
   * @return true if the spice node is voltage type
   */
  bool is_voltage_node(unsigned int n) const
  { return _node_info_array[n].type == 3; }

  /**
   * @return true if the spice node is current type
   */
  bool is_current_node(unsigned int n) const
  { return _node_info_array[n].type == 4; }

  /**
   * ckt node links to device electrode
   */
  void set_ckt_node_electrode_flag(unsigned int n)
  { _node_info_array[n].is_electrode = 1; }

  /**
   * @return the nonzero pattern for a row
   */
  unsigned int n_nonzero(unsigned int row) const
    { return _matrix_nonzero_pattern[row].first.size(); }

  /**
   * set the offset, then we can map spice matrix to petsc global matrix
   */
  void set_offset(unsigned int global, unsigned int local=invalid_uint, unsigned int array=invalid_uint)
  {
    _spice_index_offset_to_global_index = global;
    _spice_index_offset_to_local_index  = local;
    _spice_index_offset_to_array_index  = array;
  }

  /**
   * @return the offset of spice vector to global petsc vector
   */
  unsigned int spice_global_offset() const
  { return _spice_index_offset_to_global_index; }

  /**
   * @return the offset of spice vector to local petsc vector
   */
  unsigned int spice_local_offset() const
  { return _spice_index_offset_to_local_index; }

  /**
   * @return the index in array which get from global vec by VecGetArray
   */
  unsigned int spice_array_offset() const
  { return _spice_index_offset_to_array_index; }

  /**
   * @return the global offset of a spice node
   */
  unsigned int global_offset(unsigned int n) const
  { return _spice_index_offset_to_global_index + _spice_order_to_matrix_order.find(static_cast<int>(n))->second; }

  /**
   * @return the local offset of a spice node
   */
  unsigned int local_offset(unsigned int n) const
  { return _spice_index_offset_to_local_index  + _spice_order_to_matrix_order.find(static_cast<int>(n))->second; }

  /**
   * @return the index in array which get from global vec by VecGetArray
   */
  unsigned int array_offset(unsigned int n) const
  { return _spice_index_offset_to_array_index  + _spice_order_to_matrix_order.find(static_cast<int>(n))->second; }

  /**
   * give a spice row number (begin with 0), get the global cols (for petsc matrix) and their values
   */
  void ckt_matrix_row(unsigned int row, int &global_row, std::vector<int> & global_col, std::vector<double> & values) const;

  /**
   * get spice rhs
   */
  double rhs(unsigned int n) const;

  /**
   * write to spice rhs
   */
  double & rhs(unsigned int n);

  /**
   * get spice rhs_old
   */
  double rhs_old(unsigned int n) const;

  /**
   * write to spice rhs_old
   */
  double & rhs_old(unsigned int n);

  /**
   * save rhs_old to _previous_solution
   * also do sync between all the processors
   * must call in parallel
   */
  void save_solution();

  /**
   * cover rhs_old with backuped value, only the last processor can do this!
   */
  void restore_solution();

  /**
   * get solution, this function can be call on all the processors
   */
  double get_solution(unsigned int n) const
  { return _solution[n]; }

  /**
   * get the residual with global index (for petsc matrix) and their values
   */
  void ckt_residual(std::vector<int> & global_index, std::vector<double> & values) const;

  /**
   * @return the 2-norm of residual
   */
  double ckt_residual_norm2() const;

  /**
   * @return  Mode of operation of the circuit
   */
  long ckt_mode() const
  { return *_ckt_mode; }

  /**
   * set circuit mode
   */
  void set_ckt_mode(long mode, bool reset=true);

  /**
   * change the CKTMode during iteration
   */
  int change_ckt_mode(int covergence)
  { return _change_ckt_mode(covergence); }

  /**
   * set the temperature spice circuit operate at
   */
  void set_temperature(double temp);

  //---------------------------------------------------------------

  /**
   * dc operator
   */
  int init_dcop()
  { return _init_dcop(); }

  //---------------------------------------------------------------

  /**
   * dc trace curve
   */
  int init_dctrcurv()
  { return _init_dctrcurv(); }

  /**
   * @return all the voltage sources in ckt
   */
  void get_voltage_sources(std::vector<std::string> & vsrcs);

  /**
   * @return all the current sources in ckt
   */
  void get_current_sources(std::vector<std::string> & isrcs);

  /**
   * @return true when voltage source can be found in ckt
   */
  bool is_ckt_voltage_source_exist(const std::string & component);

  /**
   * @return true when current source can be found in ckt
   */
  bool is_ckt_current_source_exist(const std::string & component);

  /**
   * @return true when resistor can be found in ckt
   */
  bool is_ckt_resistor_exist(const std::string & component);

  /**
   * get voltage source from spice ckt for later usage
   */
  void get_ckt_voltage_source(const std::string & component);

  /**
   * get current source from spice ckt for later usage
   */
  void get_ckt_current_source(const std::string & component);

  /**
   * get resistor from spice ckt for later usage
   */
  void get_ckt_resistor(const std::string & component);

  /**
   * get dc voltage
   */
  double get_voltage_from(const std::string & component);

  /**
   * set dc voltage
   */
  void set_voltage_to(const std::string & component, double v);

  /**
   * get dc current
   */
  double get_current_from(const std::string & component);

  /**
   * set dc current
   */
  void set_current_to(const std::string & component, double i);

  /**
   * set resistance
   */
  void set_resistance_to(const std::string & component, double r);

  //---------------------------------------------------------------

  /**
   * transient
   */
  int init_dctran(long TRANmode)
  { return  _init_dctran(TRANmode); }

  /**
   * tell spice to do modeinittran
   */
  void set_modeinittran(double delta);

  /**
   * tell spice current time step
   */
  void set_delta(double delta);

  /**
   * tell spice current time
   */
  void set_time(double time);

  /**
   * tell spice current time integrate order
   */
  void set_time_order(int order, int max_order=2);

  /**
   * tell spice which integrate method should be used
   */
  void set_integrate_method(int method);

  /**
   * for first step of transient simulation, set state vector
   */
  void prepare_ckt_state_first_time();


  /**
   * set UIC flag to true
   */
  void set_uic(bool);

  /**
   * @return true if UIC required
   */
  long uic();


  void exchange_rhs()
  {
    double *tmp;
    tmp =  *_p_rhs_old;
    *_p_rhs_old = *_p_rhs;
    *_p_rhs = tmp;
  }

  //---------------------------------------------------------------

  /**
   * ask spice to rotate the CKTstates vectors
   */
  void rotate_state_vectors()
  { _rotate_state_vectors(); }

 /**
  * load circuit, and reorder the matrix
  */
  int circuit_load()
  {
    return _circuit_load();
  }

  bool is_reordered() const
  { return matrix_reorder; }

  /**
   * reorder the matrix, must call it in parallel
   */
  int smp_preorder();


  //---------------------------------------------------------------

  /**
   * gmin used in ckt
   */
  double ckt_gmin() const;


  /**
   * gmin used in ckt
   */
  void ckt_set_gmin(double );

private:

  /**
   * hold the name of ckt file?
   */
  const std::string  _ckt_file;

  /**
   * pointer to dynamic loaded library file
   */
#ifdef CYGWIN
  HINSTANCE                   dll_file;
#else
  void                      * dll_file;
#endif

  /**
   * number of nodes in spice circuit
   */
  unsigned int _n_nodes;

  struct SPICE_NODE
  {
    /// name of the spice node, case insensitive
    std::string name;
    /// number in spice
    int  number;
    /// node type voltage or current
    int  type;
    /// true if this node connect to device electrode
    int  is_electrode;
    /// global offset in genius matrix
    unsigned int global_dof;
    /// local offset in genius matrix
    unsigned int local_dof;
  };

  /**
   * save node information
   */
  std::vector<SPICE_NODE> _node_info_array;

  /**
   * genius electrode name -> spice node
   */
  std::map<std::string, unsigned int> _electrode_to_spice_node_map;

  /**
   * map genius electrode boundary to spice node
   */
  std::map<const BoundaryCondition *, unsigned int> _bc_to_spice_node_map;


  /**
   * spice matrix nonzero pattern, as col (row, ptr), here col and row are spice index
   */
  std::vector< std::pair< std::vector<int>, std::vector<double *> > > _matrix_nonzero_pattern;

  /**
   * a flag to show if spice matrix is reordered
   */
  bool matrix_reorder;

  /**
   * we will use reorder routine to remove zero value in matrix diagonal to prevent a zero pivot
   * this map saves spice index to reordered matrix index
   */
  std::map<int, int> _spice_order_to_matrix_order;

  /**
   * set _spice_order_to_matrix_order map
   */
  void _reset_spice_to_genius_mapping();


  /**
   * converged solution
   */
  std::vector<double> _solution;

  /**
   * pointer to spice rhs
   */
  double ** _p_rhs;

  /**
   * pointer to spice rhs_old
   */
  double ** _p_rhs_old;

  /**
   * pointer to spice cktmode
   */
  long * _ckt_mode;

  /**
   * the global_index in matrix/vec = spice_index + this_offset
   */
  unsigned int _spice_index_offset_to_global_index;

  /**
   * the local_index in matrix/vec = spice_index + this_offset
   */
  unsigned int _spice_index_offset_to_local_index;

  /**
   * the index in array which get from global vec by VecGetArray
   */
  unsigned int _spice_index_offset_to_array_index;


  std::map<const std::string, void *>  _ckt_components;

  /**
   * init D.C. Operating point analysis
   */
  int (*_init_dcop)();

  /**
   * init D.C. Transfer curve analysis
   */
  int (*_init_dctrcurv)();

  /**
   * init Transient analysis
   */
  int (*_init_dctran)(long);

  /**
   * ask spice to rotate the CKTstates vectors
   */
  void (*_rotate_state_vectors)();

  /**
   * load spice matrix and rhs
   */
  int (*_circuit_load)();

  /**
   * reorder spice matrix to remove zero diagonal value
   */
  int (*_smp_preorder)();

  /**
   * change the CKTMode during iteration
   */
  int (*_change_ckt_mode)(int);

};

#endif
