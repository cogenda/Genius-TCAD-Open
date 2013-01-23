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
#ifdef WINDOWS
  class HINSTANCE__; // Forward or never
  typedef HINSTANCE__* HINSTANCE;
#endif

class BoundaryCondition;
typedef struct sCKTnode CKTnode;

#include "schur_solver.h"


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
   * clear genius electrode info
   */
  void clear_electrode()
  {
    _spice_node_links_to_bc.clear();
    _bc_to_spice_node_map.clear();
  }

  /**
   * link genius electrode to spice node
   */
  void link_electrode(unsigned int n, const BoundaryCondition * bc)
  {
    _spice_node_links_to_bc.push_back(n);
    _bc_to_spice_node_map[bc] = n;
  }

  /**
   * after electrode settings, we can init schur solver now
   */
  void build_schur_solver();

  /**
   * @return the number of electrodes
   */
  unsigned int n_electrode() const
  { return  _bc_to_spice_node_map.size(); }

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
  unsigned int n_nonzero(unsigned int row) const;

  /**
   * @return the max nonzeros in row of the spice matrix
   */
  unsigned int max_row_nonzeros() const;

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
   * @return the global offset of spice solution
   */
  unsigned int global_offset_x(unsigned int n) const;

  /**
   * @return the local offset of spice solution
   */
  unsigned int local_offset_x(unsigned int n) const;

  /**
   * @return the index in array which get from global vec by VecGetArray
   */
  unsigned int array_offset_x(unsigned int n) const;

  /**
   * @return the global offset of spice function
   */
  unsigned int global_offset_f(unsigned int n) const;

  /**
   * @return the local offset of spice function
   */
  unsigned int local_offset_f(unsigned int n) const;

  /**
   * @return the index in array which get from global vec by VecGetArray
   */
  unsigned int array_offset_f(unsigned int n) const;

  /**
   * give a spice row number (begin with 0), get the global row index, global cols (for petsc matrix) and their values
   */
  void ckt_matrix_row(unsigned int row, int &global_row, std::vector<int> & global_col, std::vector<double> & values) const;

  /**
   * give a spice row number (begin with 0), get the local cols (begin with 0) and their values
   */
  void ckt_matrix_row(unsigned int row, std::vector<int> & col, std::vector<double> & values) const;

  /**
   * output ckt matrix
   */
  void print_ckt_matrix() const;

  /**
   * get dense schur matrix
   */
  void ckt_schur_matrix(std::vector<double> & mat);

  /**
   * get spice rhs
   */
  double rhs(unsigned int n) const;

  /**
   * get spice rhs_old
   */
  double rhs_old(unsigned int n) const;

  /**
   * write solution to spice rhs_old, which is used by ckt_load
   */
  void update_rhs_old(const std::vector<double> &rhs);

  /**
   * write solution to spice rhs_old, which is used by ckt_load
   */
  void update_rhs_old_schur(const std::vector<double> &x);

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
   * give a spice row number (begin with 0), get the residual with global index (for petsc matrix) and its value
   */
  void ckt_residual(unsigned int n, int & global_index, double & value) const;

  /**
   * give a spice row number (begin with 0), get the residual
   */
  void ckt_residual(unsigned int n, double & value) const;

  /**
   * get ckt res
   */
  void ckt_schur_residual(std::vector<double> & vec);

  /**
   * @return the 2-norm of residual
   */
  double ckt_residual_norm2() const;

  /**
   * @return  Mode of operation of the circuit
   */
  long ckt_mode() const  { return *_ckt_mode; }

  /**
   * set circuit mode
   */
  void set_ckt_mode(long mode, bool reset=true);

  /**
   * print current ckt mode, debug only
   */
  void print_ckt_mode() const;

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
   * @return true when voltage source can be found in ckt, parallel version
   */
  bool is_ckt_voltage_source_exist_sync(const std::string & component);

  /**
   * @return true when current source can be found in ckt, parallel version
   */
  bool is_ckt_current_source_exist_sync(const std::string & component);

  /**
   * get dc voltage by all the processor
   */
  double get_voltage_from_sync(const std::string & component);

  /**
   * get dc current by all the processor
   */
  double get_current_from_sync(const std::string & component);

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
   * get the ith state vector
   */
  void get_state_vector(int i, std::vector<double> &) const;

  /**
   * load circuit, and reorder the matrix
   */
  int circuit_load()
  { return _circuit_load(); }


  /**
   * load circuit, build schur matrix
   */
  int circuit_load_schur();

  /**
   * load nodeset as initial guess
   */
  void do_node_set(bool);

  /**
   * load ic as initial guess
   */
  void do_ic(bool);

  //---------------------------------------------------------------

  /**
   * gmin used in ckt
   */
  double ckt_gmin() const;


  /**
   * gmin used in ckt
   */
  void ckt_set_gmin(double );

  //---------------------------------------------------------------

  /**
   * print spice node location in spice matrix, debug only
   */
  void print_spice_node_location_in_matrix() const;


private:

  /**
   * hold the name of ckt file?
   */
  const std::string  _ckt_file;

  /**
   * pointer to dynamic loaded library file
   */
#ifdef WINDOWS
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
    /// value of the initial condition
    double ic;
    /// value of the .nodeset option
    double nodeset;
    /// flag ic given
    bool ic_given;
    /// flag nodeset given
    bool nodeset_given;
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
   * spice node pointer
   */
  std::vector<CKTnode *> _node_ptrs;

  /**
   * genius electrode name -> spice node
   */
  std::map<std::string, unsigned int> _electrode_to_spice_node_map;

  /**
   * spice node which links to device electrode
   * the order keeps the same as electrode in genius
   */
  std::vector<unsigned int> _spice_node_links_to_bc;

  /**
   * map genius electrode boundary to spice node
   */
  std::map<const BoundaryCondition *, unsigned int> _bc_to_spice_node_map;


  /**
   * spice matrix nonzero pattern, as col (row, ptr), here col and row are spice index
   */
  std::vector< std::pair< std::vector<int>, std::vector<double *> > > _matrix_nonzero_pattern;

  /**
   * we will use reorder routine to remove zero value in matrix diagonal to prevent a zero pivot
   * here saves row index to reordered row index map
   */
  std::vector<int> _row_to_matrix_row;

  std::vector<int> _col_to_matrix_col;

  /**
   * converged solution
   */
  std::vector<double> _solution;

  /**
   * spice schur solver
   */
  SchurSolver _schur_solver;

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

  /**
   * map to ckt component (voltage sources, current sources, resistors, etc)
   */
  std::map<std::string, void *>  _ckt_components;

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
