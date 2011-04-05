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

#include <cmath>
#include <cassert>

#include "genius_env.h"
#include "genius_common.h"

#ifdef CYGWIN
  #include <Windows.h>
  #undef max
  #undef min
  #define LDFUN GetProcAddress
#else
  #include <dlfcn.h>
  #define LDFUN dlsym
#endif

#include "log.h"
#include "spice_ckt.h"
#include "ngspice_interface.h"
#include "spdefs.h"
#include "spice_ckt_define.h"
#include "parallel.h"



SPICE_CKT::SPICE_CKT(const std::string & ckt_file)
    : _ckt_file(ckt_file)
{
  std::string filename =  Genius::genius_dir() + "/lib/spice.so";

#ifdef CYGWIN
  dll_file = LoadLibrary(filename.c_str());
#else
  // NOTE: Because genius exports all its symbol with --export-dynamic flag
  // It seems libspice has symbol conflict when Intel MKL is used (fblaslapack is ok).
  // Here we use RTLD_DEEPBIND to prevent symbol confilict with genius.
  // With this flag, a self-contained library will use its
  // own symbols in preference to global symbols with the same name
  // contained in libraries that have already been loaded.
  // However, RTLD_DEEPBIND flag is only valid for GLIBC >= 2.3.4
#ifdef RTLD_DEEPBIND
  dll_file = dlopen(filename.c_str(), RTLD_LAZY | RTLD_DEEPBIND);
#else
  dll_file = dlopen(filename.c_str(), RTLD_LAZY );
#endif
#endif

  if(dll_file==NULL)
  { MESSAGE<<"Load spice engine error." << '\n'; RECORD(); genius_error();}

  //read the spice netlist
  int (*init_netlist) (char * );
  init_netlist = (int (*)(char * ))LDFUN(dll_file,"ngspice_init_netlist");
  assert(init_netlist);
  if( init_netlist(const_cast<char *>(_ckt_file.c_str())) )
  { MESSAGE<<"Init spice netlist error.\n"; RECORD(); genius_error();}


  //build spice circuit
  int (*init_ckt) ();
  init_ckt = (int (*)())LDFUN(dll_file,"ngspice_init_ckt");
  assert(init_ckt);
  if( init_ckt() )
  { MESSAGE<<"Init spice circuit error.\n"; RECORD(); genius_error();}

  //how many nodes in spice circuit?
  unsigned int (*n_nodes)();
  n_nodes = (unsigned int (*)())LDFUN(dll_file,"ngspice_n_nodes");
  assert(n_nodes);
  _n_nodes = n_nodes();

  // read spice node information
  void (*get_nodal_info)(CKTnode **);
  get_nodal_info = (void (*)(CKTnode **))LDFUN(dll_file,"ngspice_get_nodal_info");
  assert(get_nodal_info);
  std::vector<CKTnode *> node_ptrs(_n_nodes);
  get_nodal_info(&node_ptrs[0]);
  for(unsigned int n=0; n<_n_nodes; ++n)
  {
    SPICE_NODE node;
    node.name   = node_ptrs[n]->name;
    node.number = node_ptrs[n]->number;
    node.type   = node_ptrs[n]->type;
    node.is_electrode = 0;
    _node_info_array.push_back(node);
  }
  _solution.resize(_n_nodes, 0.0);

  // get ndev electrode
  unsigned int (*n_electrode)();
  n_electrode = (unsigned int (*)())LDFUN(dll_file,"ngspice_n_electrode");
  assert(n_electrode);

  unsigned int _n_electrode = n_electrode();
  if(_n_electrode)
  {
    std::vector<char *> electrodes(_n_electrode);
    std::vector<CKTnode *> spice_nodes(_n_electrode);

    void (*find_electrode)(char **, CKTnode ** );
    find_electrode = (void (*)(char **, CKTnode **))LDFUN(dll_file,"ngspice_find_electrode");
    assert(find_electrode);
    find_electrode(&electrodes[0], &spice_nodes[0]);
    for(unsigned int n=0; n<_n_electrode; ++n)
      _electrode_to_spice_node_map.insert(std::make_pair(std::string(electrodes[n]), spice_nodes[n]->number));
  }

  // read pattern of spice matrix
  MatrixPtr (*get_matrix)();
  get_matrix = (MatrixPtr (*)())LDFUN(dll_file,"ngspice_get_matrix");
  assert(get_matrix);
  MatrixPtr matrix_ptr = get_matrix();

  double * (*get_matrix_entry)(int , int );
  get_matrix_entry = (double * (*)(int , int ))LDFUN(dll_file,"ngspice_get_matrix_entry");
  assert(get_matrix_entry);

  // consider 0 node
  assert(matrix_ptr->Size+1 == static_cast<int>(_n_nodes));
  _matrix_nonzero_pattern.resize(matrix_ptr->Size+1);
  _matrix_nonzero_pattern[0].first.push_back(0);
  _matrix_nonzero_pattern[0].second.push_back((double*)NULL);
  _spice_order_to_matrix_order[0]=0;
  //matrix pattern
  for(int col=1; col<=matrix_ptr->Size; ++col)
    for(ElementPtr elem_ptr=matrix_ptr->FirstInCol[col]; elem_ptr; elem_ptr = elem_ptr->NextInCol)
    {
      int row = elem_ptr->Row;

      int ext_col = matrix_ptr->IntToExtColMap[col];
      int ext_row = matrix_ptr->IntToExtRowMap[row];

      double * entity = get_matrix_entry(ext_row, ext_col);
      _matrix_nonzero_pattern[ext_row].first.push_back(ext_col);
      _matrix_nonzero_pattern[ext_row].second.push_back(entity);
      _spice_order_to_matrix_order[ext_col] = col;
    }

  // get spice rhs
  double ** (*get_rhs)();
  get_rhs = (double ** (*)())LDFUN(dll_file,"ngspice_get_rhs");
  assert(get_rhs);
  _p_rhs = get_rhs();

  // get spice rhs_old
  double ** (*get_rhs_old)();
  get_rhs_old = (double ** (*)())LDFUN(dll_file,"ngspice_get_rhs_old");
  assert(get_rhs_old);
  _p_rhs_old = get_rhs_old();

  long * (*get_ckt_mode)();
  get_ckt_mode = (long * (*)())LDFUN(dll_file,"ngspice_get_ckt_mode");
  assert(get_ckt_mode);
  _ckt_mode = get_ckt_mode();


  _change_ckt_mode = (int (*)(int))LDFUN(dll_file,"ngspice_change_ckt_mode");
  assert(_change_ckt_mode);

  //set interface functions
  _init_dcop = (int (*)())LDFUN(dll_file,"ngspice_init_dcop");
  assert(_init_dcop);

  _init_dctrcurv = (int (*)())LDFUN(dll_file,"ngspice_init_dctrcurv");
  assert(_init_dctrcurv);

  _init_dctran = (int (*)(long))LDFUN(dll_file,"ngspice_init_dctran");
  assert(_init_dctran);

  _rotate_state_vectors = (void (*)())LDFUN(dll_file,"ngspice_rotate_state_vectors");
  assert(_rotate_state_vectors);

  _circuit_load = (int (*)())LDFUN(dll_file,"ngspice_circuit_load");
  assert(_circuit_load);

  _smp_preorder = (int (*)())LDFUN(dll_file,"ngspice_smp_preorder");
  assert(_smp_preorder);

  _spice_index_offset_to_global_index = invalid_uint;
  _spice_index_offset_to_local_index  = invalid_uint;
  _spice_index_offset_to_array_index  = invalid_uint;

  matrix_reorder = false;
}




SPICE_CKT::~SPICE_CKT()
{
  if(dll_file)
  {
    int (*destroy_ckt)();
    destroy_ckt = (int (*)())LDFUN(dll_file,"ngspice_destroy_ckt");
    assert(destroy_ckt);
    destroy_ckt();

#ifdef CYGWIN
    FreeLibrary(dll_file);
#else
    dlclose( dll_file );
#endif

  }
}


void SPICE_CKT::sync()
{
  if(Genius::n_processors()==1) return;

  unsigned int root_id = Genius::last_processor_id();

  // broadcast spice node info
  Parallel::broadcast(_n_nodes, root_id);
  if(Genius::processor_id()!=root_id)
  {
    _node_info_array.resize(_n_nodes);
    _matrix_nonzero_pattern.resize(_n_nodes);
  }
  for(unsigned int n=0; n<_n_nodes; ++n)
  {
    Parallel::broadcast(_node_info_array[n].name,          root_id);
    Parallel::broadcast(_node_info_array[n].number,        root_id);
    Parallel::broadcast(_node_info_array[n].type,          root_id);
    Parallel::broadcast(_node_info_array[n].is_electrode,  root_id);
  }

  // broadcast genius electrode to spice node map
  {
    std::vector<std::string> key;
    std::vector<unsigned int> value;
    if(Genius::processor_id()==root_id)
    {
      std::map<std::string, unsigned int>::iterator it = _electrode_to_spice_node_map.begin();
      for(; it!=_electrode_to_spice_node_map.end(); ++it)
      {
        key.push_back(it->first);
        value.push_back(it->second);
      }
    }
    Parallel::broadcast(key,root_id);
    Parallel::broadcast(value,root_id);
    assert(key.size() == value.size() );
    if(Genius::processor_id()!=root_id)
    {
       for(unsigned int n=0; n<key.size(); ++n)
         _electrode_to_spice_node_map[key[n]]=value[n];
    }
  }

  // broadcast matrix pattern
  for(unsigned int n=0; n<_n_nodes; ++n)
  {
    Parallel::broadcast(_matrix_nonzero_pattern[n].first,  root_id);
  }
  Parallel::broadcast(_spice_order_to_matrix_order,  root_id);
}



void SPICE_CKT::set_ckt_mode(long mode, bool reset)
{
  if(reset)
    *_ckt_mode = (*_ckt_mode & MODEUIC) | mode;
  else
    *_ckt_mode = mode;
}


double SPICE_CKT::rhs(unsigned int n) const
{
  assert(Genius::is_last_processor());
  return  (*_p_rhs)[n];
}


double & SPICE_CKT::rhs(unsigned int n)
{
  assert(Genius::is_last_processor());
  return  (*_p_rhs)[n];
}


double SPICE_CKT::rhs_old(unsigned int n) const
{
  assert(Genius::is_last_processor());
  return  (*_p_rhs_old)[n];
}


double & SPICE_CKT::rhs_old(unsigned int n)
{
  assert(Genius::is_last_processor());
  return  (*_p_rhs_old)[n];
}


void SPICE_CKT::set_uic(bool uic)
{
  if(uic)
    (*_ckt_mode) |= 0x100000;
  else
    (*_ckt_mode) &= ~0x100000;
}

long SPICE_CKT::uic()
{ return (*_ckt_mode) & MODEUIC; }


void SPICE_CKT::set_temperature(double temp)
{
  void (*set_temp)(double);

  set_temp = (void (*)(double))LDFUN(dll_file, "ngspice_set_temperature");
  assert(set_temp);

  set_temp(temp);
}


void SPICE_CKT::save_solution()
{
  _solution.clear();

  if(Genius::is_last_processor())
  {
    for(unsigned int n=0; n<_n_nodes; ++n)
      _solution.push_back((*_p_rhs_old)[n]);
  }

  Parallel::broadcast(_solution, Genius::last_processor_id());
}


void SPICE_CKT::restore_solution()
{
  assert(_solution.size()==_n_nodes);
  assert(Genius::is_last_processor());
  for(unsigned int n=0; n<_n_nodes; ++n)
    (*_p_rhs_old)[n] = _solution[n];
}


void SPICE_CKT::get_voltage_sources(std::vector<std::string> & vsrcs)
{
  int (*find_n_voltage_source)();
  find_n_voltage_source = (int (*)())LDFUN(dll_file,"ngspice_find_n_voltage_source");
  assert(find_n_voltage_source);

  int n_vsrcs = find_n_voltage_source();
  if( n_vsrcs == 0 ) return;

  std::vector<char *> vsrc_names(n_vsrcs);

  void (*find_voltage_sources)( char **);
  find_voltage_sources = (void (*)(char **))LDFUN(dll_file,"ngspice_find_voltage_sources");
  assert(find_voltage_sources);

  find_voltage_sources(&vsrc_names[0]);
  for(unsigned int n=0; n<vsrc_names.size(); ++n)
    vsrcs.push_back(vsrc_names[n]);
}


void SPICE_CKT::get_current_sources(std::vector<std::string> & isrcs)
{
  int (*find_n_current_source)();
  find_n_current_source = (int (*)())LDFUN(dll_file,"ngspice_find_n_current_source");
  assert(find_n_current_source);

  int n_isrcs = find_n_current_source();
  if( n_isrcs == 0 ) return;

  std::vector<char *> isrc_names(n_isrcs);

  void (*find_current_sources)( char **);
  find_current_sources = (void (*)(char **))LDFUN(dll_file,"ngspice_find_current_sources");
  assert(find_current_sources);

  find_current_sources(&isrc_names[0]);
  for(unsigned int n=0; n<isrc_names.size(); ++n)
    isrcs.push_back(isrc_names[n]);
}


bool SPICE_CKT::is_ckt_voltage_source_exist(const std::string & component)
{
  void * (*find_voltage_source)(const char * );

  find_voltage_source = (void * (*)(const char *))LDFUN(dll_file,"ngspice_find_voltage_source");
  assert(find_voltage_source);

  return find_voltage_source(component.c_str())!=NULL;
}


bool SPICE_CKT::is_ckt_current_source_exist(const std::string & component)
{
  void * (*find_current_source)(const char * );

  find_current_source = (void * (*)(const char *))LDFUN(dll_file,"ngspice_find_current_source");
  assert(find_current_source);

  return find_current_source(component.c_str())!=NULL;
}


bool SPICE_CKT::is_ckt_resistor_exist(const std::string & component)
{
  void * (*find_resistor)(const char * );

  find_resistor = (void * (*)(const char *))LDFUN(dll_file,"ngspice_find_resistor");
  assert(find_resistor);

  return find_resistor(component.c_str())!=NULL;
}


void SPICE_CKT::get_ckt_voltage_source(const std::string & component)
{
  void * (*find_voltage_source)(const char * );

  find_voltage_source = (void * (*)(const char *))LDFUN(dll_file,"ngspice_find_voltage_source");
  assert(find_voltage_source);

  _ckt_components[component] = find_voltage_source(component.c_str());
}


void SPICE_CKT::get_ckt_current_source(const std::string & component)
{
  void * (*find_current_source)(const char * );

  find_current_source = (void * (*)(const char *))LDFUN(dll_file,"ngspice_find_current_source");
  assert(find_current_source);

  _ckt_components[component] = find_current_source(component.c_str());
}



void SPICE_CKT::get_ckt_resistor(const std::string & component)
{
  void * (*find_resistor)(const char * );

  find_resistor = (void * (*)(const char *))LDFUN(dll_file,"ngspice_find_resistor");
  assert(find_resistor);

  _ckt_components[component] = find_resistor(component.c_str());
}


double SPICE_CKT::get_voltage_from(const std::string & component)
{
  void (*get_voltage_source)(void *, double *);
  get_voltage_source = (void (*)(void *, double *))LDFUN(dll_file,"ngspice_get_voltage_source");
  assert(get_voltage_source);

  if(_ckt_components.find(component) == _ckt_components.end() )
    get_ckt_voltage_source(component);

  double v=0;
  void * vsrc = _ckt_components.find(component)->second;
  if(vsrc) get_voltage_source(vsrc, &v);

  return v;
}


void SPICE_CKT::set_voltage_to(const std::string & component, double v)
{
  void (*set_voltage_source)(void *, double);
  set_voltage_source = (void (*)(void *, double))LDFUN(dll_file,"ngspice_set_voltage_source");
  assert(set_voltage_source);

  if(_ckt_components.find(component) == _ckt_components.end() )
    get_ckt_voltage_source(component);

  void * vsrc = _ckt_components.find(component)->second;
  if(vsrc) set_voltage_source(vsrc, v);
}


double SPICE_CKT::get_current_from(const std::string & component)
{
  void (*get_current_source)(void *, double *);
  get_current_source = (void (*)(void *, double *))LDFUN(dll_file,"ngspice_get_current_source");
  assert(get_current_source);

  if(_ckt_components.find(component) == _ckt_components.end() )
    get_ckt_current_source(component);

  double i=0;
  void * isrc = _ckt_components.find(component)->second;
  if(isrc) get_current_source(isrc, &i);

  return i;
}


void SPICE_CKT::set_current_to(const std::string & component, double i)
{
  void (*set_current_source)(void *, double);
  set_current_source = (void (*)(void *, double))LDFUN(dll_file,"ngspice_set_current_source");
  assert(set_current_source);

  if(_ckt_components.find(component) == _ckt_components.end() )
    get_ckt_current_source(component);

  void * isrc = _ckt_components.find(component)->second;
  if(isrc) set_current_source(isrc, i);
}



void SPICE_CKT::set_resistance_to(const std::string & component, double r)
{
  void (*set_resistance)(void *, double);
  set_resistance = (void (*)(void *, double))LDFUN(dll_file,"ngspice_set_resistance");
  assert(set_resistance);

  if(_ckt_components.find(component) == _ckt_components.end() )
    get_ckt_resistor(component);

  void * res = _ckt_components.find(component)->second;
  set_resistance(res, r);
}



void SPICE_CKT::set_modeinittran(double delta)
{
  void (*set_inittran)(double);
  set_inittran = (void (*)(double))LDFUN(dll_file,"ngspice_modeinittran");
  assert(set_inittran);

  set_inittran(delta);
}



void SPICE_CKT::set_delta(double delta)
{
  void (*set_dt)(double);
  set_dt = (void (*)(double))LDFUN(dll_file,"ngspice_set_delta");
  assert(set_dt);

  set_dt(delta);
}



void SPICE_CKT::set_time(double time)
{
  void (*set_t)(double);
  set_t = (void (*)(double))LDFUN(dll_file,"ngspice_set_time");
  assert(set_t);

  set_t(time);
}



void SPICE_CKT::set_time_order(int order, int max_order)
{
  void (*set_order)(int, int);
  set_order = (void (*)(int, int))LDFUN(dll_file,"ngspice_set_time_order");
  assert(set_order);

  set_order(order, max_order);
}



void SPICE_CKT::set_integrate_method(int method)
{
  void (*set_method)(int);
  set_method = (void (*)(int))LDFUN(dll_file,"ngspice_set_integrate_method");
  assert(set_method);

  set_method(method);
}



void SPICE_CKT::prepare_ckt_state_first_time()
{
  void (*set_state)();
  set_state = (void (*)())LDFUN(dll_file,"ngspice_prepare_ckt_state_first_time");
  assert(set_state);

  set_state();
}



int SPICE_CKT::smp_preorder()
{
  if(!matrix_reorder)
  {
    if(Genius::is_last_processor()) _smp_preorder();
    _reset_spice_to_genius_mapping();
    matrix_reorder = true;
  }
  return 0;
}



void SPICE_CKT::_reset_spice_to_genius_mapping()
{
  if(Genius::is_last_processor())
  {

    MatrixPtr (*get_matrix)();
    get_matrix = (MatrixPtr (*)())LDFUN(dll_file,"ngspice_get_matrix");
    assert(get_matrix);
    MatrixPtr matrix_ptr = get_matrix();

    // consider 0 node
    _spice_order_to_matrix_order[0] = 0;

    // matrix pattern
    for(int col=1; col<=matrix_ptr->Size; ++col)
    {
      int ext_col = matrix_ptr->IntToExtColMap[col];
      _spice_order_to_matrix_order[ext_col] = col;
    }
  }

  Parallel::broadcast(_spice_order_to_matrix_order,  Genius::last_processor_id());

}


void SPICE_CKT::ckt_matrix_row(unsigned int row, int &global_row, std::vector<int> & global_col, std::vector<double> & values) const
{
  global_row = row + _spice_index_offset_to_global_index;

  if(row==0)
  {
    global_col.push_back(0+_spice_index_offset_to_global_index);
    values.push_back(1.0);
    return;
  }

  const std::vector<int> & col = _matrix_nonzero_pattern[row].first;
  const std::vector<double *> & col_value = _matrix_nonzero_pattern[row].second;

  for(unsigned int n=0; n<col.size(); ++n)
  {
    global_col.push_back(this->global_offset(col[n]));
    values.push_back(*col_value[n]);
  }

  // if we have diagonal item, that's all
  for(unsigned int n=0; n<col.size(); ++n)
  { if( int(this->global_offset(col[n])) == global_row) return; }

  // else fill 1e-20 to diagonal to avoid zero pivot?
  global_col.push_back(global_row);
  values.push_back(1e-20);

}



void SPICE_CKT::ckt_residual(std::vector<int> & global_index, std::vector<double> & values) const
{
  for(unsigned int n=0; n<_n_nodes; ++n)
  {
    global_index.push_back(n + _spice_index_offset_to_global_index);

    if(n==0)
      values.push_back((*_p_rhs_old)[0]);
    else
    {
      const std::vector<int> & col = _matrix_nonzero_pattern[n].first;
      const std::vector<double *> & col_value = _matrix_nonzero_pattern[n].second;

      double v = 0.0;
      for(unsigned int i=0; i<col.size(); ++i)
      {
        unsigned int col_index = col[i];
        v += (*col_value[i])*(*_p_rhs_old)[col_index];
      }
      v -= (*_p_rhs)[n];

      values.push_back(v);
    }
  }
}


double SPICE_CKT::ckt_residual_norm2() const
{
  double norm2 = 0.0;

  for(unsigned int n=1; n<_n_nodes; ++n)
  {
    if(!_node_info_array[n].is_electrode)
    {
      const std::vector<int> & col = _matrix_nonzero_pattern[n].first;
      const std::vector<double *> & col_value = _matrix_nonzero_pattern[n].second;

      double v = 0.0;
      for(unsigned int i=0; i<col.size(); ++i)
      {
        unsigned int col_index = col[i];
        v += (*col_value[i])*(*_p_rhs_old)[col_index];
      }
      v -= (*_p_rhs)[n];

      norm2 += v*v;
    }
  }

  return sqrt(norm2);
}



double SPICE_CKT::ckt_gmin() const
{
  double (*gmin) ();
  gmin = (double (*)())LDFUN(dll_file,"ngspice_gmin");
  assert(gmin);
  return gmin();
}


void SPICE_CKT::ckt_set_gmin(double gmin)
{
  void (*set_gmin) (double);
  set_gmin = (void (*)(double))LDFUN(dll_file,"ngspice_set_gmin");
  assert(set_gmin);
  set_gmin(gmin);
}

