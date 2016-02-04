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
#include <algorithm>

#include "genius_env.h"
#include "genius_common.h"

#ifdef WINDOWS
  #include <Windows.h>
  #undef max
  #undef min
  #define LDFUN GetProcAddress
#else
  #include <dlfcn.h>
  #define LDFUN dlsym
#endif

#include "log.h"
#include "ngspice_interface.h"
#include "spdefs.h"
#include "spice_ckt_define.h"
#include "spice_ckt.h"
#include "parallel.h"



SPICE_CKT::SPICE_CKT(const std::string & ckt_file)
    : _ckt_file(ckt_file)
{
  std::string filename =  Genius::genius_dir() + "/lib/spice.so";

#ifdef WINDOWS
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

  _node_ptrs.resize(_n_nodes);
  get_nodal_info(&(_node_ptrs[0]));
  for(unsigned int n=0; n<_n_nodes; ++n)
  {
    SPICE_NODE node;
    node.name          = _node_ptrs[n]->name;
    node.number        = _node_ptrs[n]->number;
    node.type          = _node_ptrs[n]->type;
    node.ic            = _node_ptrs[n]->ic;
    node.nodeset       = _node_ptrs[n]->nodeset;
    node.ic_given      = _node_ptrs[n]->icGiven;
    node.nodeset_given = _node_ptrs[n]->nsGiven;
    node.is_electrode  = 0;
    
    _spice_node_name_to_spice_node_map.insert(std::make_pair(node.name, _node_info_array.size()));
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
    {
      _electrode_to_spice_node_map.insert(std::make_pair(std::string(electrodes[n]), spice_nodes[n]->number));
    }
  }

  // read pattern of spice matrix
  MatrixPtr (*get_matrix)();
  get_matrix = (MatrixPtr (*)())LDFUN(dll_file,"ngspice_get_matrix");
  assert(get_matrix);
  MatrixPtr matrix_ptr = get_matrix();

  double * (*get_matrix_entry)(int , int );
  get_matrix_entry = (double * (*)(int , int ))LDFUN(dll_file,"ngspice_get_matrix_entry");
  assert(get_matrix_entry);

  // load spice matrix
  _circuit_load = (int (*)())LDFUN(dll_file,"ngspice_circuit_load");
  assert(_circuit_load);
  _circuit_load();

  // do matrix reorder
  //_smp_preorder = (int (*)())LDFUN(dll_file,"ngspice_smp_preorder");
  //assert(_smp_preorder);
  //_smp_preorder();

  // get matrix pattern
  assert(matrix_ptr->Size+1 == static_cast<int>(_n_nodes));
  _matrix_nonzero_pattern.resize(matrix_ptr->Size+1);
  _matrix_nonzero_pattern[0].first.push_back(0);
  _matrix_nonzero_pattern[0].second.push_back((double*)NULL);
  _row_to_matrix_row.resize(matrix_ptr->Size+1, 0);
  _col_to_matrix_col.resize(matrix_ptr->Size+1, 0);
  for(int col=0; col<=matrix_ptr->Size; ++col)
    for(ElementPtr elem_ptr=matrix_ptr->FirstInCol[col]; elem_ptr; elem_ptr = elem_ptr->NextInCol)
    {
      int row = elem_ptr->Row;

      int ext_col = matrix_ptr->IntToExtColMap[col];
      int ext_row = matrix_ptr->IntToExtRowMap[row];

      double * entity = get_matrix_entry(ext_row, ext_col);
      _matrix_nonzero_pattern[ext_row].first.push_back(ext_col);
      _matrix_nonzero_pattern[ext_row].second.push_back(entity);
      _row_to_matrix_row[ext_row] = row;
      _col_to_matrix_col[ext_col] = col;
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


  _spice_index_offset_to_global_index = invalid_uint;
  _spice_index_offset_to_local_index  = invalid_uint;
  _spice_index_offset_to_array_index  = invalid_uint;

}




SPICE_CKT::~SPICE_CKT()
{
  if(dll_file)
  {
    int (*destroy_ckt)();
    destroy_ckt = (int (*)())LDFUN(dll_file,"ngspice_destroy_ckt");
    assert(destroy_ckt);
    destroy_ckt();

#ifdef WINDOWS
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
    _solution.resize(_n_nodes, 0.0);
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
      std::map<std::string, unsigned int>::const_iterator it = _electrode_to_spice_node_map.begin();
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
  Parallel::broadcast(_row_to_matrix_row,  root_id);
  Parallel::broadcast(_col_to_matrix_col,  root_id);
}


unsigned int SPICE_CKT::ckt_node_index(const std::string &name)
{
  std::string name_lower_cast = name;
  std::transform(name_lower_cast.begin(), name_lower_cast.end(), name_lower_cast.begin(), ::tolower);
  
  if(_spice_node_name_to_spice_node_map.find(name_lower_cast) != _spice_node_name_to_spice_node_map.end())
    return  _spice_node_name_to_spice_node_map.find(name_lower_cast)->second;
  return invalid_uint;
}


void SPICE_CKT::set_ckt_mode(long mode, bool reset)
{
  if(!reset)
    *_ckt_mode = (*_ckt_mode & MODEUIC) | mode;
  else
    *_ckt_mode = mode;
}


void SPICE_CKT::print_ckt_mode() const
{
  std::cout<<"CKTMode: ";
  if(*_ckt_mode & MODETRAN) std::cout<< "MODETRAN ";
  if(*_ckt_mode & MODEDC) std::cout<< "MODEDC ";


  if(*_ckt_mode & MODETRANOP) std::cout<< "MODETRANOP ";

  if(*_ckt_mode & MODEINITFLOAT) std::cout<< "MODEINITFLOAT ";
  if(*_ckt_mode & MODEINITJCT) std::cout<< "MODEINITJCT ";
  if(*_ckt_mode & MODEINITFIX) std::cout<< "MODEINITFIX ";
  if(*_ckt_mode & MODEINITTRAN) std::cout<< "MODEINITTRAN ";

  if(*_ckt_mode & MODEUIC) std::cout<< "MODEUIC ";

  std::cout<<std::endl;
}


void SPICE_CKT::export_solution(const std::string &file) const
{
  if(Genius::is_last_processor())
  {
    std::ofstream fout(file.c_str());
    // export nodal value
    fout<<"* nodal voltage and branch current"<<std::endl;
    fout<<this->n_ckt_nodes()<<std::endl;
    for(unsigned int n=0; n<this->n_ckt_nodes(); ++n)
    {
      fout << std::left << std::setw(30) << this->ckt_node_name(n)  << std::setw(15) << this->rhs_old(n) << std::endl;
    }
    // export state 
    fout<<"* state"<<std::endl;
    fout<<n_state()<<std::endl;
    
    std::vector<double> state0;
    get_state_vector(0, state0);
    
    std::vector<double> state1;
    get_state_vector(1, state1);
    
    std::vector<double> state2;
    get_state_vector(2, state2);
    
    for(unsigned int n=0; n<state0.size(); ++n)
    {
      fout << std::setw(15) << state0[n] 
           << std::setw(15) << state1[n] 
           << std::setw(15) << state2[n] 
           << std::endl;
    }
    
    fout.close();
  }
}


void SPICE_CKT::import_solution(const std::string &file)
{
  if(Genius::is_last_processor())
  {
    void (*set_nodal_value)(char *, double);
    set_nodal_value = (void (*)(char *, double))LDFUN(dll_file,"ngspice_set_nodal_value");
    assert(set_nodal_value);
    
    std::ifstream fin(file.c_str());
    if(!fin.good()) return;
    
    std::stringstream ss; 
    
    while (!fin.eof())
    {
      std::string line;
      std::getline(fin, line);
      
      if(line.empty()) continue;
    
      size_t pos = line.find('*');
      if( pos != std::string::npos )
        line = line.substr(0, pos);
      
      ss << line << ' ';
    }
   
    fin.close();
    
    std::map<std::string, double> nodal_value;
    std::vector<double> state0,state1,state2;
    
    int solution_num;
    ss >> solution_num;
    for(int i=0; i<solution_num; i++)
    {
      std::string node_name;
      double value;
      ss >> node_name >> value;
      nodal_value[node_name] = value;
    }
    
    int state_num;
    ss >> state_num;
    for(int i=0; i<state_num; i++)
    {
      double s0, s1, s2;
      ss >> s0 >> s1 >> s2;
      state0.push_back(s0);
      state1.push_back(s1);
      state2.push_back(s2);
    }

    
    for(std::map<std::string, double>::const_iterator it=nodal_value.begin(); it!=nodal_value.end(); it++)
    {
      std::string node_name = it->first;
      double value = it->second;
      
      if(_spice_node_name_to_spice_node_map.find(node_name) != _spice_node_name_to_spice_node_map.end())
      {
    
        unsigned int index = _spice_node_name_to_spice_node_map.find(node_name)->second;

        SPICE_NODE & node = _node_info_array[index]; 
        node.nodeset_given = true;
        node.nodeset = value;
        (*_p_rhs_old)[index] = value;
        
        CKTnode * ckt_node = _node_ptrs[index];
        ckt_node->nsGiven = 1 ;
        ckt_node->nodeset = value;
      }
    }
    
    if(n_state() == state_num)
    {
      set_state_vector(0, state0);
      set_state_vector(1, state1);
      set_state_vector(2, state2);
    }
  }
}



void SPICE_CKT::update_rhs_old(const std::vector<double> &rhs)
{
  assert(Genius::is_last_processor());
  for(unsigned int n=0; n<rhs.size(); ++n)
    (*_p_rhs_old)[n] = rhs[n];
}


double SPICE_CKT::rhs(unsigned int n) const
{
  assert(Genius::is_last_processor());
  return  (*_p_rhs)[n];
}


double SPICE_CKT::rhs_old(unsigned int n) const
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




bool SPICE_CKT::is_ckt_voltage_source_exist_sync(const std::string & component)
{
  int flag;
  if(Genius::is_last_processor())
    flag = is_ckt_voltage_source_exist(component);
  Parallel::broadcast(flag, Genius::last_processor_id());
  return static_cast<bool>(flag);
}


bool SPICE_CKT::is_ckt_current_source_exist_sync(const std::string & component)
{
  int flag;
  if(Genius::is_last_processor())
    flag = is_ckt_current_source_exist(component);
  Parallel::broadcast(flag, Genius::last_processor_id());
  return static_cast<bool>(flag);
}



double SPICE_CKT::get_current_from_sync(const std::string & component)
{
  parallel_only();
  double i = 0;
  if(Genius::is_last_processor())
    i = get_current_from(component);
  Parallel::broadcast(i, Genius::last_processor_id());
  return i;
}



double SPICE_CKT::get_voltage_from_sync(const std::string & component)
{
  parallel_only();
  double v = 0;
  if(Genius::is_last_processor())
    v = get_voltage_from(component);
  Parallel::broadcast(v, Genius::last_processor_id());
  return v;
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


unsigned int SPICE_CKT::n_state() const
{
  unsigned int (*n_state)();
  n_state = (unsigned int (*)())LDFUN(dll_file,"ngspice_n_state");
  assert(n_state);
  return n_state();
}


void SPICE_CKT::get_state_vector(int i, std::vector<double> &state) const
{
  state.resize(n_state());

  double ** (*get_state)(int);
  get_state = (double ** (*)(int))LDFUN(dll_file,"ngspice_get_state");
  assert(get_state);

  double * state_array = *(get_state(i));
  for(unsigned int n=0; n<state.size(); ++n)
    state[n] = state_array[n];
}


void SPICE_CKT::set_state_vector(int i, const std::vector<double> &state)
{
  assert(state.size() == n_state());

  double ** (*get_state)(int);
  get_state = (double ** (*)(int))LDFUN(dll_file,"ngspice_get_state");
  assert(get_state);

  double * state_array = *(get_state(i));
  for(unsigned int n=0; n<state.size(); ++n)
    state_array[n] = state[n];
}



void SPICE_CKT::prepare_ckt_state_first_time()
{
  void (*set_state)();
  set_state = (void (*)())LDFUN(dll_file,"ngspice_prepare_ckt_state_first_time");
  assert(set_state);

  set_state();
}


void SPICE_CKT::do_node_set(bool flag)
{
  for(unsigned int n=0; n<_n_nodes; ++n)
  {
    const SPICE_NODE & node = _node_info_array[n];
    CKTnode * ckt_node = _node_ptrs[n];
    if(flag)
    {
      ckt_node->nsGiven = node.nodeset_given ? 1:0;
      ckt_node->nodeset = node.nodeset;
      if(node.nodeset_given)
        (*_p_rhs_old)[n] = node.nodeset;
    }
  }
}


void SPICE_CKT::do_ic(bool flag)
{
  for(unsigned int n=0; n<_n_nodes; ++n)
  {
    const SPICE_NODE & node = _node_info_array[n];
    CKTnode * ckt_node = _node_ptrs[n];
    if(flag)
    {
      ckt_node->icGiven = node.ic_given ? 1:0;
      ckt_node->ic = node.ic;
      if(node.ic_given)
        (*_p_rhs_old)[n] = node.ic;
    }
  }
}



unsigned int SPICE_CKT::global_offset_x(unsigned int n) const
{ return _spice_index_offset_to_global_index + _col_to_matrix_col[n]; }


unsigned int SPICE_CKT::local_offset_x(unsigned int n) const
{ return _spice_index_offset_to_local_index  + _col_to_matrix_col[n]; }


unsigned int SPICE_CKT::array_offset_x(unsigned int n) const
{ return _spice_index_offset_to_array_index  + _col_to_matrix_col[n]; }


unsigned int SPICE_CKT::global_offset_f(unsigned int n) const
{ return _spice_index_offset_to_global_index + _row_to_matrix_row[n]; }


unsigned int SPICE_CKT::local_offset_f(unsigned int n) const
{ return _spice_index_offset_to_local_index  + _row_to_matrix_row[n]; }


unsigned int SPICE_CKT::array_offset_f(unsigned int n) const
{ return _spice_index_offset_to_array_index  + _row_to_matrix_row[n]; }


unsigned int SPICE_CKT::n_nonzero(unsigned int row) const
{ return _matrix_nonzero_pattern[row].first.size(); }


unsigned int SPICE_CKT::max_row_nonzeros() const
{
  unsigned int n_nonzeros=0;
  for(unsigned int n=0; n<_n_nodes; ++n)
  {
    n_nonzeros = std::max(n_nonzeros, n_nonzero(n));
  }
  return n_nonzeros;
}


void SPICE_CKT::build_schur_solver()
{
  std::vector<unsigned int> schur_block;
  std::vector<unsigned int> non_schur_block;
  
  for(unsigned int n=0; n<_node_info_array.size(); ++n)
  {
    if(_node_info_array[n].is_electrode)
      schur_block.push_back(n);
    else  
      non_schur_block.push_back(n);
  }
  assert( schur_block.size() + non_schur_block.size() == _n_nodes );
  _schur_solver.build(_n_nodes,schur_block,non_schur_block);
 
}


int SPICE_CKT::circuit_load()
{ 
  int res = _circuit_load(); 
  //print_ckt_matrix();
  return res;
}
  
  

int SPICE_CKT::circuit_load_schur()
{
  _circuit_load();

  //print_ckt_matrix();

  _schur_solver.MatZero();
  _schur_solver.RHSZero();

  for(unsigned int n=0; n<_n_nodes; n++)
  {
    std::vector<int>  col;
    std::vector<double> values;
    ckt_matrix_row(n, col, values);
    for(unsigned int c=0; c<col.size(); ++c)
    {
      _schur_solver.MatSetValue(n, col[c], values[c]);
    }
    double r = ckt_residual(n);
    _schur_solver.RHSSetValue(n, r);
  }

  _schur_solver.SchurSolve();
  return 0;
}



void SPICE_CKT::ckt_schur_matrix(std::vector<double> & mat)
{ _schur_solver.SchurMatrix(mat); }



void SPICE_CKT::ckt_schur_residual(std::vector<double> & vec)
{ _schur_solver.SchueRHS(vec); }


void SPICE_CKT::update_rhs_old_schur(const std::vector<double> &x)
{
  _schur_solver.SchueSolveX(x);
  std::vector<double> rhs_x;
  _schur_solver.XGetValue(rhs_x);


  // do damping here?
  std::vector<double> rhs(_n_nodes, 0.0);
   
  
  for(unsigned int n=1; n<_n_nodes; ++n)
  {
    double f = 1.0;
    /*
    if(is_voltage_node(n))
    {
      if( std::abs(rhs_x[n-1]) > 1e-6 )
        f = log(1+std::abs(rhs_x[n-1])/1.0)/(std::abs(rhs_x[n-1])/1.0);
    }

    if(is_current_node(n))
    {
      if( std::abs(rhs_x[n-1]) > 1.0 )
        f = 1.0/std::abs(rhs_x[n-1]);
    }
    */
    
    rhs[n] = (*_p_rhs_old)[n] - f*rhs_x[n];
  }

  update_rhs_old(rhs);
}


void SPICE_CKT::ckt_matrix_row(unsigned int n, int &global_row, std::vector<int> & global_col, std::vector<double> & values) const
{
  global_row = this->global_offset_f(n);
  if(n==0)
  {
    global_col.push_back(0+_spice_index_offset_to_global_index);
    values.push_back(1.0);
    return;
  }

  const std::vector<int> & col = _matrix_nonzero_pattern[n].first;
  const std::vector<double *> & col_value = _matrix_nonzero_pattern[n].second;

  for(unsigned int c=0; c<col.size(); ++c)
  {
    global_col.push_back(_col_to_matrix_col[col[c]] + _spice_index_offset_to_global_index);
    values.push_back(*col_value[c]);
  }


  // if we have diagonal item, that's all
  for(unsigned int c=0; c<global_col.size(); ++c)
  { if(global_row==global_col[c]) return; }

  //  fill 1e-20 to diagonal to avoid zero pivot?
  global_col.push_back(global_row);
  values.push_back(1e-20);

}


void SPICE_CKT::ckt_matrix_row(unsigned int row, std::vector<int> & col, std::vector<double> & values) const
{
  if(row==0)
  {
    col.push_back(0);
    values.push_back(1.0);
    return;
  }

  const std::vector<int> & matrix_col = _matrix_nonzero_pattern[row].first;
  const std::vector<double *> & matrix_col_value = _matrix_nonzero_pattern[row].second;

  for(unsigned int c=0; c<matrix_col.size(); ++c)
  {
    col.push_back(matrix_col[c]);
    values.push_back(*matrix_col_value[c]);
  }

}


void SPICE_CKT::print_ckt_matrix() const
{
  std::cout<< std::left << std::setw(32) << "Spice node" << std::setw(20) << "solution" << std::setw(20) << "residual" << std::endl;
  for(unsigned int n=1; n<_n_nodes; n++)
  {
    double r = ckt_residual(n);
    double x = (*_p_rhs_old)[n];
    
    std::cout<< std::left << std::setw(32) << _node_info_array[n].name << std::setw(20) << x << std::setw(20) << r << std::endl;
  }
  std::cout<<std::endl;
  
  
  std::cout<<"spice M"<<std::endl;
  std::cout<<"(0, 0, 1.0)" << std::endl;
  for(unsigned int n=1; n<_n_nodes; n++)
  {
    const std::vector<int> & matrix_col = _matrix_nonzero_pattern[n].first;
    const std::vector<double *> & matrix_col_value = _matrix_nonzero_pattern[n].second;

    for(unsigned int c=0; c<matrix_col.size(); ++c)
    {
      std::cout<<"("<<n<<", "<<matrix_col[c]<<", "<<*matrix_col_value[c]<<") ";
    }
    std::cout<<std::endl;
  }
}




void SPICE_CKT::ckt_residual(unsigned int n, int & global_index, double & value) const
{
  global_index = this->global_offset_f(n);

  if(n==0)
     value = 0.0;
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
    value = v;
  }
}


double SPICE_CKT::ckt_residual(unsigned int n) const
{
  if(n==0)
    return 0.0;
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
    return v;
  }
  return 0.0;
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

