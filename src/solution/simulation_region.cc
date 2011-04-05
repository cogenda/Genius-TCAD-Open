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

#include "elem.h"
#include "simulation_region.h"
#include "material.h"
#include "parallel.h"

// static member
std::map<unsigned int,  SimulationRegion *>  SimulationRegion::_subdomain_id_to_region_map;



SimulationRegion::SimulationRegion(const std::string &name, const std::string &material, const PetscScalar T)
    :_region_name(name), _region_material(material), _T_external(T)
{}


SimulationRegion::~SimulationRegion()
{
  this->clear();
}


PetscScalar SimulationRegion::T_external() const
  { return _T_external; }


unsigned int SimulationRegion::n_on_processor_cell() const
{
  unsigned int ncell = 0;
  std::vector<const Elem *>:: const_iterator it = _region_cell.begin();
  for(; it != _region_cell.end(); ++it)
    if( (*it)->processor_id() == Genius::processor_id() ) ncell++;
  return ncell;
}

unsigned int SimulationRegion::n_on_processor_node() const
{
  unsigned int n_node = 0;
  std::map<unsigned int, FVM_Node *>::const_iterator it = _region_node.begin();
  for(; it != _region_node.end(); ++it)
    if( (*it).second->on_processor() ) n_node++;
  return n_node;
}


void SimulationRegion::region_node(std::vector<unsigned int> & nodes) const
{
  parallel_only();

  for(unsigned int n=0; n<_region_processor_node.size(); ++n)
  {
    const FVM_Node * fvm_node = _region_processor_node[n];
    nodes.push_back(fvm_node->root_node()->id());
  }

  Parallel::allgather( nodes );
  std::sort(nodes.begin(), nodes.end());
}


void SimulationRegion::clear()
{
  _region_cell.clear();

  for(unsigned int n=0; n<_region_cell_data.size(); ++n)
    delete _region_cell_data[n];
  _region_cell_data.clear();

  std::map<unsigned int, FVM_Node *>::iterator it = _region_node.begin();

  for( ; it != _region_node.end(); it++ )
  {
    delete (*it).second;
    (*it).second = NULL;
  }

  _region_node.clear();
  _region_local_node.clear();
  _region_processor_node.clear();

  _cell_data_storage.clear();
  _node_data_storage.clear();

  _region_edges.clear();
  _region_elem_edge_in_edges_index.clear();
  _region_neighbors.clear();
  _region_bounding_box = std::make_pair(Point(), Point());

  _hanging_node_on_elem_side.clear();
  _hanging_node_on_elem_edge.clear();
}


void SimulationRegion::reserve_data_block(unsigned int n_cell_data, unsigned int n_node_data)
{
  _cell_data_storage.reserve(n_cell_data);
  _node_data_storage.reserve(n_node_data);
}

void SimulationRegion::rebuild_region_fvm_node_list()
{
  _region_local_node.clear();
  _region_processor_node.clear();
  // fill on_local and on_processor node vector
  for(std::map<unsigned int, FVM_Node *>::iterator nodes_it = _region_node.begin(); nodes_it != _region_node.end(); nodes_it++)
  {
    FVM_Node * fvm_node = (*nodes_it).second;
    if( fvm_node->on_local() )
    {
      genius_assert(fvm_node->node_data());
      _region_local_node.push_back(fvm_node);
    }
    if( fvm_node->on_processor() ) _region_processor_node.push_back(fvm_node);
  }
}


void SimulationRegion::prepare_for_use()
{
  // first, we set std::map< const Node *, FVM_Node * > _node_neighbor for FVM_Node
  std::map<unsigned int, FVM_Node *>::iterator nodes_it = _region_node.begin();
  for(; nodes_it != _region_node.end(); ++nodes_it)
  {
    FVM_Node * fvm_node = (*nodes_it).second;

    // skip nonlocal fvm_node
    if( !fvm_node->on_local() ) continue;

    FVM_Node::fvm_neighbor_node_iterator  nb_fvm_node_it = fvm_node->neighbor_node_begin();
    for(; nb_fvm_node_it!=fvm_node->neighbor_node_end(); ++nb_fvm_node_it)
      fvm_node->set_node_neighbor( (*nb_fvm_node_it).first, region_fvm_node((*nb_fvm_node_it).first) );
  }

  // for efficient reason, let each element hold pointer to corresponding FVM_Node
  // since _region_cell is <const Elem *>, we do const_cast here
  element_iterator elem_it = elements_begin();
  for(; elem_it != elements_end(); elem_it++)
  {
    Elem * e = const_cast<Elem *> (*elem_it);
    for(unsigned int i=0; i<e->n_nodes(); i++)
      e->hold_fvm_node( i, region_fvm_node( e->get_node(i)) );
  }

  // fix extra on local cells
  // neighbor elem of an on processor element is marked as on local previously
  // they are set to hold FVM_Node
  // now, we can remove them from region cells if they don't have on processor node
  {
    std::vector<const Elem *> cells = _region_cell;
    _region_cell.clear();
    for(unsigned int n=0; n<cells.size(); ++n)
    {
      const Elem * elem =  cells[n];
      if(elem->on_processor())
      {_region_cell.push_back(elem); continue;}
      for( unsigned int m=0; m<elem->n_nodes(); m++ )
        if( elem->get_node(m)->on_processor() )
        {
          _region_cell.push_back(elem);
          break;
        }
    }
  }

  // build these two vector for fast iteration
  rebuild_region_fvm_node_list();


  // build region edges
  {
    typedef std::map<std::pair<unsigned int, unsigned int>, std::vector<std::pair<const Elem*, unsigned int> > > EdgeCellMap;
    EdgeCellMap region_edge_map;
    for(elem_it = elements_begin(); elem_it != elements_end(); elem_it++)
    {
      const Elem * elem = *elem_it; // elem are on local
      for(unsigned int n=0; n<elem->n_edges(); ++n)
      {
        std::pair<unsigned int, unsigned int> local_edge_nodes;
        elem->nodes_on_edge (n, local_edge_nodes);

        unsigned int node1_id = elem->get_node(local_edge_nodes.first)->id();
        unsigned int node2_id = elem->get_node(local_edge_nodes.second)->id();
        std::pair<unsigned int, unsigned int> edge_nodes = std::make_pair(node1_id, node2_id);
        if( edge_nodes.first > edge_nodes.second )
          std::swap(edge_nodes.first, edge_nodes.second);
        region_edge_map[edge_nodes].push_back(std::make_pair(elem, n));
      }
    }
    EdgeCellMap::const_iterator  edge_it = region_edge_map.begin();
    for(; edge_it != region_edge_map.end(); ++edge_it)
    {
      const std::pair<unsigned int, unsigned int> & edge_nodes = edge_it->first;
      unsigned int edge_index =  _region_edges.size();
      _region_edges.push_back( std::make_pair(region_fvm_node(edge_nodes.first), region_fvm_node(edge_nodes.second)) );

      const std::vector<std::pair<const Elem*, unsigned int> > & elem_shares_edge = edge_it->second;
      for(unsigned int n=0; n<elem_shares_edge.size(); ++n)
      {
        const Elem * elem = elem_shares_edge[n].first;
        unsigned int local_edge_index = elem_shares_edge[n].second;

        if(_region_elem_edge_in_edges_index.find(elem) == _region_elem_edge_in_edges_index.end())
        {
          _region_elem_edge_in_edges_index[elem].resize(elem->n_edges());
          _region_elem_edge_in_edges_index[elem][local_edge_index] = edge_index;
        }
        else
          _region_elem_edge_in_edges_index[elem][local_edge_index] = edge_index;
      }
    }
  }


  // find region neighbors
  {
    std::set<unsigned int> neighbor_region_id;
    for(elem_it = elements_begin(); elem_it != elements_end(); elem_it++)
    {
      const Elem * elem = *elem_it; // elem are on local
      if( !elem->on_processor() ) continue; //

      for(unsigned int n=0; n<elem->n_sides(); ++n)
      {
        const Elem * neighbor_elem = elem->neighbor(n);
        if(neighbor_elem && neighbor_elem->subdomain_id() != _subdomain_id)
          neighbor_region_id.insert( neighbor_elem->subdomain_id() );
      }
    }
    Parallel::allgather(neighbor_region_id);
    std::set<unsigned int>::const_iterator it = neighbor_region_id.begin();
    for( ; it != neighbor_region_id.end(); ++it)
      _region_neighbors.push_back( _subdomain_id_to_region_map.find(*it)->second );
  }


  // build bounding box of the region
  {
    std::vector<Real> min(3, 1.e30);
    std::vector<Real> max(3, -1.e30);

    for(nodes_it = _region_node.begin(); nodes_it != _region_node.end(); ++nodes_it)
    {
      const FVM_Node * fvm_node = (*nodes_it).second;
      const Node * node = fvm_node->root_node();
      for (unsigned int i=0; i<3; i++)
      {
        min[i] = std::min(min[i], (*node)(i));
        max[i] = std::max(max[i], (*node)(i));
      }
    }
    Parallel::min(min);
    Parallel::max(max);
    _region_bounding_box = std::make_pair(Point(&min[0]), Point(&max[0]));
  }

  // set fvm_node volumn for all the FVM_Node
  {
    std::map<unsigned int, Real> fvm_node_volumn_map;
    for(nodes_it = _region_node.begin(); nodes_it != _region_node.end(); ++nodes_it)
    {
      const FVM_Node * fvm_node = (*nodes_it).second;
      if( !fvm_node->on_processor() ) continue;
      fvm_node_volumn_map.insert( std::make_pair(fvm_node->root_node()->id(), fvm_node->volume()) );
    }

    Parallel::allgather(fvm_node_volumn_map);

    for(nodes_it = _region_node.begin(); nodes_it != _region_node.end(); ++nodes_it)
    {
      FVM_Node * fvm_node = (*nodes_it).second;
      fvm_node->set_control_volume( fvm_node_volumn_map.find(fvm_node->root_node()->id())->second );
    }
  }
}


bool SimulationRegion::is_neighbor(const SimulationRegion *r) const
{
  for(unsigned int n=0; n<_region_neighbors.size(); ++n)
    if(_region_neighbors[n] == const_cast<SimulationRegion *>(r) )
      return true;
  return false;
}


void SimulationRegion::remove_remote_object()
{
  std::vector<unsigned int> remote_nodes;

  std::map< unsigned int, FVM_Node * >::iterator it = _region_node.begin();
  for( ; it != _region_node.end(); ++it)
  {
    FVM_Node * fvm_node = it->second;
    if( !fvm_node->on_local() )
    {
      delete fvm_node;
      remote_nodes.push_back(it->first);
    }
  }

  for(unsigned int n=0; n<remote_nodes.size(); ++n)
    _region_node.erase( remote_nodes[n] );
}




unsigned int SimulationRegion::add_variable(const SimulationVariable &v)
{
  std::string var_name = v.variable_name;
  unsigned int var_index = v.variable_index;

  if( v.variable_data_location == POINT_CENTER )
  {
    // when variable not exist, add it to variable map
    if(_region_point_variables.find(var_name) == _region_point_variables.end() )
      _region_point_variables[var_name] = v;

    // when the variable already exist, check memory allocation
    std::map<std::string, SimulationVariable>::iterator it = _region_point_variables.find(var_name);

    switch( v.variable_data_type )
    {
        case SCALAR  :  it->second.variable_index = _node_data_storage.add_scalar_variable(var_index, it->second.variable_valid); break;
        case COMPLEX :  it->second.variable_index = _node_data_storage.add_complex_variable(var_index, it->second.variable_valid); break;
        case VECTOR  :  it->second.variable_index = _node_data_storage.add_vector_variable(var_index, it->second.variable_valid); break;
        case TENSOR  :  it->second.variable_index = _node_data_storage.add_tensor_variable(var_index, it->second.variable_valid); break;
    }

    genius_assert(it->second.variable_index!=invalid_uint);
    return it->second.variable_index;
  }

  if( v.variable_data_location == CELL_CENTER )
  {
    // when variable not exist, add it to variable map
    if(_region_cell_variables.find(var_name) == _region_cell_variables.end() )
      _region_cell_variables[var_name] = v;

    // when the variable already exist, check memory allocation
    std::map<std::string, SimulationVariable>::iterator it = _region_cell_variables.find(var_name);

    switch( v.variable_data_type )
    {
        case SCALAR  :  it->second.variable_index = _cell_data_storage.add_scalar_variable(var_index, it->second.variable_valid); break;
        case COMPLEX :  it->second.variable_index = _cell_data_storage.add_complex_variable(var_index, it->second.variable_valid); break;
        case VECTOR  :  it->second.variable_index = _cell_data_storage.add_vector_variable(var_index, it->second.variable_valid); break;
        case TENSOR  :  it->second.variable_index = _cell_data_storage.add_tensor_variable(var_index, it->second.variable_valid); break;
    }

    genius_assert(it->second.variable_index!=invalid_uint);
    return it->second.variable_index;
  }

  return invalid_uint;
}


bool SimulationRegion::add_variable(const std::string &var_name, DataLocation location)
{
  if( location == POINT_CENTER )
  {
    // when variable not exist, add it to variable map
    if(_region_point_variables.find(var_name) == _region_point_variables.end() )
      return false;

    // when the variable already exist, check memory allocation
    std::map<std::string, SimulationVariable>::iterator it = _region_point_variables.find(var_name);
    if(!it->second.variable_valid) it->second.variable_valid = true;
    switch( it->second.variable_data_type )
    {
        case SCALAR  :  it->second.variable_index = _node_data_storage.add_scalar_variable(it->second.variable_index, it->second.variable_valid); break;
        case COMPLEX :  it->second.variable_index = _node_data_storage.add_complex_variable(it->second.variable_index, it->second.variable_valid); break;
        case VECTOR  :  it->second.variable_index = _node_data_storage.add_vector_variable(it->second.variable_index, it->second.variable_valid); break;
        case TENSOR  :  it->second.variable_index = _node_data_storage.add_tensor_variable(it->second.variable_index, it->second.variable_valid); break;
    }
  }

  if( location == CELL_CENTER )
  {
    // when variable not exist, add it to variable map
    if(_region_cell_variables.find(var_name) == _region_cell_variables.end() )
      return false;

    // when the variable already exist, check memory allocation
    std::map<std::string, SimulationVariable>::iterator it = _region_cell_variables.find(var_name);
    if(!it->second.variable_valid) it->second.variable_valid = true;
    switch( it->second.variable_data_type )
    {
        case SCALAR  :  it->second.variable_index = _cell_data_storage.add_scalar_variable(it->second.variable_index, it->second.variable_valid); break;
        case COMPLEX :  it->second.variable_index = _cell_data_storage.add_complex_variable(it->second.variable_index, it->second.variable_valid); break;
        case VECTOR  :  it->second.variable_index = _cell_data_storage.add_vector_variable(it->second.variable_index, it->second.variable_valid); break;
        case TENSOR  :  it->second.variable_index = _cell_data_storage.add_tensor_variable(it->second.variable_index, it->second.variable_valid); break;
    }
  }

  return true;
}



bool SimulationRegion::has_variable(const std::string &var_name, DataLocation location) const
{
  if( location == POINT_CENTER )
  {
    if( _region_point_variables.find(var_name) == _region_point_variables.end() ) return false;
  }

  if( location == CELL_CENTER )
  {
    if( _region_cell_variables.find(var_name) == _region_cell_variables.end() ) return false;
  }

  return true;
}


const SimulationVariable & SimulationRegion::get_variable(const std::string &var_name, DataLocation location) const
{
  if( location == POINT_CENTER )
  {
    genius_assert( _region_point_variables.find(var_name) != _region_point_variables.end() );
    return _region_point_variables.find(var_name)->second;
  }

  if( location == CELL_CENTER )
  {
    genius_assert( _region_cell_variables.find(var_name) != _region_cell_variables.end() );
    return _region_cell_variables.find(var_name)->second;
  }

  // prevent compiler warning
  genius_assert(0);
  return _region_point_variables.find("")->second;
}

bool SimulationRegion::get_variable(const std::string &var_name, DataLocation location, SimulationVariable & v) const
{
  if( location == POINT_CENTER )
  {
    if( _region_point_variables.find(var_name) == _region_point_variables.end() ) return false;
    v = _region_point_variables.find(var_name)->second;
  }

  if( location == CELL_CENTER )
  {
    if( _region_cell_variables.find(var_name) == _region_cell_variables.end() ) return false;
    v = _region_cell_variables.find(var_name)->second;
  }

  return true;
}


void SimulationRegion::get_user_defined_variable(DataLocation location, DataType type, std::vector<SimulationVariable> & v) const
{
  if( location == POINT_CENTER )
  {
    std::map<std::string, SimulationVariable>::const_iterator it = _region_point_variables.begin();
    for( ; it != _region_point_variables.end(); ++it)
      if( it->second.variable_user_defined )
        v.push_back(it->second);
  }

  if( location == CELL_CENTER )
  {
    std::map<std::string, SimulationVariable>::const_iterator it = _region_cell_variables.begin();
    for( ; it != _region_cell_variables.end(); ++it)
      if( it->second.variable_user_defined )
        v.push_back(it->second);
  }
}



template <typename T>
bool SimulationRegion::get_variable_data(const std::string &var_name, DataLocation location, std::vector<T> &sv) const
{
  parallel_only();

  if( location == POINT_CENTER )
  {
    if( _region_point_variables.find(var_name) == _region_point_variables.end() ) return false;

    const SimulationVariable & variable = _region_point_variables.find(var_name)->second;
    unsigned int variable_index = variable.variable_index;
    Real unit = variable.variable_unit;

    std::map<unsigned int, T> value;
    for(unsigned int n=0; n<_region_processor_node.size(); ++n)
    {
      const FVM_Node * fvm_node = _region_processor_node[n];
      unsigned int offset = fvm_node->node_data()->offset();
      value.insert( std::make_pair(fvm_node->root_node()->id(), _node_data_storage.data<T>(variable_index, offset)/unit ) );
    }
    Parallel::allgather(value);

    sv.reserve(value.size());
    typename std::map<unsigned int, T>::const_iterator it = value.begin();
    for( ; it != value.end(); ++it)
      sv.push_back( it->second );
  }

  if( location == CELL_CENTER )
  {
    if( _region_cell_variables.find(var_name) == _region_cell_variables.end() ) return false;

    const SimulationVariable & variable = _region_cell_variables.find(var_name)->second;
    unsigned int variable_index = variable.variable_index;
    Real unit = variable.variable_unit;

    std::map<unsigned int, T> value;
    for(unsigned int n=0; n<_region_cell.size(); ++n)
    {
      const Elem * elem = _region_cell[n];
      if( elem->on_processor() )
      {
        unsigned int offset = _region_cell_data[n]->offset();
        value.insert( std::make_pair(elem->id(), _cell_data_storage.data<T>(variable_index, offset)/unit ) );
      }
    }
    Parallel::allgather(value);

    sv.reserve(value.size());
    typename std::map<unsigned int, T>::const_iterator it = value.begin();
    for( ; it != value.end(); ++it)
      sv.push_back( it->second );
  }

  return true;

}


void SimulationRegion::set_pmi(const std::string &type, const std::string &model_name, std::vector<Parser::Parameter> & pmi_parameters)
{
  get_material_base()->set_pmi(type,model_name,pmi_parameters);

  local_node_iterator it = on_local_nodes_begin();
  for ( ; it!=on_local_nodes_end(); ++it)
  {
    FVM_Node * node = (*it);
    FVM_NodeData * node_data = (*it)->node_data();
    genius_assert(node_data!=NULL);
    get_material_base()->init_node(type, (*it)->root_node(), node_data);
  }
}


std::string SimulationRegion::get_pmi_info(const std::string& type, const int verbosity)
{
  std::stringstream output;

  output << "-------------------------------------------------PHYSICAL MODEL REPORT--------------------------------------------------"
  << std::endl;
  output << "  Region Name   : " <<  name() << std::endl;
  output << "  Material Name : " << material() << std::endl;
  output << "------------------------------------------------------------------------------------------------------------------------"
  << std::endl;
  output << get_material_base()->get_pmi_info(type, verbosity);
  output << "------------------------------------------------------------------------------------------------------------------------";
  return output.str();
}


bool SimulationRegion::has_2d_hanging_node() const
{
  bool flag = ( _hanging_node_on_elem_side.size() != 0 && _hanging_node_on_elem_edge.size() == 0 );
  Parallel::max(flag);
  return flag;
}

bool SimulationRegion::has_3d_hanging_node() const
{
  bool flag = ( _hanging_node_on_elem_edge.size() != 0 );
  Parallel::max(flag);
  return flag;
}

//explicit instantiation
template
bool SimulationRegion::get_variable_data<Real>(const std::string &var_name, DataLocation location, std::vector<Real> &sv) const;

template
bool SimulationRegion::get_variable_data<Complex>(const std::string &var_name, DataLocation location, std::vector<Complex> &sv) const;

template
bool SimulationRegion::get_variable_data< VectorValue<Real> >(const std::string &var_name, DataLocation location, std::vector< VectorValue<Real> > &sv) const;

template
bool SimulationRegion::get_variable_data< TensorValue<Real> >(const std::string &var_name, DataLocation location, std::vector< TensorValue<Real> > &sv) const;

