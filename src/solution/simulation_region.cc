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
#include "boundary_condition.h"
#include "material.h"
#include "parallel.h"

// static member
std::map<unsigned int,  SimulationRegion *>  SimulationRegion::_subdomain_id_to_region_map;



SimulationRegion::SimulationRegion(const std::string &name, const std::string &material, const double T, unsigned int dim, const double z)
  :_region_name(name), _region_material(material), _T_external(T), _mesh_dim(dim), _z_width(z)
{}


SimulationRegion::~SimulationRegion()
{
  this->clear();
}


double SimulationRegion::T_external() const
{ return _T_external; }



double SimulationRegion::z_width() const
{ return _z_width; }


void SimulationRegion::add_boundary(BoundaryCondition * bc)
{
  _region_boundaries[bc->boundary_id()] = bc;
}


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
  return _region_processor_node.size();
}


void SimulationRegion::region_node(std::vector<unsigned int> & nodes) const
{
  parallel_only();

  nodes.clear();

  for(unsigned int n=0; n<_region_processor_node.size(); ++n)
  {
    const FVM_Node * fvm_node = _region_processor_node[n];
    nodes.push_back(fvm_node->root_node()->id());
  }

  Parallel::allgather( nodes );
  std::sort(nodes.begin(), nodes.end());
}


void SimulationRegion::region_node_vol(std::vector<double> & vols) const
{
  parallel_only();

  vols.clear();

  std::map<unsigned int, double> node_vol_map;
  for(unsigned int n=0; n<_region_processor_node.size(); ++n)
  {
    const FVM_Node * fvm_node = _region_processor_node[n];
    node_vol_map[fvm_node->root_node()->id()]= fvm_node->volume();
  }

  Parallel::allgather( node_vol_map );

  std::map<unsigned int, double>::const_iterator it=node_vol_map.begin();
  for( ; it!=node_vol_map.end(); ++it)
    vols.push_back(it->second);
}


FVM_Node * SimulationRegion::region_fvm_node(const Node* node) const
{
  std::map<unsigned int, FVM_Node *>::const_iterator it = _region_node.find( node->id() );
  if( it!=_region_node.end() )
    return (*it).second;
  return NULL;
}


FVM_Node * SimulationRegion::region_fvm_node(unsigned int id) const
{
  std::map<unsigned int, FVM_Node *>::const_iterator it = _region_node.find( id );
  if( it!=_region_node.end() )
    return (*it).second;
  return NULL;
}


FVM_NodeData * SimulationRegion::region_node_data(const Node* node) const
{
  std::map<unsigned int, FVM_Node *>::const_iterator it = _region_node.find( node->id() );
  if( it!=_region_node.end() )
    return (*it).second->node_data();
  return NULL;
}


FVM_NodeData * SimulationRegion::region_node_data(unsigned int id) const
{
  std::map<unsigned int, FVM_Node *>::const_iterator it = _region_node.find( id );
  if( it!=_region_node.end() )
    return (*it).second->node_data();
  return NULL;
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
  _region_ghost_node.clear();
  _region_image_node.clear();


  _cell_data_storage.clear();
  _node_data_storage.clear();

  _region_edges.clear();
  _region_elem_edge_in_edges_index.clear();
  _region_neighbors.clear();
  _region_boundaries.clear();
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
  _region_ghost_node.clear();
  _region_image_node.clear();


  // fill on_local and on_processor node vector
  for(std::map<unsigned int, FVM_Node *>::iterator nodes_it = _region_node.begin(); nodes_it != _region_node.end(); nodes_it++)
  {
    FVM_Node * fvm_node = (*nodes_it).second;
    if( fvm_node->on_local() )
    {
      genius_assert(fvm_node->node_data());
      _region_local_node.push_back(fvm_node);
    }
    if( fvm_node->on_processor() )
      _region_processor_node.push_back(fvm_node);
    if( fvm_node->on_local() && !fvm_node->on_processor())
      _region_ghost_node.push_back(fvm_node);
  }

  std::set<unsigned int> ghost_nodes;
  for(unsigned int n=0; n<_region_ghost_node.size(); ++n)
    ghost_nodes.insert( _region_ghost_node[n]->root_node()->id() );
  Parallel::allgather(ghost_nodes);

  std::set<unsigned int>::const_iterator it = ghost_nodes.begin();
  for( ; it != ghost_nodes.end(); ++it)
  {
    unsigned int id = *it;
    if( _region_node.find(id) == _region_node.end() ) continue;

    FVM_Node * fvm_node = _region_node.find(id)->second;
    if( fvm_node->on_processor() )
      _region_image_node.push_back(fvm_node);
  }

}


void SimulationRegion::prepare_for_use()
{
  START_LOG("prepare_for_use()", "SimulationRegion");

  std::map<unsigned int, FVM_Node *>::iterator nodes_it = _region_node.begin();
  for(; nodes_it != _region_node.end(); ++nodes_it)
  {
    FVM_Node * fvm_node = (*nodes_it).second;
    fvm_node->prepare_for_use();

    // skip not on processor fvm_node
    if( !fvm_node->on_processor() ) continue;

     // fix FVM_Node if laplace operator < 0.0
    if( fvm_node->laplace_unit() < 0.0 )
    {
      fvm_node->truncate_cv_surface_area();
    }

    // fix FVM_Node if it has a negative surface area to its neighbor
    fvm_node->posotive_cv_surface_area_to_each_neighbor(true);
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



  // build region edges
  {
    typedef std::map<std::pair<unsigned int, unsigned int>, std::vector<std::pair<const Elem*, unsigned int> > > EdgeCellMap;
    EdgeCellMap region_edge_map;
    for(element_iterator elem_it = elements_begin(); elem_it != elements_end(); elem_it++)
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

  STOP_LOG("prepare_for_use()", "SimulationRegion");
}


void SimulationRegion::prepare_for_use_parallel()
{
  START_LOG("prepare_for_use_parallel()", "SimulationRegion");

  // build these two vector for fast iteration
  rebuild_region_fvm_node_list();

  // statistic all the region nodes
  _n_region_node = this->n_on_processor_node();
  Parallel::sum(_n_region_node);

  // build bounding box of the region
  {
    std::vector<Real> min(3, 1.e30);
    std::vector<Real> max(3, -1.e30);

    processor_node_iterator nodes_it = on_processor_nodes_begin();
    for(; nodes_it != on_processor_nodes_end(); ++nodes_it)
    {
      const FVM_Node * fvm_node = *nodes_it;
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

  // find region neighbors
  {
    std::set<unsigned int> neighbor_region_id;
    for(element_iterator elem_it = elements_begin(); elem_it != elements_end(); elem_it++)
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

  sync_fvm_node_volume();

  // hanging node flag
  {
    std::vector<int> hanging_node_flags;
    hanging_node_flags.push_back(static_cast<int>( _hanging_node_on_elem_side.size() != 0 && _hanging_node_on_elem_edge.size() == 0 ));
    hanging_node_flags.push_back(static_cast<int>(  _hanging_node_on_elem_edge.size() != 0 ));
    Parallel::max(hanging_node_flags);
    _hanging_node_on_elem_side_flag = static_cast<bool>(hanging_node_flags[0]);
    _hanging_node_on_elem_edge_flag = static_cast<bool>(hanging_node_flags[1]);
  }

  STOP_LOG("prepare_for_use_parallel()", "SimulationRegion");
}



void SimulationRegion::sync_fvm_node_volume()
{
  // reset fvm_node volume for all the FVM_Node (also sync ghost nodes)
  //NOTE zero/negative fvm_node volume due to bad mesh will cause simulation fail.
  // here we force all the fvm_node volume to be positive
  {
    std::map<unsigned int, Real> fvm_node_volume_map;
    std::map<unsigned int, FVM_Node *>::iterator nodes_it = _region_node.begin();
    for(; nodes_it != _region_node.end(); ++nodes_it)
    {
      const FVM_Node * fvm_node = (*nodes_it).second;
      if( !fvm_node->on_processor() ) continue;
      fvm_node_volume_map.insert( std::make_pair(fvm_node->root_node()->id(), fvm_node->volume()) );
    }

    Parallel::allgather(fvm_node_volume_map);
    const double min_fvm_volume = std::pow(0.01*PhysicalUnit::nm, (double)_mesh_dim);
    for(nodes_it = _region_node.begin(); nodes_it != _region_node.end(); ++nodes_it)
    {
      FVM_Node * fvm_node = (*nodes_it).second;
      fvm_node->set_control_volume( std::max(min_fvm_volume, std::abs(fvm_node_volume_map.find(fvm_node->root_node()->id())->second)) );
    }
  }

}



Real SimulationRegion::fvm_cell_quality() const
{
#if 0
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    if( fvm_node->total_cv_surface_area() < 0.0 )
    {
      std::cout<< *fvm_node << std::endl;
    }
  }
#endif
  return 1.0;
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
    PetscScalar unit = variable.variable_unit;

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
    PetscScalar unit = variable.variable_unit;

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


template <typename T>
bool SimulationRegion::set_variable_data(const std::string &var_name, DataLocation location, const T val)
{
  parallel_only();

  if( location == POINT_CENTER )
  {
    if( _region_point_variables.find(var_name) == _region_point_variables.end() ) return false;

    const SimulationVariable & variable = _region_point_variables.find(var_name)->second;
    unsigned int variable_index = variable.variable_index;
    PetscScalar unit = variable.variable_unit;

    for(unsigned int n=0; n<_region_processor_node.size(); ++n)
    {
      const FVM_Node * fvm_node = _region_processor_node[n];
      unsigned int offset = fvm_node->node_data()->offset();
      _node_data_storage.data<T>(variable_index, offset) = val*unit;
    }
  }

  if( location == CELL_CENTER )
  {
    if( _region_cell_variables.find(var_name) == _region_cell_variables.end() ) return false;

    const SimulationVariable & variable = _region_cell_variables.find(var_name)->second;
    unsigned int variable_index = variable.variable_index;
    PetscScalar unit = variable.variable_unit;

    std::map<unsigned int, T> value;
    for(unsigned int n=0; n<_region_cell.size(); ++n)
    {
      const Elem * elem = _region_cell[n];
      if( elem->on_processor() )
      {
        unsigned int offset = _region_cell_data[n]->offset();
        _cell_data_storage.data<T>(variable_index, offset) = val*unit;
      }
    }
  }

  return true;

}


template <typename T>
bool SimulationRegion::set_variable_data(const std::string &var_name, DataLocation location, const T val, const double unit)
{
  parallel_only();

  if( location == POINT_CENTER )
  {
    if( _region_point_variables.find(var_name) == _region_point_variables.end() ) return false;

    const SimulationVariable & variable = _region_point_variables.find(var_name)->second;
    unsigned int variable_index = variable.variable_index;

    for(unsigned int n=0; n<_region_processor_node.size(); ++n)
    {
      const FVM_Node * fvm_node = _region_processor_node[n];
      unsigned int offset = fvm_node->node_data()->offset();
      _node_data_storage.data<T>(variable_index, offset) = val*(PetscScalar)unit;
    }
  }

  if( location == CELL_CENTER )
  {
    if( _region_cell_variables.find(var_name) == _region_cell_variables.end() ) return false;

    const SimulationVariable & variable = _region_cell_variables.find(var_name)->second;
    unsigned int variable_index = variable.variable_index;

    std::map<unsigned int, T> value;
    for(unsigned int n=0; n<_region_cell.size(); ++n)
    {
      const Elem * elem = _region_cell[n];
      if( elem->on_processor() )
      {
        unsigned int offset = _region_cell_data[n]->offset();
        _cell_data_storage.data<T>(variable_index, offset) = val*(PetscScalar)unit;
      }
    }
  }

  return true;

}



template <typename T>
bool SimulationRegion::sync_point_variable(const std::string &var_name)
{
  parallel_only();

  if( _region_point_variables.find(var_name) == _region_point_variables.end() ) return false;

  const SimulationVariable & variable = _region_point_variables.find(var_name)->second;
  unsigned int variable_index = variable.variable_index;

  std::map<unsigned int, T> ghost_values;
  for(unsigned int n=0; n<_region_image_node.size(); ++n)
  {
    const FVM_Node * fvm_node = _region_image_node[n];
    unsigned int offset = fvm_node->node_data()->offset();
    ghost_values.insert( std::make_pair(fvm_node->root_node()->id(), _node_data_storage.data<T>(variable_index, offset)) );
  }
  Parallel::allgather(ghost_values);

  for(unsigned int n=0; n<_region_ghost_node.size(); ++n)
  {
    const FVM_Node * fvm_node = _region_ghost_node[n];
    T value = ghost_values.find(fvm_node->root_node()->id())->second;
    unsigned int offset = fvm_node->node_data()->offset();
    _node_data_storage.data<T>(variable_index, offset) = value;
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



BoundaryCondition * SimulationRegion::floating_metal() const
{
  std::map<short int, BoundaryCondition *>::const_iterator it = _region_boundaries.begin();
  std::map<short int, BoundaryCondition *>::const_iterator it_end = _region_boundaries.end();
  for(; it!=it_end; it++)
  {
    BoundaryCondition * bc = it->second;
    if( bc->bc_type() == ChargedContact )
      return bc->inter_connect_hub();
  }
  return 0;
}



void SimulationRegion::add_hanging_node_on_side(const Node * node, const Elem * elem, unsigned int s)
{
  std::map<unsigned int, FVM_Node *>::iterator it = _region_node.find(node->id());
  genius_assert( it!=_region_node.end() );

  const FVM_Node * fvm_node = (*it).second;
  _hanging_node_on_elem_side[fvm_node] = std::pair<const Elem *, unsigned int>(elem, s);
}


void SimulationRegion::add_hanging_node_on_edge(const Node * node, const Elem * elem, unsigned int e)
{
  std::map<unsigned int, FVM_Node *>::iterator it = _region_node.find(node->id());
  genius_assert( it!=_region_node.end() );

  const FVM_Node * fvm_node = (*it).second;
  _hanging_node_on_elem_edge[fvm_node] = std::pair<const Elem *, unsigned int>(elem, e);

}


size_t SimulationRegion::memory_size() const
{
  size_t counter = sizeof(*this);

  counter += _region_neighbors.capacity()*sizeof(SimulationRegion *);
  // FIXME std::map<short int, BoundaryCondition *> _region_boundaries;
  counter += _region_cell.capacity()*sizeof(const Elem *);
  counter += _region_cell_data.capacity()*sizeof(FVM_CellData *);
  counter += _cell_data_storage.memory_size();
  std::map<unsigned int, FVM_Node *>::const_iterator it = _region_node.begin();
  for( ; it != _region_node.end(); it++ )
  {
    counter += it->second->memory_size();
    counter += sizeof(it->first);
  }

  counter += _region_local_node.capacity()*sizeof(FVM_Node *);
  counter += _region_processor_node.capacity()*sizeof(FVM_Node *);
  counter += _region_ghost_node.capacity()*sizeof(FVM_Node *);
  counter += _region_image_node.capacity()*sizeof(FVM_Node *);
  counter +=  _node_data_storage.memory_size();
  counter += _region_edges.capacity()*sizeof(std::pair<FVM_Node *, FVM_Node *>);

  return counter;
}




//explicit instantiation
template
bool SimulationRegion::get_variable_data<PetscScalar>(const std::string &var_name, DataLocation location, std::vector<PetscScalar> &sv) const;

template
bool SimulationRegion::get_variable_data< std::complex<PetscScalar> >(const std::string &var_name, DataLocation location, std::vector< std::complex<PetscScalar> > &sv) const;

template
bool SimulationRegion::get_variable_data< VectorValue<PetscScalar> >(const std::string &var_name, DataLocation location, std::vector< VectorValue<PetscScalar> > &sv) const;

template
bool SimulationRegion::get_variable_data< TensorValue<PetscScalar> >(const std::string &var_name, DataLocation location, std::vector< TensorValue<PetscScalar> > &sv) const;

//explicit instantiation
template
bool SimulationRegion::set_variable_data<PetscScalar>(const std::string &var_name, DataLocation location, const PetscScalar sv);

template
bool SimulationRegion::set_variable_data< std::complex<PetscScalar> >(const std::string &var_name, DataLocation location, const std::complex<PetscScalar> sv);

template
bool SimulationRegion::set_variable_data< VectorValue<PetscScalar> >(const std::string &var_name, DataLocation location, const VectorValue<PetscScalar> sv);

template
bool SimulationRegion::set_variable_data< TensorValue<PetscScalar> >(const std::string &var_name, DataLocation location, const TensorValue<PetscScalar> sv);

template
bool SimulationRegion::set_variable_data<PetscScalar>(const std::string &var_name, DataLocation location, const PetscScalar sv, const double unit);

template
bool SimulationRegion::set_variable_data< std::complex<PetscScalar> >(const std::string &var_name, DataLocation location, const std::complex<PetscScalar> sv, const double unit);

template
bool SimulationRegion::set_variable_data< VectorValue<PetscScalar> >(const std::string &var_name, DataLocation location, const VectorValue<PetscScalar> sv, const double unit);

template
bool SimulationRegion::set_variable_data< TensorValue<PetscScalar> >(const std::string &var_name, DataLocation location, const TensorValue<PetscScalar> sv, const double unit);



template
bool SimulationRegion::sync_point_variable<PetscScalar>(const std::string &var_name);

