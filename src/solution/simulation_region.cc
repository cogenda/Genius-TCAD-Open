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
#include "simulation_system.h"
#include "simulation_region.h"
#include "material.h"
#include "parallel.h"

SimulationRegion::SimulationRegion(const std::string &name, const std::string &material, SimulationSystem & system)
  :_region_name(name), _region_material(material), _system(system)
{}


SimulationRegion::~SimulationRegion()
{
  this->clear();
}


PetscScalar SimulationRegion::T_external() const
{ return _system.T_external(); }


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
  const_node_iterator it = _region_node.begin();
  for(; it != _region_node.end(); ++it)
    if( (*it).second->root_node()->processor_id() == Genius::processor_id() ) n_node++;
  return n_node;
}

void SimulationRegion::clear()
{
  _region_cell.clear();

  for(unsigned int n=0; n<_region_cell_data.size(); ++n)
    delete _region_cell_data[n];
  _region_cell_data.clear();

  node_iterator it = _region_node.begin();

  for( ; it != _region_node.end(); it++ )
  {
    delete (*it).second;
    (*it).second = NULL;
  }

  _region_node.clear();

  _hanging_node_on_elem_side.clear();

  _hanging_node_on_elem_edge.clear();
}


void SimulationRegion::prepare_for_use()
{

    // first, we set std::map< const Node *, FVM_Node * > _node_neighbor for FVM_Node
  node_iterator nodes_it = nodes_begin();
  for(; nodes_it != nodes_end(); nodes_it++)
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
}

void SimulationRegion::set_pmi(const std::string &type, const std::string &model_name,
                               const std::vector<Parser::Parameter> & pmi_parameters)
{
  get_material_base()->set_pmi(type,model_name,pmi_parameters);
  node_iterator it = _region_node.begin();

  for ( ; it!=_region_node.end(); ++it)
  {
      // when node_data is not NULL, this node belongs to current process
      // or at least it is a ghost node
    if((*it).second->node_data() == NULL) continue;

    FVM_NodeData * node_data = (*it).second->node_data();
    get_material_base()->init_node(type,(*it).second->root_node(), node_data);
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


bool SimulationRegion::is_2d_hanging_node() const
{
  bool flag = ( _hanging_node_on_elem_side.size() != 0 && _hanging_node_on_elem_edge.size() == 0 );
  Parallel::max(flag);
  return flag;
}

bool SimulationRegion::is_3d_hanging_node() const
{
  bool flag = ( _hanging_node_on_elem_edge.size() != 0 );
  Parallel::max(flag);
  return flag;
}

