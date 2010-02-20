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

//  $Id: fvm_node_info.cc,v 1.8 2008/07/09 05:58:16 gdiso Exp $

#include "elem.h"
#include "fvm_node_info.h"
#include "boundary_info.h"


FVM_Node::FVM_Node(const Node *n)
    : _node(n),
    _node_data(0),
    _elem_has_this_node(0),
    _node_neighbor(0),
    _cv_surface_area(0),
    _ghost_nodes(0),
    _volume(0),
    _boundary_id(BoundaryInfo::invalid_id),
    _global_offset(invalid_uint),
    _local_offset(invalid_uint)
{}



void FVM_Node::PDE_node_pattern(std::vector<std::pair<unsigned int, unsigned int> > & v_region_nodes, bool elem_based) const
{
  //only consider neighbor nodes, link this node by an edge
  if( elem_based==false )
  {

    // the PDE node pattern in thie region. n_node_neighbors()+1 means "this" node is also considered.
    v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(subdomain_id (), n_node_neighbors()+1) );

    // consider ghost nodes in other regions
    if( _ghost_nodes!=NULL && !_ghost_nodes->empty() )
      for(fvm_ghost_node_iterator  git = ghost_node_begin(); git!=ghost_node_end(); ++git)
      {
        FVM_Node *ghost_node = (*git).first;
        if(ghost_node)
          v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(ghost_node->subdomain_id (), ghost_node->n_node_neighbors()+1) );
      }
    return;
  }
  // consider all the nodes belongs to neighbor elements
  else
  {
    std::map<unsigned int, std::set<const FVM_Node *> > _region_node_map;

    // search for all the neighbor elements in this region
    fvm_element_iterator element_it = elem_begin();
    for( ; element_it != elem_end(); ++element_it)
    {
      const Elem * e = (*element_it).first;
      for(unsigned int v=0; v<e->n_vertices(); v++)
      {
        const FVM_Node * fvm_node = e->get_fvm_node(v) ;
        _region_node_map[fvm_node->subdomain_id()].insert(fvm_node);
      }
    }

    // consider ghost nodes in other regions
    if( _ghost_nodes!=NULL && !_ghost_nodes->empty() )
    {
      for(fvm_ghost_node_iterator  git = ghost_node_begin(); git != ghost_node_end(); ++git)
      {
        FVM_Node *ghost_node = (*git).first;
        if(ghost_node)
          for(element_it = ghost_node->elem_begin(); element_it != ghost_node->elem_end(); ++element_it)
          {
            const Elem * e = (*element_it).first;
            for(unsigned int v=0; v<e->n_vertices(); v++)
            {
              const FVM_Node * fvm_node = e->get_fvm_node(v) ;
              _region_node_map[fvm_node->subdomain_id()].insert(fvm_node);
            }
          }
      }
    }

    // statistic all the PDE pattern infomation
    std::map<unsigned int, std::set<const FVM_Node *> >::iterator it= _region_node_map.begin();
    for(; it!= _region_node_map.end(); ++it)
      v_region_nodes.push_back(std::pair<unsigned int, unsigned int>((*it).first, (*it).second.size()));
    return;
  }

}



void FVM_Node::PDE_off_processor_node_pattern(std::vector<std::pair<unsigned int, unsigned int> > & v_region_nodes, bool elem_based) const
{
  //only consider neighbor nodes, link this node by an edge
  if( elem_based==false )
  {
    unsigned int nodes = 0;
    //neighbor not on this processor
    fvm_neighbor_node_iterator it= neighbor_node_begin();
    for(; it!= neighbor_node_end(); ++it)
      if( (*it).first->on_local() && (*it).first->processor_id() != Genius::processor_id() )   nodes++;
    v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(subdomain_id (), nodes) );

    // we should consider ghost node, since it has the same root node
    if( _ghost_nodes!=NULL && !_ghost_nodes->empty() )
      for(fvm_ghost_node_iterator  git = ghost_node_begin(); git!=ghost_node_end(); ++git)
      {
        FVM_Node *ghost_node = (*git).first;
        if( ghost_node == NULL) continue;
        unsigned int nodes = 0;
        fvm_neighbor_node_iterator gnit = ghost_node->neighbor_node_begin();
        for(; gnit!= ghost_node->neighbor_node_end(); ++gnit)
          if( (*gnit).first->on_local() && (*gnit).first->processor_id() != Genius::processor_id() )   nodes++;
        v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(ghost_node->subdomain_id (), nodes) );
      }

  }
  // consider all the nodes belongs to neighbor elements
  else
  {
    std::map<unsigned int, std::set<const FVM_Node *> > _region_node_map;

    // search for all the neighbor elements in this region
    fvm_element_iterator element_it = elem_begin();
    for( ; element_it != elem_end(); ++element_it)
    {
      const Elem * e = (*element_it).first;
      for(unsigned int v=0; v<e->n_vertices(); v++)
      {
        const FVM_Node * fvm_node = e->get_fvm_node(v) ;
        if( fvm_node->root_node()->on_local() && fvm_node->root_node()->processor_id() != Genius::processor_id() )
          _region_node_map[fvm_node->subdomain_id()].insert(fvm_node);
      }
    }

    // consider ghost nodes in other regions
    if( _ghost_nodes!=NULL && !_ghost_nodes->empty() )
    {
      for(fvm_ghost_node_iterator  git = ghost_node_begin(); git != ghost_node_end(); ++git)
      {
        FVM_Node *ghost_node = (*git).first;
        if(ghost_node)
          for(element_it = ghost_node->elem_begin(); element_it != ghost_node->elem_end(); ++element_it)
          {
            const Elem * e = (*element_it).first;
            for(unsigned int v=0; v<e->n_vertices(); v++)
            {
              const FVM_Node * fvm_node = e->get_fvm_node(v) ;
              if( fvm_node->root_node()->on_local() && fvm_node->root_node()->processor_id() != Genius::processor_id() )
                _region_node_map[fvm_node->subdomain_id()].insert(fvm_node);
            }
          }
      }
    }

    // statistic all the off-processor PDE pattern infomation
    std::map<unsigned int, std::set<const FVM_Node *> >::iterator it= _region_node_map.begin();
    for(; it!= _region_node_map.end(); ++it)
      v_region_nodes.push_back(std::pair<unsigned int, unsigned int>((*it).first, (*it).second.size()));
    return;
  }

}



void FVM_Node::operator += (const FVM_Node &other_node)
{
  // check if has the same root node
  genius_assert( this->_node == other_node.root_node() );

  // skip nonlocal fvm_node
  if( !this->on_local() ) return;

  // combine element
  fvm_element_iterator elem_it = other_node.elem_begin();
  for( ; elem_it!= other_node.elem_end(); ++elem_it)
    _elem_has_this_node->push_back( *elem_it );

  // combine neighbor node
  fvm_neighbor_node_iterator node_it = other_node.neighbor_node_begin();
  for( ; node_it!= other_node.neighbor_node_end(); ++node_it)
    _node_neighbor->insert( *node_it );

  // add volume
  _volume += other_node.volume();

  // add cv surface
  std::map< const Node *, Real >::const_iterator cv_it = other_node. _cv_surface_area->begin();
  for(; cv_it != other_node._cv_surface_area->end(); ++cv_it)
  {
    if(_cv_surface_area->find((*cv_it).first) != _cv_surface_area->end())
      (*_cv_surface_area->find((*cv_it).first)).second += (*cv_it).second;
    else
      // insert a new item
      (*_cv_surface_area)[(*cv_it).first] = (*cv_it).second;
  }

}
