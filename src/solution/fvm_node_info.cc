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

#include <algorithm>

#include "elem.h"
#include "fvm_node_info.h"
#include "boundary_info.h"

// TNT matrix-vector library
#include <TNT/tnt.h>
#include <TNT/jama_lu.h>


unsigned int FVM_Node::_solver_index=0;


FVM_Node::FVM_Node(const Node *n)
    : _node(n),
    _node_data(0),
    _ghost_nodes(0),
    _volume(0),
    _boundary_id(BoundaryInfo::invalid_id),
    _bc_type(INVALID_BC_TYPE),
    _subdomain_id(invalid_uint)
{
  for(unsigned int n=0; n<4; ++n)
  {
    _global_offset[n] = invalid_uint;
    _local_offset[n] = invalid_uint;
  }
}


FVM_Node::~FVM_Node()
{
  _node = 0;
  delete _node_data;
  delete _ghost_nodes;
}


void FVM_Node::add_elem_it_belongs(const Elem * el, unsigned int n)
{
  _elem_has_this_node.push_back( std::pair<const Elem *,unsigned int>(el,n) );
}


void FVM_Node::add_fvm_node_neighbor(FVM_Node *nb, Real s)
{
  for(unsigned int n=0; n<_fvm_node_neighbor.size(); ++n)
  {
    if( _fvm_node_neighbor[n].first->root_node() == nb->root_node() )
    {
      _fvm_node_neighbor[n].second.first  += s;
      _fvm_node_neighbor[n].second.second += std::abs(s);
      return;
    }
  }

  _fvm_node_neighbor.push_back( std::make_pair(nb, std::make_pair(s, std::abs(s))) );
}


void FVM_Node::set_fvm_node_neighbor(FVM_Node *nb, Real s, Real abs)
{
  for(unsigned int n=0; n<_fvm_node_neighbor.size(); ++n)
  {
    if( _fvm_node_neighbor[n].first->root_node() == nb->root_node() )
    {
      _fvm_node_neighbor[n].second.first  = s;
      _fvm_node_neighbor[n].second.second = abs;
      return;
    }
  }

  _fvm_node_neighbor.push_back( std::make_pair(nb, std::make_pair(s, abs)) );
}


void FVM_Node::set_ghost_node(FVM_Node * fn, unsigned int sub_id, Real area)
{
  if( _ghost_nodes == NULL)
    _ghost_nodes = new std::map< FVM_Node *, std::pair<unsigned int, Real> >;

  std::pair<unsigned int, Real> gf(sub_id,area);
  _ghost_nodes->insert( std::pair< FVM_Node *, std::pair<unsigned int, Real> >(fn, gf) );
}


void FVM_Node::set_ghost_node_area(unsigned int sub_id, Real area)
{
  // the sub_id of boundary face equal to this FVM_Node and _ghost_nodes are empty
  // this is a boundary face, not interface face
  if( sub_id == _subdomain_id && _ghost_nodes==NULL )
  {
    _ghost_nodes = new std::map< FVM_Node *, std::pair<unsigned int, Real> >;
    std::pair<unsigned int, Real> gf(invalid_uint, area);
    _ghost_nodes->insert( std::pair< FVM_Node *, std::pair<unsigned int, Real> >((FVM_Node *)NULL, gf) );
    return;
  }

  // else we find in ghost nodes which matches sub_id
  genius_assert(_ghost_nodes);
  std::map< FVM_Node *, std::pair<unsigned int, Real> >::iterator it = _ghost_nodes->begin();
  for(; it!=_ghost_nodes->end(); ++it)
    if( (*it).second.first ==  sub_id )
    {
      (*it).second.second = area;
    }
}


FVM_Node * FVM_Node::ghost_fvm_node(unsigned int i) const
{
  fvm_ghost_node_iterator  it = ghost_node_begin();
  for(unsigned int l=0; l<i; ++i, ++it);
  return (*it).first;
}



Real FVM_Node::cv_surface_area(const FVM_Node * neighbor) const
{
  unsigned int size = _fvm_node_neighbor.size();
  if( size < 5 )
  {
    for(unsigned int n=0; n<size; ++n)
      if(_fvm_node_neighbor[n].first == neighbor)
        return _fvm_node_neighbor[n].second.first;
  }
  else
  {
    unsigned int left  = 0;
    unsigned int right = size-1;
    while(left <= right)
    {
      unsigned int middle = (left+right)/2;
      if (neighbor < _fvm_node_neighbor[middle].first)
        right = middle - 1;
      else if (neighbor > _fvm_node_neighbor[middle].first)
        left = middle + 1;
      else
        return _fvm_node_neighbor[middle].second.first;
    }
  }
  return 0.0;
}


Real FVM_Node::total_cv_surface_area() const
{
  Real total_surface_area = 0.0;
  for(unsigned int n=0; n<_fvm_node_neighbor.size(); ++n)
  {
    total_surface_area += _fvm_node_neighbor[n].second.first;
  }
  return total_surface_area;
}


Real FVM_Node::cv_abs_surface_area(const FVM_Node * neighbor) const
{
  unsigned int size = _fvm_node_neighbor.size();
  if( size < 5 )
  {
    for(unsigned int n=0; n<size; ++n)
      if(_fvm_node_neighbor[n].first == neighbor)
        return _fvm_node_neighbor[n].second.second;
  }
  else
  {
    unsigned int left  = 0;
    unsigned int right = size-1;
    while(left <= right)
    {
      unsigned int middle = (left+right)/2;
      if (neighbor < _fvm_node_neighbor[middle].first)
        right = middle - 1;
      else if (neighbor > _fvm_node_neighbor[middle].first)
        left = middle + 1;
      else
        return _fvm_node_neighbor[middle].second.second;
    }
  }
  return 0.0;
}


Real FVM_Node::total_abs_cv_surface_area() const
{
  Real total_surface_area = 0.0;
  for(unsigned int n=0; n<_fvm_node_neighbor.size(); ++n)
  {
    total_surface_area += _fvm_node_neighbor[n].second.second;
  }
  return total_surface_area;
}


bool FVM_Node::posotive_cv_surface_area_to_each_neighbor(bool fix)
{
  std::map<const Node *, Real> cv_surface_area_map;

  for( unsigned int n=0; n<_fvm_node_neighbor.size(); ++n)
  {
    const Node * node = _fvm_node_neighbor[n].first->root_node();
    Real  cv = _fvm_node_neighbor[n].second.first;
    cv_surface_area_map[node] += cv;
  }

  if(_ghost_nodes)
  {
    std::map< FVM_Node *, std::pair<unsigned int, Real> >::const_iterator it = _ghost_nodes->begin();
    for(; it != _ghost_nodes->end(); ++it)
    {
      const FVM_Node * ghost_fvm_node = it->first;
      if(!ghost_fvm_node) continue;
      for( unsigned int n=0; n<ghost_fvm_node->_fvm_node_neighbor.size(); ++n)
      {
        const Node * node = ghost_fvm_node->_fvm_node_neighbor[n].first->root_node();
        Real  cv = ghost_fvm_node->_fvm_node_neighbor[n].second.first;
        cv_surface_area_map[node] += cv;
      }
    }
  }


  bool positive = true;
  std::map<const Node *, Real>::const_iterator it = cv_surface_area_map.begin();
  for(; it != cv_surface_area_map.end(); ++it)
    if( it->second < 0.0 ) positive = false;
  if(positive) return true;

  // force to positive
  if(fix)
  {
    for( unsigned int n=0; n<_fvm_node_neighbor.size(); ++n)
    {
      const Node * node = _fvm_node_neighbor[n].first->root_node();
      if( cv_surface_area_map[node] < 0.0 )
      {
        Real cv     = std::abs(_fvm_node_neighbor[n].second.first);
        Real cv_abs = _fvm_node_neighbor[n].second.second;
        // force to positive
        _fvm_node_neighbor[n].second.first  = cv;
        // also update my neighbor
        FVM_Node * neighbor_fvm_node = _fvm_node_neighbor[n].first;
        neighbor_fvm_node->set_fvm_node_neighbor(this, cv, cv_abs);
      }
    }

    if(_ghost_nodes)
    {
      std::map< FVM_Node *, std::pair<unsigned int, Real> >::const_iterator it = _ghost_nodes->begin();
      for(; it != _ghost_nodes->end(); ++it)
      {
        FVM_Node * ghost_fvm_node = it->first;
        if(!ghost_fvm_node) continue;
        for( unsigned int n=0; n<ghost_fvm_node->_fvm_node_neighbor.size(); ++n)
        {
          const Node * node = ghost_fvm_node->_fvm_node_neighbor[n].first->root_node();
          if( cv_surface_area_map[node] < 0.0 )
          {
            Real cv     = std::abs(ghost_fvm_node->_fvm_node_neighbor[n].second.first);
            Real cv_abs = ghost_fvm_node->_fvm_node_neighbor[n].second.second;
            // force to positive
            ghost_fvm_node->_fvm_node_neighbor[n].second.first  = cv;
            // also update my neighbor
            FVM_Node * neighbor_fvm_node = ghost_fvm_node->_fvm_node_neighbor[n].first;
            neighbor_fvm_node->set_fvm_node_neighbor(ghost_fvm_node, cv, cv_abs);
          }
        }
      }
    }
  }

  return false;
}



void FVM_Node::truncate_cv_surface_area()
{
  for(unsigned int n=0; n<_fvm_node_neighbor.size(); ++n)
  {
    if(_fvm_node_neighbor[n].second.first < 0.0)
    {
      Real cv     = std::abs(_fvm_node_neighbor[n].second.first);
      Real cv_abs = _fvm_node_neighbor[n].second.second;
      // force to positive
      _fvm_node_neighbor[n].second.first  = cv;

      // also update my neighbor
      FVM_Node * neighbor_fvm_node = _fvm_node_neighbor[n].first;

      neighbor_fvm_node->set_fvm_node_neighbor(this, cv, cv_abs);
    }
  }
}


bool FVM_Node::on_boundary() const
{
  genius_error();
  return false;
}


std::vector<unsigned int> FVM_Node::subdomains() const
{
  genius_assert(_ghost_nodes);

  std::set<unsigned int> subdomains_set;
  subdomains_set.insert(_subdomain_id);
  std::map< FVM_Node *, std::pair<unsigned int, Real> >::const_iterator it = _ghost_nodes->begin();
  for( ; it != _ghost_nodes->end(); ++it)
  {
    if( !it->first ) continue;
    subdomains_set.insert(  it->second.first );
  }

  std::vector<unsigned int> subdomains_vector(subdomains_set.begin(), subdomains_set.end());
  return subdomains_vector;
}


Real FVM_Node::outside_boundary_surface_area() const
{
  genius_assert(_ghost_nodes);
  Real area = 0.0;
  fvm_ghost_node_iterator  it = ghost_node_begin();
  for(; it!=ghost_node_end(); ++it )
    area += (*it).second.second;
  return area;
}


Real FVM_Node::laplace_unit() const
{
  Real r = 0.0;

  for(unsigned int n=0; n<_fvm_node_neighbor.size(); ++n)
  {
    r+= _fvm_node_neighbor[n].second.first/this->distance(_fvm_node_neighbor[n].first);
  }

  return r/this->volume();
}


Real FVM_Node::sup_laplace_unit() const
{
  Real r = 0.0;

  for(unsigned int n=0; n<_fvm_node_neighbor.size(); ++n)
  {
    r+= std::abs(_fvm_node_neighbor[n].second.second)/this->distance(_fvm_node_neighbor[n].first);
  }

  return r/this->volume();
}


unsigned int FVM_Node::n_pure_ghost_node() const
{
  // sun NULL ghost node
  if( _ghost_nodes->find(NULL)!=_ghost_nodes->end() )
    return _ghost_nodes->size() -1 ;
  return _ghost_nodes->size();
}


void FVM_Node::PDE_node_pattern(std::vector<std::pair<unsigned int, unsigned int> > & v_region_nodes, bool elem_based) const
{
  //only consider neighbor nodes, link this node by an edge
  if( elem_based==false )
  {

    // the PDE node pattern in thie region. fvm_node_neighbors()+1 means "this" node is also considered.
    v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(subdomain_id (), fvm_node_neighbors()+1) );

    // consider ghost nodes in other regions
    if( _ghost_nodes!=NULL && !_ghost_nodes->empty() )
      for(fvm_ghost_node_iterator  git = ghost_node_begin(); git!=ghost_node_end(); ++git)
      {
        FVM_Node *ghost_node = (*git).first;
        if(ghost_node)
          v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(ghost_node->subdomain_id (), ghost_node->fvm_node_neighbors()+1) );
      }
    return;
  }
  // consider all the nodes belongs to neighbor elements
  else
  {
    std::set<const Elem *> elems;
    std::map<unsigned int, std::set<const FVM_Node *> > region_node_map;

    // search for all the neighbor elements in this region
    fvm_element_iterator element_it = elem_begin();
    for( ; element_it != elem_end(); ++element_it)
    {
      const Elem * e = (*element_it).first;
      elems.insert(e);
      for(unsigned int n=0; n<e->n_sides(); ++n)
        if( e->neighbor(n) ) elems.insert(e->neighbor(n));
    }

    // consider ghost nodes in other regions
    if( _ghost_nodes!=NULL && !_ghost_nodes->empty() )
    {
      for(fvm_ghost_node_iterator  git = ghost_node_begin(); git != ghost_node_end(); ++git)
      {
        FVM_Node *ghost_node = (*git).first;
        if(!ghost_node) continue;

        for(element_it = ghost_node->elem_begin(); element_it != ghost_node->elem_end(); ++element_it)
        {
          const Elem * e = (*element_it).first;
          elems.insert(e);
          for(unsigned int n=0; n<e->n_sides(); ++n)
            if( e->neighbor(n) ) elems.insert(e->neighbor(n));
        }
      }
    }

    std::set<const Elem *>::const_iterator elem_it = elems.begin();
    for( ; elem_it != elems.end(); ++elem_it)
    {
      const Elem * e = *elem_it;
      for(unsigned int v=0; v<e->n_vertices(); v++)
      {
        const FVM_Node * fvm_node = e->get_fvm_node(v) ;
        region_node_map[fvm_node->subdomain_id()].insert(fvm_node);
      }
    }

    // statistic all the PDE pattern infomation
    std::map<unsigned int, std::set<const FVM_Node *> >::iterator it= region_node_map.begin();
    for(; it!= region_node_map.end(); ++it)
      v_region_nodes.push_back(std::pair<unsigned int, unsigned int>((*it).first, (*it).second.size()));
    return;
  }

}


void FVM_Node::PDE_node_pattern(std::vector<std::pair<unsigned int, unsigned int> > & v_region_nodes,
                                std::vector<std::pair<unsigned int, unsigned int> > & v_off_region_nodes,
                                bool elem_based) const
{
  //only consider neighbor nodes, link this node by an edge
  if( elem_based==false )
  {
    // the PDE node pattern in thie region. fvm_node_neighbors()+1 means "this" node is also considered.
    v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(subdomain_id (), fvm_node_neighbors()+1) );

    //neighbor not on this processor
    unsigned int off_nodes = 0;
    for(fvm_neighbor_node_iterator it= neighbor_node_begin(); it!= neighbor_node_end(); ++it)
      if( (*it).first->on_local() && (*it).first->root_node()->processor_id() != Genius::processor_id() )   off_nodes++;
    v_off_region_nodes.push_back( std::pair<unsigned int, unsigned int>(subdomain_id (), off_nodes) );

    // consider ghost nodes in other regions
    if( _ghost_nodes!=NULL && !_ghost_nodes->empty() )
    {
      for(fvm_ghost_node_iterator  git = ghost_node_begin(); git!=ghost_node_end(); ++git)
      {
        FVM_Node *ghost_node = (*git).first;
        if( ghost_node == NULL) continue;

        v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(ghost_node->subdomain_id (), ghost_node->fvm_node_neighbors()+1) );

        unsigned int off_nodes = 0;
        fvm_neighbor_node_iterator gnit = ghost_node->neighbor_node_begin();
        for(; gnit!= ghost_node->neighbor_node_end(); ++gnit)
          if( (*gnit).first->on_local() && (*gnit).first->root_node()->processor_id() != Genius::processor_id() )   off_nodes++;
        v_off_region_nodes.push_back( std::pair<unsigned int, unsigned int>(ghost_node->subdomain_id (), off_nodes) );
      }
    }
    return;
  }
  // consider all the nodes belongs to neighbor elements
  else
  {
    std::set<const Elem *> elems;
    std::map<unsigned int, std::set<const FVM_Node *> > region_node_map;
    std::map<unsigned int, std::set<const FVM_Node *> > region_off_node_map;

    // search for all the neighbor elements in this region
    fvm_element_iterator element_it = elem_begin();
    for( ; element_it != elem_end(); ++element_it)
    {
      const Elem * e = (*element_it).first;
      elems.insert(e);
      for(unsigned int n=0; n<e->n_sides(); ++n)
        if( e->neighbor(n) ) elems.insert(e->neighbor(n));
    }

    // consider ghost nodes in other regions
    if( _ghost_nodes!=NULL && !_ghost_nodes->empty() )
    {
      for(fvm_ghost_node_iterator  git = ghost_node_begin(); git != ghost_node_end(); ++git)
      {
        FVM_Node *ghost_node = (*git).first;
        if(!ghost_node) continue;

        for(element_it = ghost_node->elem_begin(); element_it != ghost_node->elem_end(); ++element_it)
        {
          const Elem * e = (*element_it).first;
          elems.insert(e);
          for(unsigned int n=0; n<e->n_sides(); ++n)
            if( e->neighbor(n) ) elems.insert(e->neighbor(n));
        }
      }
    }

    std::set<const Elem *>::const_iterator elem_it = elems.begin();
    for( ; elem_it != elems.end(); ++elem_it)
    {
      const Elem * e = *elem_it;
      for(unsigned int v=0; v<e->n_vertices(); v++)
      {
        const FVM_Node * fvm_node = e->get_fvm_node(v) ;
        region_node_map[fvm_node->subdomain_id()].insert(fvm_node);

        if( fvm_node->on_local() && !fvm_node->on_processor() )
          region_off_node_map[fvm_node->subdomain_id()].insert(fvm_node);
      }
    }

    // statistic all the PDE pattern infomation
    std::map<unsigned int, std::set<const FVM_Node *> >::iterator it;
    for(it= region_node_map.begin(); it!= region_node_map.end(); ++it)
      v_region_nodes.push_back(std::pair<unsigned int, unsigned int>((*it).first, (*it).second.size()));

    // statistic all the off-processor PDE pattern infomation
    for(it= region_off_node_map.begin(); it!= region_off_node_map.end(); ++it)
      v_off_region_nodes.push_back(std::pair<unsigned int, unsigned int>((*it).first, (*it).second.size()));

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
      if( (*it).first->on_local() && (*it).first->root_node()->processor_id() != Genius::processor_id() )   nodes++;
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
          if( (*gnit).first->on_local() && (*gnit).first->root_node()->processor_id() != Genius::processor_id() )   nodes++;
        v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(ghost_node->subdomain_id (), nodes) );
      }

  }
  // consider all the nodes belongs to neighbor elements
  else
  {
    std::set<const Elem *> elems;
    std::map<unsigned int, std::set<const FVM_Node *> > region_node_map;

    // search for all the neighbor elements in this region
    fvm_element_iterator element_it = elem_begin();
    for( ; element_it != elem_end(); ++element_it)
    {
      const Elem * e = (*element_it).first;
      elems.insert(e);
      for(unsigned int n=0; n<e->n_sides(); ++n)
        if( e->neighbor(n) ) elems.insert(e->neighbor(n));
    }

    // consider ghost nodes in other regions
    if( _ghost_nodes!=NULL && !_ghost_nodes->empty() )
    {
      for(fvm_ghost_node_iterator  git = ghost_node_begin(); git != ghost_node_end(); ++git)
      {
        FVM_Node *ghost_node = (*git).first;
        if(!ghost_node) continue;

        for(element_it = ghost_node->elem_begin(); element_it != ghost_node->elem_end(); ++element_it)
        {
          const Elem * e = (*element_it).first;
          elems.insert(e);
          for(unsigned int n=0; n<e->n_sides(); ++n)
            if( e->neighbor(n) ) elems.insert(e->neighbor(n));
        }
      }
    }

    std::set<const Elem *>::const_iterator elem_it = elems.begin();
    for( ; elem_it != elems.end(); ++elem_it)
    {
      const Elem * e = *elem_it;
      for(unsigned int v=0; v<e->n_vertices(); v++)
      {
        const FVM_Node * fvm_node = e->get_fvm_node(v) ;
        if( fvm_node->on_local() && !fvm_node->on_processor() )
          region_node_map[fvm_node->subdomain_id()].insert(fvm_node);
      }
    }

    // statistic all the off-processor PDE pattern infomation
    std::map<unsigned int, std::set<const FVM_Node *> >::iterator it= region_node_map.begin();
    for(; it!= region_node_map.end(); ++it)
      v_region_nodes.push_back(std::pair<unsigned int, unsigned int>((*it).first, (*it).second.size()));
    return;
  }

}



PetscScalar FVM_Node::variable(SolutionVariable var) const
{
  genius_assert( _node_data );
  return _node_data->get_variable_real(var);
}


VectorValue<PetscScalar> FVM_Node::gradient(SolutionVariable var, bool ghost) const
{
  PetscScalar a=0,b=0,c=0,d=0,e=0,f=0;
  PetscScalar r1=0,r2=0,r3=0;

  fvm_neighbor_node_iterator neighbor_it = this->neighbor_node_begin();
  fvm_neighbor_node_iterator neighbor_it_end =  this->neighbor_node_end();
  for(;neighbor_it!=neighbor_it_end; ++neighbor_it)
  {
    const FVM_Node * neighbor_node = neighbor_it->first;
    genius_assert(neighbor_node->on_local());

    PetscScalar w = 1.0/ this->distance(neighbor_node);
    PetscScalar dx = neighbor_node->root_node()->x() - this->root_node()->x();
    PetscScalar dy = neighbor_node->root_node()->y() - this->root_node()->y();
    PetscScalar dz = neighbor_node->root_node()->z() - this->root_node()->z();
    PetscScalar dphi = neighbor_node->variable(var) - this->variable(var);
    a += w*dx*dx;
    b += w*dx*dy;
    c += w*dx*dz;
    d += w*dy*dy;
    e += w*dy*dz;
    f += w*dz*dz;

    r1 += w*dx*dphi;
    r2 += w*dy*dphi;
    r3 += w*dz*dphi;
  }

  if(ghost)
  {
    std::map< FVM_Node *, std::pair<unsigned int, Real> >::const_iterator g_it = _ghost_nodes->begin();
    for( ; g_it != _ghost_nodes->end(); ++ g_it)
    {
      const FVM_Node *ghost_node = g_it->first;
      if(!ghost_node) continue;

      fvm_neighbor_node_iterator neighbor_it = ghost_node->neighbor_node_begin();
      fvm_neighbor_node_iterator neighbor_it_end =  ghost_node->neighbor_node_end();
      for(;neighbor_it!=neighbor_it_end; ++neighbor_it)
      {
        const FVM_Node * neighbor_node = neighbor_it->first;
        genius_assert(neighbor_node->on_local());

        PetscScalar w = 1.0/ ghost_node->distance(neighbor_node);
        PetscScalar dx = neighbor_node->root_node()->x() - ghost_node->root_node()->x();
        PetscScalar dy = neighbor_node->root_node()->y() - ghost_node->root_node()->y();
        PetscScalar dz = neighbor_node->root_node()->z() - ghost_node->root_node()->z();
        PetscScalar dphi = neighbor_node->variable(var) - ghost_node->variable(var);
        a += w*dx*dx;
        b += w*dx*dy;
        c += w*dx*dz;
        d += w*dy*dy;
        e += w*dy*dz;
        f += w*dz*dz;

        r1 += w*dx*dphi;
        r2 += w*dy*dphi;
        r3 += w*dz*dphi;
      }
    }
  }

  // for 2D case, f is 0, which breaks LU solver
  // so we set it to 1.0
  if(fabs(f)<1e-6) f=1.0;

  TNT::Array2D<Real> A (3, 3);
  TNT::Array1D<Real> r (3);
  A[0][0] = a; A[0][1] = b; A[0][2] = c;
  A[1][0] = b; A[1][1] = d; A[1][2] = e;
  A[2][0] = c; A[2][1] = e; A[2][2] = f;

  r[0] = r1;
  r[1] = r2;
  r[2] = r3;

  JAMA::LU<Real> solver(A);
  TNT::Array1D<Real> dphi = solver.solve(r);

  return VectorValue<PetscScalar>(dphi[0], dphi[1], dphi[2]);
}


void FVM_Node::operator += (const FVM_Node &other_node)
{
  // check if has the same root node
  genius_assert( this->_node == other_node.root_node() );

  // skip nonlocal fvm_node
  if( !this->on_local() ) return;

  // combine element
  _elem_has_this_node.insert(_elem_has_this_node.end(),  other_node.elem_begin(),  other_node.elem_end());

  // combine neighbor node
  for(unsigned int n=0; n<other_node._fvm_node_neighbor.size(); ++n)
    this->add_fvm_node_neighbor(other_node._fvm_node_neighbor[n].first, other_node._fvm_node_neighbor[n].second.first);

  // add volume
  _volume += other_node._volume;

}



void FVM_Node::prepare_for_use()
{
  //std::vector< std::pair<const Elem *, unsigned int> >(_elem_has_this_node).swap(_elem_has_this_node);
  //std::vector< std::pair<FVM_Node *, Real> >(_fvm_node_neighbor).swap(_fvm_node_neighbor);
  FNLess less;
  std::sort( _fvm_node_neighbor.begin(), _fvm_node_neighbor.end(), less );
}



size_t FVM_Node::memory_size() const
{
  size_t counter = sizeof(*this);

  counter += sizeof(FVM_NodeData);
  counter += _elem_has_this_node.capacity()*sizeof(std::pair<const Elem *, unsigned int>);
  counter += _fvm_node_neighbor.capacity()*sizeof(std::pair<FVM_Node *, std::pair<Real, Real> >);

  if(_ghost_nodes)
  {
    std::map< FVM_Node *, std::pair<unsigned int, Real> >::const_iterator it = _ghost_nodes->begin();
    for( ; it != _ghost_nodes->end(); ++it )
      counter += sizeof(it->first) + sizeof(it->second);
  }

  return counter;
}


void FVM_Node::print(std::ostream& os) const
{
  os << "FVM_Node id " << this->root_node()->id() << ", region " << _subdomain_id << ", volume " << this->volume() << std::endl;
  for(unsigned int n=0; n<_fvm_node_neighbor.size(); ++n)
  {
    os << "  Neighbor " << _fvm_node_neighbor[n].first->root_node()->id() << ", surface area " << _fvm_node_neighbor[n].second.first << std::endl;
  }
  os << "Total surface area " << this->total_cv_surface_area() << std::endl;
}

