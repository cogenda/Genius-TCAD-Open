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

// TNT matrix-vector library
#include <TNT/tnt.h>
#include <TNT/jama_lu.h>


unsigned int FVM_Node::_solver_index=0;


FVM_Node::FVM_Node(const Node *n)
    : _node(n),
    _node_data(0),
    _ghost_nodes(0),
    _volume(0),
    _boundary_id(BoundaryInfo::invalid_id)
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
    // the PDE node pattern in thie region. n_node_neighbors()+1 means "this" node is also considered.
    v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(subdomain_id (), n_node_neighbors()+1) );

    //neighbor not on this processor
    unsigned int off_nodes = 0;
    for(fvm_neighbor_node_iterator it= neighbor_node_begin(); it!= neighbor_node_end(); ++it)
      if( (*it).first->on_local() && (*it).first->processor_id() != Genius::processor_id() )   off_nodes++;
    v_off_region_nodes.push_back( std::pair<unsigned int, unsigned int>(subdomain_id (), off_nodes) );

    // consider ghost nodes in other regions
    if( _ghost_nodes!=NULL && !_ghost_nodes->empty() )
    {
      for(fvm_ghost_node_iterator  git = ghost_node_begin(); git!=ghost_node_end(); ++git)
      {
        FVM_Node *ghost_node = (*git).first;
        if( ghost_node == NULL) continue;

        v_region_nodes.push_back( std::pair<unsigned int, unsigned int>(ghost_node->subdomain_id (), ghost_node->n_node_neighbors()+1) );

        unsigned int off_nodes = 0;
        fvm_neighbor_node_iterator gnit = ghost_node->neighbor_node_begin();
        for(; gnit!= ghost_node->neighbor_node_end(); ++gnit)
          if( (*gnit).first->on_local() && (*gnit).first->processor_id() != Genius::processor_id() )   off_nodes++;
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
    const FVM_Node * neighbor_node = neighbor_it->second;
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
        const FVM_Node * neighbor_node = neighbor_it->second;
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
  _node_neighbor.insert(other_node.neighbor_node_begin(),  other_node.neighbor_node_end());

  // add volume
  _volume += other_node._volume;

  // add cv surface
  std::map< const Node *, Real >::const_iterator cv_it = other_node._cv_surface_area.begin();
  for(; cv_it != other_node._cv_surface_area.end(); ++cv_it)
  {
    if(_cv_surface_area.find(cv_it->first) != _cv_surface_area.end())
      (*_cv_surface_area.find(cv_it->first)).second += cv_it->second;
    else
      // insert a new item
      _cv_surface_area[(*cv_it).first] = cv_it->second;
  }

}
