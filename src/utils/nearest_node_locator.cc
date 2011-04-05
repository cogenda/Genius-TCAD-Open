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

#include <set>
#include "mesh_base.h"
#include "nearest_node_locator.h"


inline Real get_location( const Node *n, unsigned int k ) { return n->coord(k); }



NearestNodeLocator::NearestNodeLocator(const MeshBase& mesh)
  :_mesh(mesh)
{
  std::vector< std::set<const Node *> > subdomain_nodes(mesh.n_subdomains());

  MeshBase::const_element_iterator       el  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();
  for ( ; el != end; ++el)
  {
    const Elem * elem = *el;
    unsigned int subdomain = elem->subdomain_id();
    std::set<const Node *> & node_set = subdomain_nodes[subdomain];
    for (unsigned int n=0; n<elem->n_nodes(); ++n)
      node_set.insert(elem->get_node(n));
  }

  for(unsigned int n=0; n<subdomain_nodes.size(); ++n)
  {
    const std::set<const Node *> & node_set = subdomain_nodes[n];
    kdtree_type * kd_tree = new kdtree_type(std::ptr_fun(get_location));
    std::set<const Node *>::const_iterator node_it =node_set.begin();
    for( ; node_it!=node_set.end(); ++node_it)
      kd_tree->insert(*node_it);

    _kdtrees.push_back(kd_tree);
  }
}


NearestNodeLocator::~NearestNodeLocator()
{
  for(unsigned int n=0; n<_kdtrees.size(); ++n)
    delete _kdtrees[n];
}


Real NearestNodeLocator::distance_to_nearest_node(const Point &p, unsigned int subdomain) const
{
  const kdtree_type * kd_tree = _kdtrees[subdomain];
  Node source_node(p);

  std::pair<kdtree_type::const_iterator,  kdtree_type::distance_type> pItr = kd_tree->find_nearest(&source_node, 1e30);
  assert(pItr.first != kd_tree->end() );

  return pItr.second;
}


std::vector<const Node * > NearestNodeLocator::nearest_nodes(const Point &p, Real radius, unsigned int subdomain) const
{
  const kdtree_type * kd_tree = _kdtrees[subdomain];
  Node source_node(p);

  std::vector< const Node * > nn;
  kd_tree->find_within_range(&source_node, radius, std::back_inserter(nn));

  return nn;
}


std::vector<const Node * > NearestNodeLocator::nearest_nodes(const Point &p1, const Point &p2, Real radius, unsigned int subdomain) const
{
  std::vector<const Node *> nn;
  std::set<const Node *> nn_set;
  this->nearest_nodes(p1, p2, radius, subdomain, nn_set);
  nn.insert(nn.end(), nn_set.begin(), nn_set.end());
  return nn;
}


void NearestNodeLocator::nearest_nodes(const Point &p1, const Point &p2, Real radius, unsigned int subdomain, std::set<const Node * >& nns) const
{
  unsigned int nn_num = nns.size();
  std::vector< const Node * > nn1 = this->nearest_nodes(p1, radius, subdomain);
  std::vector< const Node * > nn2 = this->nearest_nodes(p2, radius, subdomain);
  nns.insert(nn1.begin(), nn1.end());
  nns.insert(nn2.begin(), nn2.end());

  if( nns.size() > nn_num )
  {
    this->nearest_nodes( p1, 0.5*(p1+p2), radius, subdomain, nns );
    this->nearest_nodes( 0.5*(p1+p2), p2, radius, subdomain, nns );
  }
}



