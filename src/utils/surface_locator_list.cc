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


// C++ includes

// Local Includes
#include "mesh_base.h"
#include "mesh_tools.h"
#include "boundary_info.h"
#include "surface_locator_list.h"




//------------------------------------------------------------------
// SurfaceLocator methods
SurfaceLocatorList::SurfaceLocatorList (const MeshBase& mesh, const unsigned int subdomain) :
    SurfaceLocatorBase (mesh, subdomain)
{
  this->init();
}


SurfaceLocatorList::SurfaceLocatorList (const MeshBase& mesh, const short int boundary) :
    SurfaceLocatorBase (mesh, boundary)
{
  this->init();
}


SurfaceLocatorList::~SurfaceLocatorList ()
{
  this->clear ();
}




void SurfaceLocatorList::clear ()
{
  for(unsigned int n=0; n<_surface_element_list.size(); ++n)
    delete _surface_element_list[n];
  _surface_element_list.clear();
}



void SurfaceLocatorList::init()
{
  std::vector<unsigned int>       elems;
  std::vector<unsigned short int> sides;
  std::vector<short int>          bds;
  _mesh.boundary_info->build_active_side_list (elems, sides, bds);

  for(unsigned int n=0; n<elems.size(); ++n)
  {
    const Elem * elem = _mesh.elem(elems[n]);
    if(subdomain_locator() && elem->subdomain_id() != _subdomain) continue;
    if(boundary_locator() && bds[n] != _boundary) continue;

    const Elem * surface_elem = elem->build_side(sides[n], false).release();
    _surface_element_list.push_back(surface_elem);
    _surface_to_volume_element_map[surface_elem] = std::make_pair(elem, sides[n]);
  }

  genius_assert( !_surface_element_list.empty() );

  // ready for take-off
  this->_initialized = true;
}



std::pair<const Elem*, unsigned int> SurfaceLocatorList::operator() (const Point& p, Point & project_point, const Real dist) const
{
  const Elem * nearest_elem = 0;
  Point nearest_point;

  Real  voting_dist = 1e30;
  for(unsigned int n=0; n<_surface_element_list.size(); ++n)
  {
    const Elem * surface_elem = _surface_element_list[n];

    Real d;
    Point np = surface_elem->nearest_point(p, &d);
    if( d < voting_dist )
    {
      nearest_elem = surface_elem;
      nearest_point = np;
      voting_dist = d;
    }
  }
  if( !nearest_elem ) return std::make_pair((const Elem*)0, invalid_uint);

  project_point = nearest_point;
  return _surface_to_volume_element_map.find(nearest_elem)->second;;
}



