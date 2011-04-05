// $Id: surface_locator_base.cc,v 1.1 2008/05/22 14:13:34 gdiso Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// C++ includes


// Local Includes
#include "mesh_base.h"
#include "boundary_info.h"

#include "surface_locator_hub.h"
#include "surface_locator_base.h"




//------------------------------------------------------------------
// PointLocatorBase methods


SurfaceLocatorHub::SurfaceLocatorHub (const MeshBase& mesh, const SurfaceLocatorType t)
  :_mesh(mesh), _type(t)
{
  _subdomain_surface_locators.resize(mesh.n_subdomains(), NULL);
}


void SurfaceLocatorHub::clear()
{
  for(unsigned int n=0; n<_subdomain_surface_locators.size(); ++n)
    delete _subdomain_surface_locators[n];

  std::map<short int, const SurfaceLocatorBase *>::iterator it = _boundary_surface_locators.begin();
  for(; it != _boundary_surface_locators.end(); ++it)
    delete it->second;
}


std::pair<const Elem*, unsigned int> SurfaceLocatorHub::operator() (const Point& p, const unsigned int subdomain, Point & project_point, const Real dist)
{
  const SurfaceLocatorBase * & subdomain_surface_locator = _subdomain_surface_locators[subdomain];
  if( subdomain_surface_locator == NULL )
    subdomain_surface_locator = SurfaceLocatorBase::build(_type, _mesh, subdomain).release();
  return (*subdomain_surface_locator)(p, project_point, dist);
}


std::pair<const Elem*, unsigned int> SurfaceLocatorHub::operator() (const Point& p, const short int boundary, Point & project_point, const Real dist)
{
  if( _boundary_surface_locators.find(boundary) == _boundary_surface_locators.end() )
  {
    const SurfaceLocatorBase * _locator = SurfaceLocatorBase::build(_type, _mesh, boundary).release();
    _boundary_surface_locators.insert( std::make_pair(boundary,_locator ) );
  }

  return (*_boundary_surface_locators.find(boundary)->second)(p, project_point, dist);
}
