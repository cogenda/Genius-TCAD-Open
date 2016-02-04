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
#include "elem.h"
#include "surface_locator_sphere.h"

#define DEBUG

using namespace SurfaceLocator;




//------------------------------------------------------------------
// SurfaceLocator methods
SurfaceLocatorSphere::SurfaceLocatorSphere (const MeshBase& mesh, const unsigned int subdomain) :
    SurfaceLocatorBase (mesh, subdomain)
{
  this->init();
}


SurfaceLocatorSphere::SurfaceLocatorSphere (const MeshBase& mesh, const short int boundary) :
    SurfaceLocatorBase (mesh, boundary)
{
  this->init();
}


SurfaceLocatorSphere::~SurfaceLocatorSphere ()
{
  this->clear ();
}




void SurfaceLocatorSphere::clear ()
{
  std::map<unsigned int, SurfaceLocator::Voxel *>::const_iterator vol_it =  _space_division.begin();
  for( ; vol_it !=  _space_division.end(); vol_it++ )
    delete vol_it->second;

  delete _exact_dist;

  for(unsigned int n=0; n<_surface_element_list.size(); ++n)
    delete _surface_element_list[n];
  _surface_element_list.clear();
}





void SurfaceLocatorSphere::init ()
{
  std::vector<unsigned int>       elems;
  std::vector<unsigned short int> sides;
  std::vector<short int>          bds;
  _mesh.boundary_info->build_active_side_list (elems, sides, bds);


  // build the bounding box of element
  // find the minimal size of element, set as an reference of radius of bounding sphere
  Point min(1.e30,   1.e30,  1.e30);
  Point max(-1.e30, -1.e30, -1.e30);
  for(unsigned int e=0; e<elems.size(); ++e)
  {
    const Elem * elem = _mesh.elem(elems[e]);
    if( subdomain_locator() && elem->subdomain_id() != _subdomain ) continue;
    if( boundary_locator() && bds[e] != _boundary ) continue;

    for (unsigned int n=0; n<elem->n_nodes(); n++)
      for (unsigned int i=0; i<3; i++)
      {
        min(i) = std::min(min(i), elem->point(n)(i));
        max(i) = std::max(max(i), elem->point(n)(i));
      }
  }
  _bounding_box = std::make_pair(min, max);

  // set radius_size, we must limit the size to about 100*100*100
  _x_size = (_bounding_box.second.x() - _bounding_box.first.x())/100.0+1e-30;
  _y_size = (_bounding_box.second.y() - _bounding_box.first.y())/100.0+1e-30;
  _z_size = (_bounding_box.second.z() - _bounding_box.first.z())/100.0+1e-30;
  _radius_size =  0.5*sqrt(_x_size*_x_size + _y_size*_y_size + _z_size*_z_size)*(1.0+1e-10);


  _dim_x = static_cast<unsigned int>((_bounding_box.second.x() - _bounding_box.first.x())/_x_size) + 1;
  _dim_y = static_cast<unsigned int>((_bounding_box.second.y() - _bounding_box.first.y())/_y_size) + 1;
  _dim_z = static_cast<unsigned int>((_bounding_box.second.z() - _bounding_box.first.z())/_z_size) + 1;
  unsigned int voxel_size = _dim_x*_dim_y*_dim_z;


  for(unsigned int n=0; n<elems.size(); ++n)
  {
    const Elem * elem = _mesh.elem(elems[n]);
    if( subdomain_locator() && elem->subdomain_id() != _subdomain ) continue;
    if( boundary_locator() && bds[n] != _boundary ) continue;

    const Elem * surface_elem = elem->build_side(sides[n], false).release();
    _surface_element_list.push_back(surface_elem);
    _surface_to_volume_element_map[surface_elem] = std::make_pair(elem, sides[n]);

    Trinity low, high;
    this->_get_elem_range(surface_elem, low, high);
    for(unsigned int nx=low.nx; nx<=high.nx; ++nx)
      for(unsigned int ny=low.ny; ny<=high.ny; ++ny)
        for(unsigned int nz=low.nz; nz<=high.nz; ++nz)
        {
          Voxel * voxel = _get_voxel(nx, ny, nz);
          if( _is_elem_in_voxel(voxel, surface_elem) )
            voxel->voxel_elements.push_back(surface_elem);
        }
  }


  _exact_dist = new tree_type(std::ptr_fun(get_location));
  //for all the non-empty voxel
  std::map<unsigned int, SurfaceLocator::Voxel *>::const_iterator vol_it =  _space_division.begin();
  for( ; vol_it !=  _space_division.end(); vol_it++ )
  {
    const Voxel * voxel = vol_it->second;
    if(voxel->voxel_elements.empty()) continue;
    _exact_dist->insert(voxel);
  }
  
  // ready for take-off
  this->_initialized = true;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}


void SurfaceLocatorSphere::_get_elem_range(const Elem * elem, Trinity &low, Trinity &high) const
{
  Point min(1.e30,   1.e30,  1.e30);
  Point max(-1.e30, -1.e30, -1.e30);

  for (unsigned int n=0; n<elem->n_nodes(); n++)
    for (unsigned int i=0; i<3; i++)
    {
      min(i) = std::min(min(i), elem->point(n)(i));
      max(i) = std::max(max(i), elem->point(n)(i));
    }

  low.nx =  static_cast<unsigned int>( (min.x() - _bounding_box.first.x()) / _x_size );
  low.ny =  static_cast<unsigned int>( (min.y() - _bounding_box.first.y()) / _y_size );
  low.nz =  static_cast<unsigned int>( (min.z() - _bounding_box.first.z()) / _z_size );

  high.nx =  static_cast<unsigned int>( (max.x() - _bounding_box.first.x()) / _x_size ) + 1;
  high.ny =  static_cast<unsigned int>( (max.y() - _bounding_box.first.y()) / _y_size ) + 1;
  high.nz =  static_cast<unsigned int>( (max.z() - _bounding_box.first.z()) / _z_size ) + 1;

}



bool SurfaceLocatorSphere::_is_elem_in_voxel(const Voxel * voxel, const Elem * elem ) const
{
  for (unsigned int n=0; n<elem->n_nodes(); n++)
  {
    if( (voxel->center - elem->point(n)).size() < voxel->r  ) return true;
  }

  return false;
}


#if 0
std::pair<const Elem*, unsigned int> SurfaceLocatorSphere::operator() (const Point& p, Point & project_point, const Real dist) const
{
  // find the nearest voxel which is not empty by kdtree
  Voxel temp;
  temp.center = p;
  std::pair<tree_type::const_iterator,  tree_type::distance_type> pItr = _exact_dist->find_nearest((&temp), dist + 2*_radius_size);

  // no find within given distance
  if(pItr.first == _exact_dist->end() )
    return std::make_pair((const Elem*)0, invalid_uint);


  // extend by 2*_radius_size
  Real d = pItr.second + 2*_radius_size;

  // and find all the voxel in this dist
  std::vector< const Voxel * > voting_voxel;
  _exact_dist->find_within_range((&temp), d, std::back_inserter(voting_voxel));
  genius_assert(!voting_voxel.empty());

  std::set<const Elem *> voting_elem;
  for(unsigned int n=0 ; n<voting_voxel.size() ; ++n)
  {
    const Voxel * voxel = voting_voxel[n];
    voting_elem.insert(voxel->voxel_elements.begin(), voxel->voxel_elements.end());
  }
  genius_assert(!voting_elem.empty());


  // search at these element
  const Elem * nearest_elem = 0;
  Point nearest_point;
  Real  voting_dist = 1e30;

  std::set<const Elem *>::const_iterator elem_it = voting_elem.begin();
  for( ; elem_it != voting_elem.end(); ++elem_it)
  {
    const Elem * surface_elem = *elem_it;

    Real d;
    Point np = surface_elem->nearest_point(p, &d);
    if( d < voting_dist )
    {
      nearest_elem = surface_elem;
      nearest_point = np;
      voting_dist = d;
    }
  }

  project_point = nearest_point;
  return _surface_to_volume_element_map.find(nearest_elem)->second;;

}
#endif


std::pair<const Elem*, unsigned int> SurfaceLocatorSphere::operator() (const Point& p, Point & project_point, const Real dist) const
{
  // find the nearest voxel which is not empty by kdtree
  Voxel temp;
  temp.center = p;
  std::pair<tree_type::const_iterator,  tree_type::distance_type> pItr = _exact_dist->find_nearest((&temp), dist + 2*_radius_size);

  // no find within given distance
  if(pItr.first == _exact_dist->end() )
    return std::make_pair((const Elem*)0, invalid_uint);

  const std::vector<const Elem *> & voting_elem = (*pItr.first)->voxel_elements;

  // search at these element
  const Elem * nearest_elem = 0;
  Point nearest_point;
  Real  voting_dist = 1e30;

  std::vector<const Elem *>::const_iterator elem_it = voting_elem.begin();
  for( ; elem_it != voting_elem.end(); ++elem_it)
  {
    const Elem * surface_elem = *elem_it;

    Real d;
    Point np = surface_elem->nearest_point(p, &d);
    if( d < voting_dist )
    {
      nearest_elem = surface_elem;
      nearest_point = np;
      voting_dist = d;
    }
  }

  project_point = nearest_point;
  return _surface_to_volume_element_map.find(nearest_elem)->second;;

}
