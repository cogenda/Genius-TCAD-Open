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


#ifndef __surface_locator_sphere_h__
#define __surface_locator_sphere_h__

#include "surface_locator_base.h"

#include <functional>
#include <kdtree.hpp>

// Forward Declarations
class MeshBase;
class Point;
class Elem;



namespace SurfaceLocator
{
  /**
   * a voxel is a sphere in space, it also contains element which in/intersect this sphere
   */
  struct Voxel
  {
    Point center;
    Real  r;
    std::vector<const Elem *> voxel_elements;
  };

  inline Real get_location( const Voxel *v, unsigned int k ) { return v->center[k]; }

}




/**
 * This is the base class for surface element locators.
 * They locate surface element of a specified region in space
 * given a start point they return the element
 * which nearest to given point in specified region
 */

// ------------------------------------------------------------
// SurfaceLocatorSphere class definition
class SurfaceLocatorSphere : public SurfaceLocatorBase
{
public:

  /**
   * Constructor.
   */
  SurfaceLocatorSphere (const MeshBase& mesh, const unsigned int subdomain);

  /**
   * Constructor.
   */
  SurfaceLocatorSphere (const MeshBase& mesh, const short int boundary);

  /**
   * Destructor.
   */
  ~SurfaceLocatorSphere ();

  /**
   * Clears the \p SurfaceLocator.
   */
  void clear();

  /**
   * Initializes the surface locator, so that the \p operator() methods can
   * be used.  Pure virtual.
   */
  void init();

  /**
   * Locates the element with specified subdomain which is nearest to given point p
   */
  std::pair<const Elem*, unsigned int> operator() (const Point& p, Point & project_point, const Real dist=1e30) const ;


private:

  typedef std::pair<Point, Point> BoundingBox;

  /**
   * the bounding box of mesh
   */
  BoundingBox _bounding_box;

  /**
   * x size of Voxel
   */
  Real _x_size;

  /**
   * y size of Voxel
   */
  Real _y_size;

  /**
   * z size of Voxel
   */
  Real _z_size;

  /**
   * radius size of Voxel sphere
   */
  Real _radius_size;

  /**
   * dim size in x
   */
  unsigned int _dim_x;

  /**
   * dim size in y
   */
  unsigned int _dim_y;

  /**
   * dim size in z
   */
  unsigned int _dim_z;


  struct Trinity
  {
    unsigned int nx, ny, nz;
  };

  /**
   * arrays of voxel
   */
  std::map<unsigned int, SurfaceLocator::Voxel *>  _space_division;

  /**
   * get a voxel by index
   */
  SurfaceLocator::Voxel * _get_voxel(unsigned int nx, unsigned int ny, unsigned int nz) 
  {
    unsigned int key = nx*(_dim_y*_dim_z) + ny*_dim_z + nz;
    if( _space_division.find(key) == _space_division.end() )
    {
      SurfaceLocator::Voxel * voxel = new SurfaceLocator::Voxel;
      voxel->center = _bounding_box.first + Point(_x_size*nx, _y_size*ny, _z_size*nz);
      voxel->r = _radius_size;
      _space_division[key] = voxel;
      return voxel;
    }
    return _space_division[key];
  }


  void _get_elem_range(const Elem * elem, Trinity &low, Trinity &high) const;


  bool _is_elem_in_voxel( const SurfaceLocator::Voxel * voxel, const Elem * elem ) const;


  typedef KDTree::KDTree<3, const SurfaceLocator::Voxel *, std::pointer_to_binary_function<const SurfaceLocator::Voxel *, unsigned int, Real> > tree_type;
  tree_type *_exact_dist;

  /**
   * store surface elem here, delete them when exit
   */
  std::vector<const Elem *> _surface_element_list;


  std::map<const Elem *, std::pair<const Elem *, unsigned int> > _surface_to_volume_element_map;

};


#endif
