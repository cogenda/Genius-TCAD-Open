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


#ifndef __surface_locator_base_h__
#define __surface_locator_base_h__

// Local Includes
#include "genius_env.h"
#include "genius_common.h"
#include "auto_ptr.h"
#include "enum_surface_locator_type.h"

/**
 * This is the base class for surface element locators.
 * They locate surface element of a specified region in space
 * given a start point they return the element
 * which nearest to given point in specified region
 */

// ------------------------------------------------------------
// SurfaceLocatorBase class definition
class SurfaceLocatorBase
{
protected:

  /**
   * Constructor.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  SurfaceLocatorBase (const MeshBase& mesh, const unsigned int subdomain):
    _mesh(mesh),
    _subdomain(subdomain),
    _boundary(-1234),
    _initialized(false)
  {}


  /**
   * Constructor.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  SurfaceLocatorBase (const MeshBase& mesh, const short int boundary):
    _mesh(mesh),
    _subdomain(invalid_uint),
    _boundary(boundary),
    _initialized(false)
    {}

public:

  /**
   * Destructor.
   */
  virtual ~SurfaceLocatorBase () {}

  /**
   * Clears the \p SurfaceLocator.
   */
  virtual void clear() = 0;

  /**
   * Initializes the surface locator, so that the \p operator() methods can
   * be used.  Pure virtual.
   */
  virtual void init() = 0;


  static AutoPtr<SurfaceLocatorBase> build (const SurfaceLocatorType t, const MeshBase& mesh, unsigned int subdomain);

  static AutoPtr<SurfaceLocatorBase> build (const SurfaceLocatorType t, const MeshBase& mesh, short int boundary);


  /**
   * Locates the surface element with specified subdomain which is nearest to given point p
   */
  virtual std::pair<const Elem*, unsigned int> operator() (const Point& p, Point & project_point, Real dist=1e30) const = 0;

  /**
   * @returns \p true when this object is properly initialized
   * and ready for use, \p false otherwise.
   */
  bool initialized () const { return _initialized; }


protected:

  bool subdomain_locator() const { return _subdomain!=invalid_uint; }

  bool boundary_locator() const { return _boundary!=-1234; }

  /**
   * constant reference to the mesh
   */
  const MeshBase& _mesh;

  /**
   * which region surface to be find
   */
  const unsigned int _subdomain;

  /**
   * which boundary surface to be find
   */
  const short int _boundary;

  /**
     * \p true when properly initialized, \p false otherwise.
   */
  bool _initialized;

};


#endif
