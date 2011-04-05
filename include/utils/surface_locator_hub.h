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


#ifndef __surface_locator_hub_h__
#define __surface_locator_hub_h__

// Local Includes
#include "genius_env.h"
#include "genius_common.h"
#include "enum_surface_locator_type.h"

class MeshBase;
class Elem;
class Point;
class SurfaceLocatorBase;

/**
 * This is the base class for surface element locators.
 * They locate surface element of a specified region in space
 * given a start point they return the element
 * which nearest to given point in specified region
 */

// ------------------------------------------------------------
// SurfaceLocatorHub class definition
class SurfaceLocatorHub
{
public:

  /**
   * Constructor.
   */
  SurfaceLocatorHub (const MeshBase& mesh, const SurfaceLocatorType t=SurfaceLocator_SPHERE);

  /**
   * Destructor.
   */
  virtual ~SurfaceLocatorHub ()
  { this->clear(); }


  /**
   * Clears the \p SurfaceLocator.
   */
  void clear();


  /**
   * Locates the surface element with specified subdomain which is nearest to given point p
   */
  std::pair<const Elem*, unsigned int> operator() (const Point& p, const unsigned int subdomain, Point & project_point, const Real dist=1e30);

  /**
   * Locates the surface element with specified boundary which is nearest to given point p
   */
  std::pair<const Elem*, unsigned int> operator() (const Point& p, const short int boundary, Point & project_point, const Real dist=1e30);

private:

  /**
   * constant reference to the mesh
   */
  const MeshBase& _mesh;

  SurfaceLocatorType _type;

  std::vector<const SurfaceLocatorBase *> _subdomain_surface_locators;

  std::map<short int, const SurfaceLocatorBase *> _boundary_surface_locators;

};


#endif
