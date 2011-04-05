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


#ifndef __surface_locator_list_h__
#define __surface_locator_list_h__

#include "surface_locator_base.h"



// Forward Declarations
class MeshBase;
class Point;
class Elem;

/**
 * class for surface element locators
 * just iterate all the possible surface elements.
 * They locate surface element of a specified region in space
 * given a start point they return the element
 * which nearest to given point in specified region
 */

// ------------------------------------------------------------
// SurfaceLocatorList class definition
class SurfaceLocatorList : public SurfaceLocatorBase
{
public:

  /**
   * Constructor.
   */
  SurfaceLocatorList (const MeshBase& mesh, const unsigned int subdomain);


  /**
   * Constructor.
   */
  SurfaceLocatorList (const MeshBase& mesh, const short int boundary);

  /**
   * Destructor.
   */
  ~SurfaceLocatorList ();

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

  std::vector<const Elem *> _surface_element_list;

  std::map<const Elem *, std::pair<const Elem *, unsigned int> > _surface_to_volume_element_map;
};


#endif
