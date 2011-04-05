// $Id: enum_surface_locator_type.h,v 1.3 2008/05/22 14:14:47 gdiso Exp $

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



#ifndef __enum_surface_locator_type_h__
#define __enum_surface_locator_type_h__

// ------------------------------------------------------------
// enum SurfaceLocatorType definition
namespace MeshEnums {

  /**
 * \enum MeshEnums::SurfaceLocatorType defines an \p enum for the types
   * of surface locators (given a point with global coordinates,
   * locate the corresponding nearest surface element in space) available in libMesh.
   */
  enum SurfaceLocatorType {
    SurfaceLocator_SPHERE = 0,
    SurfaceLocator_LIST,
    INVALID_SurfaceLocator};
}

using namespace MeshEnums;

#endif // #ifndef __enum_surface_locator_type_h__




