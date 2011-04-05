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
#include "surface_locator_base.h"
#include "surface_locator_list.h"
#include "surface_locator_sphere.h"



//------------------------------------------------------------------
// SurfaceLocatorBase methods


AutoPtr<SurfaceLocatorBase> SurfaceLocatorBase::build (const SurfaceLocatorType t, const MeshBase& mesh, unsigned int subdomain)
{
  switch (t)
  {
      case SurfaceLocator_SPHERE:
      {
        AutoPtr<SurfaceLocatorBase> ap(new SurfaceLocatorSphere(mesh, subdomain));
        return ap;
      }

      case SurfaceLocator_LIST:
      {
        AutoPtr<SurfaceLocatorBase> ap(new SurfaceLocatorList(mesh, subdomain));
        return ap;
      }

      default:
      {
        std::cerr << "ERROR: Bad SurfaceLocatorType = " << t << std::endl;
        genius_error();
      }
  }

  genius_error();
  AutoPtr<SurfaceLocatorBase> ap(NULL);
  return ap;
}


AutoPtr<SurfaceLocatorBase> SurfaceLocatorBase::build (const SurfaceLocatorType t, const MeshBase& mesh, short int boundary)
{
  switch (t)
  {
    case SurfaceLocator_SPHERE:
    {
      AutoPtr<SurfaceLocatorBase> ap(new SurfaceLocatorSphere(mesh, boundary));
      return ap;
    }

    case SurfaceLocator_LIST:
    {
      AutoPtr<SurfaceLocatorBase> ap(new SurfaceLocatorList(mesh, boundary));
      return ap;
    }

    default:
    {
      std::cerr << "ERROR: Bad SurfaceLocatorType = " << t << std::endl;
      genius_error();
    }
  }

  genius_error();
  AutoPtr<SurfaceLocatorBase> ap(NULL);
  return ap;
}

