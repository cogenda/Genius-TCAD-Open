// $Id: dof_object.cc,v 1.2 2008/05/22 04:53:44 gdiso Exp $

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

// Local includes
#include "dof_object.h"



// ------------------------------------------------------------
// DofObject class static member initialization
const unsigned int       DofObject::invalid_id           = invalid_uint;
const unsigned short int DofObject::invalid_processor_id = static_cast<unsigned short int>(-1);



// ------------------------------------------------------------
// DofObject class members
// Copy Constructor
DofObject::DofObject (const DofObject& dof_obj):
  _id            (dof_obj._id),
  _processor_id  (dof_obj._processor_id),
  _on_local      (dof_obj._on_local)
{

}





