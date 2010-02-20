// $Id: petsc_macro.h,v 1.2 2007-10-21 20:48:43 benkirk Exp $

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

#ifndef __petsc_macro_h__
#define __petsc_macro_h__

// Local includes
#include "config.h"

#ifdef HAVE_PETSC


#include "petsc.h"
#include "petscversion.h"


// A set of convenient macro for comparing PETSc versions.

#define PETSC_VERSION_LT(major,minor,subminor)			                    \
  ((PETSC_VERSION_MAJOR < (major) ||						    \
    (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR < (minor) ||	    \
				  (PETSC_VERSION_MINOR == (minor) &&		    \
				   PETSC_VERSION_SUBMINOR < (subminor))))) ? 1 : 0)

#define PETSC_VERSION_LE(major,minor,subminor)			                    \
  ((PETSC_VERSION_MAJOR < (major) ||						    \
    (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR < (minor) ||	    \
				  (PETSC_VERSION_MINOR == (minor) &&		    \
				   PETSC_VERSION_SUBMINOR <= (subminor))))) ? 1 : 0)

#define PETSC_VERSION_EQ(major,minor,subminor)  \
   ((PETSC_VERSION_MAJOR == (major) && \
     PETSC_VERSION_MINOR == (minor)  && \
     PETSC_VERSION_SUBMINOR == (subminor)) ? 1 : 0 )

#define PETSC_VERSION_GT(major,minor,subminor)			                    \
  ((PETSC_VERSION_MAJOR > (major) ||						    \
    (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR > (minor) ||	    \
				  (PETSC_VERSION_MINOR == (minor) &&		    \
				   PETSC_VERSION_SUBMINOR > (subminor))))) ? 1 : 0)

#define PETSC_VERSION_GE(major,minor,subminor)			                    \
  ((PETSC_VERSION_MAJOR == (major) ||						    \
    (PETSC_VERSION_MAJOR == (major) && (PETSC_VERSION_MINOR == (minor) ||	    \
				  (PETSC_VERSION_MINOR == (minor) &&		    \
				   PETSC_VERSION_SUBMINOR >= (subminor))))) ? 1 : 0)

#endif // HAVE_PETSC

#endif // __petsc_macro_h__
