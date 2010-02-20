// $Id: enum_elem_type.h,v 1.4 2008/04/17 13:25:47 gdiso Exp $

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



#ifndef __enum_elem_type_h__
#define __enum_elem_type_h__

/**
 * The \p libMeshEnums namespace is the namespace all \p enum definitions
 * should be put into.
 */

// ------------------------------------------------------------
// enum ElemType definition
namespace libMeshEnums {

  /**
   * Defines an \p enum for geometric element types.
   */
  enum ElemType {
         EDGE2=0,    // 0
         EDGE2_FVM,  // 1
         EDGE3,      // 2
         EDGE4,      // 3

         TRI3,       // 4
         TRI3_FVM,   // 5
         TRI6,       // 6

         QUAD4,      // 7
         QUAD4_FVM,  // 8
         QUAD8,      // 9
         QUAD9,      // 10

         TET4,       // 11
         TET4_FVM,   // 12
         TET10,      // 13

         HEX8,       // 14
         HEX8_FVM,   // 15
         HEX20,      // 16
         HEX27,      // 17

         PRISM6,     // 18
         PRISM6_FVM, // 19
         PRISM15,    // 20
         PRISM18,    // 21

         PYRAMID5,    // 22
         PYRAMID5_FVM,// 23

         INFEDGE2,   // 24

         INFQUAD4,   // 25
         INFQUAD6,   // 26

         INFHEX8,    // 27
         INFHEX16,   // 28
         INFHEX18,   // 29

         INFPRISM6,  // 30
         INFPRISM12, // 31

         NODEELEM,   // 32

         REMOTEELEM,   // 33

         INVALID_ELEM};  // 34 - should always be last
}

using namespace libMeshEnums;

#endif // #define __enum_elem_type_h__




