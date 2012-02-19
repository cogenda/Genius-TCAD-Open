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
         EDGE2=0,     // 0
         EDGE2_FVM,   // 1
         EDGE3,       // 2
         EDGE4,       // 3

         TRI3,        // 4
         TRI3_FVM,    // 5
         TRI3_CY_FVM, // 6
         TRI6,        // 7

         QUAD4,       // 8
         QUAD4_FVM,   // 9
         QUAD4_CY_FVM,// 10
         QUAD8,       // 11
         QUAD9,       // 12

         TET4,        // 13
         TET4_FVM,    // 14
         TET10,       // 15

         HEX8,        // 16
         HEX8_FVM,    // 17
         HEX20,       // 18
         HEX27,       // 19

         PRISM6,      // 20
         PRISM6_FVM,  // 21
         PRISM15,     // 22
         PRISM18,     // 23

         PYRAMID5,    // 24
         PYRAMID5_FVM,// 25

         INFEDGE2,    // 26

         INFQUAD4,    // 27
         INFQUAD6,    // 28

         INFHEX8,     // 29
         INFHEX16,    // 30
         INFHEX18,    // 31

         INFPRISM6,   // 32
         INFPRISM12,  // 33

         NODEELEM,    // 34

         REMOTEELEM,  // 35

         INVALID_ELEM};  // 36 - should always be last
}

using namespace libMeshEnums;

#endif // #define __enum_elem_type_h__




