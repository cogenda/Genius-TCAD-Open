
// $Id: location_maps.h,v 1.1 2008/05/17 11:09:54 gdiso Exp $

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



#ifndef __location_maps_h__
#define __location_maps_h__

#include "genius_common.h"

// C++ Includes   -----------------------------------
#if defined(HAVE_TR1_UNORDERED_MAP)
# include <tr1/unordered_map>
#elif defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER) || defined(HAVE_UNORDERED_MAP)
# include <unordered_map>
#else
# include <map>
#endif



// Forward Declarations -----------------------------
class Elem;
class MeshBase;
class Node;


  /**
   * Data structures that enable location-based lookups
   * The key is a hash of the Point location.
   * For efficiency we will use a hashed multimap if it is
   * available, otherwise a regular multimap.
   */
template <typename T>
class LocationMap
{
#if defined(HAVE_UNORDERED_MAP)
  typedef std::unordered_multimap<unsigned int, T*> map_type;
#elif defined(HAVE_TR1_UNORDERED_MAP) || defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER)
  typedef std::tr1::unordered_multimap<unsigned int, T*> map_type;
#else
  typedef std::multimap<unsigned int, T*>  map_type;
#endif

public:
  void init(MeshBase&);

  void clear() { _map.clear(); }

  void insert(T&);

  bool empty() const { return _map.empty(); }

  T* find(const Point&,
	  const Real tol = TOLERANCE);

  Point point_of(const T&) const;

protected:
  unsigned int key(const Point&);

  void fill(MeshBase&);

private:
  map_type          _map;
  std::vector<Real> _lower_bound;
  std::vector<Real> _upper_bound;
};


#endif // #define __location_maps_h__
