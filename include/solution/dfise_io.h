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



#ifndef __df_ise_io_h__
#define __df_ise_io_h__

// C++ includes
#include <vector>
#include <fstream>

// Local includes
#include "config.h"
#include "auto_ptr.h"
#include "field_input.h"
#include "field_output.h"
#include "enum_elem_type.h"

#if defined(HAVE_TR1_UNORDERED_MAP)
#include <tr1/unordered_map>
#elif defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER) || defined(HAVE_UNORDERED_MAP)
#include <unordered_map>
#else
#include <map>
#endif


// Forward declarations
class SimulationSystem;
class Node;
class Elem;

namespace DFISE {
  class DFISE_MESH;
  class Element;
  class INFO;
  class DATASET;
}

/**
 * This class implements reading and writing solution data in the DF-ISE format.
 * The native grid/data file used by Synopsys sentaurus
 */

// ------------------------------------------------------------
// DFISEIO class definition
class DFISEIO : public FieldInput<SimulationSystem>,
      public FieldOutput<SimulationSystem>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  DFISEIO (SimulationSystem& system);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  DFISEIO (const SimulationSystem& system);


  /**
   * This method implements reading a mesh from a specified file
   * in DF-ISE format.
   */
  virtual void read (const std::string& );


  /**
   * This method implements writing a mesh to a specified DF-ISE file.
   */
  virtual void write (const std::string& );

private:

  DFISE::DFISE_MESH *ise_reader;

  // map node * to dfise node index
  std::map<Node *, unsigned int>       _node_to_dfise_node_index_map;

  void _set_mesh_element(const DFISE::INFO & grid_info, const DFISE::Element &);

  void _set_boundary(const DFISE::INFO & grid_info);

  // hold the boundary elems here
  std::vector<std::pair<Elem *, int> > _boundary_elem_regions;

  std::map<unsigned int, unsigned int> _node_id_to_dfise_node_index_map;

  void _set_region();


private:

  DFISE::DFISE_MESH *ise_writer;

#if defined(HAVE_UNORDERED_MAP)
  typedef std::unordered_multimap<unsigned int, const Elem *> _multimap_type;
#elif defined(HAVE_TR1_UNORDERED_MAP) || defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER)
  typedef std::tr1::unordered_multimap<unsigned int, const Elem *> _multimap_type;
#else
  typedef std::multimap<unsigned int, const Elem *>  _multimap_type;
#endif

  // private members to write df-ise

  void _build_element();


  /**
   * the '<' operator used in set/map key
   */
  struct _lt_pair_int
  {
    bool operator()(const std::pair<int, int> &e1, const std::pair<int, int> &e2) const
    {
      int e1_min = std::min(e1.first, e1.second);
      int e1_max = std::max(e1.first, e1.second);
      int e2_min = std::min(e2.first, e2.second);
      int e2_max = std::max(e2.first, e2.second);

      if( e1_min != e2_min ) return e1_min < e2_min;
      return e1_max < e2_max;
    }
  };

  /**
   * edge map
   */
  std::map< std::pair<int, int>, int, _lt_pair_int> _edges_map;

  /**
   * all the mesh edges
   */
  std::vector< std::pair<int, int> > _edges;

  /**
   * function to add an edge to _edges_map
   */
  bool _add_edge( int n1, int n2 );

  /**
   * build the edge index of a given face
   */
  std::vector<int> _edge_index(const Elem *face) const;



  /**
   * face map
   */
  _multimap_type _faces_map;

  /**
   * all the mesh faces
   */
  std::vector<const Elem *> _faces;

  /**
   * function to add an edge to _faces_map and _faces
   */
  bool _add_face( Elem * face );

  /**
   * build the face index of a given cell
   */
  std::vector<int> _face_index(const Elem *cell) const;



  // location
  std::vector<char> _location;

  // pre-computed element faces
  std::map<const Elem *, std::vector<int> > _elem_faces;


  void _build_region();

  /**
   * all the elements
   */
  std::vector<const Elem *> _elems;

  /**
   * boundary elements, should delete them
   */
  std::vector<const Elem *> _boundary_elems;

  /**
   * region + electrode boundary
   */
  std::vector<std::string> _regions;

  /**
   * region material + electrode boundary material "Contact"
   */
  std::vector<std::string> _materials;

  /**
   * region + electrode boundary element list
   */
  std::vector< std::vector<unsigned int> > _region_elements;

  /**
   * @return DFISE element code
   */
  int _elem_code(ElemType ) const;

  /**
   * build dfise data set by given variable
   */
  DFISE::DATASET * _build_dataset( const std::string & variable ,  const std::string & function);

  /**
   * free everything
   */
  void _write_on_exit();
};



// ------------------------------------------------------------
// DFISEIO inline members
inline
DFISEIO::DFISEIO (SimulationSystem& system) :
    FieldInput<SimulationSystem> (system),
    FieldOutput<SimulationSystem> (system)
{
  ise_reader = 0;
  ise_writer = 0;
}



inline
DFISEIO::DFISEIO (const SimulationSystem& system) :
    FieldOutput<SimulationSystem>(system)
{
  ise_reader = 0;
  ise_writer = 0;
}



#endif // #define __df_ise_io_h__
