// $Id: unv_io.h 4951 2011-11-10 21:14:30Z roystgnr $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __unv_io_h__
#define __unv_io_h__


// C++ inludes
#include <vector>
#include <map>
#include <string>

// Local includes
#include "genius_common.h"
#include "field_input.h"
#include "field_output.h"
#include "simulation_system.h"


// Forward declarations
class MeshBase;
class Elem;

/**
 * The \p UNVIO class implements the Ideas \p UNV universal
 * file format.  This class enables both reading and writing
 * \p UNV files.
 */

// ------------------------------------------------------------
// UNVIO class definition
class UNVIO : public FieldInput<SimulationSystem>,
              public FieldOutput<SimulationSystem>
{

 public:
  
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  UNVIO (SimulationSystem& system);

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  UNVIO (const SimulationSystem& system);
  
  /**
   * Destructor.
   */
  virtual ~UNVIO ();
  
  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string& );
  
  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string& )  { genius_error(); }

  /**
   * Set the flag indicationg if we should be verbose.
   */
  bool & verbose ();

  
 private:
  

  /**
   * The actual implementation of the read function.
   * The public read interface simply decides which
   * type of stream to pass the implementation.
   */
  void read_implementation (std::istream& in_stream);
  
  
  void read_post_process();
  
  
  void build_boundary();
  
  /**
   * Clears the data structures to a pristine
   * state.
   */
  void clear();


  //-------------------------------------------------------------
  // read support methods
  /**
   * When reading, counting the nodes first
   * helps pre-allocation.  Also determine
   * whether we need to convert from "D" to "e".
   */
  void count_nodes (std::istream& in_file);
  
  /**
   * Method reads nodes from \p in_file and stores them in
   * vector<Node*> \p nodes in the order they come in.
   * The original node labels are being stored in
   * \p _assign_nodes in order to assign the elements to
   * the correct nodes later.  In addition, provided it is 
   * active, the \p MeshData gets to know the node id from
   * the Universal file, too.
   */
  void node_in (std::istream& in_file);

  /**
   * When reading, counting the elements first
   * helps pre-allocation.
   */
  void count_elements (std::istream& in_file);

  /**
   * Method reads elements and stores them in
   * \p std::vector<Elem*> \p _elements in the same order as they
   * come in. Within \p UNVIO, element labels are
   * ignored, but \p MeshData takes care of such things
   * (if active).
   */
  void element_in (std::istream& in_file);
  
  /**
   * read group info in the unv file
   */
  void read_groups (std::istream& in_file);

  /**
   * @returns \p false when error occured, \p true otherwise.
   * Adjusts the \p in_stream to the beginning of the
   * dataset \p ds_name.
   */
  bool beginning_of_dataset (std::istream& in_file, 
			     const std::string& ds_name) const;

  /**
   * Method for converting exponential notation
   * from "D" to "e", for example
   * \p 3.141592654D+00 \p --> \p 3.141592654e+00
   * in order to make it readable for C++.
   */
  Real D_to_e (std::string& number) const;


  //-------------------------------------------------------------
  // local data

  /**
   * should be be verbose?
   */
  bool _verbose;

  /**
   * maps node id's from UNV to internal.  Used when reading.
   */
  std::vector<unsigned int> _assign_nodes;

  /**
   * stores positions of relevant datasets in the file, should
   * help to re-read the data faster.  Used when reading.
   */
  std::map<std::string,std::streampos> _ds_position;
  
  /**
   * max dim of mesh element
   */
  unsigned int _n_dim;

  /**
   * total number of nodes, determined through \p count_nodes().
   * Primarily used when reading.
   */
  unsigned int _n_nodes;

  /**
   * total number of elements, determined through 
   * \p count_elements().  Primarily used when reading.
   */
  unsigned int _n_elements;
  
  
  /**
   * number of groups in the mesh
   * group may be: elem with dim (region), elem with dim-1 (boundary)
   */
  unsigned int _n_groups;
  
  /**
   * group name
   */
  std::vector<std::string> _groups;
 
  /**
   * map element id to group index
   */
  std::map<unsigned int, unsigned int> _element_in_group;

  /// groupd id -> region index
  std::map<unsigned int, unsigned int> _regions;
  
  // group id -> boundary index
  //std::set<unsigned int> _boundaries;
  

  /**
   * label for the node dataset
   */
  static const std::string _label_dataset_nodes;

  /**
   * label for the element dataset
   */
  static const std::string _label_dataset_elements;

  /**
   * label for the group dataset
   */
  static const std::string _label_dataset_groups;  

  /**
   * whether we need to convert notation of exponentials.
   * Used when reading.
   */
  bool _need_D_to_e;


};



// ------------------------------------------------------------
// MeshIO inline members
inline
UNVIO::UNVIO (SimulationSystem& system) :
    FieldInput<SimulationSystem> (system),
    FieldOutput<SimulationSystem> (system)
{
  _n_dim       = 0;
  _n_nodes     = 0;
  _n_elements  = 0;
  _n_groups    = 0;
  _need_D_to_e = true;
}



inline
UNVIO::UNVIO (const SimulationSystem& system) :
    FieldOutput<SimulationSystem>(system)
{
  _n_dim       = 0; 
  _n_nodes     = 0;
  _n_elements  = 0;
  _n_groups    = 0;
  _need_D_to_e = true;
}



inline
UNVIO::~UNVIO ()
{
  this->clear ();
}



inline
bool & UNVIO::verbose ()
{
  return _verbose;
}






#endif // #define __unv_io_h__
