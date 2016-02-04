// $Id: gmsh_io.h 4951 2011-11-10 21:14:30Z roystgnr $

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



#ifndef __gmsh_io_h__
#define __gmsh_io_h__


#include "genius_common.h"
#include "field_input.h"
#include "field_output.h"
#include "simulation_system.h"

// Forward declarations
class MeshBase;
class Elem;


/**
 * This class implements writing geometry info to
 * Geometry Description Markup Language (GDML) file
 * This file format is used by GEANT4 for build its detector.
 */
class GmshIO : public FieldInput<SimulationSystem>,
               public FieldOutput<SimulationSystem>
{
 public:

  /**
   * Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements via the read() command.
   */
  GmshIO (SimulationSystem& system);
  
  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  GmshIO (const SimulationSystem& mesh);

  /**
   * Reads in a mesh in the Gmsh *.msh format
   * from the ASCII file given by name.
   *
   * Note that for this method to work (in 2d and 3d) you have to
   * explicitly set the mesh dimension prior to calling GmshIO::read()
   * and that Mesh::prepare_for_use() must be called after reading the
   * mesh and before using it.
   */
  virtual void read (const std::string& name);

  /**
   * This method implements writing a mesh to a specified file
   * in the Gmsh *.msh format.
   */
  virtual void write (const std::string& name)
  { genius_error(); }

private:
  void read_info (const std::string& name);
  
  // physical, name-material
  std::map< int, std::pair<std::string, std::string> > region_info;
  // physical, name
  std::map< int, std::string > boundary_info;
   
  
  void read_mesh (const std::string& name);
};



// ------------------------------------------------------------
// GmshIO inline members
inline
GmshIO::GmshIO (SimulationSystem& system) :
    FieldInput<SimulationSystem> (system),
    FieldOutput<SimulationSystem> (system)
{
}



inline
GmshIO::GmshIO (const SimulationSystem& system) :
    FieldOutput<SimulationSystem>(system)
{
}


#endif // #define __gmsh_io_h__
