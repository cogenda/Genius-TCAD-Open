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



#ifndef __vtk2_io_h__
#define __vtk2_io_h__

// C++ includes
#include <set>
#include <map>


// Local includes
#include "genius_common.h"
#include "field_input.h"
#include "field_output.h"
#include "simulation_system.h"
#include "simulation_region.h"

// Forward declarations
class MeshBase;


#ifdef HAVE_VTK
class vtkUnstructuredGrid;
#endif


/**
 * This class implements reading and writing solution data in the VTK format.
 * Format description:
 * cf. <a href="http://www.vtk.org/">VTK home page</a>.
 */

// ------------------------------------------------------------
// VTK2IO class definition
class VTK2IO : public FieldInput<SimulationSystem>,
               public FieldOutput<SimulationSystem>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  VTK2IO (SimulationSystem& system);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  VTK2IO (const SimulationSystem& system, const std::vector<std::string> & variables);


  /**
   * This method implements reading a mesh from a specified file
   * in VTK format.
   */
  virtual void read (const std::string& ) {}

  /**
   * This method implements writing a mesh to a specified fil
   * in VTK format
   */
  virtual void write (const std::string& );

private:

  std::set<std::string> _variables;

  // boundary info
  std::vector<unsigned int>       _el;
  std::vector<unsigned short int> _sl;
  std::vector<short int>          _il;

  // < <region, id>, id >
  std::map< std::pair<unsigned int, unsigned int>, unsigned int > _region_node_id_map;

private:

#ifdef HAVE_VTK
  /**
   * write the nodes from the mesh into a vtkUnstructuredGrid
   */
  void nodes_to_vtk(const MeshBase& mesh);

  /**
   * write the cells from the mesh into a vtkUnstructuredGrid
   */
  void cells_to_vtk(const MeshBase& mesh);

  /**
   * write out mesh data to the VTK file, this might come in handy to display
   * material data and partition information
   */
  void meshinfo_to_vtk(const MeshBase& mesh);

  /**
   * write mesh relative solution data to vtk file
   */
  void solution_to_vtk();

  /**
   * write extra region name/material information and boundary information into xml vtk file
   */
  std::string export_extra_info();

  /**
   * determin solutions to output
   */
  void build_export_solutions();

  /**
   * aux function to write a solution to vtkUnstructuredGrid
   */
  void write_node_solution(const std::string & sol_name);

  /**
   * aux function to write a scaler solution to vtkUnstructuredGrid
   */
  void write_node_solution_scalar(const std::string & sol_name, const std::string & sol_unit);

  /**
   * aux function to write a complex solution to vtkUnstructuredGrid
   */
  void write_node_solution_complex(const std::string & sol_name, const std::string & sol_unit);

  /**
   * aux function to write a vector solution to vtkUnstructuredGrid
   */
  void write_node_solution_vector(const std::string & sol_name, const std::string & sol_unit);

  /**
   * aux function to write a tensor solution to vtkUnstructuredGrid
   */
  void write_node_solution_tensor(const std::string & sol_name, const std::string & sol_unit);


  /**
   * aux function to write a scaler solution to vtkUnstructuredGrid
   */
  void write_node_solution_scalar(const std::string & sol_name, std::vector<float> &sol_value);

  /**
   * pointer to the VTK grid
   */
  vtkUnstructuredGrid * _vtk_grid;

  class XMLUnstructuredGridWriter;

#endif

};



// ------------------------------------------------------------
// VTKIO inline members
inline
VTK2IO::VTK2IO (SimulationSystem& system) :
    FieldInput<SimulationSystem> (system),
    FieldOutput<SimulationSystem> (system)
{
#ifdef HAVE_VTK
  _vtk_grid = 0;
#endif
}



inline
VTK2IO::VTK2IO (const SimulationSystem& system, const std::vector<std::string> & variables) :
    FieldOutput<SimulationSystem>(system)
{
  for(unsigned int n=0; n<variables.size(); ++n)
    _variables.insert(variables[n]);
#ifdef HAVE_VTK
  _vtk_grid = 0;
#endif
}



#endif // #define __vtk2_io_h__


