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



#ifndef __vtk_io_h__
#define __vtk_io_h__

// C++ includes
#include <map>
#include <fstream>

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
// VTKIO class definition
class VTKIO : public FieldInput<SimulationSystem>,
          public FieldOutput<SimulationSystem>
{
public:
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  VTKIO (SimulationSystem& system);

  /**
   * Constructor.  Takes a read-only reference to a mesh object.
   * This is the constructor required to write a mesh.
   */
  VTKIO (const SimulationSystem& system);


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

  // boundary info
  std::vector<unsigned int>       _el;
  std::vector<unsigned short int> _sl;
  std::vector<short int>          _il;


  // < <region, id>, id >
  std::map< std::pair<unsigned int, unsigned int>, unsigned int > _region_node_id_map;

private:

  /**
   * sometimes data may have numerical error. when the data range is very small,
   * the numerical error will affect plot effect. here we try to smooth the numerical error
   * when max-min \< tol, set all the data to 0.5(max+min)
   */
  void _smooth_numerical_error(std::map<unsigned int, float> &, float tol);

#ifdef HAVE_VTK
  /**
   * write the nodes from the mesh into a vtkUnstructuredGrid
   */
  void nodes_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* grid);

  /**
   * write the cells from the mesh into a vtkUnstructuredGrid
   */
  void cells_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* grid);

  /**
   * write out mesh data to the VTK file, this might come in handy to display
   * material data and partition information
   */
  void meshinfo_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* grid);

  /**
   * write mesh relative solution data to vtk file
   */
  void solution_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* grid);


  /**
   * write extra region name/material information and boundary information into xml vtk file
   */
  std::string export_extra_info();

  /**
   * aux function to write a scaler solution to vtkUnstructuredGrid
   */
  void write_node_scaler_solution(const std::vector< std::vector<unsigned int> >& ,
                                  const std::vector< std::vector<float> > &,
                                  const std::string & sol_name, vtkUnstructuredGrid* grid);

  /**
   * aux function to write a complex solution to vtkUnstructuredGrid
   */
  void write_node_complex_solution(const std::vector< std::vector<unsigned int> >& ,
                                   const std::vector< std::vector<std::complex<float> > > & ,
                                   const std::string & sol_name, vtkUnstructuredGrid* grid);

  /**
   * aux function to write a vector solution to vtkUnstructuredGrid
   */
  void write_node_vector_solution (const std::vector< std::vector<unsigned int> >& ,
                                   const std::vector< std::vector<float> > &,
                                   const std::vector< std::vector<float> > &,
                                   const std::vector< std::vector<float> > &,
                                   const std::string & sol_name, vtkUnstructuredGrid* grid);

  /**
   * aux function to write a cell based scaler solution to vtkUnstructuredGrid
   */
  void write_cell_scaler_solution(const std::vector<unsigned int> & order, 
                                  std::vector<float> &sol, 
                                  const std::string & sol_name, 
                                  vtkUnstructuredGrid* grid);

  /**
   * aux function to write a cell based vector solution to vtkUnstructuredGrid
   */
  void write_cell_vector_solution (const std::vector<unsigned int> & order,
                                   std::vector<float > & sol_x,
                                   std::vector<float > & sol_y,
                                   std::vector<float > & sol_z,
                                   const std::string & sol_name, vtkUnstructuredGrid* grid);

  /**
   * pointer to the VTK grid
   */
  vtkUnstructuredGrid* _vtk_grid;

  class XMLUnstructuredGridWriter;
#endif


};



// ------------------------------------------------------------
// VTKIO inline members
inline
VTKIO::VTKIO (SimulationSystem& system) :
    FieldInput<SimulationSystem> (system),
    FieldOutput<SimulationSystem> (system)
{
#ifdef HAVE_VTK
  _vtk_grid = NULL;
#endif
}



inline
VTKIO::VTKIO (const SimulationSystem& system) :
    FieldOutput<SimulationSystem>(system)
{
#ifdef HAVE_VTK
  _vtk_grid = NULL;
#endif
}



#endif // #define __vtk_io_h__
