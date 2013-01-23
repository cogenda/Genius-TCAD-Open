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

//  $Id: simulation_system.h,v 1.31 2008/07/09 05:58:16 gdiso Exp $

#ifndef __simulation_system_h__
#define __simulation_system_h__


#include "vector_value.h"
#include "enum_solution.h"
#include "enum_solver_specify.h"
#include "error_vector.h"
#include "physical_unit.h"
#include "interpolation_base.h"

namespace Parser {
class InputParser;
class Card;
}

class MeshBase;
class SimulationRegion;
class BoundaryConditionCollector;
class ElectricalSource;
class FieldSource;
class SPICE_CKT;

/**
 * @brief the main structure for mesh and solution data storage
 */
class SimulationSystem
{
public:

  /**
   * @brief empty constructor
   * @param mesh    Mesh writable reference
   *
   * @note set PhysicalUnit structure and T_external with default value
   */
  SimulationSystem(MeshBase & mesh);

  /**
   * @brief constructor
   * @param mesh    Mesh writable reference
   * @param decks   Parser::InputParser writable reference
   *
   * @note also build PhysicalUnit structure and T_external
   */
  SimulationSystem(MeshBase & mesh, Parser::InputParser & decks);


  /**
   * @brief destructor
   */
  virtual ~SimulationSystem();

  /**
   * clear system data, including mesh, region and bc data structure
   */
  void clear(bool clear_mesh=true);

  /**
   * @return true iff we run in the resistive metal mode
   */
  bool resistive_metal_mode() const
  { return _resistive_metal_mode; }

  /**
   * set resistive metal mode
   */
  void set_resistive_metal_mode(bool flag)
  { _resistive_metal_mode = flag; }

  /**
   * @brief build the simulation system from mesh and mesh boundary
   */
  void build_simulation_system();

  /**
   * @brief  data initialization for each region
   * @note   call it after build_simulation_system() when doping concentration and mole fraction have been set.
   */
  void init_region();

  /**
   * @brief  data re-initialization for each region
   * @note   call it after data file (i.e. cgns or tif)  import
   */
  void reinit_region_after_import();

  /**
   * post process after init_region / reinit_region_after_import
   */
  void init_region_post_process();

  /**
   * @brief build boundary condition array
   */
  void build_boundary_conditions();

  /**
   * @return the number of region (subdomain)
   */
  unsigned int n_regions() const;

  /**
   * @return the cluster of region (subdomain) for partition
   */
  std::vector< std::vector<unsigned int > > build_subdomain_cluster();

  /**
   * the enveriment temperature
   */
  double T_external() const
  { return  _T_external;}

  /**
   * @brief the main saving / loading mechanism, from cgnsfile
   */
  void import_cgns(const std::string& filename);

  /**
   * @brief load system from xml vtk file, only for debug reason
   */
  void import_vtk(const std::string& filename);

  /**
   * @brief load device information from Silvaco file, which is used by Silvaco ATLAS and other tools.
   */
  void import_silvaco(const std::string& filename);

  /**
   * @brief load device information from TIF file, which is used by synopsys Medici and Tsuprem.
   */
  void import_tif(const std::string& filename);

  /**
   * @brief load device information from TIF3D file, which is used by our GDSII to 3D model tool.
   */
  void import_tif3d(const std::string& filename);

  /**
   * @brief load device information from DF-ISE file, which is used by synopsys sentaurus.
   */
  void import_ise(const std::string& filename);


  /**
   * @brief save complete system information to cgns file, which can be loaded again
   */
  void export_cgns(const std::string& filename) const;

  /**
   * @brief save mesh and doping to df-ise file
   */
  void export_ise(const std::string& filename) const;

  /**
   * @brief export solution to vtk file
   */
  void export_vtk(const std::string& filename, bool ascii) const;

  /**
   * @brief write geometry and material info to gdml file
   */
  void export_gdml_surface(const std::string& filename) const;

  /**
   * @brief write geometry and material info to gdml file
   */
  void export_gdml_body(const std::string& filename) const;

  /**
   * @brief save node location to file
   */
  void export_node_location(const std::string& filename, const PetscScalar unit, const bool number=true) const;

  /**
   * @return true if the system is empty
   */
  bool empty() const { return _simulation_regions.empty(); }

  /**
   * @return the pointer to nth region
   */
  SimulationRegion * region(unsigned int n);

  /**
   * @return the const pointer to nth region
   */
  const SimulationRegion * region(unsigned int n) const;

  /**
   * @return the pointer to region with label region_name
   */
  SimulationRegion * region(const std::string & region_label);

  /**
   * @return the const pointer to region with label region_name
   */
  const SimulationRegion * region(const std::string & region_label) const;

  /**
   * @return true when compound semiconductor masterial exist
   */
  bool has_single_compound_semiconductor_region() const;

  /**
   * @return true when compound semiconductor masterial exist
   */
  bool has_complex_compound_semiconductor_region() const;

  /**
   * @return writable reference to mesh
   */
  MeshBase & mesh()
  { return _mesh; }

  /**
   * @return const reference to mesh
   */
  const MeshBase & mesh() const
  { return _mesh; }

  /**
   * @return dim of the mesh
   */
  unsigned int dim() const;

  /**
   * @return Boundary Condition Collector pointer
   */
  BoundaryConditionCollector  * get_bcs()
  { return _bcs;}

  /**
   * @return const Boundary Condition Collector pointer
   */
  const BoundaryConditionCollector  * get_bcs() const
  { return _bcs;}

  /**
   * @return pointer to ElectricalSource
   */
  ElectricalSource  * get_electrical_source()
  { return _electrical_source;}

  /**
   * @return const pointer to ElectricalSource
   */
  const ElectricalSource  * get_electrical_source() const
  { return _electrical_source;}

  /**
   * @return pointer to FieldSource
   */
  FieldSource * get_field_source()
  { return _field_source; }

  /**
   * @return const pointer to FieldSource
   */
  const FieldSource * get_field_source() const
  { return _field_source; }

  /**
   * @return magnetic field
   */
  const VectorValue<double> & get_magnetic_field() const
  { return _magnetic_field; }

  /**
   * @return circuit data
   */
  SPICE_CKT * get_circuit()
  { return _spice_ckt; }

  /**
   * @return circuit data
   */
  const SPICE_CKT * get_circuit() const
  { return _spice_ckt; }

  /**
   * @brief output debug information
   */
  void print_info (std::ostream& os=std::cout) const;

  /**
   * @brief Synchronized output debug information
   */
  void sync_print_info () const;

  /**
   * @return flag indicate that we have global z.width parameter
   */
  bool global_z_width() const { return _global_z_width; }

  /**
   * @return the dimension in Z direction.
   * @note only used for 2D mesh which lies on xy plane
   */
  double z_width() const;

  /**
   * @return ture when system is self consistent (semiconductor region satisfy poisson's equation)
   */
  bool self_consistent() const;

  /**
   * compute the error for each
   * cell and place it in the "error_per_cell" vector.
   */
  void estimate_error (const Parser::Card &c, ErrorVector & error_per_cell) const;

  /**
   * fill system data (mesh and node data) into interpolator for later usage
   * type can be -- linear interpolation
   *             -- signed log interpolation
   *             -- asinh interpolation
   */
  void fill_interpolator(InterpolationBase *, const std::string &, InterpolationBase::InterpolationType /* type */) const;

  /**
   * get data from interpolator after mesh refinement
   */
  void do_interpolation(const InterpolationBase *, const std::string &);

  /**
   * set unique solver name to _solver_active_history
   */
  void record_active_solver(SolverSpecify::SolverType solver_type)
  { _solver_active_history.push_back(solver_type); }

  /**
   * @return the sequence of solves which have performed on this system.
   */
  const std::vector<SolverSpecify::SolverType> & solve_history() const
  { return _solver_active_history; }


  /**
   * approx memory usage
   */
  virtual size_t memory_size() const;

private:

  /**
   * the enveriment temperature
   */
  double _T_external;

  /**
   * the reference to mesh
   */
  MeshBase & _mesh;

  /**
   * flag for cylindrically symmetric mesh
   */
  bool _cylindrical_mesh;
  
  /**
   * enable distributed mesh 
   */
  bool _distributed_mesh;

  /**
   * the array for all the simulation region.
   * each region contains mesh and mesh data
   */
  std::vector<SimulationRegion *>  _simulation_regions;

  /**
   * create resistive metal regions instead of electrode region
   */
  bool _resistive_metal_mode;

  /**
   * partition resistive metal region into same block
   */
  bool _block_partition;

  /**
   * data structure for fvm solver
   * only build nodes which belongs to local processor
   */
  void build_region_fvm_mesh();

  /**
   * all the boundary conditions
   */
  BoundaryConditionCollector    * _bcs;

  /**
   * all the predefined electrical source
   */
  ElectricalSource              * _electrical_source;

  /**
   * Field source: light or particle source
   */
  FieldSource                   * _field_source;


  /**
   * static magnetic field
   */
  VectorValue<double>            _magnetic_field;


  /**
   * mixed simulation with spice
   */
  SPICE_CKT                     * _spice_ckt;

  /**
   * flag to indicate that we have a global z.width
   */
  bool    _global_z_width;

  /**
   * when 2D mesh is used, we may need the dimension in z direction
   */
  double  _z_width;

  /**
   * each solver should record itself in this _solver_active_history vector when active
   * we can determine the solve sequence by this vector
   */
  std::vector<SolverSpecify::SolverType> _solver_active_history;

};





#endif

