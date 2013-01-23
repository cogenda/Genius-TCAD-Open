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

#include <cstdlib>
#include <limits>
#include <iomanip>

#include "mesh_base.h"
#include "boundary_info.h"
#include "semiconductor_region.h"
#include "solver_base.h"
#include "mob_monitor_hook.h"
#include "MXMLUtil.h"
#include "parallel.h"

#ifdef HAVE_VTK

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkConfigure.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"

#endif //HAVE_VTK


using PhysicalUnit::s;
using PhysicalUnit::eV;
using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::K;

/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
MobMonitorHook::MobMonitorHook ( SolverBase & solver, const std::string & name, void * param)
  : Hook ( solver, name ),  mesh(solver.get_system().mesh())
{
  this->_ddm_solver = false;
  this->solution_count=0;

  _t_start = 0;
  _t_end   = std::numeric_limits<double>::infinity();
  _surface = false;

  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for(std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
      parm_it != parm_list.end(); parm_it++)
  {
    if(parm_it->name() == "tstart" && parm_it->type() == Parser::REAL)
      _t_start=parm_it->get_real() * PhysicalUnit::s;
    if(parm_it->name() == "tstop" && parm_it->type() == Parser::REAL)
      _t_end=parm_it->get_real() * PhysicalUnit::s;
    if(parm_it->name() == "surface" && parm_it->type() == Parser::BOOL)
      _surface=parm_it->get_bool();
  }
}


/*----------------------------------------------------------------------
 * destructor, close file
 */
MobMonitorHook::~MobMonitorHook()
{}



/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void MobMonitorHook::on_init()
{
  this->_ddm_solver = false;

  if( _solver.solver_type() == SolverSpecify::DDML1 ||
      _solver.solver_type() == SolverSpecify::DDML2 ||
      _solver.solver_type() == SolverSpecify::EBML3 ||
      _solver.solver_type() == SolverSpecify::DDML1MIXA ||
      _solver.solver_type() == SolverSpecify::DDML2MIXA ||
      _solver.solver_type() == SolverSpecify::EBML3MIXA ||
      _solver.solver_type() == SolverSpecify::HALLDDML1
    )
  {
    _ddm_solver = true;
  }
#ifdef HAVE_VTK
  // buffer mesh node/edges, must run in parallel
  mesh.pack_nodes(points);

  std::vector< std::pair<unsigned int, unsigned int> > edge_array;

  if(_surface)
    mesh.pack_boundary_egeds(edge_array);
  else
    mesh.pack_egeds(edge_array);
  for(unsigned int n=0; n<edge_array.size(); ++n)
  {
    edges.insert( std::make_pair( edge_array[n], n ) );
  }
#endif
}


/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void MobMonitorHook::post_solve()
{
  // only monitor ddm solver
  if(!this->_ddm_solver) return;

  if(SolverSpecify::Type == SolverSpecify::TRANSIENT)
  {
    if( SolverSpecify::clock < _t_start || SolverSpecify::clock > _t_end ) return;
  }

  const SimulationSystem &system = _solver.get_system();

#ifdef HAVE_VTK
  std::string vtk_prefix = SolverSpecify::out_prefix+".mob.monitor";
  std::ostringstream vtk_filename;
  vtk_filename << vtk_prefix << '.' << this->solution_count << ".vtu";

  // use vtk library routine to export mesh and solution as base64 binary file
  grid = vtkUnstructuredGrid::New();

  nodes_to_vtk();
  edges_to_vtk();
  node_id_to_vtk();
  solution_to_vtk();

  // only processor 0 write VTK file
  if(Genius::processor_id() == 0)
  {
    std::cout<<"MobMonitorHook dump " << vtk_filename.str()<<std::endl;
    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetInput(grid);

    writer->SetFileName(vtk_filename.str().c_str());
    writer->Write();
    writer->Delete();
  }
  //clean up
  grid->Delete();

#endif

  this->solution_count++;
}


/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void MobMonitorHook::on_close()
{}




//---------------------------------------------------------------------



// private functions
#ifdef HAVE_VTK

void MobMonitorHook::nodes_to_vtk()
{
  if(Genius::processor_id() == 0)
  {
    vtkPoints* vtk_points = vtkPoints::New();

    // write out nodal data
    unsigned int n_nodes = mesh.n_nodes();
    for ( unsigned int n=0; n<n_nodes; ++n )
    {
      float tuple[3];
      {
        tuple[0] =  points[3*n+0]/um; //scale to um
        tuple[1] =  points[3*n+1]/um; //scale to um
        tuple[2] =  points[3*n+2]/um; //scale to um
      }
      vtk_points->InsertPoint(n, tuple);
    }

    grid->SetPoints(vtk_points);

    //clean up
    vtk_points->Delete();
  }
}


void MobMonitorHook::edges_to_vtk()
{
  if(Genius::processor_id() == 0)
  {
    grid->Allocate(edges.size());

    // insert edges to vtkUnstructuredGrid
    std::map< std::pair<unsigned int, unsigned int>, unsigned int >::const_iterator it = edges.begin();
    for(; it != edges.end(); ++it)
    {
      vtkIdList *pts = vtkIdList::New();
      pts->SetNumberOfIds(2);
      pts->SetId(0, it->first.first);
      pts->SetId(1, it->first.second);
      grid->InsertNextCell(VTK_LINE, pts);
      pts->Delete();
    }
  }
}


void MobMonitorHook::node_id_to_vtk()
{
  if (Genius::processor_id() == 0)
  {
    unsigned int n_nodes = mesh.n_nodes();
    //create vtk data array
    vtkIntArray *vtk_node_id_array = vtkIntArray::New();
    vtk_node_id_array->SetName("node");
    vtk_node_id_array->SetNumberOfValues(n_nodes);

    // save data into vtk data array
    for(unsigned int n=0; n<n_nodes; ++n)
    {
      // set scalar value to vtkfloatArray
      vtk_node_id_array->InsertValue(n, n);
    }

    grid->GetPointData()->AddArray(vtk_node_id_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_node_id_array->Delete();
  }
}

void MobMonitorHook::solution_to_vtk()
{
  std::vector<float> elec_mob(edges.size(), 0.0);
  std::vector<float> hole_mob(edges.size(), 0.0);
  std::vector<float> surface(edges.size(), 0.0);

  const SimulationSystem & system = _solver.get_system();
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = system.region(r);
    if(region->type() != SemiconductorRegion) continue;

    std::vector< std::pair<unsigned int, unsigned int> > edge;
    std::vector< std::pair<double, double> > mob;
    std::vector< double > weight;
    region->Mob_Evaluation(edge, mob, weight);

    for(unsigned int n=0; n<edge.size(); ++n)
    {
      if(edges.find(edge[n]) == edges.end()) continue;

      unsigned int index = edges.find(edge[n])->second;
      elec_mob[index] = mob[n].first/(cm*cm/V/s);
      hole_mob[index] = mob[n].second/(cm*cm/V/s);

      const FVM_Node * f1 = region->region_fvm_node(edge[n].first);
      const FVM_Node * f2 = region->region_fvm_node(edge[n].second);
      surface[index] = f1->cv_surface_area(f2)/(um*um);
    }
  }

  write_cell_scaler_solution(elec_mob, "elec_mob");
  write_cell_scaler_solution(hole_mob, "hole_mob");
  write_cell_scaler_solution(surface, "cv_area");
}


void MobMonitorHook::write_cell_scaler_solution( const std::vector<float> &sol, const std::string & sol_name)
{
  if (Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();
    vtk_sol_array->SetName(sol_name.c_str());
    vtk_sol_array->SetNumberOfValues(sol.size());

    // save data into vtk data array
    for(unsigned int n=0; n<sol.size(); ++n)
    {
      // set scalar value to vtkfloatArray
      vtk_sol_array->InsertValue(n, sol[n]);
    }

    grid->GetCellData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }
}


#endif

#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new MobMonitorHook ( solver, name, fun_data );
  }

}

#endif


