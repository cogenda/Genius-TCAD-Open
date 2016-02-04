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

#include "jflux1.h"

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
using PhysicalUnit::kb;
using PhysicalUnit::e;

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
  const SimulationSystem & system = _solver.get_system();

  double z_width = system.z_width();
  double T = system.T_external();
  double Vt = kb*T/e;
  double cmc = std::pow(cm,-3);

  std::vector<float> elec_mob_array(edges.size(), 0.0);
  std::vector<float> hole_mob_array(edges.size(), 0.0);
  std::vector<float> In_array(edges.size(), 0.0);
  std::vector<float> Ip_array(edges.size(), 0.0);
  std::vector<float> I_array(edges.size(), 0.0);
  std::vector<float> surface_array(edges.size(), 0.0);


  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = system.region(r);
    if(region->type() != SemiconductorRegion) continue;
    const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
    Material::MaterialSemiconductor * mt = semi_region->material();

    std::vector< std::pair<unsigned int, unsigned int> > edge;
    std::vector< std::pair<double, double> > mob;
    std::vector< double > weight;
    semi_region->Mob_Evaluation(edge, mob, weight);

    for(unsigned int n=0; n<edge.size(); ++n)
    {
      if(edges.find(edge[n]) == edges.end()) continue;

      unsigned int index = edges.find(edge[n])->second;
      double mun = mob[n].first;
      double mup = mob[n].second;
      elec_mob_array[index] = mun/(cm*cm/V/s);
      hole_mob_array[index] = mup/(cm*cm/V/s);

      const FVM_Node * f1 = region->region_fvm_node(edge[n].first);
      const FVM_Node * f2 = region->region_fvm_node(edge[n].second);

      const FVM_NodeData * n1_data = f1->node_data();
      const FVM_NodeData * n2_data = f2->node_data();

      double length = f1->distance(f2);

      mt->mapping(f1->root_node(), n1_data, SolverSpecify::clock);
      double V1 = n1_data->psi();
      double n1 = n1_data->n();
      double p1 = n1_data->p();
      double Ec1 =  -(e*V1 + n1_data->affinity() + mt->band->EgNarrowToEc(p1, n1, T) + kb*T*log(n1_data->Nc()));
      double Ev1 =  -(e*V1 + n1_data->affinity() - mt->band->EgNarrowToEv(p1, n1, T) - kb*T*log(n1_data->Nv()) + mt->band->Eg(T));

      mt->mapping(f2->root_node(), n2_data, SolverSpecify::clock);
      double V2 = n2_data->psi();
      double n2 = n2_data->n();
      double p2 = n2_data->p();
      double Ec2 =  -(e*V2 + n2_data->affinity() + mt->band->EgNarrowToEc(p2, n2, T) + kb*T*log(n2_data->Nc()));
      double Ev2 =  -(e*V2 + n2_data->affinity() - mt->band->EgNarrowToEv(p2, n2, T) - kb*T*log(n2_data->Nv()) + mt->band->Eg(T));

      double area = f1->cv_surface_area(f2)*z_width;
      double In = mun*In_dd(Vt,(Ec2-Ec1)/e,n1,n2,length)*area;
      double Ip = mup*Ip_dd(Vt,(Ev2-Ev1)/e,p1,p2,length)*area;

      double I  = In+Ip;

      surface_array[index] = area/(um*um);


      In_array[index] = In/(A);
      Ip_array[index] = Ip/(A);
      I_array[index] = I/(A);
    }
  }

  write_cell_scaler_solution(elec_mob_array, "elec_mob");
  write_cell_scaler_solution(hole_mob_array, "hole_mob");
  write_cell_scaler_solution(In_array, "In");
  write_cell_scaler_solution(Ip_array, "Ip");
  write_cell_scaler_solution(I_array, "current");
  write_cell_scaler_solution(surface_array, "cv_area");
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


