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
#include "fvm_flex_nonlinear_solver.h"
#include "ddm_monitor_hook.h"
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
using PhysicalUnit::C;
using PhysicalUnit::W;


#ifdef HAVE_VTK
class DDMMonitorHook::XMLUnstructuredGridWriter : public vtkXMLUnstructuredGridWriter
{
protected:
  XMLUnstructuredGridWriter() :vtkXMLUnstructuredGridWriter()
  {}

  virtual int WriteHeader()
  {
    ostream& os = *(this->Stream);
    os << _header;

    return vtkXMLUnstructuredGridWriter::WriteHeader();
  }
public:
  static XMLUnstructuredGridWriter* New()
  {
    return new XMLUnstructuredGridWriter;
  }

  void setExtraHeader(const std::string &header)
  {
    _header = header;
  }
private:
  std::string _header;
};

inline vtkIdType elem_type_vtk(const Elem * elem)
{
  vtkIdType celltype = VTK_EMPTY_CELL; // initialize to something to avoid compiler warning

  switch(elem->type())
  {
  case EDGE2:
  case EDGE2_FVM:
    celltype = VTK_LINE;
    break;
  case EDGE3:
    celltype = VTK_QUADRATIC_EDGE;
    break;// 1
  case TRI3:
  case TRI3_FVM:
  case TRI3_CY_FVM:
    celltype = VTK_TRIANGLE;
    break;// 3
  case TRI6:
    celltype = VTK_QUADRATIC_TRIANGLE;
    break;// 4
  case QUAD4:
  case QUAD4_FVM:
  case QUAD4_CY_FVM:
    celltype = VTK_QUAD;
    break;// 5
  case QUAD8:
    celltype = VTK_QUADRATIC_QUAD;
    break;// 6
  case TET4:
  case TET4_FVM:
    celltype = VTK_TETRA;
    break;// 8
  case TET10:
    celltype = VTK_QUADRATIC_TETRA;
    break;// 9
  case HEX8:
  case HEX8_FVM:
    celltype = VTK_HEXAHEDRON;
    break;// 10
  case HEX20:
    celltype = VTK_QUADRATIC_HEXAHEDRON;
    break;// 12
  case PRISM6:
  case PRISM6_FVM:
    celltype = VTK_WEDGE;
    break;// 13
  case PRISM15:
    celltype = VTK_HIGHER_ORDER_WEDGE;
    break;// 14
  case PRISM18:
    break;// 15
  case PYRAMID5:
  case PYRAMID5_FVM:
    celltype = VTK_PYRAMID;
    break;// 16
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
  case QUAD9:
    celltype = VTK_BIQUADRATIC_QUAD;
    break;
#else
  case QUAD9:
#endif
  case EDGE4:
  case HEX27:
  case NODEELEM:
  case INVALID_ELEM:
  default:
    {
      std::cerr<<"element type "<<elem->type()<<" not implemented"<<std::endl;
      genius_error();
    }
  }

  return celltype;
}


#endif


/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
DDMMonitorHook::DDMMonitorHook ( SolverBase & solver, const std::string & name, void * param)
    : Hook ( solver, name ),  mesh(solver.get_system().mesh())
{
  this->_poisson_solver = false;
  this->_ddm_solver = false;
  this->_ric_solver = false;
  this->solution_count=0;
  this->iteration_count=0;

  _t_start = 0;
  _t_end   = std::numeric_limits<double>::infinity();

  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for(std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
      parm_it != parm_list.end(); parm_it++)
  {
    if(parm_it->name() == "tstart" && parm_it->type() == Parser::REAL)
      _t_start=parm_it->get_real() * PhysicalUnit::s;
    if(parm_it->name() == "tstop" && parm_it->type() == Parser::REAL)
      _t_end=parm_it->get_real() * PhysicalUnit::s;
  }

  mesh.boundary_info->build_on_processor_side_list (el, sl, il);
}


/*----------------------------------------------------------------------
 * destructor, close file
 */
DDMMonitorHook::~DDMMonitorHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void DDMMonitorHook::on_init()
{
  this->_poisson_solver = false;
  this->_ddm_solver = false;
  this->_tid_drift_solver = false;
  this->_ric_solver = false;
  this->_dictat_solver = false;

  if(_solver.solver_type() == SolverSpecify::POISSON)
  {
    _poisson_solver = true;
  }

  if( _solver.solver_type() == SolverSpecify::DDML1 ||
      _solver.solver_type() == SolverSpecify::DDML2 ||
      _solver.solver_type() == SolverSpecify::EBML3 ||
      _solver.solver_type() == SolverSpecify::DDML1MIXA ||
      _solver.solver_type() == SolverSpecify::DDML2MIXA ||
      _solver.solver_type() == SolverSpecify::EBML3MIXA ||
      _solver.solver_type() == SolverSpecify::DDML1MIX  ||
      _solver.solver_type() == SolverSpecify::DDML2MIX  ||
      _solver.solver_type() == SolverSpecify::EBML3MIX  ||
      _solver.solver_type() == SolverSpecify::HALLDDML1 ||
      _solver.solver_type() == SolverSpecify::DENSITY_GRADIENT
    )
  {
    _ddm_solver = true;
  }
  
  if(_solver.solver_type() == SolverSpecify::TID_DRIFT)
  {
    _tid_drift_solver = true;
  }

  if(_solver.solver_type() == SolverSpecify::RIC )
  {
    _ric_solver = true;
  }

  if(_solver.solver_type() == SolverSpecify::DICTAT)
  {
    _dictat_solver = true;
  }
}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void DDMMonitorHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void DDMMonitorHook::post_solve()
{
  this->solution_count++;
  this->iteration_count=0;
}


/*----------------------------------------------------------------------
 *  This is executed before each (nonlinear) iteration
 */
void DDMMonitorHook::pre_iteration()
{

}


/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void DDMMonitorHook::post_check(void * f, void * x, void * dx, void * w, bool & , bool &)
{
  // only monitor ddm solver / poisson / ric /tid solver
  if(!this->_ddm_solver && 
     !this->_poisson_solver && 
     !this->_tid_drift_solver && 
     !this->_ric_solver && 
     !this->_dictat_solver) return;

  if(SolverSpecify::Type == SolverSpecify::TRANSIENT)
  {
    if( SolverSpecify::clock < _t_start || SolverSpecify::clock > _t_end ) return;
  }

  const SimulationSystem &system = _solver.get_system();

#ifdef HAVE_VTK
  std::string vtk_prefix = SolverSpecify::out_prefix+".monitor";
  std::ostringstream vtk_filename;
  vtk_filename << vtk_prefix << '.' << this->solution_count<< '.' << this->iteration_count << ".vtu";

  // use vtk library routine to export mesh and solution as base64 binary file
  grid = vtkUnstructuredGrid::New();

  nodes_to_vtk();
  cells_to_vtk();
  meshinfo_to_vtk();

  if(this->_ddm_solver)
    solution_to_vtk_ddm(f, x, dx, w);

  if(this->_poisson_solver)
    solution_to_vtk_poisson(f, x, dx, w);
  
  if(this->_tid_drift_solver)
    solution_to_vtk_tid_drift(f, x, dx, w);

  if(this->_ric_solver)
    solution_to_vtk_ric(f, x, dx, w);

  if(this->_dictat_solver)
    solution_to_vtk_dictat(f, x, dx, w);

  // only processor 0 write VTK file
  if(Genius::processor_id() == 0)
  {
    std::cout<<"DDMMonitorHook dump " << vtk_filename.str()<<std::endl;
    XMLUnstructuredGridWriter* writer = XMLUnstructuredGridWriter::New();
    writer->SetInput(grid);
    writer->setExtraHeader(this->export_extra_info());

    writer->SetFileName(vtk_filename.str().c_str());
    writer->Write();
    writer->Delete();
  }
  //clean up
  grid->Delete();

#endif

#if 0
  const FVM_FlexNonlinearSolver & nonlinear_solver = dynamic_cast<FVM_FlexNonlinearSolver &>(_solver);

  std::string matrix_prefix = SolverSpecify::out_prefix+".monitor";
  std::ostringstream matrix_filename;
  matrix_filename << matrix_prefix << '.' << this->solution_count<< '.' << this->iteration_count << ".mat";


  nonlinear_solver.dump_matrix_asc(nonlinear_solver.jacobian_matrix(), matrix_filename.str());


  std::string vector_prefix = SolverSpecify::out_prefix+".monitor";
  std::ostringstream vector_filename;
  vector_filename << vector_prefix << '.' << this->solution_count<< '.' << this->iteration_count << ".vec";

  nonlinear_solver.dump_vector_petsc(nonlinear_solver.rhs_vector(), vector_filename.str());

#endif

  this->iteration_count++;

}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void DDMMonitorHook::on_close()
{}


//---------------------------------------------------------------------

// private functions
#ifdef HAVE_VTK

void DDMMonitorHook::nodes_to_vtk()
{
  const SimulationSystem & system = _solver.get_system();

  unsigned int region_n_nodes = 0;
  region_node_id_map.clear();

  std::vector<unsigned int> system_nodes;
  for( unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);

    std::vector<unsigned int> nodes;
    region->region_node(nodes);

    for(unsigned int n=0; n<nodes.size(); ++n)
    {
      unsigned int node_id = nodes[n];
      region_node_id_map.insert( std::make_pair(std::make_pair(r,node_id), region_n_nodes+n));
    }
    region_n_nodes += nodes.size();
    system_nodes.insert(system_nodes.end(), nodes.begin(), nodes.end());
  }


  if(Genius::processor_id() == 0)
  {
    vtkPoints* points = vtkPoints::New();

    // write out nodal data
    for ( unsigned int n=0; n<system_nodes.size(); ++n )
    {
      float tuple[3];
      {
        const Point & p = mesh.point(system_nodes[n]);
        tuple[0] =  p[0]/um; //scale to um
        tuple[1] =  p[1]/um; //scale to um
        tuple[2] =  p[2]/um; //scale to um
      }
      points->InsertPoint(n, tuple);
    }

    grid->SetPoints(points);

    //clean up
    points->Delete();
  }

}


void DDMMonitorHook::cells_to_vtk()
{
  const SimulationSystem & system = _solver.get_system();
  const MeshBase & mesh = system.mesh();
  
  // must run in parallel
  std::vector<int> vtk_cell_conns;
  std::vector<unsigned int> vtk_cell_ids;
  std::vector<unsigned int> vtk_cell_regions;
  {
    MeshBase::const_element_iterator       it  = mesh.active_this_pid_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_this_pid_elements_end();
    for ( ; it != end; ++it)
    {
      const Elem *elem  = (*it);
      vtkIdType celltype = elem_type_vtk(elem);

      // get the connectivity for this element
      std::vector<unsigned int> conn;
      elem->connectivity(0,VTK,conn);

      vtk_cell_ids.push_back(elem->id());
      vtk_cell_regions.push_back(elem->subdomain_id());
      vtk_cell_conns.push_back(celltype);
      vtk_cell_conns.push_back(conn.size());
      for(unsigned int i=0;i<conn.size();++i)
        vtk_cell_conns.push_back(conn[i]);

    } // end loop over active elements

    Parallel::gather(0, vtk_cell_ids);
    Parallel::gather(0, vtk_cell_regions);
    Parallel::gather(0, vtk_cell_conns);
  }

  std::vector<int> vtk_boundary_cell_conns;
  std::vector<unsigned int> vtk_boundary_cell_ids(el);
  std::vector<unsigned int> vtk_boundary_cell_regions;
  std::vector<unsigned short int> vtk_boundary_cell_sides(sl);
  {
    for(unsigned int n=0; n<il.size(); ++n)
    {
      const Elem * elem = mesh.elem(el[n]);
      AutoPtr<Elem> boundary_elem =  elem->build_side(sl[n]);

      vtkIdType celltype = elem_type_vtk(boundary_elem.get());

      vtk_boundary_cell_regions.push_back(elem->subdomain_id());
      // get the connectivity for this element
      std::vector<unsigned int> conn;
      boundary_elem->connectivity(0,VTK,conn);

      vtk_boundary_cell_conns.push_back(celltype);
      vtk_boundary_cell_conns.push_back(conn.size());
      for(unsigned int i=0;i<conn.size();++i)
        vtk_boundary_cell_conns.push_back(conn[i]);
    }

    Parallel::gather(0, vtk_boundary_cell_conns);
    Parallel::gather(0, vtk_boundary_cell_ids);
    Parallel::gather(0, vtk_boundary_cell_regions);
    Parallel::gather(0, vtk_boundary_cell_sides);
  }

  if(Genius::processor_id() == 0)
  {
    grid->Allocate(vtk_cell_ids.size()+vtk_boundary_cell_ids.size());

    // reorder cell by its id
    std::map<unsigned int, std::pair<vtkIdType, vtkIdList *> > reorder_vtk_cells;
    {
      unsigned int cnt = 0;
      for(unsigned int n=0; n<vtk_cell_ids.size(); ++n)
      {
        vtkIdType celltype = vtk_cell_conns[cnt++];
        unsigned int region = vtk_cell_regions[n];
        vtkIdList *pts = vtkIdList::New();
        unsigned int n_nodes = vtk_cell_conns[cnt++];
        pts->SetNumberOfIds(n_nodes);
        for(unsigned int i=0;i<n_nodes;++i)
        {
          unsigned int node = static_cast<unsigned int>(vtk_cell_conns[cnt++]);
          unsigned int node_mapped = region_node_id_map.find(std::make_pair(region, node))->second;
          pts->SetId(i, node_mapped);
        }
        reorder_vtk_cells.insert( std::make_pair(vtk_cell_ids[n], std::make_pair(celltype, pts)));
      }
    }

    // insert reordered cell to vtkUnstructuredGrid
    std::map<unsigned int, std::pair<vtkIdType, vtkIdList *> >::const_iterator vtk_cells_it = reorder_vtk_cells.begin();
    for(; vtk_cells_it != reorder_vtk_cells.end(); ++vtk_cells_it)
    {
      vtkIdType celltype = vtk_cells_it->second.first;
      vtkIdList *pts = vtk_cells_it->second.second;
      grid->InsertNextCell(celltype, pts);
      pts->Delete();
    }


    typedef std::pair<unsigned int, unsigned short int> boundary_elem_key;
    std::map< boundary_elem_key, std::pair<vtkIdType, vtkIdList *> > reorder_vtk_boundary_cells;
    {
      unsigned int cnt = 0;
      for(unsigned int n=0; n<vtk_boundary_cell_ids.size(); ++n)
      {
        boundary_elem_key key = std::make_pair(vtk_boundary_cell_ids[n], vtk_boundary_cell_sides[n]);
        unsigned int region = vtk_boundary_cell_regions[n];
        vtkIdType celltype = vtk_boundary_cell_conns[cnt++];
        vtkIdList *pts = vtkIdList::New();
        unsigned int n_nodes = vtk_boundary_cell_conns[cnt++];
        pts->SetNumberOfIds(n_nodes);
        for(unsigned int i=0;i<n_nodes;++i)
        {
          unsigned int node = static_cast<unsigned int>(vtk_boundary_cell_conns[cnt++]);
          unsigned int node_mapped = region_node_id_map.find(std::make_pair(region, node))->second;
          pts->SetId(i, node_mapped);
        }
        reorder_vtk_boundary_cells.insert( std::make_pair(key, std::make_pair(celltype, pts)));
      }
    }

    // insert reordered cell to vtkUnstructuredGrid
    std::map< boundary_elem_key, std::pair<vtkIdType, vtkIdList *> >::const_iterator vtk_boundary_cells_it;
    for(vtk_boundary_cells_it=reorder_vtk_boundary_cells.begin();
        vtk_boundary_cells_it != reorder_vtk_boundary_cells.end(); ++vtk_boundary_cells_it)
    {
      vtkIdType celltype = vtk_boundary_cells_it->second.first;
      vtkIdList *pts = vtk_boundary_cells_it->second.second;
      grid->InsertNextCell(celltype, pts);
      pts->Delete();
    }

  }


}






void DDMMonitorHook::meshinfo_to_vtk()
{
  const SimulationSystem & system = _solver.get_system();
  const MeshBase & mesh = system.mesh();
  
  //write cell based region and partition info to vtk

  // collect info of mesh elements
  std::vector<int> elem_ids;
  std::vector<int> elem_region;
  std::vector<int> elem_boundary;
  std::vector<int> elem_partition;
  {
    MeshBase::const_element_iterator       it  = mesh.active_this_pid_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_this_pid_elements_end();
    for (int i=0; it != end; ++it, ++i)
    {
      elem_ids.push_back((*it)->id());
      elem_region.push_back((*it)->subdomain_id());
      elem_boundary.push_back(0);
      elem_partition.push_back((*it)->processor_id());
    }
  }
  Parallel::gather(0, elem_ids);
  Parallel::gather(0, elem_region);
  Parallel::gather(0, elem_boundary);
  Parallel::gather(0, elem_partition);


  // collect info of boundary elements
  std::vector<unsigned int> boundary_elem_ids(el);
  std::vector<unsigned short int> boundary_elem_sides(sl);
  std::vector<int> boundary_elem_region;
  std::vector<int> boundary_elem_boundary;
  std::vector<int> boundary_elem_partition;
  {
    for(unsigned int n=0; n<il.size(); ++n)
    {
      const Elem * elem = mesh.elem(el[n]);
      boundary_elem_region.push_back(elem->subdomain_id());
      boundary_elem_boundary.push_back(il[n]);
      boundary_elem_partition.push_back(elem->processor_id());
    }
  }
  Parallel::gather(0, boundary_elem_ids);
  Parallel::gather(0, boundary_elem_sides);
  Parallel::gather(0, boundary_elem_region);
  Parallel::gather(0, boundary_elem_boundary);
  Parallel::gather(0, boundary_elem_partition);

  if(Genius::processor_id() == 0)
  {
    const unsigned int n_elems = elem_ids.size() + boundary_elem_ids.size();

    vtkIntArray *region_info    = vtkIntArray::New();
    vtkIntArray *boundary_info  = vtkIntArray::New();
    vtkIntArray *partition_info  = vtkIntArray::New();

    region_info->SetName("region");
    region_info->SetNumberOfValues(n_elems);

    boundary_info->SetName("boundary");
    boundary_info->SetNumberOfValues(n_elems);

    partition_info->SetName("partition");
    partition_info->SetNumberOfValues(n_elems);

    for (unsigned int n=0; n<elem_ids.size(); ++n)
    {
      int id = elem_ids[n];
      region_info->SetValue(id, elem_region[n]);
      boundary_info->SetValue(id, elem_boundary[n]);
      partition_info->SetValue(id, elem_partition[n]);
    }


    typedef std::pair<unsigned int, unsigned short int> boundary_elem_key;
    std::map<boundary_elem_key, unsigned int> boundary_face_order;
    for(unsigned int n=0; n<boundary_elem_ids.size(); ++n)
    {
      boundary_elem_key key = std::make_pair(boundary_elem_ids[n], boundary_elem_sides[n]);
      boundary_face_order.insert(std::make_pair(key, n));
    }

    std::map<boundary_elem_key, unsigned int>::const_iterator boundary_face_it = boundary_face_order.begin();
    for(int loc=elem_ids.size(); boundary_face_it != boundary_face_order.end(); ++boundary_face_it, ++loc)
    {
      int n = boundary_face_it->second;
      region_info->SetValue(loc, boundary_elem_region[n]);
      boundary_info->SetValue(loc, boundary_elem_boundary[n]);
      partition_info->SetValue(loc, boundary_elem_partition[n]);
    }

    grid->GetCellData()->AddArray(region_info);
    grid->GetCellData()->AddArray(boundary_info);
    grid->GetCellData()->AddArray(partition_info);

    region_info->Delete();
    boundary_info->Delete();
    partition_info->Delete();
  }

}


std::string DDMMonitorHook::export_extra_info()
{
  std::stringstream os;
  const SimulationSystem & system = _solver.get_system();

  os << "<Genius xmlns=\"http://www.cogenda.com/xmlns/Genius/vtu/Genius/1.0\">" << std::endl;

  os << "  <mesh dimension=\"" << system.mesh().mesh_dimension() << "\" />" << std::endl;

  // region info
  {
    os << "  <regions>" << std::endl;
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      os << "    <region id=\"" << r << "\" "
      << "name=\"" << region->name() << "\" "
      << "material=\"" << region->material()  << "\" "
      <<"/>"  <<std::endl;
    }
    os << "  </regions>" << std::endl;
  }

  // boundary info
  {
    os << "  <boundaries>" << std::endl;
    const BoundaryConditionCollector * bcs = system.get_bcs();
    for( unsigned int b=0; b<bcs->n_bcs(); b++)
    {
      const BoundaryCondition * bc = bcs->get_bc(b);
      if( bc->boundary_type() == INTER_CONNECT ) continue;

      os << "    <boundary id=\"" << bcs->get_bd_id_by_bc_index(b) << "\" "
      << "name=\"" << bc->label() << "\" "
      << "bc=\"" << bc->bc_type_name() << "\" "
      << "/>"  << std::endl;
    }
    os << "  </boundaries>" << std::endl;
  }

  // time info
  {
    if (SolverSpecify::Type==SolverSpecify::TRANSIENT)
    {
      os << "  <time value=\"" << SolverSpecify::clock/PhysicalUnit::s << "\" />"
      << std::endl;
    }
  }

  // contact info
  {
    os << "  <contacts>" << std::endl;

    if (SolverSpecify::Type==SolverSpecify::DCSWEEP || SolverSpecify::Type==SolverSpecify::TRANSIENT )
    {
      {
        const BoundaryConditionCollector * bcs = system.get_bcs();

        for(unsigned int i=0; i<bcs->n_bcs(); i++)
        {
          const BoundaryCondition * bc = bcs->get_bc(i);
          // skip bc which is not electrode
          if( !bc->is_electrode() ) continue;
          //
          std::string bc_label = bc->label();
          if(!bc->electrode_label().empty())
            bc_label = bc->electrode_label();

          os << "    <contact name=\"" << bc_label << "\" "
          << "voltage=\"" << bc->ext_circuit()->Vapp()/PhysicalUnit::V << "\" "
          << "potential=\"" << bc->ext_circuit()->potential()/PhysicalUnit::V << "\" "
          << "current=\"" << bc->ext_circuit()->current()/PhysicalUnit::A << "\" "
          << "/>" << std::endl;
        }
      }
    }
    os << "  </contacts>" << std::endl;
  }

  os << "</Genius>" << std::endl;
  return os.str();
}


void DDMMonitorHook::solution_to_vtk_ddm(void * _f, void * _x, void * _dx, void * _w)
{

  bool temperature = (_solver.solver_type() == SolverSpecify::DDML2 || _solver.solver_type() == SolverSpecify::DDML2MIXA);
  bool quantum     = (_solver.solver_type() == SolverSpecify::DENSITY_GRADIENT);

  const SimulationSystem & system = _solver.get_system();
  Vec f  = Vec(_f);  // previous function residual
  Vec x  = Vec(_x);  // previous iterate value
  Vec dx = Vec(_dx); // new search direction and length
  Vec w  = Vec(_w);  // current candidate iterate

  PetscScalar    *ff;
  PetscScalar    *xx;
  PetscScalar    *dxx;
  PetscScalar    *ww;

  VecGetArray(f, &ff);
  VecGetArray(x, &xx);
  VecGetArray(dx, &dxx);
  VecGetArray(w, &ww);

  unsigned int nr = system.n_regions();
  std::vector< std::vector<unsigned int> > region_node_ids(nr);

  // gather data from all the processor
  std::vector< std::vector<float> > Na(nr), Nd(nr), net_doping(nr);
  std::vector< std::vector<float> > fpsi(nr), fn(nr), fp(nr), fT(nr);    // residual at this iteration
  std::vector< std::vector<float> > psi(nr), n(nr), p(nr), T(nr);    // solution at this iteration
  std::vector< std::vector<float> > dpsi(nr), dn(nr), dp(nr), dT(nr); // newton update
  std::vector< std::vector<float> > psi_(nr), n_(nr), p_(nr), T_(nr); // next solution
  std::vector< std::vector<float> > R(nr), G(nr); // recombination and generation
  std::vector< std::vector<float> > dndt1(nr), dpdt1(nr), dndt2(nr), dpdt2(nr); // for time derivative

  std::vector< std::vector<float> > fqc(nr), fqv(nr);    // residual at this iteration
  std::vector< std::vector<float> > Eqc(nr), Eqv(nr);    // solution at this iteration
  std::vector< std::vector<float> > dEqc(nr), dEqv(nr);  // newton update

  const double concentration_scale = pow(cm, -3);

  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);
    if(_solver.solver_type() == SolverSpecify::EBML3 && region->get_advanced_model()->enable_Tl()) temperature=true;

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = *it;
      const FVM_NodeData * node_data = fvm_node->node_data();
      unsigned int id = fvm_node->root_node()->id();

      region_node_ids[r].push_back(id);

      switch( system.region(fvm_node->subdomain_id())->type() )
      {
          case SemiconductorRegion:
          {

            Na[r].push_back( node_data->Total_Na()/concentration_scale );
            Nd[r].push_back( node_data->Total_Nd()/concentration_scale );
            net_doping[r].push_back( node_data->Net_doping()/concentration_scale );

            fpsi[r].push_back( ff[fvm_node->local_offset()+0]/C );
            fn[r].push_back  ( ff[fvm_node->local_offset()+1]/A );
            fp[r].push_back  ( ff[fvm_node->local_offset()+2]/A );

            psi[r].push_back( xx[fvm_node->local_offset()+0]/V );
            n[r].push_back  ( xx[fvm_node->local_offset()+1]/concentration_scale );
            p[r].push_back  ( xx[fvm_node->local_offset()+2]/concentration_scale );

            dpsi[r].push_back( -dxx[fvm_node->local_offset()+0]/V );
            dn[r].push_back  ( -dxx[fvm_node->local_offset()+1]/concentration_scale );
            dp[r].push_back  ( -dxx[fvm_node->local_offset()+2]/concentration_scale );

            psi_[r].push_back( ww[fvm_node->local_offset()+0]/V );
            n_[r].push_back  ( ww[fvm_node->local_offset()+1]/concentration_scale );
            p_[r].push_back  ( ww[fvm_node->local_offset()+2]/concentration_scale );

            if(temperature)
            {
              fT[r].push_back( ff[fvm_node->local_offset()+3]/W );
              T[r].push_back( xx[fvm_node->local_offset()+3]/K );
              dT[r].push_back( -dxx[fvm_node->local_offset()+3]/K );
              T_[r].push_back( ww[fvm_node->local_offset()+3]/K );
            }

            if(quantum)
            {
              unsigned int eqc_index = region->dg_variable_offset(EQC);
              unsigned int eqv_index = region->dg_variable_offset(EQV);
              if(eqc_index!=invalid_uint)
              {
                fqc[r].push_back( ff[fvm_node->local_offset()+eqc_index]/eV );
                Eqc[r].push_back( xx[fvm_node->local_offset()+eqc_index]/eV );
                dEqc[r].push_back( -dxx[fvm_node->local_offset()+eqc_index]/eV );
              }
              else
              {
                fqc[r].push_back( 0.0 );
                Eqc[r].push_back( 0.0 );
                dEqc[r].push_back( 0.0 );
              }

              if(eqv_index!=invalid_uint)
              {
                fqv[r].push_back( ff[fvm_node->local_offset()+eqv_index]/eV );
                Eqv[r].push_back( xx[fvm_node->local_offset()+eqv_index]/eV );
                dEqv[r].push_back( -dxx[fvm_node->local_offset()+eqv_index]/eV );
              }
              else
              {
                fqv[r].push_back( 0.0 );
                Eqv[r].push_back( 0.0 );
                dEqv[r].push_back( 0.0 );
              }
            }

            if(SolverSpecify::TimeDependent == true)
            {
              dndt1[r].push_back( (xx[fvm_node->local_offset()+1] - node_data->n())/SolverSpecify::dt*fvm_node->volume() / (concentration_scale/s) );
              dpdt1[r].push_back( (xx[fvm_node->local_offset()+2] - node_data->p())/SolverSpecify::dt*fvm_node->volume() / (concentration_scale/s) );

              if(SolverSpecify::BDF2_LowerOrder==false)
              {
                const double r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
                dndt2[r].push_back( ((2-r)/(1-r)*xx[fvm_node->local_offset()+1] - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                                 / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume() / (concentration_scale/s) );
                dpdt2[r].push_back( ((2-r)/(1-r)*xx[fvm_node->local_offset()+2] - 1.0/(r*(1-r))*node_data->p() + (1-r)/r*node_data->p_last())
                                 / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume() / (concentration_scale/s) );
              }
              else
              {
                dndt2[r].push_back( (xx[fvm_node->local_offset()+1] - node_data->n())/SolverSpecify::dt*fvm_node->volume() / (concentration_scale/s) );
                dpdt2[r].push_back( (xx[fvm_node->local_offset()+2] - node_data->p())/SolverSpecify::dt*fvm_node->volume() / (concentration_scale/s) );
              }
            }

            const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *>(system.region(fvm_node->subdomain_id()));
            semiconductor_region->material()->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);      // map this node and its data to material database
            R[r].push_back( semiconductor_region->material()->band->Recomb(xx[fvm_node->local_offset()+2], xx[fvm_node->local_offset()+1], node_data->T()) / (concentration_scale/s));
            G[r].push_back(node_data->Field_G()/ (concentration_scale/s));
            break;
          }
          case InsulatorRegion :
          case ElectrodeRegion :
          case MetalRegion :
          {
            Na[r].push_back( 0 );
            Nd[r].push_back( 0 );
            net_doping[r].push_back( 0 );

            fpsi[r].push_back( ff[fvm_node->local_offset()+0]/C );
            fn[r].push_back  ( 0 );
            fp[r].push_back  ( 0 );

            psi[r].push_back( xx[fvm_node->local_offset()+0]/V );
            n[r].push_back  ( 0 );
            p[r].push_back  ( 0 );

            dpsi[r].push_back( -dxx[fvm_node->local_offset()+0]/V );
            dn[r].push_back  ( 0 );
            dp[r].push_back  ( 0 );

            psi_[r].push_back( ww[fvm_node->local_offset()+0]/V );
            n_[r].push_back  ( 0 );
            p_[r].push_back  ( 0 );

            if(temperature)
            {
              fT[r].push_back( ff[fvm_node->local_offset()+1]/W );
              T[r].push_back( xx[fvm_node->local_offset()+1]/K );
              dT[r].push_back( -dxx[fvm_node->local_offset()+1]/K );
              T_[r].push_back( ww[fvm_node->local_offset()+1]/K );
            }

            if(quantum)
            {
              fqc[r].push_back( 0.0 );
              Eqc[r].push_back( 0.0 );
              dEqc[r].push_back( 0.0 );

              fqv[r].push_back( 0.0 );
              Eqv[r].push_back( 0.0 );
              dEqv[r].push_back( 0.0 );
             }

            if(SolverSpecify::TimeDependent == true)
            {
              dndt1[r].push_back( 0 );
              dpdt1[r].push_back( 0 );
              dndt2[r].push_back( 0 );
              dpdt2[r].push_back( 0 );
            }
            R[r].push_back(0);
            G[r].push_back(0);
            break;
          }
          case VacuumRegion:
          break;
          default:
          genius_error(); //we should never reach here
      }

    }
  }


  write_node_scaler_solution(region_node_ids, Na, "Na[cm-3]");
  write_node_scaler_solution(region_node_ids, Nd, "Nd[cm-3]");
  write_node_scaler_solution(region_node_ids, net_doping, "Net Doping[cm-3]");

  write_node_scaler_solution(region_node_ids, fpsi, "fphi[C]");
  write_node_scaler_solution(region_node_ids, fn, "fn[A]");
  write_node_scaler_solution(region_node_ids, fp, "fp[A]");

  write_node_scaler_solution(region_node_ids, psi, "phi[V]");
  write_node_scaler_solution(region_node_ids, n, "n[cm-3]");
  write_node_scaler_solution(region_node_ids, p, "p[cm-3]");

  write_node_scaler_solution(region_node_ids, dpsi, "phi update[V]");
  write_node_scaler_solution(region_node_ids, dn, "n update[cm-3]");
  write_node_scaler_solution(region_node_ids, dp, "p update[cm-3]");

  write_node_scaler_solution(region_node_ids, psi_, "phi next[V]");
  write_node_scaler_solution(region_node_ids, n_, "n next[cm-3]");
  write_node_scaler_solution(region_node_ids, p_, "p next[cm-3]");


  if(temperature)
  {
    write_node_scaler_solution(region_node_ids, fT, "fT[W]");
    write_node_scaler_solution(region_node_ids, T, "T[K]");
    write_node_scaler_solution(region_node_ids, dT, "T update[K]");
    write_node_scaler_solution(region_node_ids, T_, "T next[K]");
  }

  if(quantum)
  {
    write_node_scaler_solution(region_node_ids, fqc, "fqc[eV]");
    write_node_scaler_solution(region_node_ids, Eqc, "Eqc[eV]");
    write_node_scaler_solution(region_node_ids, dEqc, "Eqc update[eV]");

    write_node_scaler_solution(region_node_ids, fqv, "fqv[eV]");
    write_node_scaler_solution(region_node_ids, Eqv, "Eqv[eV]");
    write_node_scaler_solution(region_node_ids, dEqv, "Eqv update[eV]");
  }

  if(SolverSpecify::TimeDependent == true)
  {
    write_node_scaler_solution(region_node_ids, dndt1, "dndt BDF1[A]");
    write_node_scaler_solution(region_node_ids, dndt2, "dndt BDF2[A]");
    write_node_scaler_solution(region_node_ids, dpdt1, "dpdt BDF1[A]");
    write_node_scaler_solution(region_node_ids, dpdt2, "dpdt BDF2[A]");
  }
  write_node_scaler_solution(region_node_ids, R, "recombination[cm-3/s]");
  write_node_scaler_solution(region_node_ids, G, "generation[cm-3/s]");

  VecRestoreArray(f, &ff);
  VecRestoreArray(x, &xx);
  VecRestoreArray(dx, &dxx);
  VecRestoreArray(w, &ww);
}



void DDMMonitorHook::solution_to_vtk_poisson(void * _f, void * _x, void * _dx, void * _w)
{

  const SimulationSystem & system = _solver.get_system();

  Vec f  = Vec(_f);  // previous function residual
  Vec x  = Vec(_x);  // previous iterate value
  Vec dx = Vec(_dx); // new search direction and length
  Vec w  = Vec(_w);  // current candidate iterate

  PetscScalar    *ff;
  PetscScalar    *xx;
  PetscScalar    *dxx;
  PetscScalar    *ww;

  VecGetArray(f, &ff);
  VecGetArray(x, &xx);
  VecGetArray(dx, &dxx);
  VecGetArray(w, &ww);

  unsigned int nr = system.n_regions();
  std::vector< std::vector<unsigned int> > region_node_ids(nr);

  // gather data from all the processor
  std::vector< std::vector<float> > Na(nr), Nd(nr), net_doping(nr);
  std::vector< std::vector<float> > fpsi(nr);   // residual at this iteration
  std::vector< std::vector<float> > psi(nr);    // solution at this iteration
  std::vector< std::vector<float> > dpsi(nr);   // newton update

  const double concentration_scale = pow(cm, -3);

  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = system.region(r);

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = *it;
      unsigned int id = fvm_node->root_node()->id();
      region_node_ids[r].push_back(id);

      switch( system.region(fvm_node->subdomain_id())->type() )
      {
          case SemiconductorRegion:
          {
            Na[r].push_back( fvm_node->node_data()->Total_Na()/concentration_scale );
            Nd[r].push_back( fvm_node->node_data()->Total_Nd()/concentration_scale );
            net_doping[r].push_back( fvm_node->node_data()->Net_doping()/concentration_scale );

            fpsi[r].push_back( ff[fvm_node->local_offset()] );
            psi[r].push_back( xx[fvm_node->local_offset()]/V );
            dpsi[r].push_back( -dxx[fvm_node->local_offset()]/V );
            break;
          }
          case InsulatorRegion :
          case ElectrodeRegion :
          case MetalRegion :
          {
            Na[r].push_back( 0 );
            Nd[r].push_back( 0 );
            net_doping[r].push_back( 0 );

            fpsi[r].push_back( ff[fvm_node->local_offset()+0] );
            psi[r].push_back( xx[fvm_node->local_offset()+0]/V );
            dpsi[r].push_back( -dxx[fvm_node->local_offset()+0]/V );
            break;
          }
          case VacuumRegion:
          break;
          default:
          genius_error(); //we should never reach here
      }

    }
  }

  write_node_scaler_solution(region_node_ids, Na, "Na[cm-3]");
  write_node_scaler_solution(region_node_ids, Nd, "Nd[cm-3]");
  write_node_scaler_solution(region_node_ids, net_doping, "Net Doping[cm-3]");

  write_node_scaler_solution(region_node_ids, fpsi, "fphi[C]");
  write_node_scaler_solution(region_node_ids, psi, "phi[V]");
  write_node_scaler_solution(region_node_ids, dpsi, "phi update[V]");

  VecRestoreArray(f, &ff);
  VecRestoreArray(x, &xx);
  VecRestoreArray(dx, &dxx);
  VecRestoreArray(w, &ww);
}


void DDMMonitorHook::solution_to_vtk_tid_drift(void * _f, void * _x, void * _dx, void * _w)
{
  const SimulationSystem & system = _solver.get_system();
  Vec f  = Vec(_f);  // previous function residual
  Vec x  = Vec(_x);  // previous iterate value
  Vec dx = Vec(_dx); // new search direction and length
  Vec w  = Vec(_w);  // current candidate iterate

  PetscScalar    *ff;
  PetscScalar    *xx;
  PetscScalar    *dxx;
  PetscScalar    *ww;

  VecGetArray(f, &ff);
  VecGetArray(x, &xx);
  VecGetArray(dx, &dxx);
  VecGetArray(w, &ww);

  unsigned int nr = system.n_regions();
  std::vector< std::vector<unsigned int> > region_node_ids(nr);

  // gather data from all the processor
  std::vector< std::vector<float> > Na(nr), Nd(nr), net_doping(nr);
  std::vector< std::vector<float> > fpsi(nr), fn(nr), fp(nr), fpt(nr);    // residual at this iteration
  std::vector< std::vector<float> > psi(nr),  n(nr), p(nr), pt(nr);    // solution at this iteration
  std::vector< std::vector<float> > dpsi(nr), dn(nr), dp(nr), dpt(nr); // newton update
  std::vector< std::vector<float> > psi_(nr), n_(nr), p_(nr), pt_(nr); // next solution
  std::vector< std::vector<float> > dndt1(nr), dpdt1(nr), dptdt1(nr); // for time derivative


  const double concentration_scale = pow(cm, -3);

  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = *it;
      const FVM_NodeData * node_data = fvm_node->node_data();
      unsigned int id = fvm_node->root_node()->id();

      region_node_ids[r].push_back(id);

      switch( system.region(fvm_node->subdomain_id())->type() )
      {
        case InsulatorRegion :
        {
          Na[r].push_back( 0 );
          Nd[r].push_back( 0 );
          net_doping[r].push_back( 0 );

          fpsi[r].push_back( ff[fvm_node->local_offset()+0]/C );
          fn[r].push_back  ( ff[fvm_node->local_offset()+1]/A );
          fp[r].push_back  ( ff[fvm_node->local_offset()+2]/A );
          fpt[r].push_back  ( ff[fvm_node->local_offset()+3]/A );

          psi[r].push_back( xx[fvm_node->local_offset()+0]/V );
          n[r].push_back  ( xx[fvm_node->local_offset()+1]/concentration_scale );
          p[r].push_back  ( xx[fvm_node->local_offset()+2]/concentration_scale );
          pt[r].push_back  ( xx[fvm_node->local_offset()+3]/concentration_scale );

          dpsi[r].push_back( -dxx[fvm_node->local_offset()+0]/V );
          dn[r].push_back  ( -dxx[fvm_node->local_offset()+1]/concentration_scale );
          dp[r].push_back  ( -dxx[fvm_node->local_offset()+2]/concentration_scale );
          dpt[r].push_back  ( -dxx[fvm_node->local_offset()+3]/concentration_scale );

          psi_[r].push_back( ww[fvm_node->local_offset()+0]/V );
          n_[r].push_back  ( ww[fvm_node->local_offset()+1]/concentration_scale );
          p_[r].push_back  ( ww[fvm_node->local_offset()+2]/concentration_scale );
          pt_[r].push_back  ( ww[fvm_node->local_offset()+3]/concentration_scale );


          if(SolverSpecify::TimeDependent == true)
          {
            dndt1[r].push_back( (xx[fvm_node->local_offset()+1] - node_data->n())/SolverSpecify::dt*fvm_node->volume() / (concentration_scale/s) );
            dpdt1[r].push_back( (xx[fvm_node->local_offset()+2] - node_data->p())/SolverSpecify::dt*fvm_node->volume() / (concentration_scale/s) );
            dptdt1[r].push_back( (xx[fvm_node->local_offset()+3] - node_data->trap_a())/SolverSpecify::dt*fvm_node->volume() / (concentration_scale/s) );
          }

          break;
        }
        case SemiconductorRegion :
        {
          Na[r].push_back( node_data->Total_Na()/concentration_scale );
          Nd[r].push_back( node_data->Total_Nd()/concentration_scale );
          net_doping[r].push_back( node_data->Net_doping()/concentration_scale );
          
          fpsi[r].push_back( 0 );
          fn[r].push_back  ( 0 );
          fp[r].push_back  ( 0 );
          fpt[r].push_back  ( 0 );

          psi[r].push_back( node_data->psi() );
          n[r].push_back  ( 0 );
          p[r].push_back  ( 0 );
          pt[r].push_back  ( 0 );

          dpsi[r].push_back( 0 );
          dn[r].push_back  ( 0 );
          dp[r].push_back  ( 0 );
          dpt[r].push_back  ( 0 );

          psi_[r].push_back( node_data->psi() );
          n_[r].push_back  ( 0 );
          p_[r].push_back  ( 0 );
          pt_[r].push_back  ( 0 );

          if(SolverSpecify::TimeDependent == true)
          {
            dndt1[r].push_back( 0 );
            dpdt1[r].push_back( 0 );
            dptdt1[r].push_back( 0 );
          }
          break;
        }
        case ElectrodeRegion :
        case MetalRegion :
        case VacuumRegion:
        {
          Na[r].push_back( 0 );
          Nd[r].push_back( 0 );
          net_doping[r].push_back( 0 );
          
          fpsi[r].push_back( 0 );
          fn[r].push_back  ( 0 );
          fp[r].push_back  ( 0 );
          fpt[r].push_back  ( 0 );

          psi[r].push_back( node_data->psi() );
          n[r].push_back  ( 0 );
          p[r].push_back  ( 0 );
          pt[r].push_back  ( 0 );

          dpsi[r].push_back( 0 );
          dn[r].push_back  ( 0 );
          dp[r].push_back  ( 0 );
          dpt[r].push_back  ( 0 );

          psi_[r].push_back( node_data->psi() );
          n_[r].push_back  ( 0 );
          p_[r].push_back  ( 0 );
          pt_[r].push_back  ( 0 );

          if(SolverSpecify::TimeDependent == true)
          {
            dndt1[r].push_back( 0 );
            dpdt1[r].push_back( 0 );
            dptdt1[r].push_back( 0 );
          }
          break;
        }
        default:
          genius_error(); //we should never reach here
      }

    }
  }

  write_node_scaler_solution(region_node_ids, Na, "Na[cm-3]");
  write_node_scaler_solution(region_node_ids, Nd, "Nd[cm-3]");
  write_node_scaler_solution(region_node_ids, net_doping, "Net Doping[cm-3]");
  
  write_node_scaler_solution(region_node_ids, fpsi, "fphi[C]");

  write_node_scaler_solution(region_node_ids, fn, "fn[A]");
  write_node_scaler_solution(region_node_ids, fp, "fp[A]");
  write_node_scaler_solution(region_node_ids, fpt, "fpt[A]");

  write_node_scaler_solution(region_node_ids, psi, "phi[V]");
  write_node_scaler_solution(region_node_ids, n, "n[cm-3]");
  write_node_scaler_solution(region_node_ids, p, "p[cm-3]");
  write_node_scaler_solution(region_node_ids, pt, "pt[cm-3]");

  write_node_scaler_solution(region_node_ids, dpsi, "phi update[V]");
  write_node_scaler_solution(region_node_ids, dn, "n update[cm-3]");
  write_node_scaler_solution(region_node_ids, dp, "p update[cm-3]");
  write_node_scaler_solution(region_node_ids, dpt, "pt update[cm-3]");

  write_node_scaler_solution(region_node_ids, psi_, "phi next[V]");
  write_node_scaler_solution(region_node_ids, n_, "n next[cm-3]");
  write_node_scaler_solution(region_node_ids, p_, "p next[cm-3]");
  write_node_scaler_solution(region_node_ids, pt_, "pt next[cm-3]");


  if(SolverSpecify::TimeDependent == true)
  {
    write_node_scaler_solution(region_node_ids, dndt1, "dndt BDF1[A]");
    write_node_scaler_solution(region_node_ids, dpdt1, "dpdt BDF1[A]");
    write_node_scaler_solution(region_node_ids, dptdt1, "dptdt BDF1[A]");
  }

  
  VecRestoreArray(f, &ff);
  VecRestoreArray(x, &xx);
  VecRestoreArray(dx, &dxx);
  VecRestoreArray(w, &ww);

}


void DDMMonitorHook::solution_to_vtk_ric(void * _f, void * _x, void * _dx, void * _w)
{
  const SimulationSystem & system = _solver.get_system();
  Vec f  = Vec(_f);  // previous function residual
  Vec x  = Vec(_x);  // previous iterate value
  Vec dx = Vec(_dx); // new search direction and length
  Vec w  = Vec(_w);  // current candidate iterate

  PetscScalar    *ff;
  PetscScalar    *xx;
  PetscScalar    *dxx;
  PetscScalar    *ww;

  VecGetArray(f, &ff);
  VecGetArray(x, &xx);
  VecGetArray(dx, &dxx);
  VecGetArray(w, &ww);

  unsigned int nr = system.n_regions();
  std::vector< std::vector<unsigned int> > region_node_ids(nr);

  // gather data from all the processor
  std::vector< std::vector<float> >  fpsi(nr), fn(nr);    // residual at this iteration
  std::vector< std::vector<float> >  psi(nr), n(nr);    // solution at this iteration
  std::vector< std::vector<float> >  dpsi(nr), dn(nr); // newton update
  std::vector< std::vector<float> >  psi_(nr), n_(nr); // next solution
  std::vector< std::vector<float> >  dndt1(nr), dndt2(nr); // for time derivative


  const double concentration_scale = pow(cm, -3);

  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = *it;
      unsigned int id = fvm_node->root_node()->id();
      const FVM_NodeData * node_data = fvm_node->node_data();

      region_node_ids[r].push_back(id);

      switch( system.region(fvm_node->subdomain_id())->type() )
      {
        case InsulatorRegion :
        {

          fpsi[r].push_back( ff[fvm_node->local_offset()+0]/C );
          fn[r].push_back  ( ff[fvm_node->local_offset()+1]/A );

          psi[r].push_back( xx[fvm_node->local_offset()+0]/V );
          n[r].push_back  ( xx[fvm_node->local_offset()+1]/concentration_scale );

          dpsi[r].push_back( -dxx[fvm_node->local_offset()+0]/V );
          dn[r].push_back  ( -dxx[fvm_node->local_offset()+1]/concentration_scale );

          psi_[r].push_back( ww[fvm_node->local_offset()+0]/V );
          n_[r].push_back  ( ww[fvm_node->local_offset()+1]/concentration_scale );


          if(SolverSpecify::TimeDependent == true)
          {
            dndt1[r].push_back( (xx[fvm_node->local_offset()+1] - node_data->n())/SolverSpecify::dt*fvm_node->volume() / (concentration_scale/s) );

            if(SolverSpecify::BDF2_LowerOrder==false)
            {
              const double r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
              dndt2[r].push_back( ((2-r)/(1-r)*xx[fvm_node->local_offset()+1] - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                                 / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume() / (concentration_scale/s) );
            }
            else
            {
              dndt2[r].push_back( (xx[fvm_node->local_offset()+1] - node_data->n())/SolverSpecify::dt*fvm_node->volume() / (concentration_scale/s) );
            }
          }

          break;
        }
        case ElectrodeRegion :
        case MetalRegion :
        case VacuumRegion:
        {
          fpsi[r].push_back( ff[fvm_node->local_offset()+0]/W );
          fn[r].push_back  ( 0 );

          psi[r].push_back( xx[fvm_node->local_offset()+0]/V );
          n[r].push_back  ( 0 );

          dpsi[r].push_back( -dxx[fvm_node->local_offset()+0]/V );
          dn[r].push_back  ( 0 );

          psi_[r].push_back( ww[fvm_node->local_offset()+0]/V );
          n_[r].push_back  ( 0 );

          if(SolverSpecify::TimeDependent == true)
          {
            dndt1[r].push_back( 0 );
            dndt2[r].push_back( 0 );
          }
          break;
        }
        default:
          genius_error(); //we should never reach here
      }

    }
  }

  write_node_scaler_solution(region_node_ids, fpsi, "fphi[C]");
  write_node_scaler_solution(region_node_ids, fn, "fn[A]");

  write_node_scaler_solution(region_node_ids, psi, "phi[V]");
  write_node_scaler_solution(region_node_ids, n, "n[cm-3]");

  write_node_scaler_solution(region_node_ids, dpsi, "phi update[V]");
  write_node_scaler_solution(region_node_ids, dn, "n update[cm-3]");

  write_node_scaler_solution(region_node_ids, psi_, "phi next[V]");
  write_node_scaler_solution(region_node_ids, n_, "n next[cm-3]");


  if(SolverSpecify::TimeDependent == true)
  {
    write_node_scaler_solution(region_node_ids, dndt1, "dndt BDF1[A]");
    write_node_scaler_solution(region_node_ids, dndt2, "dndt BDF2[A]");
  }

  VecRestoreArray(f, &ff);
  VecRestoreArray(x, &xx);
  VecRestoreArray(dx, &dxx);
  VecRestoreArray(w, &ww);
}


void DDMMonitorHook::solution_to_vtk_dictat(void * _f, void * _x, void * _dx, void * _w)
{
  const SimulationSystem & system = _solver.get_system();
  Vec f  = Vec(_f);  // previous function residual
  Vec x  = Vec(_x);  // previous iterate value
  Vec dx = Vec(_dx); // new search direction and length
  Vec w  = Vec(_w);  // current candidate iterate

  PetscScalar    *ff;
  PetscScalar    *xx;
  PetscScalar    *dxx;
  PetscScalar    *ww;

  VecGetArray(f, &ff);
  VecGetArray(x, &xx);
  VecGetArray(dx, &dxx);
  VecGetArray(w, &ww);

  unsigned int nr = system.n_regions();
  std::vector< std::vector<unsigned int> > region_node_ids(nr);

  // gather data from all the processor
  std::vector< std::vector<float> >  fpsi(nr);    // residual at this iteration
  std::vector< std::vector<float> >  psi(nr);    // solution at this iteration
  std::vector< std::vector<float> >  dpsi(nr); // newton update
  std::vector< std::vector<float> >  psi_(nr); // next solution



  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = *it;
      unsigned int id = fvm_node->root_node()->id();
      const FVM_NodeData * node_data = fvm_node->node_data();

      region_node_ids[r].push_back(id);

      switch( system.region(fvm_node->subdomain_id())->type() )
      {
        case InsulatorRegion :
        {

          fpsi[r].push_back( ff[fvm_node->local_offset()+0]/C );

          psi[r].push_back( xx[fvm_node->local_offset()+0]/V );

          dpsi[r].push_back( -dxx[fvm_node->local_offset()+0]/V );

          psi_[r].push_back( ww[fvm_node->local_offset()+0]/V );


          break;
        }
        case ElectrodeRegion :
        case MetalRegion :
        case VacuumRegion:
        {
          fpsi[r].push_back( ff[fvm_node->local_offset()+0]/W );

          psi[r].push_back( xx[fvm_node->local_offset()+0]/V );

          dpsi[r].push_back( -dxx[fvm_node->local_offset()+0]/V );

          psi_[r].push_back( ww[fvm_node->local_offset()+0]/V );

          break;
        }
        default:
          genius_error(); //we should never reach here
      }

    }
  }

  write_node_scaler_solution(region_node_ids, fpsi, "fphi[C]");

  write_node_scaler_solution(region_node_ids, psi, "phi[V]");

  write_node_scaler_solution(region_node_ids, dpsi, "phi update[V]");

  write_node_scaler_solution(region_node_ids, psi_, "phi next[V]");

  VecRestoreArray(f, &ff);
  VecRestoreArray(x, &xx);
  VecRestoreArray(dx, &dxx);
  VecRestoreArray(w, &ww);
}




void DDMMonitorHook::write_node_scaler_solution(const std::vector< std::vector<unsigned int> >& region_order,
                                       const std::vector< std::vector<float> > &region_sol,
                                       const std::string & sol_name)
{
  // this should run on parallel for all the processor
  std::vector<float> solution;

  for(unsigned int r=0; r<region_order.size(); ++r)
  {
    const std::vector<unsigned int> & order = region_order[r];
    const std::vector<float> & sol = region_sol[r];

    std::map<unsigned int, float> region_variable;
    for(unsigned int i=0; i<order.size(); ++i)
    {
      region_variable.insert(std::make_pair(order[i], sol[i]));
    }

    Parallel::allgather(region_variable);

    std::map<unsigned int, float>::const_iterator it = region_variable.begin();
    for( ; it != region_variable.end(); it++)
      solution.push_back(it->second);
  }

  if (Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();
    vtk_sol_array->SetName(sol_name.c_str());
    vtk_sol_array->SetNumberOfValues(solution.size());

    // save data into vtk data array
    for(unsigned int n=0; n<solution.size(); ++n)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array->InsertValue(n, solution[n]);
    }

    grid->GetPointData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }
}






void DDMMonitorHook::write_cell_scaler_solution(const std::vector<unsigned int> & order, std::vector<float> &sol,
                                       const std::string & sol_name)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if (Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();
    vtk_sol_array->SetName(sol_name.c_str());
    vtk_sol_array->SetNumberOfValues(order.size());

    std::vector<float> sol_re_order(order.size());
    for(unsigned int n=0; n<order.size(); ++n)
      sol_re_order[order[n]] = sol[n];

    // save data into vtk data array
    for(unsigned int n=0; n<order.size(); ++n)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array->InsertValue(n, sol_re_order[n]);
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
    return new DDMMonitorHook ( solver, name, fun_data );
  }

}

#endif

