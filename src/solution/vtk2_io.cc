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


// C++ includes
#include <fstream>
#include <sstream>

// Local includes
#include "vtk2_io.h"

#include "elem.h"
#include "mesh_base.h"
#include "boundary_info.h"
#include "material_define.h"
#include "parallel.h"
#include "boundary_condition_collector.h"
#include "field_source.h"
#include "ngspice_interface.h"
#include "spice_ckt.h"


#ifdef HAVE_VTK

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkConfigure.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"

#endif //HAVE_VTK


using PhysicalUnit::s;
using PhysicalUnit::eV;
using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::K;
using PhysicalUnit::J;
using PhysicalUnit::kg;


#ifdef HAVE_VTK
class VTK2IO::XMLUnstructuredGridWriter : public vtkXMLUnstructuredGridWriter
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


void VTK2IO::nodes_to_vtk(const MeshBase& mesh)
{
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();

  unsigned int region_n_nodes = 0;
  _region_node_id_map.clear();

  std::vector<unsigned int> system_nodes;
  for( unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);

    std::vector<unsigned int> nodes;
    region->region_node(nodes);

    for(unsigned int n=0; n<nodes.size(); ++n)
    {
      unsigned int node_id = nodes[n];
      _region_node_id_map.insert( std::make_pair(std::make_pair(r,node_id), region_n_nodes+n));
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

    _vtk_grid->SetPoints(points);

    //clean up
    points->Delete();
  }

}


void VTK2IO::cells_to_vtk(const MeshBase& mesh)
{
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
  std::vector<unsigned int> vtk_boundary_cell_ids(_el);
  std::vector<unsigned int> vtk_boundary_cell_regions;
  std::vector<unsigned short int> vtk_boundary_cell_sides(_sl);
  {
    for(unsigned int n=0; n<_il.size(); ++n)
    {
      const Elem * elem = mesh.elem(_el[n]);
      AutoPtr<Elem> boundary_elem =  elem->build_side(_sl[n]);

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
    _vtk_grid->Allocate(vtk_cell_ids.size()+vtk_boundary_cell_ids.size());

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
          unsigned int node_mapped = _region_node_id_map.find(std::make_pair(region, node))->second;
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
      _vtk_grid->InsertNextCell(celltype, pts);
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
          unsigned int node_mapped = _region_node_id_map.find(std::make_pair(region, node))->second;
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
      _vtk_grid->InsertNextCell(celltype, pts);
      pts->Delete();
    }

  }

}


void VTK2IO::meshinfo_to_vtk(const MeshBase& mesh)
{
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
  std::vector<unsigned int> boundary_elem_ids(_el);
  std::vector<unsigned short int> boundary_elem_sides(_sl);
  std::vector<int> boundary_elem_region;
  std::vector<int> boundary_elem_boundary;
  std::vector<int> boundary_elem_partition;
  {
    for(unsigned int n=0; n<_il.size(); ++n)
    {
      const Elem * elem = mesh.elem(_el[n]);
      boundary_elem_region.push_back(elem->subdomain_id());
      boundary_elem_boundary.push_back(_il[n]);
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

    _vtk_grid->GetCellData()->AddArray(region_info);
    _vtk_grid->GetCellData()->AddArray(boundary_info);
    _vtk_grid->GetCellData()->AddArray(partition_info);

    region_info->Delete();
    boundary_info->Delete();
    partition_info->Delete();
  }

}


std::string VTK2IO::export_extra_info()
{
  std::stringstream os;
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();

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
      if(SolverSpecify::Solver == SolverSpecify::DDML1MIXA ||
          SolverSpecify::Solver == SolverSpecify::DDML2MIXA ||
          SolverSpecify::Solver == SolverSpecify::EBML3MIXA
        )
      {
        // in MixedA mode

        const SPICE_CKT * spice_ckt = system.get_circuit();
        const std::map<std::string, unsigned int> & elec_node_map = spice_ckt->get_electrode_info();

        for(std::map<std::string, unsigned int>::const_iterator it = elec_node_map.begin();
            it != elec_node_map.end(); it++)
        {
          const std::string &label = it->first;
          const int &node = it->second;

          if (spice_ckt->is_voltage_node(node))
          {
            os << "    <contact name=\"" << label << "\" "
            << "voltage=\"" << spice_ckt->get_solution(node) << "\" "
            << "/>" << std::endl;
          }
        }
      }
      else
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



void VTK2IO::build_export_solutions()
{
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();

  // check which data should be export

  bool semiconductor_material = false;
  for(unsigned int n=0; n<system.n_regions(); ++n)
  {
    const SimulationRegion * region = system.region(n);
    if(Material::IsSemiconductor(region->material()))
    { semiconductor_material=true; break; }
  }


  bool ebm_solution_data     = false;
  bool dg_solution_data      = false;
  bool optical_generation    = false;
  bool particle_generation   = false;
  bool field_generation      = false;
  bool ddm_ac_data           = false;
  bool optical_complex_field = false;
  bool tid_data              = false;

  const std::vector<SolverSpecify::SolverType> & solve_history = system.solve_history();

  for(unsigned int n=0; n<solve_history.size(); ++n)
  {
    switch  (solve_history[n])
    {
    case SolverSpecify::EBML3      :
    case SolverSpecify::EBML3MIX   :
      ebm_solution_data = true;
      break;
    case SolverSpecify::DENSITY_GRADIENT :
      dg_solution_data = true;
      break;
    case SolverSpecify::EM_FEM_2D  :
    case SolverSpecify::EM_FEM_3D  :
      optical_complex_field = true;
      optical_generation = true;
      break;
    case SolverSpecify::RAY_TRACE  :
      optical_generation = true;
      break;
    case SolverSpecify::DDMAC      :
      ddm_ac_data = true;
      break;
    case SolverSpecify::RIC  :
      field_generation = true;
      break;
    case SolverSpecify::TID_TRAP   :
    case SolverSpecify::TID_DRIFT  :
    case SolverSpecify::TID_FULL   :
    case SolverSpecify::TIDOP      :
      tid_data = true;
      break;
    default : break;
    }
  }

  if(system.get_field_source()->is_light_source_exist())
    optical_generation = true;

  if(system.get_field_source()->is_particle_source_exist())
    particle_generation = true;

  _variables.insert("potential");
  _variables.insert("electron");
  _variables.insert("hole");
  _variables.insert("temperature");
  _variables.insert("Na");
  _variables.insert("Nd");
  _variables.insert("mole_x");
  _variables.insert("mole_y");
  _variables.insert("recombination");
  _variables.insert("recombination_dir");
  _variables.insert("recombination_srh");
  _variables.insert("recombination_auger");
  _variables.insert("impact_ionization");
  _variables.insert("mun");
  _variables.insert("mup");
  _variables.insert("Ec");
  _variables.insert("Ev");
  _variables.insert("qfn");
  _variables.insert("qfp");

  if(ebm_solution_data)
  {
    _variables.insert("elec_temperature");
    _variables.insert("hole_temperature");
  }

  if(dg_solution_data)
  {
    _variables.insert("Eqc");
    _variables.insert("Eqv");
  }

  if(optical_generation)
  {
    _variables.insert("optical_generation");
  }
  if(particle_generation)
  {
    _variables.insert("particle_generation");
  }

  if(field_generation)
  {
    _variables.insert("field_generation");
  }

  if(ddm_ac_data)
  {
    _variables.insert("potential.ac");
    _variables.insert("electron.ac");
    _variables.insert("hole.ac");
  }

  if(tid_data)
  {
    _variables.insert("trap.p");
  }
}


void VTK2IO::solution_to_vtk()
{
  build_export_solutions();

  std::set<std::string>::const_iterator it = _variables.begin();
  for( ;  it != _variables.end(); it++)
    write_node_solution(*it);

}



void VTK2IO::write_node_solution(const std::string & sol_name)
{
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();

  std::string unit;
  DataType data_type = INVALID_DATATYPE;
  for( unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);
    if(region->has_variable(sol_name, POINT_CENTER))
    {
      const SimulationVariable & var = region->get_variable(sol_name, POINT_CENTER);

      if(data_type == INVALID_DATATYPE)
        data_type = var.variable_data_type;
      else
        genius_assert(data_type == var.variable_data_type);

      if(unit.empty())
        unit=var.variable_unit_string;
      else
        genius_assert(unit==var.variable_unit_string);
    }
  }

  std::string sol_unit = '['+unit+']';

  switch(data_type)
  {
    case SCALAR  : write_node_solution_scalar(sol_name, sol_unit);  break;
    case COMPLEX : write_node_solution_complex(sol_name, sol_unit); break;
    case VECTOR  : write_node_solution_vector(sol_name, sol_unit); break;
    case TENSOR  : write_node_solution_tensor(sol_name, sol_unit);  break;
    default: break;
  }

}


void VTK2IO::write_node_solution_scalar(const std::string & sol_name, const std::string & sol_unit)
{
  std::vector<float> sol;

  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();

  for( unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);

    if(region->has_variable(sol_name, POINT_CENTER))
    {
      std::vector<PetscScalar> region_sol;
      region->get_variable_data<PetscScalar>(sol_name, POINT_CENTER, region_sol);
      sol.insert(sol.end(), region_sol.begin(), region_sol.end());
    }
    else
    {
      std::vector<float> region_sol(region->n_node(), 0.0);
      sol.insert(sol.end(), region_sol.begin(), region_sol.end());
    }
  }

  if (Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();
    vtk_sol_array->SetName((sol_name+sol_unit).c_str());
    vtk_sol_array->SetNumberOfValues(sol.size());

    // save data into vtk data array
    for(unsigned int n=0; n<sol.size(); ++n)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array->InsertValue(n, (sol[n]));
    }

    _vtk_grid->GetPointData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }
}



void VTK2IO::write_node_solution_complex(const std::string & sol_name, const std::string & sol_unit)
{
  std::vector< std::complex<float> > sol;

  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();

  for( unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);

    if(region->has_variable(sol_name, POINT_CENTER))
    {
      std::vector< std::complex<PetscScalar> > region_sol;
      region->get_variable_data< std::complex<PetscScalar> >(sol_name, POINT_CENTER, region_sol);
      for(size_t i=0; i<region_sol.size(); i++)
        sol.push_back(std::complex<float>(region_sol[i].real(), region_sol[i].imag()));
    }
    else
    {
      std::vector< std::complex<float> > region_sol(region->n_node(), 0.0);
      for(size_t i=0; i<region_sol.size(); i++)
        sol.push_back(std::complex<float>(region_sol[i].real(), region_sol[i].imag()));
    }
  }

  if (Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array_real = vtkFloatArray::New();
    vtk_sol_array_real->SetName((sol_name+".real"+sol_unit).c_str());
    vtk_sol_array_real->SetNumberOfValues(sol.size());

    vtkFloatArray *vtk_sol_array_imag = vtkFloatArray::New();
    vtk_sol_array_imag->SetName((sol_name+".imag"+sol_unit).c_str());
    vtk_sol_array_imag->SetNumberOfValues(sol.size());

    // save data into vtk data array
    for(unsigned int n=0; n<sol.size(); ++n)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array_real->InsertValue(n, (sol[n].real()));
      vtk_sol_array_imag->InsertValue(n, (sol[n].imag()));
    }

    _vtk_grid->GetPointData()->AddArray(vtk_sol_array_real);
    _vtk_grid->GetPointData()->AddArray(vtk_sol_array_imag);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array_real->Delete();
    vtk_sol_array_imag->Delete();
  }

}


void VTK2IO::write_node_solution_vector(const std::string & sol_name, const std::string & sol_unit)
{
  std::vector< VectorValue<float> > sol;

  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();

  for( unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);

    if(region->has_variable(sol_name, POINT_CENTER))
    {
      std::vector< VectorValue<PetscScalar> > region_sol;
      region->get_variable_data< VectorValue<PetscScalar> >(sol_name, POINT_CENTER, region_sol);
      for(unsigned int n=0; n<region_sol.size(); ++n)
      {
        const VectorValue<PetscScalar> & v = region_sol[n];
        sol.push_back(VectorValue<float>(v[0],v[1],v[2]));
      }
    }
    else
    {
      std::vector< VectorValue<float> > region_sol(region->n_node(), 0.0);
      sol.insert(sol.end(), region_sol.begin(), region_sol.end());
    }
  }


  if (Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();

    vtk_sol_array->SetNumberOfComponents(3);
    vtk_sol_array->SetNumberOfTuples(sol.size());
    vtk_sol_array->SetName((sol_name+sol_unit).c_str());

    // save data into vtk data array
    for(unsigned int n=0; n<sol.size(); ++n)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array->InsertTuple(n, &(sol[n])[0] );
    }

    _vtk_grid->GetPointData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }

}



void VTK2IO::write_node_solution_tensor(const std::string & sol_name, const std::string & sol_unit)
{
  std::vector< TensorValue<float> > sol;

  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();

  for( unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);

    if(region->has_variable(sol_name, POINT_CENTER))
    {
      std::vector< TensorValue<PetscScalar> > region_sol;
      region->get_variable_data< TensorValue<PetscScalar> >(sol_name, POINT_CENTER, region_sol);
      for(unsigned int n=0; n<region_sol.size(); ++n)
      {
        const TensorValue<PetscScalar> & t = region_sol[n];
        sol.push_back(TensorValue<float>(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8]));
      }
    }
    else
    {
      std::vector< TensorValue<float> > region_sol(region->n_node(), 0.0);
      sol.insert(sol.end(), region_sol.begin(), region_sol.end());
    }
  }


  if (Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();

    vtk_sol_array->SetNumberOfComponents(9);
    vtk_sol_array->SetNumberOfTuples(sol.size());
    vtk_sol_array->SetName((sol_name+sol_unit).c_str());

    // save data into vtk data array
    for(unsigned int n=0; n<sol.size(); ++n)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array->InsertTuple(n, &(sol[n])[0] );
    }

    _vtk_grid->GetPointData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }

}



void VTK2IO::write_node_solution_scalar(const std::string & sol_name, std::vector<float> &sol_value)
{
  if (Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();
    vtk_sol_array->SetName(sol_name.c_str());
    vtk_sol_array->SetNumberOfValues(sol_value.size());

    // save data into vtk data array
    for(unsigned int n=0; n<sol_value.size(); ++n)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array->InsertValue(n, (sol_value[n]));
    }

    _vtk_grid->GetPointData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }

}



#endif

// ------------------------------------------------------------
// vtkIO class members
//

/**
 * This method implements writing to a .vtu (VTK Unstructured Grid) file.
 * This is one of the new style XML dataformats, binary output is used to keep
 * the file size down.
 */
void VTK2IO::write (const std::string& name)
{

  const MeshBase& mesh = FieldOutput<SimulationSystem>::system().mesh();
  mesh.boundary_info->build_on_processor_side_list (_el, _sl, _il);

  // vtk file extension have a ".vtu" format?
  if(name.rfind(".vtu") < name.size())
  {
#ifdef HAVE_VTK
    // use vtk library routine to export mesh and solution as base64 binary file
    _vtk_grid = vtkUnstructuredGrid::New();


    nodes_to_vtk(mesh);
    cells_to_vtk(mesh);
    meshinfo_to_vtk(mesh);
    solution_to_vtk();

    // only processor 0 write VTK file
    if(Genius::processor_id() == 0)
    {
      XMLUnstructuredGridWriter* writer = XMLUnstructuredGridWriter::New();
      writer->SetInput(_vtk_grid);
      writer->setExtraHeader(this->export_extra_info());

      writer->SetFileName(name.c_str());
      writer->Write();
      writer->Delete();

    }
    //clean up
    _vtk_grid->Delete();
#endif

  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  // it seems vtk may generate FE_INVALID flag. clear it here
  feclearexcept(FE_INVALID);
#endif


}



