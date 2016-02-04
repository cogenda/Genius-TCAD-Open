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
#include "vtk_io.h"
#include "elem.h"
#include "mesh_base.h"
#include "boundary_info.h"
#include "material_define.h"
#include "parallel.h"
#include "mesh_communication.h"
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "boundary_condition_collector.h"
#include "field_source.h"
#include "ngspice_interface.h"
#include "spice_ckt.h"
#include "material.h"
#include "solver_specify.h"

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
using PhysicalUnit::J;
using PhysicalUnit::kg;
using PhysicalUnit::Pa;

#ifdef HAVE_VTK
class VTKIO::XMLUnstructuredGridWriter : public vtkXMLUnstructuredGridWriter
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
#endif

// private functions
#ifdef HAVE_VTK


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


void VTKIO::nodes_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* grid)
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


void VTKIO::cells_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* grid)
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


void VTKIO::meshinfo_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* grid)
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


std::string VTKIO::export_extra_info()
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


void VTKIO::solution_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* grid)
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


  // write node based data
  {
    unsigned int nr = system.n_regions();

    std::vector< std::vector<unsigned int> > region_node_ids(nr);


    // gather data from all the processor
    std::vector< std::vector<float> > dmin(nr);
    std::vector< std::vector<float> > Na(nr), Nd(nr), net_doping(nr), net_charge(nr);
    std::vector< std::vector<float> > psi(nr), Ec(nr), Ev(nr), qFn(nr), qFp(nr);
    std::vector< std::vector<float> > n(nr), p(nr), T(nr), Tn(nr), Tp(nr), Eqc(nr), Eqv(nr), Qcn(nr), Qcp(nr); //scalar value
    std::vector< std::vector<float> > mole_x(nr), mole_y(nr);
    std::vector< std::vector<float> > OptG(nr), PatG(nr), FG(nr), DoseRate(nr);
    std::vector< std::vector<float> > mun(nr), mup(nr);
    std::vector< std::vector<float> > R(nr), Rdir(nr), Rsrh(nr), Rauger(nr), Taun(nr), Taup(nr);
    std::vector< std::vector<float> > II(nr);
    std::vector< std::vector<float> > Emag(nr);
    std::vector< std::vector<float> > stress(nr), strain(nr);
    std::vector< std::vector<float> > dEcStrain(nr), dEvStrain(nr);

    std::vector< std::vector<float> > trap_a(nr),trap_b(nr),interface_charge(nr);

    std::vector< std::vector<std::complex<float> > > psi_ac(nr), n_ac(nr), p_ac(nr), T_ac(nr), Tn_ac(nr), Tp_ac(nr); //complex value
    std::vector< std::vector<std::complex<float> > > OptE_complex(nr), OptH_complex(nr);


    // scale back to normal unit
    double concentration_scale = std::pow(cm, -3);


    // search all the node belongs to current processor
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      SimulationRegion::const_processor_node_iterator on_processor_nodes_it = system.region(r)->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator on_processor_nodes_it_end = system.region(r)->on_processor_nodes_end();
      for(; on_processor_nodes_it!=on_processor_nodes_it_end; ++on_processor_nodes_it)
      {
        const FVM_Node * fvm_node = *on_processor_nodes_it;
        const FVM_NodeData * node_data = fvm_node->node_data();
        const SimulationRegion * region = system.region(r);
        const unsigned int id = fvm_node->root_node()->id();


        assert ( node_data != NULL );
        region_node_ids[r].push_back(id);

        {
          dmin[r].push_back(static_cast<float>(node_data->dmin()/um));
          if(semiconductor_material)
          {
            Na[r].push_back(static_cast<float>(node_data->Total_Na()/concentration_scale));
            Nd[r].push_back(static_cast<float>(node_data->Total_Nd()/concentration_scale));
            net_doping[r].push_back(static_cast<float>(node_data->Net_doping()/concentration_scale));
            net_charge[r].push_back(static_cast<float>(node_data->Net_charge()/concentration_scale));
            mole_x[r].push_back(static_cast<float>(node_data->mole_x()));
            mole_y[r].push_back(static_cast<float>(node_data->mole_y()));
            R[r].push_back(static_cast<float>(node_data->Recomb()/(concentration_scale/s)));
            Rdir[r].push_back(static_cast<float>(node_data->Recomb_Dir()/(concentration_scale/s)));
            Rsrh[r].push_back(static_cast<float>(node_data->Recomb_SRH()/(concentration_scale/s)));
            Rauger[r].push_back(static_cast<float>(node_data->Recomb_Auger()/(concentration_scale/s)));
            II[r].push_back(static_cast<float>(node_data->ImpactIonization()/(concentration_scale/s)));
            mun[r].push_back(static_cast<float>(node_data->mun()/(cm*cm/V/s)));
            mup[r].push_back(static_cast<float>(node_data->mup()/(cm*cm/V/s)));
            Ec[r].push_back(static_cast<float>(node_data->Ec()/eV));
            Ev[r].push_back(static_cast<float>(node_data->Ev()/eV));
            qFn[r].push_back(static_cast<float>(node_data->qFn()/eV));
            qFp[r].push_back(static_cast<float>(node_data->qFp()/eV));
            stress[r].push_back(static_cast<float>(node_data->stress().size()/Pa));
            strain[r].push_back(static_cast<float>(node_data->strain().size()));
            dEcStrain[r].push_back(static_cast<float>(node_data->dEcStrain())/eV);
            dEvStrain[r].push_back(static_cast<float>(node_data->dEvStrain())/eV);
          }

          psi[r].push_back(static_cast<float>(node_data->psi()/V));
          Emag[r].push_back(static_cast<float>(node_data->E().size()/(V/cm)));
          n[r].push_back(static_cast<float>(node_data->n()/concentration_scale));
          p[r].push_back(static_cast<float>(node_data->p()/concentration_scale));
          T[r].push_back(static_cast<float>(node_data->T()/K));

          if(ebm_solution_data)
          {
            Tn[r].push_back(static_cast<float>(node_data->Tn()/K));
            Tp[r].push_back(static_cast<float>(node_data->Tp()/K));
          }

          if(dg_solution_data)
          {
            Eqc[r].push_back(static_cast<float>(node_data->Eqc()/eV));
            Eqv[r].push_back(static_cast<float>(node_data->Eqv()/eV));
            Qcn[r].push_back(static_cast<float>((node_data->Eqc() - node_data->Ec())/eV));
            Qcp[r].push_back(static_cast<float>((node_data->Eqv() - node_data->Ev())/eV));
          }

          if(optical_generation)
          {
            OptG[r].push_back(static_cast<float>(node_data->OptG()/(concentration_scale/s)));
          }
          if(particle_generation)
          {
            PatG[r].push_back(static_cast<float>(node_data->PatG()/(concentration_scale/s)));
          }

          if(field_generation)
          {
            FG[r].push_back(static_cast<float>(node_data->Field_G()/(concentration_scale/s)));
            DoseRate[r].push_back(static_cast<float>(node_data->DoseRate()/(0.01*J/kg/s)));
          }

          if(ddm_ac_data)
          {
            psi_ac[r].push_back(std::complex<float>(node_data->psi_ac().real()/V, node_data->psi_ac().imag()/V));
            n_ac[r].push_back(  std::complex<float>(node_data->n_ac().real()/concentration_scale, node_data->n_ac().imag()/concentration_scale));
            p_ac[r].push_back(  std::complex<float>(node_data->p_ac().real()/concentration_scale, node_data->p_ac().imag()/concentration_scale));
            T_ac[r].push_back(  std::complex<float>(node_data->T_ac().real()/K, node_data->T_ac().imag()/K));
            Tn_ac[r].push_back(  std::complex<float>(node_data->Tn_ac().real()/K, node_data->Tn_ac().imag()/K));
            Tp_ac[r].push_back(  std::complex<float>(node_data->Tp_ac().real()/K, node_data->Tp_ac().imag()/K));
          }

          if(optical_complex_field)
          {
            OptE_complex[r].push_back(std::complex<float>(node_data->OptE_complex().real()/(V/cm), node_data->OptE_complex().imag()/(V/cm)));
            OptH_complex[r].push_back(std::complex<float>(node_data->OptH_complex().real()/(A/cm), node_data->OptH_complex().imag()/(A/cm)));
          }

          if(tid_data)
          {
            trap_a[r].push_back(static_cast<float>(node_data->trap_a()/concentration_scale));
            trap_b[r].push_back(static_cast<float>(node_data->trap_b()/concentration_scale));
            interface_charge[r].push_back(static_cast<float>(node_data->interface_charge()/std::pow(cm, -2)));
          }

          /*
          if( region->type() == SemiconductorRegion )
          {
            const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
            semiconductor_region->material()->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);
            Taun.push_back(static_cast<float>(semiconductor_region->material()->band->TAUN(node_data->T())/s));
            Taup.push_back(static_cast<float>(semiconductor_region->material()->band->TAUP(node_data->T())/s));
          }
          else
          {
            Taun.push_back(static_cast<float>(0.0));
            Taup.push_back(static_cast<float>(0.0));
          }
          */
        }
      }
    }

    write_node_scaler_solution(region_node_ids, dmin, "dmin[um]", grid);
    write_node_scaler_solution(region_node_ids, psi,  "potential[V]", grid);
    write_node_scaler_solution(region_node_ids, Emag, "E[V/cm]", grid);
    write_node_scaler_solution(region_node_ids, n,    "elec_density[cm-3]", grid);
    write_node_scaler_solution(region_node_ids, p,    "hole_density[cm-3]", grid);

    if(semiconductor_material)
    {
      write_node_scaler_solution(region_node_ids, Na,  "Na[cm-3]", grid);
      write_node_scaler_solution(region_node_ids, Nd,  "Nd[cm-3]", grid);

      write_node_scaler_solution(region_node_ids, net_doping,  "net_doping[cm-3]", grid);
      write_node_scaler_solution(region_node_ids, net_charge,  "net_charge[cm-3]", grid);
      write_node_scaler_solution(region_node_ids, mole_x, "mole_x", grid);
      write_node_scaler_solution(region_node_ids, mole_y, "mole_y", grid);
      write_node_scaler_solution(region_node_ids, R, "recombination[cm-3/s]", grid);
      write_node_scaler_solution(region_node_ids, Rdir, "recombination_dir[cm-3/s]", grid);
      write_node_scaler_solution(region_node_ids, Rsrh, "recombination_srh[cm-3/s]", grid);
      write_node_scaler_solution(region_node_ids, Rauger, "recombination_auger[cm-3/s]", grid);
      write_node_scaler_solution(region_node_ids, mun, "elec_mobility[cm2/V/s]", grid);
      write_node_scaler_solution(region_node_ids, mup, "hole_mobility[cm2/V/s]", grid);
      write_node_scaler_solution(region_node_ids, II, "impact_ionization[cm-3/s]", grid);
      write_node_scaler_solution(region_node_ids, Ec,  "Ec[eV]",  grid);
      write_node_scaler_solution(region_node_ids, Ev,  "Ev[eV]",  grid);
      write_node_scaler_solution(region_node_ids, qFn, "elec_quasi_Fermi[eV]", grid);
      write_node_scaler_solution(region_node_ids, qFp, "hole_quasi_Fermi[eV]", grid);
      write_node_scaler_solution(region_node_ids, stress, "stress[Pa]", grid);
      write_node_scaler_solution(region_node_ids, strain, "strain", grid);
      write_node_scaler_solution(region_node_ids, dEcStrain, "dEc_strain[eV]", grid);
      write_node_scaler_solution(region_node_ids, dEvStrain, "dEv_strain[eV]", grid);
    }

    write_node_scaler_solution(region_node_ids, T,   "temperature[K]", grid);

    if(ebm_solution_data)
    {
      write_node_scaler_solution(region_node_ids, Tn,  "elec_temperature[K]", grid);
      write_node_scaler_solution(region_node_ids, Tp,  "hole_temperature[K]", grid);
    }

    if(dg_solution_data)
    {
      write_node_scaler_solution(region_node_ids, Eqc,  "quantum_conduction_band[eV]", grid);
      write_node_scaler_solution(region_node_ids, Eqv,  "quantum_valence_band[eV]", grid);
      write_node_scaler_solution(region_node_ids, Qcn,  "quantum_potential_elec[V]", grid);
      write_node_scaler_solution(region_node_ids, Qcp,  "quantum_potential_hole[V]", grid);
    }

    if(ddm_ac_data)
    {
      write_node_complex_solution(region_node_ids, psi_ac, "potential_AC[V]", grid);
      write_node_complex_solution(region_node_ids, n_ac,   "elec_AC[cm-3]", grid);
      write_node_complex_solution(region_node_ids, p_ac,   "hole_AC[cm-3]", grid);
      write_node_complex_solution(region_node_ids, T_ac,   "temperature_AC[K]", grid);
      write_node_complex_solution(region_node_ids, Tn_ac,  "elec_temperature_AC[K]", grid);
      write_node_complex_solution(region_node_ids, Tp_ac,  "hole_temperature_AC[K]", grid);
    }

    if(optical_complex_field)
    {
      write_node_complex_solution(region_node_ids, OptE_complex, "opt_E[V/cm]", grid);
      write_node_complex_solution(region_node_ids, OptH_complex, "opt_H[A/cm]", grid);
    }

    if(optical_generation)
    {
      write_node_scaler_solution(region_node_ids, OptG,   "optical_generation[cm-3/s]", grid);
    }
    if(particle_generation)
    {
      write_node_scaler_solution(region_node_ids, PatG,   "radiation_generation[cm-3/s]", grid);
    }
    if(field_generation)
    {
      write_node_scaler_solution(region_node_ids, FG,       "field_generation[cm-3/s]", grid);
      write_node_scaler_solution(region_node_ids, DoseRate, "dose_rate[rad/s]", grid);
    }

    if(tid_data)
    {
      write_node_scaler_solution(region_node_ids, trap_a,   "trap_A_density[cm-3]", grid);
      write_node_scaler_solution(region_node_ids, trap_b,   "trap_B_density[cm-3]", grid);
      write_node_scaler_solution(region_node_ids, interface_charge,   "interface_charge[cm-2]", grid);
    }

    //write_node_scaler_solution(region_node_ids, Taun,   "Taun", grid);
    //write_node_scaler_solution(region_node_ids, Taup,   "Taup", grid);

  }

  // write cell based data
  {
    std::map<const Elem *, const FVM_CellData *> elem_to_elem_data_map;
    std::vector<unsigned int> order;

    //std::vector<float> mos_channel_flag;

    std::vector<float> Ex,  Ey,  Ez;
    std::vector<float> Jnx,  Jny,  Jnz;
    std::vector<float> Jpx,  Jpy,  Jpz;

    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      for(unsigned int n=0; n<region->n_cell(); ++n)
      {
        const Elem * elem = region->get_region_elem(n);
        if( elem->processor_id() != Genius::processor_id() ) continue;
        const FVM_CellData * elem_data = region->get_region_elem_data(n);
        elem_to_elem_data_map.insert( std::make_pair(elem, elem_data) );
        order.push_back(elem->id());
        /*
        if( region->type()==SemiconductorRegion)
        {
          const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
          mos_channel_flag.push_back(semiconductor_region->is_elem_in_mos_channel(elem));
        }
        else
          mos_channel_flag.push_back(0.0);
        */

        Ex.push_back(static_cast<float>(elem_data->E()(0)/(V/cm)));
        Ey.push_back(static_cast<float>(elem_data->E()(1)/(V/cm)));
        Ez.push_back(static_cast<float>(elem_data->E()(2)/(V/cm)));

        Jnx.push_back(static_cast<float>(elem_data->Jn()(0)/(A/(cm*cm))));
        Jny.push_back(static_cast<float>(elem_data->Jn()(1)/(A/(cm*cm))));
        Jnz.push_back(static_cast<float>(elem_data->Jn()(2)/(A/(cm*cm))));

        Jpx.push_back(static_cast<float>(elem_data->Jp()(0)/(A/(cm*cm))));
        Jpy.push_back(static_cast<float>(elem_data->Jp()(1)/(A/(cm*cm))));
        Jpz.push_back(static_cast<float>(elem_data->Jp()(2)/(A/(cm*cm))));

      }
    }

    // write solution for boundary elements, just keep the same as mesh elems
    {
      std::vector<unsigned int> boundary_elem_ids(_el);
      std::vector<unsigned short int> boundary_elem_sides(_sl);
      Parallel::gather(0, boundary_elem_ids);
      Parallel::gather(0, boundary_elem_sides);

      typedef std::pair<unsigned int, unsigned short int> boundary_elem_key;
      std::map<boundary_elem_key, unsigned int> boundary_face_order;
      for(unsigned int n=0; n<boundary_elem_ids.size(); ++n)
      {
        boundary_elem_key key = std::make_pair(boundary_elem_ids[n], boundary_elem_sides[n]);
        boundary_face_order.insert(std::make_pair(key, 0));
      }
      std::map< boundary_elem_key, unsigned int >::iterator boundary_face_it=boundary_face_order.begin();
      for(unsigned int i=0;boundary_face_it != boundary_face_order.end(); ++boundary_face_it, ++i)
      {
        boundary_face_it->second = i;
      }

      unsigned int id_offset = mesh.n_elem();
      for(unsigned int n=0; n<_il.size(); ++n)
      {
        const Elem * elem = mesh.elem(_el[n]);
        if( elem->processor_id() != Genius::processor_id() ) continue;
        genius_assert(elem_to_elem_data_map.find(elem) != elem_to_elem_data_map.end());
        const FVM_CellData * elem_data = elem_to_elem_data_map.find(elem)->second;

        boundary_elem_key key = std::make_pair(_el[n], _sl[n]);
        order.push_back(id_offset + boundary_face_order.find(key)->second);

        //std::cout<<id_offset + boundary_face_order.find(key)->second<<std::endl;
        //mos_channel_flag.push_back(0.5);


#if 1
        Ex.push_back(0);
        Ey.push_back(0);
        Ez.push_back(0);

        Jnx.push_back(0);
        Jny.push_back(0);
        Jnz.push_back(0);

        Jpx.push_back(0);
        Jpy.push_back(0);
        Jpz.push_back(0);
#endif
      }
    }

    Parallel::gather(0, order);
    //write_cell_scaler_solution(order, mos_channel_flag,  "mos channel", grid);
    write_cell_vector_solution(order, Ex,  Ey,  Ez,  "electrical_field[V/cm]", grid);
    write_cell_vector_solution(order, Jnx, Jny, Jnz, "elec_current[A/cm^2]", grid);
    write_cell_vector_solution(order, Jpx, Jpy, Jpz, "hole_current[A/cm^2]", grid);
  }
}







void VTKIO::write_node_scaler_solution(const std::vector< std::vector<unsigned int> >& region_order,
                                       const std::vector< std::vector<float> > &region_sol,
                                       const std::string & sol_name, vtkUnstructuredGrid* grid)
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


void VTKIO::write_node_complex_solution(const std::vector< std::vector<unsigned int> >& region_order,
                                        const std::vector< std::vector<std::complex<float> > > & region_sol,
                                        const std::string & sol_name, vtkUnstructuredGrid* grid)
{
  // this should run on parallel for all the processor
  std::vector<float> solution_real;
  std::vector<float> solution_imag;

  for(unsigned int r=0; r<region_order.size(); ++r)
  {
    const std::vector<unsigned int> & order = region_order[r];
    const std::vector<std::complex<float> > & sol = region_sol[r];

    std::map<unsigned int, float> region_variable_real;
    std::map<unsigned int, float> region_variable_imag;
    for(unsigned int i=0; i<order.size(); ++i)
    {
      region_variable_real.insert(std::make_pair(order[i], sol[i].real()));
      region_variable_imag.insert(std::make_pair(order[i], sol[i].imag()));
    }

    Parallel::allgather(region_variable_real);
    Parallel::allgather(region_variable_imag);

    std::map<unsigned int, float>::const_iterator real_it = region_variable_real.begin();
    for( ; real_it != region_variable_real.end(); real_it++)
      solution_real.push_back(real_it->second);

    std::map<unsigned int, float>::const_iterator imag_it = region_variable_imag.begin();
    for( ; imag_it != region_variable_imag.end(); imag_it++)
      solution_imag.push_back(imag_it->second);
  }

  genius_assert(solution_real.size() == solution_imag.size());

  if ( Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array_magnitude = vtkFloatArray::New();
    vtkFloatArray *vtk_sol_array_angle     = vtkFloatArray::New();

    vtk_sol_array_magnitude->SetName((sol_name+" abs").c_str());
    vtk_sol_array_magnitude->SetNumberOfValues(solution_real.size());

    vtk_sol_array_angle->SetName((sol_name+" angle").c_str());
    vtk_sol_array_angle->SetNumberOfValues(solution_real.size());

    // save data into vtk data array
    for(unsigned int n=0; n<solution_real.size(); ++n)
    {
      std::complex<float> data(solution_real[n], solution_imag[n]);
      // set scalar value to vtkFloatArray
      vtk_sol_array_magnitude->InsertValue(n, std::abs(data));
      vtk_sol_array_angle->InsertValue(n, std::arg(data));
    }

    grid->GetPointData()->AddArray(vtk_sol_array_magnitude);
    grid->GetPointData()->AddArray(vtk_sol_array_angle);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array_magnitude->Delete();
    vtk_sol_array_angle->Delete();
  }

}


void VTKIO::write_node_vector_solution(const std::vector< std::vector<unsigned int> >& region_order,
                                       const std::vector< std::vector<float> > & region_x,
                                       const std::vector< std::vector<float> > & region_y,
                                       const std::vector< std::vector<float> > & region_z,
                                       const std::string & sol_name, vtkUnstructuredGrid* grid)
{
  // this should run on parallel for all the processor
  std::vector<float> solution_x;
  std::vector<float> solution_y;
  std::vector<float> solution_z;

  for(unsigned int r=0; r<region_order.size(); ++r)
  {
    const std::vector<unsigned int> & order = region_order[r];
    const std::vector<float> & vec_x = region_x[r];
    const std::vector<float> & vec_y = region_y[r];
    const std::vector<float> & vec_z = region_z[r];

    std::map<unsigned int, float> region_variable_x;
    std::map<unsigned int, float> region_variable_y;
    std::map<unsigned int, float> region_variable_z;

    for(unsigned int i=0; i<order.size(); ++i)
    {
      region_variable_x.insert(std::make_pair(order[i], vec_x[i]));
      region_variable_y.insert(std::make_pair(order[i], vec_y[i]));
      region_variable_z.insert(std::make_pair(order[i], vec_z[i]));
    }

    Parallel::allgather(region_variable_x);
    Parallel::allgather(region_variable_y);
    Parallel::allgather(region_variable_z);

    std::map<unsigned int, float>::const_iterator x_it = region_variable_x.begin();
    for( ; x_it != region_variable_x.end(); x_it++)
      solution_x.push_back(x_it->second);

    std::map<unsigned int, float>::const_iterator y_it = region_variable_y.begin();
    for( ; y_it != region_variable_y.end(); y_it++)
      solution_y.push_back(y_it->second);

    std::map<unsigned int, float>::const_iterator z_it = region_variable_z.begin();
    for( ; z_it != region_variable_z.end(); z_it++)
      solution_z.push_back(z_it->second);
  }

  if ( Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();

    vtk_sol_array->SetNumberOfComponents(3);
    vtk_sol_array->SetNumberOfTuples(solution_x.size());
    vtk_sol_array->SetName(sol_name.c_str());

    // save data into vtk data array
    for(unsigned int n=0; n<solution_x.size(); ++n)
    {
      float sol[3]={solution_x[n], solution_y[n], solution_z[n]};
      vtk_sol_array->InsertTuple(n, &sol[0] );
    }

    grid->GetPointData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }
}



void VTKIO::write_cell_scaler_solution(const std::vector<unsigned int> & order, std::vector<float> &sol,
                                       const std::string & sol_name, vtkUnstructuredGrid* grid)
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


void VTKIO::write_cell_vector_solution(const std::vector<unsigned int> & order,
                                       std::vector<float > & sol_x,
                                       std::vector<float > & sol_y,
                                       std::vector<float > & sol_z,
                                       const std::string & sol_name, vtkUnstructuredGrid* grid)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol_x);
  Parallel::gather(0, sol_y);
  Parallel::gather(0, sol_z);

  if ( Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();

    vtk_sol_array->SetNumberOfComponents(3);
    vtk_sol_array->SetNumberOfTuples(order.size());
    vtk_sol_array->SetName(sol_name.c_str());

    std::vector<float> sol_re_order_x(order.size());
    std::vector<float> sol_re_order_y(order.size());
    std::vector<float> sol_re_order_z(order.size());
    for(unsigned int n=0; n<order.size(); ++n)
    {
      sol_re_order_x[order[n]] = sol_x[n];
      sol_re_order_y[order[n]] = sol_y[n];
      sol_re_order_z[order[n]] = sol_z[n];
    }

    // save data into vtk data array
    for( unsigned int n=0; n<order.size(); n++)
    {
      float sol[3]={sol_re_order_x[n], sol_re_order_y[n], sol_re_order_z[n]};
      vtk_sol_array->InsertTuple(n, &sol[0] );
    }

    grid->GetCellData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }
}


#endif




void VTKIO::_smooth_numerical_error(std::map<unsigned int, float> &dataset, float tol)
{
  float min=1e38, max=-1e38;
  std::map<unsigned int, float>::iterator it = dataset.begin();
  for(; it!=dataset.end(); ++it)
  {
    float v = it->second;
    min = std::min(min, v);
    max = std::max(max, v);
  }

  if(max-min > tol) return;

  float smooth_v = 0.5*(min+max);
  for(it = dataset.begin(); it!=dataset.end(); ++it)
  {
    it->second = smooth_v;
  }

}


// ------------------------------------------------------------
// vtkIO class members
//

/**
 * This method implements writing to a .vtu (VTK Unstructured Grid) file.
 * This is one of the new style XML dataformats, binary output is used to keep
 * the file size down.
 */
void VTKIO::write (const std::string& name)
{

  const MeshBase& mesh = FieldOutput<SimulationSystem>::system().mesh();
  mesh.boundary_info->build_on_processor_side_list (_el, _sl, _il);

  // vtk file extension have a ".vtu" format?
  if(name.rfind(".vtu") < name.size())
  {
#ifdef HAVE_VTK
    // use vtk library routine to export mesh and solution as base64 binary file
    _vtk_grid = vtkUnstructuredGrid::New();

    nodes_to_vtk(mesh, _vtk_grid);
    cells_to_vtk(mesh, _vtk_grid);
    meshinfo_to_vtk(mesh, _vtk_grid);
    solution_to_vtk(mesh, _vtk_grid);


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






