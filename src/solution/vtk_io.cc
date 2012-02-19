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

void VTKIO::nodes_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* grid)
{
  unsigned int n_nodes = mesh.n_nodes();
  // collect node location, must run in parallel
  std::vector<Real> pts;
  mesh.pack_nodes(pts);

  if(Genius::processor_id() == 0)
  {
    vtkPoints* points = vtkPoints::New();

    // write out nodal data
    for ( unsigned int n=0; n<n_nodes; ++n )
    {
      float tuple[3];
      {
        tuple[0] =  pts[3*n+0]/um; //scale to um
        tuple[1] =  pts[3*n+1]/um; //scale to um
        tuple[2] =  pts[3*n+2]/um; //scale to um
      }
      points->InsertPoint(n, tuple);
    }

    grid->SetPoints(points);

    //clean up
    points->Delete();
  }

}


void VTKIO::cells_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* grid)
{

  // must run in parallel
  std::vector<int> vtk_cell_conns;
  std::vector<unsigned int> vtk_cell_ids;
  {
    MeshBase::const_element_iterator       it  = mesh.active_this_pid_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_this_pid_elements_end();
    for ( ; it != end; ++it)
    {
      const Elem *elem  = (*it);
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


      // get the connectivity for this element
      std::vector<unsigned int> conn;
      elem->connectivity(0,VTK,conn);

      vtk_cell_ids.push_back(elem->id());

      vtk_cell_conns.push_back(celltype);
      vtk_cell_conns.push_back(conn.size());
      for(unsigned int i=0;i<conn.size();++i)
        vtk_cell_conns.push_back(conn[i]);

    } // end loop over active elements

    Parallel::gather(0, vtk_cell_ids);
    Parallel::gather(0, vtk_cell_conns);
  }

  std::vector<int> vtk_boundary_cell_conns;
  std::vector<unsigned int> vtk_boundary_cell_ids(_el);
  std::vector<unsigned short int> vtk_boundary_cell_sides(_sl);
  {
    for(unsigned int n=0; n<_il.size(); ++n)
    {
      const Elem * elem = mesh.elem(_el[n]);
      AutoPtr<Elem> boundary_elem =  elem->build_side(_sl[n]);

      vtkIdType celltype = VTK_EMPTY_CELL; // initialize to something to avoid compiler warning

      switch(boundary_elem->type())
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
          celltype = VTK_TRIANGLE;
          break;// 3
          case TRI6:
          celltype = VTK_QUADRATIC_TRIANGLE;
          break;// 4
          case QUAD4:
          case QUAD4_FVM:
          celltype = VTK_QUAD;
          break;// 5
          case QUAD8:
          celltype = VTK_QUADRATIC_QUAD;
          break;// 6
          default:
          {
            std::cerr<<"element type "<<elem->type()<<" not implemented"<<std::endl;
            genius_error();
          }
      }

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
        vtkIdList *pts = vtkIdList::New();
        unsigned int n_nodes = vtk_cell_conns[cnt++];
        pts->SetNumberOfIds(n_nodes);
        for(unsigned int i=0;i<n_nodes;++i)
          pts->SetId(i, vtk_cell_conns[cnt++]);
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

        vtkIdType celltype = vtk_boundary_cell_conns[cnt++];
        vtkIdList *pts = vtkIdList::New();
        unsigned int n_nodes = vtk_boundary_cell_conns[cnt++];
        pts->SetNumberOfIds(n_nodes);
        for(unsigned int i=0;i<n_nodes;++i)
          pts->SetId(i, vtk_boundary_cell_conns[cnt++]);
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

    grid->GetCellData()->AddArray(region_info);
    grid->GetCellData()->AddArray(boundary_info);
    grid->GetCellData()->AddArray(partition_info);

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
      semiconductor_material=true;
  }


  bool ebm_solution_data     = false;
  bool optical_generation    = false;
  bool particle_generation   = false;
  bool ddm_ac_data           = false;
  bool optical_complex_field = false;

  const std::vector<SolverSpecify::SolverType> & solve_history = system.solve_history();

  for(unsigned int n=0; n<solve_history.size(); ++n)
  {
    switch  (solve_history[n])
    {
        case SolverSpecify::EBML3      :
        case SolverSpecify::EBML3MIX   :
        ebm_solution_data = true;
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
        default : break;
    }
  }

  if(system.get_field_source()->is_light_source_exist())
    optical_generation = true;

  if(system.get_field_source()->is_particle_source_exist())
    particle_generation = true;


  // write node based data
  {
    // scale back to normal unit
    double concentration_scale = pow(cm, -3);

    std::set<unsigned int> recorder;
    std::vector<unsigned int> order;

    // gather data from all the processor
    std::vector<float> Na, Nd, net_doping, net_charge;
    std::vector<float> psi, Ec, Ev, qFn, qFp;
    std::vector<float> n, p, T, Tn, Tp; //scalar value
    std::vector<float> mole_x, mole_y;
    std::vector<float> OptG, OptE, PatG, PatE;
    std::vector<float> R, Taun, Taup;
    std::vector<float> II;

    std::vector<std::complex<float> >psi_ac, n_ac, p_ac, T_ac, Tn_ac, Tp_ac; //complex value
    std::vector<std::complex<float> > OptE_complex, OptH_complex;

    std::vector<float> Jnx, Jny, Jnz; //vector value
    std::vector<float> Jpx, Jpy, Jpz;
    std::vector<float> Vnx, Vny, Vnz;
    std::vector<float> Vpx, Vpy, Vpz;
    std::vector<float> Ex,  Ey,  Ez;

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

        // if the fvm_node lies on the interface of two material regions,
        // we shall use the node data in the more important region.
        // it just for visualization reason
        if( fvm_node->boundary_id() != BoundaryInfo::invalid_id )
        {
          // test if we had already processed this node
          if( recorder.find(id) != recorder.end() ) continue;
          recorder.insert(id);

          unsigned int bc_index = system.get_bcs()->get_bc_index_by_bd_id(fvm_node->boundary_id());
          const BoundaryCondition * bc = system.get_bcs()->get_bc(bc_index);
          const FVM_Node * primary_fvm_node = (*bc->region_node_begin(fvm_node->root_node())).second.second;
          fvm_node = primary_fvm_node;
          node_data = fvm_node->node_data();
          region = system.region(primary_fvm_node->subdomain_id());
        }

        assert ( node_data != NULL );
        order.push_back(id);

        {
          if(semiconductor_material)
          {
            Na.push_back(static_cast<float>(node_data->Total_Na()/concentration_scale));
            Nd.push_back(static_cast<float>(node_data->Total_Nd()/concentration_scale));
            net_doping.push_back(static_cast<float>(node_data->Net_doping()/concentration_scale));
            net_charge.push_back(static_cast<float>(node_data->Net_charge()/concentration_scale));
            n.push_back(static_cast<float>(node_data->n()/concentration_scale));
            p.push_back(static_cast<float>(node_data->p()/concentration_scale));
            mole_x.push_back(static_cast<float>(node_data->mole_x()));
            mole_y.push_back(static_cast<float>(node_data->mole_y()));
            R.push_back(static_cast<float>(node_data->Recomb()/(concentration_scale/s)));
            II.push_back(static_cast<float>(node_data->ImpactIonization()/(concentration_scale/s)));
            /*
            Jnx.push_back(static_cast<float>(node_data->Jn()(0)/(A/cm)));
            Jny.push_back(static_cast<float>(node_data->Jn()(1)/(A/cm)));
            Jnz.push_back(static_cast<float>(node_data->Jn()(2)/(A/cm)));

            Jpx.push_back(static_cast<float>(node_data->Jp()(0)/(A/cm)));
            Jpy.push_back(static_cast<float>(node_data->Jp()(1)/(A/cm)));
            Jpz.push_back(static_cast<float>(node_data->Jp()(2)/(A/cm)));

            Ex.push_back(static_cast<float>(node_data->E()(0)/(V/cm)));
            Ey.push_back(static_cast<float>(node_data->E()(1)/(V/cm)));
            Ez.push_back(static_cast<float>(node_data->E()(2)/(V/cm)));
            */
          }

          psi.push_back(static_cast<float>(node_data->psi()/V));
          Ec.push_back(static_cast<float>(node_data->Ec()/eV));
          Ev.push_back(static_cast<float>(node_data->Ev()/eV));
          qFn.push_back(static_cast<float>(node_data->qFn()/eV));
          qFp.push_back(static_cast<float>(node_data->qFp()/eV));
          T.push_back(static_cast<float>(node_data->T()/K));

          if(ebm_solution_data)
          {
            Tn.push_back(static_cast<float>(node_data->Tn()/K));
            Tp.push_back(static_cast<float>(node_data->Tp()/K));
          }

          if(optical_generation)
          {
            OptE.push_back(static_cast<float>(node_data->OptE()));
            OptG.push_back(static_cast<float>(node_data->Field_G()/(concentration_scale/s)));
          }
          if(particle_generation)
          {
            PatE.push_back(static_cast<float>(node_data->PatE()/(eV)));
            PatG.push_back(static_cast<float>(node_data->Field_G()/(concentration_scale/s)));
          }

          if(ddm_ac_data)
          {
            psi_ac.push_back(static_cast<std::complex<float> >(node_data->psi_ac()/V));
            n_ac.push_back(static_cast<std::complex<float> >(node_data->n_ac()/concentration_scale));
            p_ac.push_back(static_cast<std::complex<float> >(node_data->p_ac()/concentration_scale));
            T_ac.push_back(static_cast<std::complex<float> >(node_data->T_ac()/K));
            Tn_ac.push_back(static_cast<std::complex<float> >(node_data->Tn_ac()/K));
            Tp_ac.push_back(static_cast<std::complex<float> >(node_data->Tp_ac()/K));
          }

          if(optical_complex_field)
          {
            OptE_complex.push_back(static_cast<std::complex<float> >(node_data->OptE_complex()/(V/cm)));
            OptH_complex.push_back(static_cast<std::complex<float> >(node_data->OptH_complex()/(A/cm)));
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

    Parallel::gather(0, order);
    if (Genius::processor_id() == 0)
      genius_assert( order.size() == mesh.n_nodes() );

    write_node_scaler_solution(order, psi, "potential", grid);
    write_node_scaler_solution(order, Ec,  "Ec",  grid);
    write_node_scaler_solution(order, Ev,  "Ev",  grid);
    write_node_scaler_solution(order, qFn, "elec_quasi_Fermi_level", grid);
    write_node_scaler_solution(order, qFp, "hole_quasi_Fermi_level", grid);

    if(semiconductor_material)
    {
      write_node_scaler_solution(order, Na,  "Na", grid);
      write_node_scaler_solution(order, Nd,  "Nd", grid);
      write_node_scaler_solution(order, n,   "electron_density", grid);
      write_node_scaler_solution(order, p,   "hole_density", grid);
      write_node_scaler_solution(order, net_doping,  "net_doping", grid);
      write_node_scaler_solution(order, net_charge,  "net_charge", grid);
      write_node_scaler_solution(order, mole_x, "mole_x", grid);
      write_node_scaler_solution(order, mole_y, "mole_y", grid);
      write_node_scaler_solution(order, R, "recombination", grid);
      write_node_scaler_solution(order, II, "impact_ionization", grid);
      //write_node_vector_solution(order, Vnx, Vny, Vnz, "Vn_node", grid);
      //write_node_vector_solution(order, Vpx, Vpy, Vpz, "Vp_node", grid);
      //write_node_vector_solution(order, Jnx, Jny, Jnz, "Jn_node", grid);
      //write_node_vector_solution(order, Jpx, Jpy, Jpz, "Jp_node", grid);
      //write_node_vector_solution(order, Ex, Ey, Ez, "E_node", grid);
    }

    write_node_scaler_solution(order, T,   "temperature", grid);

    if(ebm_solution_data)
    {
      write_node_scaler_solution(order, Tn,  "elec_temperature", grid);
      write_node_scaler_solution(order, Tp,  "hole_temperature", grid);
    }

    if(ddm_ac_data)
    {
      write_node_complex_solution(order, psi_ac, "potential_AC", grid);
      write_node_complex_solution(order, n_ac,   "electron_AC", grid);
      write_node_complex_solution(order, p_ac,   "hole_AC", grid);
      write_node_complex_solution(order, T_ac,   "temperature_AC", grid);
      write_node_complex_solution(order, Tn_ac,  "elec_temperature_AC", grid);
      write_node_complex_solution(order, Tp_ac,  "hole_temperature_AC", grid);
    }

    if(optical_complex_field)
    {
      write_node_complex_solution(order, OptE_complex, "Optical_E", grid);
      write_node_complex_solution(order, OptH_complex, "Optical_H", grid);
    }

    if(optical_generation)
    {
      write_node_scaler_solution(order, OptE,   "optical_energy_density", grid);
      write_node_scaler_solution(order, OptG,   "optical_generation", grid);
    }
    if(particle_generation)
    {
      write_node_scaler_solution(order, PatE,   "radiation_energy_density", grid);
      write_node_scaler_solution(order, PatG,   "radiation_generation", grid);
    }
    //write_node_scaler_solution(order, Taun,   "Taun", grid);
    //write_node_scaler_solution(order, Taup,   "Taup", grid);
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

        Jnx.push_back(static_cast<float>(elem_data->Jn()(0)/(A/cm)));
        Jny.push_back(static_cast<float>(elem_data->Jn()(1)/(A/cm)));
        Jnz.push_back(static_cast<float>(elem_data->Jn()(2)/(A/cm)));

        Jpx.push_back(static_cast<float>(elem_data->Jp()(0)/(A/cm)));
        Jpy.push_back(static_cast<float>(elem_data->Jp()(1)/(A/cm)));
        Jpz.push_back(static_cast<float>(elem_data->Jp()(2)/(A/cm)));

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
#if 0
        Ex.push_back(static_cast<float>(elem_data->E()(0)/(V/cm)));
        Ey.push_back(static_cast<float>(elem_data->E()(1)/(V/cm)));
        Ez.push_back(static_cast<float>(elem_data->E()(2)/(V/cm)));

        Jnx.push_back(static_cast<float>(elem_data->Jn()(0)/(A/cm)));
        Jny.push_back(static_cast<float>(elem_data->Jn()(1)/(A/cm)));
        Jnz.push_back(static_cast<float>(elem_data->Jn()(2)/(A/cm)));

        Jpx.push_back(static_cast<float>(elem_data->Jp()(0)/(A/cm)));
        Jpy.push_back(static_cast<float>(elem_data->Jp()(1)/(A/cm)));
        Jpz.push_back(static_cast<float>(elem_data->Jp()(2)/(A/cm)));
#endif

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
    write_cell_vector_solution(order, Ex,  Ey,  Ez,  "electrical_field", grid);
    write_cell_vector_solution(order, Jnx, Jny, Jnz, "electron_current", grid);
    write_cell_vector_solution(order, Jpx, Jpy, Jpz, "hole_current", grid);
  }
}


void VTKIO::write_node_scaler_solution(const std::vector<unsigned int> & order, std::vector<float> &sol,
                                       const std::string & sol_name, vtkUnstructuredGrid* grid)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if (Genius::processor_id() == 0)
  {
    genius_assert(order.size() == sol.size());
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

    grid->GetPointData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }
}


void VTKIO::write_node_complex_solution(const std::vector<unsigned int> & order, std::vector<std::complex<float> > & sol,
                                        const std::string & sol_name, vtkUnstructuredGrid* grid)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if ( Genius::processor_id() == 0)
  {
    genius_assert(order.size() == sol.size());

    //create vtk data array
    vtkFloatArray *vtk_sol_array_magnitude = vtkFloatArray::New();
    vtkFloatArray *vtk_sol_array_angle     = vtkFloatArray::New();

    vtk_sol_array_magnitude->SetName((sol_name+" abs").c_str());
    vtk_sol_array_magnitude->SetNumberOfValues(order.size());

    vtk_sol_array_angle->SetName((sol_name+" angle").c_str());
    vtk_sol_array_angle->SetNumberOfValues(order.size());

    std::vector<float> sol_re_order_abs(order.size());
    for(unsigned int n=0; n<order.size(); ++n)
      sol_re_order_abs[order[n]] = std::abs(sol[n]);

    std::vector<float> sol_re_order_arg(order.size());
    for(unsigned int n=0; n<order.size(); ++n)
      sol_re_order_arg[order[n]] = std::arg(sol[n]);

    // save data into vtk data array
    for(unsigned int n=0; n<order.size(); ++n)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array_magnitude->InsertValue(n, sol_re_order_abs[n]);
      vtk_sol_array_angle->InsertValue(n, sol_re_order_arg[n]);
    }

    grid->GetPointData()->AddArray(vtk_sol_array_magnitude);
    grid->GetPointData()->AddArray(vtk_sol_array_angle);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array_magnitude->Delete();
    vtk_sol_array_angle->Delete();
  }

}


void VTKIO::write_node_vector_solution(const std::vector<unsigned int> & order,
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
    for(unsigned int n=0; n<order.size(); ++n)
    {
      float sol[3]={sol_re_order_x[n], sol_re_order_y[n], sol_re_order_z[n]};
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



// private functions
void VTKIO::nodes_to_vtk(const MeshBase& mesh, std::ofstream & out)
{
  unsigned int n_nodes = mesh.n_nodes();
  // collect node location, must run in parallel
  std::vector<Real> pts;
  mesh.pack_nodes(pts);

  if(Genius::processor_id() == 0)
  {
    out << "# vtk DataFile Version 3.0"      <<'\n';
    out << "Date Generated by Genius TCAD"   <<'\n';
    out << "ASCII"                           <<'\n';
    out << "DATASET UNSTRUCTURED_GRID"       <<'\n';
    out << "POINTS " << mesh.n_nodes()       << " float" << '\n';

    // write out nodal data
    for ( unsigned int n=0; n<n_nodes; ++n )
    {
      out << pts[3*n+0]/um << " " << pts[3*n+1]/um << " " << pts[3*n+2]/um << '\n';
    }
    out << std::endl;
  }

}


void VTKIO::cells_to_vtk(const MeshBase& mesh,  std::ofstream & out)
{
  int cell_size = mesh.n_elem();

  int cell_nodes = 0;
  std::map<unsigned int, unsigned int> elem_type;
  MeshBase::const_element_iterator       elem_it  = mesh.active_this_pid_elements_begin();
  const MeshBase::const_element_iterator elem_it_end = mesh.active_this_pid_elements_end();
  for (; elem_it != elem_it_end; ++elem_it)
  {
    const Elem * elem = *elem_it;
    cell_nodes += elem->n_nodes();
    elem_type.insert( std::make_pair( elem->id(), static_cast<unsigned int>(elem->type())) );
  }
  Parallel::sum(cell_nodes);
  Parallel::gather(0, elem_type);

  std::vector<int> conn;
  mesh.pack_elems(conn);


  // only do it at root node
  if(Genius::processor_id() == 0)
  {
    out<<"CELLS "<<mesh.n_elem()<<" "<<cell_size+cell_nodes<<'\n';

    unsigned int cnt = 0;
    while (cnt < conn.size())
    {
      // Unpack the element header
      const ElemType elem_type    = static_cast<ElemType>(conn[cnt]);
      AutoPtr<Elem> elem = Elem::build(elem_type);
#ifdef ENABLE_AMR
      cnt += 10;
#else
      cnt += 4;
#endif
      out << elem->n_nodes();
      for(unsigned int i=0;i<elem->n_nodes();++i)
        out << " " << conn[cnt++];

      out << std::endl;

      cnt += elem->n_sides();
    }
    out << std::endl;

    out << "CELL_TYPES " << mesh.n_elem() << '\n';

    std::map<unsigned int, unsigned int>::const_iterator it = elem_type.begin();
    for (; it != elem_type.end(); ++it)
    {
      unsigned int celltype;
      const ElemType type    = static_cast<ElemType>(it->second);
      switch(type)
      {
          case EDGE2:
          case EDGE2_FVM:
          celltype = 3;  //VTK_LINE;
          break;
          case EDGE3:
          celltype = 21; //VTK_QUADRATIC_EDGE;
          break;// 1
          case TRI3:
          case TRI3_FVM:
          case TRI3_CY_FVM:
          celltype = 5;  //VTK_TRIANGLE;
          break;// 3
          case TRI6:
          celltype = 22; //VTK_QUADRATIC_TRIANGLE;
          break;// 4
          case QUAD4:
          case QUAD4_FVM:
          case QUAD4_CY_FVM:
          celltype = 9;  //VTK_QUAD;
          break;// 5
          case QUAD8:
          celltype = 23; //VTK_QUADRATIC_QUAD;
          break;// 6
          case TET4:
          case TET4_FVM:
          celltype = 10; //VTK_TETRA;
          break;// 8
          case TET10:
          celltype = 24; //VTK_QUADRATIC_TETRA;
          break;// 9
          case HEX8:
          case HEX8_FVM:
          celltype = 12; //VTK_HEXAHEDRON;
          break;// 10
          case HEX20:
          celltype = 25; //VTK_QUADRATIC_HEXAHEDRON;
          break;// 12
          case PRISM6:
          case PRISM6_FVM:
          celltype = 13; //VTK_WEDGE;
          break;// 13
          case PYRAMID5:
          case PYRAMID5_FVM:
          celltype = 14; //VTK_PYRAMID;
          break;// 16
          default:
          {
            std::cerr<<"element type "<<type<<" not implemented"<<std::endl;
            genius_error();
          }
      }

      out << celltype << std::endl;
    }

    out << std::endl;
  }
}


void VTKIO::meshinfo_cell_to_vtk(const MeshBase& mesh, std::ofstream & out)
{

  std::map<unsigned int, unsigned int> cell_processor_id;
  std::map<unsigned int, unsigned int> cell_subdomain_id;

  MeshBase::const_element_iterator       elem_it  = mesh.active_this_pid_elements_begin();
  const MeshBase::const_element_iterator elem_it_end = mesh.active_this_pid_elements_end();
  for (; elem_it != elem_it_end; ++elem_it)
  {
    const Elem * elem = *elem_it;
    cell_processor_id.insert( std::make_pair(elem->id(), elem->processor_id()) );
    cell_subdomain_id.insert( std::make_pair(elem->id(), elem->subdomain_id()) );
  }

  Parallel::gather(0, cell_processor_id);
  Parallel::gather(0, cell_subdomain_id);

  // only do it at root node
  if(Genius::processor_id() == 0)
  {
    //write cell based region and partition info to vtk

    // region information
    out << "SCALARS region float 1"  <<'\n';
    out << "LOOKUP_TABLE default"    <<'\n';
    {
      std::map<unsigned int, unsigned int>::const_iterator it = cell_subdomain_id.begin();
      for ( ; it != cell_subdomain_id.end(); ++it)
        out << static_cast<float>(it->second)<<'\n';
    }
    out<<std::endl;

    // partition information
    out<<"SCALARS partition float 1"  <<'\n';
    out<<"LOOKUP_TABLE default"       <<'\n';
    {
      std::map<unsigned int, unsigned int>::const_iterator it = cell_processor_id.begin();
      for ( ; it != cell_processor_id.end(); ++it)
        out << static_cast<float>(it->second)<<'\n';
    }
    out<<std::endl;
  }

}




void VTKIO::meshinfo_node_to_vtk(const MeshBase& mesh, std::ofstream & out)
{
  std::set<unsigned int> node_id;
  std::map<unsigned int, unsigned int> node_processor_id;

  MeshBase::const_node_iterator nd = mesh.this_pid_nodes_begin();
  MeshBase::const_node_iterator nd_end = mesh.this_pid_nodes_end();
  for (; nd != nd_end; ++nd)
  {
    const Node * node = *nd;
    node_id.insert(node->id());
    node_processor_id.insert( std::make_pair(node->id(), node->processor_id()) );
  }

  Parallel::gather(0, node_id);
  Parallel::gather(0, node_processor_id);


  // only do it at root node
  if(Genius::processor_id() == 0)
  {
    out<<"SCALARS node_order float 1" <<'\n';
    out<<"LOOKUP_TABLE default"       <<'\n';
    {
      std::set<unsigned int>::const_iterator it =  node_id.begin();
      for (; it !=  node_id.end(); ++it)
        out<<static_cast<float>(*it)<<'\n' ;
    }
    out<<std::endl;


    out<<"SCALARS node_partition float 1" <<'\n';
    out<<"LOOKUP_TABLE default"       <<'\n';
    {
      std::map<unsigned int, unsigned int>::const_iterator it = node_processor_id.begin();
      for (; it != node_processor_id.end(); ++it)
        out<<static_cast<float>(it->second)<<'\n' ;
    }
    out<<std::endl;
  }
}


void VTKIO::solution_to_vtk(const MeshBase& mesh, std::ofstream & out)
{
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();

  // check which data should be export

  bool semiconductor_material = false;
  for(unsigned int n=0; n<system.n_regions(); ++n)
  {
    const SimulationRegion * region = system.region(n);
    if(Material::IsSemiconductor(region->material()))
      semiconductor_material=true;
  }


  bool ebm_solution_data     = false;
  bool optical_generation    = false;
  bool particle_generation   = false;
  bool ddm_ac_data           = false;
  bool optical_complex_field = false;

  const std::vector<SolverSpecify::SolverType> & solve_history = system.solve_history();

  for(unsigned int n=0; n<solve_history.size(); ++n)
  {
    switch  (solve_history[n])
    {
        case SolverSpecify::EBML3      :
        case SolverSpecify::EBML3MIX   :
        ebm_solution_data = true;
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
        default : break;
    }
  }

  if(system.get_field_source()->is_light_source_exist())
    optical_generation = true;

  if(system.get_field_source()->is_particle_source_exist())
    particle_generation = true;


  // write cell based data
  {
    const unsigned int n_elems = mesh.n_elem();
    out << "CELL_DATA "<<n_elems     <<'\n';

    meshinfo_cell_to_vtk(mesh, out);

    std::vector<unsigned int> order;

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

        order.push_back(elem->id());

        const FVM_CellData * elem_data = region->get_region_elem_data(n);

        Ex.push_back(static_cast<float>(elem_data->E()(0)/(V/cm)));
        Ey.push_back(static_cast<float>(elem_data->E()(1)/(V/cm)));
        Ez.push_back(static_cast<float>(elem_data->E()(2)/(V/cm)));

        Jnx.push_back(static_cast<float>(elem_data->Jn()(0)/(A/cm)));
        Jny.push_back(static_cast<float>(elem_data->Jn()(1)/(A/cm)));
        Jnz.push_back(static_cast<float>(elem_data->Jn()(2)/(A/cm)));

        Jpx.push_back(static_cast<float>(elem_data->Jp()(0)/(A/cm)));
        Jpy.push_back(static_cast<float>(elem_data->Jp()(1)/(A/cm)));
        Jpz.push_back(static_cast<float>(elem_data->Jp()(2)/(A/cm)));
      }
    }
    Parallel::gather(0, order);
    write_cell_vector_solution(order, Ex,  Ey,  Ez,  "electrical_field", out);
    write_cell_vector_solution(order, Jnx, Jny, Jnz, "electron_current", out);
    write_cell_vector_solution(order, Jpx, Jpy, Jpz, "hole_current", out);
  }


  // write node based data
  {
    //write node based boundary info to vtk
    const unsigned int n_nodes = mesh.n_nodes();
    out<<"POINT_DATA "<<n_nodes      <<'\n';

    meshinfo_node_to_vtk(mesh, out);

    // scale back to normal unit
    double concentration_scale = pow(cm, -3);

    std::vector<unsigned int> order;

    // gather data from all the processor
    std::vector<float> Na, Nd, net_doping, net_charge, psi, Ec, Ev, qFn, qFp, n, p, T, Tn, Tp; //scalar value
    std::vector<float> mole_x, mole_y;
    std::vector<float> OptG, PatG;
    std::vector<float> R;

    std::vector<std::complex<float> >psi_ac, n_ac, p_ac, T_ac, Tn_ac, Tp_ac; //complex value
    std::vector<std::complex<float> > OptE_complex, OptH_complex;

    std::vector<float> Jnx, Jny, Jnz; //vector value
    std::vector<float> Jpx, Jpy, Jpz;
    std::vector<float> Vnx, Vny, Vnz;
    std::vector<float> Vpx, Vpy, Vpz;
    std::vector<float> Ex,  Ey,  Ez;

    // search all the node belongs to current processor
    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      SimulationRegion::const_processor_node_iterator on_processor_nodes_it = region->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator on_processor_nodes_it_end = region->on_processor_nodes_end();
      for(; on_processor_nodes_it!=on_processor_nodes_it_end; ++on_processor_nodes_it)
      {
        const FVM_Node * fvm_node = *on_processor_nodes_it;
        const FVM_NodeData * node_data = fvm_node->node_data();
        const unsigned int id = fvm_node->root_node()->id();

        order.push_back(id);

        // if the fvm_node lies on the interface of two material regions,
        // we shall use the node data in the more important region.
        // it just for visualization reason
        if( fvm_node->boundary_id() != BoundaryInfo::invalid_id )
        {
          unsigned int bc_index = system.get_bcs()->get_bc_index_by_bd_id(fvm_node->boundary_id());
          const BoundaryCondition * bc = system.get_bcs()->get_bc(bc_index);
          const FVM_Node * primary_fvm_node = (*bc->region_node_begin(fvm_node->root_node())).second.second;
          node_data = primary_fvm_node->node_data();
        }

        assert ( node_data != NULL );
        {
          if(semiconductor_material)
          {
            Na.push_back(static_cast<float>(node_data->Total_Na()/concentration_scale));
            Nd.push_back(static_cast<float>(node_data->Total_Nd()/concentration_scale));
            net_doping.push_back(static_cast<float>(node_data->Net_doping()/concentration_scale));
            net_charge.push_back(static_cast<float>(node_data->Net_charge()/concentration_scale));
            n.push_back(static_cast<float>(node_data->n()/concentration_scale));
            p.push_back(static_cast<float>(node_data->p()/concentration_scale));
            mole_x.push_back(static_cast<float>(node_data->mole_x()));
            mole_y.push_back(static_cast<float>(node_data->mole_y()));
            R.push_back(static_cast<float>(node_data->Recomb()/(concentration_scale/s)));

            Jnx.push_back(static_cast<float>(node_data->Jn()(0)/(A/cm)));
            Jny.push_back(static_cast<float>(node_data->Jn()(1)/(A/cm)));
            Jnz.push_back(static_cast<float>(node_data->Jn()(2)/(A/cm)));

            Jpx.push_back(static_cast<float>(node_data->Jp()(0)/(A/cm)));
            Jpy.push_back(static_cast<float>(node_data->Jp()(1)/(A/cm)));
            Jpz.push_back(static_cast<float>(node_data->Jp()(2)/(A/cm)));

            Ex.push_back(static_cast<float>(node_data->E()(0)/(V/cm)));
            Ey.push_back(static_cast<float>(node_data->E()(1)/(V/cm)));
            Ez.push_back(static_cast<float>(node_data->E()(2)/(V/cm)));
          }

          psi.push_back(static_cast<float>(node_data->psi()/V));
          Ec.push_back(static_cast<float>(node_data->Ec()/eV));
          Ev.push_back(static_cast<float>(node_data->Ev()/eV));
          qFn.push_back(static_cast<float>(node_data->qFn()/eV));
          qFp.push_back(static_cast<float>(node_data->qFp()/eV));
          T.push_back(static_cast<float>(node_data->T()/K));

          if(ebm_solution_data)
          {
            Tn.push_back(static_cast<float>(node_data->Tn()/K));
            Tp.push_back(static_cast<float>(node_data->Tp()/K));
          }

          if(optical_generation)
            OptG.push_back(static_cast<float>(node_data->OptG()/(concentration_scale/s)));

          if(particle_generation)
            PatG.push_back(static_cast<float>(node_data->PatG()/(concentration_scale/s)));

          if(ddm_ac_data)
          {
            psi_ac.push_back(static_cast<std::complex<float> >(node_data->psi_ac()/V));
            n_ac.push_back(static_cast<std::complex<float> >(node_data->n_ac()/concentration_scale));
            p_ac.push_back(static_cast<std::complex<float> >(node_data->p_ac()/concentration_scale));
            T_ac.push_back(static_cast<std::complex<float> >(node_data->T_ac()/K));
            Tn_ac.push_back(static_cast<std::complex<float> >(node_data->Tn_ac()/K));
            Tp_ac.push_back(static_cast<std::complex<float> >(node_data->Tp_ac()/K));
          }

          if(optical_complex_field)
          {
            OptE_complex.push_back(static_cast<std::complex<float> >(node_data->OptE_complex()/(V/cm)));
            OptH_complex.push_back(static_cast<std::complex<float> >(node_data->OptH_complex()/(A/cm)));
          }

        }
      }
    }

    Parallel::gather(0, order);

    write_node_scaler_solution(order, psi, "psi", out);
    write_node_scaler_solution(order, Ec,  "Ec",  out);
    write_node_scaler_solution(order, Ev,  "Ev",  out);
    write_node_scaler_solution(order, qFn, "elec_quasi_Fermi_level", out);
    write_node_scaler_solution(order, qFp, "hole_quasi_Fermi_level", out);

    if(semiconductor_material)
    {
      write_node_scaler_solution(order, Na,  "Na", out);
      write_node_scaler_solution(order, Nd,  "Nd", out);
      write_node_scaler_solution(order, n,   "electron_density", out);
      write_node_scaler_solution(order, p,   "hole_density", out);
      write_node_scaler_solution(order, net_doping,  "net_doping", out);
      write_node_scaler_solution(order, net_charge,  "net_charge", out);
      write_node_scaler_solution(order, mole_x, "mole_x", out);
      write_node_scaler_solution(order, mole_y, "mole_y", out);
      write_node_scaler_solution(order, R, "recombination", out);
      //write_node_vector_solution(order, Jnx, Jny, Jnz, "Jn", out);
      //write_node_vector_solution(order, Jpx, Jpy, Jpz, "Jp", out);
    }

    write_node_scaler_solution(order, T,   "temperature", out);

    if(ebm_solution_data)
    {
      write_node_scaler_solution(order, Tn,  "elec_temperature", out);
      write_node_scaler_solution(order, Tp,  "hole_temperature", out);
    }

    if(ddm_ac_data)
    {
      write_node_complex_solution(order, psi_ac, "psi_AC", out);
      write_node_complex_solution(order, n_ac,   "electron_AC", out);
      write_node_complex_solution(order, p_ac,   "hole_AC", out);
      write_node_complex_solution(order, T_ac,   "temperature_AC", out);
      write_node_complex_solution(order, Tn_ac,  "elec_temperature_AC", out);
      write_node_complex_solution(order, Tp_ac,  "hole_temperature_AC", out);
    }

    if(optical_complex_field)
    {
      write_node_complex_solution(order, OptE_complex, "Optical_E", out);
      write_node_complex_solution(order, OptH_complex, "Optical_H", out);
    }

    if(optical_generation)
      write_node_scaler_solution(order, OptG,   "Optical_Generation", out);

    if(particle_generation)
      write_node_scaler_solution(order, PatG,   "Radiation_Generation", out);
  }

}

void VTKIO::write_node_scaler_solution(const std::vector<unsigned int> & order, std::vector<float> & sol,
                                       const std::string & sol_name, std::ofstream & out)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);


  if ( Genius::processor_id() == 0)
  {
    out << "SCALARS "<<sol_name<<" float 1" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;

    std::vector<float> re_order(order.size());
    for(unsigned int n=0; n<order.size(); ++n)
      re_order[order[n]] = sol[n];
    for(unsigned int n=0; n<re_order.size(); ++n)
    {
      out << re_order[n] << std::endl;
    }
    out<<std::endl;
  }
}


void VTKIO::write_node_complex_solution(const std::vector<unsigned int> & order, std::vector<std::complex<float> > & sol,
                                        const std::string & sol_name, std::ofstream & out)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if ( Genius::processor_id() == 0)
  {
    // out put magnitude
    out << "SCALARS "<<sol_name<<"_abs float 1" << std::endl;
    out << "LOOKUP_TABLE default"   << std::endl;
    std::vector<float> re_order_abs(order.size());
    for(unsigned int n=0; n<order.size(); ++n)
      re_order_abs[order[n]] = std::abs(sol[n]);
    for(unsigned int n=0; n<re_order_abs.size(); ++n)
    {
      out << re_order_abs[n] << std::endl;
    }
    out<<std::endl;


    // out put phase angle
    out << "SCALARS "<<sol_name<<"_angle float 1" << std::endl;
    out << "LOOKUP_TABLE default"   << std::endl;
    std::vector<float> re_order_arg(order.size());
    for(unsigned int n=0; n<order.size(); ++n)
      re_order_arg[order[n]] = std::arg(sol[n]);
    for(unsigned int n=0; n<re_order_arg.size(); ++n)
    {
      out << re_order_arg[n] << std::endl;
    }
    out<<std::endl;
  }

}


void VTKIO::write_node_vector_solution(const std::vector<unsigned int> & order,
                                       std::vector<float > & sol_x,
                                       std::vector<float > & sol_y,
                                       std::vector<float > & sol_z,
                                       const std::string & sol_name, std::ofstream & out)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol_x);
  Parallel::gather(0, sol_y);
  Parallel::gather(0, sol_z);

  if ( Genius::processor_id() == 0)
  {
    out << "VECTORS "<<sol_name<<" float" << std::endl;

    std::vector<float> re_order_x(order.size());
    std::vector<float> re_order_y(order.size());
    std::vector<float> re_order_z(order.size());
    for(unsigned int n=0; n<order.size(); ++n)
    {
      re_order_x[order[n]] = sol_x[n];
      re_order_y[order[n]] = sol_y[n];
      re_order_z[order[n]] = sol_z[n];
    }

    for( unsigned int i=0; i<sol_x.size(); i++)
    {
      out << re_order_x[i] <<" "<< re_order_y[i] <<" "<< re_order_z[i] <<std::endl;
    }
    out<<std::endl;
  }
}


void VTKIO::write_cell_scaler_solution(const std::vector<unsigned int> & order, std::vector<float> & sol,
                                       const std::string & sol_name, std::ofstream & out)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if ( Genius::processor_id() == 0)
  {
    out << "SCALARS "<<sol_name<<" float 1" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    std::vector<float> re_order(order.size());
    for(unsigned int n=0; n<order.size(); ++n)
      re_order[order[n]] = sol[n];
    for(unsigned int n=0; n<re_order.size(); ++n)
    {
      out << re_order[n] << std::endl;
    }
    out<<std::endl;
  }
}

void VTKIO::write_cell_vector_solution(const std::vector<unsigned int> & order,
                                       std::vector<float > & sol_x,
                                       std::vector<float > & sol_y,
                                       std::vector<float > & sol_z,
                                       const std::string & sol_name, std::ofstream & out)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol_x);
  Parallel::gather(0, sol_y);
  Parallel::gather(0, sol_z);

  if ( Genius::processor_id() == 0)
  {
    out << "VECTORS "<<sol_name<<" float" << std::endl;

    std::vector<float> re_order_x(order.size());
    std::vector<float> re_order_y(order.size());
    std::vector<float> re_order_z(order.size());
    for(unsigned int n=0; n<order.size(); ++n)
    {
      re_order_x[order[n]] = sol_x[n];
      re_order_y[order[n]] = sol_y[n];
      re_order_z[order[n]] = sol_z[n];
    }
    for( unsigned int i=0; i<sol_x.size(); i++)
    {
      out << re_order_x[i] <<" "<< re_order_y[i] <<" "<< re_order_z[i] <<std::endl;
    }
    out<<std::endl;
  }
}



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

  // vtk file extension have a ".vtk" format?
  if(name.rfind(".vtk") < name.size())
  {
    // write mesh and solution as txt file
    if(Genius::processor_id() == 0)
    {
      _out.open(name.c_str(), std::ofstream::trunc);
    }

    nodes_to_vtk    (mesh, _out);
    cells_to_vtk    (mesh, _out);
    solution_to_vtk (mesh, _out);

    if(Genius::processor_id()==0)
    {
      _out.close();
    }
  }

}






void VTKIO::read (const std::string& name)
{
#ifdef HAVE_VTK

  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  // clear the system
  system.clear();

  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  if(Genius::processor_id() == 0)
  {
    vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
    reader->SetFileName( name.c_str() );

    // Force reading
    reader->Update();

    // read in the grid
    //  vtkUnstructuredGrid *grid = reader->GetOutput();
    _vtk_grid = reader->GetOutput();
    _vtk_grid->Update();

    // always numbered nicely??, so we can loop like this
    // I'm pretty sure it is numbered nicely
    for (unsigned int i=0; i < static_cast<unsigned int>(_vtk_grid->GetNumberOfPoints()); ++i)
    {
      // add to the id map
      // and add the actual point
      double * pnt = _vtk_grid->GetPoint(static_cast<vtkIdType>(i));
      Point xyz(pnt[0],pnt[1],pnt[2]);
      Node* newnode = mesh.add_point(xyz,i);
    }

    // read cells
    for (unsigned int i=0; i < static_cast<unsigned int>(_vtk_grid->GetNumberOfCells()); ++i)
    {
      vtkCell* cell = _vtk_grid->GetCell(i);
      Elem* elem = NULL;  // Initialize to avoid compiler warning
      switch(cell->GetCellType())
      {
          case VTK_TRIANGLE:
          elem = Elem::build(TRI3).release();
          break;
          case VTK_QUAD:
          elem = Elem::build(QUAD4).release();
          break;
          case VTK_TETRA:
          elem = Elem::build(TET4).release();
          break;
          case VTK_WEDGE:
          elem = Elem::build(PRISM6).release();
          break;
          case VTK_HEXAHEDRON:
          elem = Elem::build(HEX8).release();
          break;
          case VTK_PYRAMID:
          elem = Elem::build(PYRAMID5).release();
          break;
          default:
          std::cerr << "element type not implemented in vtkinterface " << cell->GetCellType() << std::endl;
          genius_error();
      }

      // get the straightforward numbering from the VTK cells
      for(unsigned int j=0;j<elem->n_nodes();++j)
      {
        elem->set_node(j) = mesh.node_ptr(cell->GetPointId(j));
      }
      // then get the connectivity
      std::vector<unsigned int> conn;
      elem->connectivity(0,VTK,conn);
      // then reshuffle the nodes according to the connectivity, this
      // two-time-assign would evade the definition of the vtk_mapping
      for(unsigned int j=0;j<conn.size();++j)
      {
        elem->set_node(j) = mesh.node_ptr(conn[j]);
      }

      elem->set_id(i);
      mesh.add_elem(elem);
    } // end loop over VTK cells

    // 2d/3d mesh?
    {
      vtkCell* cell = _vtk_grid->GetCell(0);
      switch(cell->GetCellType())
      {
          case VTK_TRIANGLE:
          case VTK_QUAD:
          mesh.magic_num() = 237;
          break;
          case VTK_TETRA:
          case VTK_WEDGE:
          case VTK_HEXAHEDRON:
          case VTK_PYRAMID:
          mesh.magic_num() = 72371;
          break;
          default:
          std::cerr << "element type not implemented in vtkinterface " << cell->GetCellType() << std::endl;
          genius_error();
      }
    }


    // read cell/point data
    vtkCellData * cell_data = _vtk_grid->GetCellData();
    vtkPointData * point_data = _vtk_grid->GetPointData();

    // process region/material information
    {
      vtkFloatArray * region_info   = dynamic_cast<vtkFloatArray *>(cell_data->GetScalars("region"));
      vtkFloatArray * material_info = dynamic_cast<vtkFloatArray *>(cell_data->GetScalars("material"));

      std::map<unsigned int, unsigned int> sub_id_to_material_id;
      for (unsigned int i=0; i < static_cast<unsigned int>(_vtk_grid->GetNumberOfCells()); ++i)
      {
        unsigned int sub_id      = static_cast<unsigned int>(region_info->GetValue(i)+0.5);
        unsigned int material_id = static_cast<unsigned int>(material_info->GetValue(i)+0.5);

        Elem* elem = mesh.elem(i);
        elem->subdomain_id() = sub_id;
        sub_id_to_material_id[sub_id] = material_id;
      }

      mesh.set_n_subdomains() = sub_id_to_material_id.size();

      std::map<unsigned int, unsigned int>::iterator it = sub_id_to_material_id.begin();
      for(; it!=sub_id_to_material_id.end(); ++it)
      {
        std::stringstream ss;
        std::string index;
        ss<<it->first;
        ss>>index;

        mesh.set_subdomain_label( it->first, "region"+index );
        mesh.set_subdomain_material( it->first, Material::get_material_by_id(it->second));
      }
    }

    {
      //NOTE: we havn't process all the interface and region external boundaries
      //build neighbor information for mesh. then elem->neighbor() is functional
      mesh.find_neighbors();

      std::map<std::string, int> bd_map;
      typedef std::map<std::string, int>::iterator Bd_It;

      const MeshBase::element_iterator el_end = mesh.elements_end();
      for (MeshBase::element_iterator el = mesh.elements_begin(); el != el_end; ++el)
      {
        const Elem* elem = *el;
        for (unsigned int ms=0; ms<elem->n_neighbors(); ms++)
        {
          //when this side is already exist, skip it
          if( mesh.boundary_info->boundary_id (elem, ms) != BoundaryInfo::invalid_id) continue;

          // region outer boundary
          if(elem->neighbor(ms)==NULL)
          {
            std::string bd_label = mesh.subdomain_label_by_id(elem->subdomain_id()) + "_Neumann" ;
            short int bd_id;
            if(bd_map.find(bd_label)!=bd_map.end())
              bd_id = bd_map[bd_label];
            else
            {
              bd_id = bd_map.size();
              bd_map[bd_label] = bd_id;
            }
            mesh.boundary_info->add_side(elem, ms, bd_id);
          }
          // region interface
          else
          {
            const Elem* neighbor = elem->neighbor(ms);
            // the element and its neighbor in diffetent subdomain?
            unsigned int sbd_id1 = elem->subdomain_id();
            unsigned int sbd_id2 = neighbor->subdomain_id();

            if (sbd_id1 == sbd_id2) continue;

            // build the label for the interface, which has the form of RegionLabel1_to_RegionLabel2,
            // the two region is alpha ordered.
            std::string bd_label;
            std::string region1 = mesh.subdomain_label_by_id(sbd_id1);
            std::string region2 = mesh.subdomain_label_by_id(sbd_id2);
            if( region1 < region2)
              bd_label = region1 + "_to_" + region2;
            else
              bd_label = region2 + "_to_" + region1;

            short int bd_id;
            if(bd_map.find(bd_label)!=bd_map.end())
              bd_id = bd_map[bd_label];
            else
            {
              bd_id = bd_map.size();
              bd_map[bd_label] = bd_id;
            }
            mesh.boundary_info->add_side(elem, ms, bd_id);
          }
        }
      }


      // after the previous remove and insert operation, the number of boundary id maybe changed. renumber it here.
      mesh.boundary_info->rebuild_ids();

      //write down bd labels
      Bd_It bd_it = bd_map.begin();
      for(; bd_it != bd_map.end(); ++bd_it)
      {
        mesh.boundary_info->set_label_to_id( (*bd_it).second, (*bd_it).first );
      }
    }

    reader->Delete();

  }

  // broadcast mesh to all the processor
  MeshCommunication mesh_comm;
  mesh_comm.broadcast(mesh);

  // build simulation system
  system.build_simulation_system();
  system.sync_print_info();

#endif // HAVE_VTK

}



