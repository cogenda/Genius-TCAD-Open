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
#include "mesh.h"
#include "boundary_info.h"
#include "parallel.h"
#include "mesh_communication.h"
#include "simulation_system.h"
#include "boundary_condition_collector.h"
#include "field_source.h"

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

#ifdef HAVE_VTK
#include "vtkObjectFactory.h"

class XMLUnstructuredGridWriter : public vtkXMLUnstructuredGridWriter
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
  static XMLUnstructuredGridWriter* New();

  void setExtraHeader(const std::string &header)
  {
    _header = header;
  }
private:
  std::string _header;
};

vtkStandardNewMacro(XMLUnstructuredGridWriter)

#endif

// private functions
#ifdef HAVE_VTK

void VTKIO::nodes_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* &grid)
{
  // only do it at root node
  if(Genius::processor_id() != 0) return;

  vtkPoints* points = vtkPoints::New();
  // FIXME change to local iterators when I've figured out how to write in parallel
  //  MeshBase::const_node_iterator nd = mesh.local_nodes_begin();
  //  MeshBase::const_node_iterator nd_end = mesh.local_nodes_end();
  MeshBase::const_node_iterator nd = mesh.active_nodes_begin();
  MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();
  // write out nodal data
  for ( ; nd!=nd_end; ++nd )
  {
    float tuple[3]; //, dsp[3];
    const Node * node = (*nd);

    {
      tuple[0] =  (*node)(0)/um; //scale to um
      tuple[1] =  (*node)(1)/um; //scale to um
      tuple[2] =  (*node)(2)/um; //scale to um
    }

    points->InsertPoint(node->id(),tuple);

  }

  grid->SetPoints(points);

  //clean up
  points->Delete();

}


void VTKIO::cells_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* &grid)
{

  // only do it at root node
  if(Genius::processor_id() != 0) return;

  //  MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
  //  const MeshBase::const_element_iterator end = mesh.active_local_elements_end();
  MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();
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

    vtkIdList *pts = vtkIdList::New();
    pts->SetNumberOfIds(elem->n_nodes());
    // get the connectivity for this element
    std::vector<unsigned int> conn;
    elem->connectivity(0,VTK,conn);
    for(unsigned int i=0;i<conn.size();++i)
    {
      pts->SetId(i,conn[i]);
    }
    grid->InsertNextCell(celltype,pts);

    pts->Delete();
  } // end loop over active elements


  // also export boundary elems (DIM-1 elem)
  {
    std::vector<unsigned int>       el;
    std::vector<unsigned short int> sl;
    std::vector<short int>          il;
    mesh.boundary_info->build_side_list (el, sl, il);
    for(unsigned int n=0; n<il.size(); ++n)
    {
      const Elem * elem = mesh.elem(el[n]);
      AutoPtr<Elem> boundary_elem =  elem->build_side(sl[n], false);

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

      vtkIdList *pts = vtkIdList::New();
      pts->SetNumberOfIds(boundary_elem->n_nodes());
      // get the connectivity for this element
      std::vector<unsigned int> conn;
      boundary_elem->connectivity(0,VTK,conn);
      for(unsigned int i=0;i<conn.size();++i)
      {
        pts->SetId(i,conn[i]);
      }
      grid->InsertNextCell(celltype, pts);

      pts->Delete();
    }
  }
}


void VTKIO::meshinfo_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* &grid)
{

  // only do it at root node
  if(Genius::processor_id() != 0) return;

  //write cell based region and partition info to vtk

  const unsigned int n_elems = mesh.n_active_elem() + mesh.boundary_info->n_boundary_conds();

  vtkIntArray *region_info    = vtkIntArray::New();
  vtkIntArray *boundary_info  = vtkIntArray::New();
  vtkIntArray *material_info  = vtkIntArray::New();
  vtkIntArray *partition_info = vtkIntArray::New();

  region_info->SetName("region");
  region_info->SetNumberOfValues(n_elems);

  boundary_info->SetName("boundary");
  boundary_info->SetNumberOfValues(n_elems);

  material_info->SetName("material");
  material_info->SetNumberOfValues(n_elems);

  std::map<unsigned int, unsigned int> sub_id_to_material_id;
  for( unsigned int n=0; n<mesh.n_subdomains(); n++)
  {
    const std::string &material = mesh.subdomain_material(n);
    sub_id_to_material_id[n] = Material::get_id_by_material(material);
  }

  partition_info->SetName("partition");
  partition_info->SetNumberOfValues(n_elems);

  // write info for mesh elements
  {
    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();
    for (int i=0; it != end; ++it, ++i)
    {
      unsigned int sub_id = (*it)->subdomain_id();
      region_info->InsertValue(i, sub_id);
      boundary_info->InsertValue(i, 0);
      material_info->InsertValue(i, sub_id_to_material_id[sub_id]);
      partition_info->InsertValue(i, (*it)->processor_id());
    }
  }

  // write info for boundary elements
  {
    std::vector<unsigned int>       el;
    std::vector<unsigned short int> sl;
    std::vector<short int>          il;
    unsigned int id = mesh.n_active_elem();
    mesh.boundary_info->build_side_list (el, sl, il);
    for(unsigned int n=0; n<il.size(); ++n, ++id)
    {
      const Elem * elem = mesh.elem(el[n]);
      AutoPtr<Elem> boundary_elem =  elem->build_side(sl[n], false);

      unsigned int sub_id = elem->subdomain_id();
      region_info->InsertValue(id, sub_id);
      boundary_info->InsertValue(id, il[n]);
      material_info->InsertValue(id, sub_id_to_material_id[sub_id]);
      partition_info->InsertValue(id, elem->processor_id());
    }
  }

  grid->GetCellData()->AddArray(region_info);
  grid->GetCellData()->AddArray(boundary_info);
  grid->GetCellData()->AddArray(material_info);
  grid->GetCellData()->AddArray(partition_info);

  region_info->Delete();
  boundary_info->Delete();
  material_info->Delete();
  partition_info->Delete();

  //write node based boundary info to vtk, does this really necessary?
#if 0
  const unsigned int n_nodes = mesh.n_nodes();

  MeshBase::const_node_iterator nd = mesh.active_nodes_begin();
  MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();

  vtkFloatArray *boundary_info = vtkFloatArray::New();
  boundary_info->SetName("boundary");
  boundary_info->SetNumberOfValues(n_nodes);

  vtkFloatArray *node_order = vtkFloatArray::New();
  node_order->SetName("node order");
  node_order->SetNumberOfValues(n_nodes);

  vtkFloatArray *node_partition = vtkFloatArray::New();
  node_partition->SetName("node partition");
  node_partition->SetNumberOfValues(n_nodes);

  for (int i=0 ; nd != nd_end; ++i, ++nd)
  {
    if (mesh.boundary_info->boundary_id(*nd) != BoundaryInfo::invalid_id )
    {
      boundary_info->InsertValue((*nd)->id(), static_cast<float>(mesh.boundary_info->boundary_id(*nd)));
    }
    else
      boundary_info->InsertValue((*nd)->id(), 0.0);

    node_order->InsertValue((*nd)->id(), static_cast<float>((*nd)->id()) );

    node_partition->InsertValue((*nd)->id(), static_cast<float>((*nd)->processor_id()) );
  }

  grid->GetPointData()->AddArray(boundary_info);
  grid->GetPointData()->AddArray(node_order);
  grid->GetPointData()->AddArray(node_partition);

  boundary_info->Delete();
  node_order->Delete();
  node_partition->Delete();
#endif

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
      os << "    <boundary id=\"" << bcs->get_bd_id_by_bc_index(b) << "\" "
      << "name=\"" << bc->label() << "\" "
      << "bc=\"" << BC_enum_to_string(bc->bc_type()) << "\" "
      << "/>"  << std::endl;
    }
    os << "  </boundaries>" << std::endl;
  }

  os << "</Genius>" << std::endl;
  return os.str();
}


void VTKIO::solution_to_vtk(const MeshBase& mesh, vtkUnstructuredGrid* &grid)
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

  const std::vector<SolverSpecify::SolverType> & solve_histroy = system.solve_histroy();

  for(unsigned int n=0; n<solve_histroy.size(); ++n)
  {
    switch  (solve_histroy[n])
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

    // gather data from all the processor
    std::map<unsigned int, float> Na, Nd, net_doping, net_charge, psi, Ec, Ev, qFn, qFp, n, p, T, Tn, Tp; //scalar value
    std::map<unsigned int, float> mole_x, mole_y;
    std::map<unsigned int, float> OptG, PatG;

    std::map<unsigned int, std::complex<float> >psi_ac, n_ac, p_ac, T_ac, Tn_ac, Tp_ac; //complex value
    std::map<unsigned int, std::complex<float> > OptE_complex, OptH_complex;

    std::map<unsigned int, float> Jnx, Jny, Jnz; //vector value
    std::map<unsigned int, float> Jpx, Jpy, Jpz;

    // search all the node belongs to current processor
    MeshBase::const_node_iterator nd = mesh.pid_nodes_begin(Genius::processor_id());
    MeshBase::const_node_iterator nd_end = mesh.pid_nodes_end(Genius::processor_id());
    for (; nd != nd_end; ++nd)
    {
      if( (*nd)->processor_id() != Genius::processor_id() ) continue;

      for( unsigned int r=0; r<system.n_regions(); r++)
      {
        const SimulationRegion * region = system.region(r);
        const FVM_Node * fvm_node = region->region_fvm_node(*nd);
        if(fvm_node==NULL) continue;

        const FVM_NodeData * node_data = region->region_node_data(*nd);

        // if the fvm_node lies on the interface of two material regions,
        // we shall use the node data in the more important region.
        // it just for visualization reason
        if( fvm_node->boundary_id() != BoundaryInfo::invalid_id )
        {
          unsigned int bc_index = system.get_bcs()->get_bc_index_by_bd_id(fvm_node->boundary_id());
          const BoundaryCondition * bc = system.get_bcs()->get_bc(bc_index);
          const FVM_Node * primary_fvm_node = (*bc->region_node_begin(fvm_node->root_node())).second.second;
          genius_assert(primary_fvm_node->root_node()==(*nd));
          node_data = primary_fvm_node->node_data();
        }


        assert ( node_data != NULL );
        {
          if(semiconductor_material)
          {
            Na        [(*nd)->id()]  = static_cast<float>(node_data->Total_Na()/concentration_scale);
            Nd        [(*nd)->id()]  = static_cast<float>(node_data->Total_Nd()/concentration_scale);
            net_doping[(*nd)->id()]  = static_cast<float>(node_data->Net_doping()/concentration_scale);
            net_charge[(*nd)->id()]  = static_cast<float>(node_data->Net_charge()/concentration_scale);
            n         [(*nd)->id()]  = static_cast<float>(node_data->n()/concentration_scale);
            p         [(*nd)->id()]  = static_cast<float>(node_data->p()/concentration_scale);
            mole_x    [(*nd)->id()]  = static_cast<float>(node_data->mole_x());
            mole_y    [(*nd)->id()]  = static_cast<float>(node_data->mole_y());

            Jnx[(*nd)->id()]  = static_cast<float>(node_data->Jn()(0)/(A/cm));
            Jny[(*nd)->id()]  = static_cast<float>(node_data->Jn()(1)/(A/cm));
            Jnz[(*nd)->id()]  = static_cast<float>(node_data->Jn()(2)/(A/cm));

            Jpx[(*nd)->id()]  = static_cast<float>(node_data->Jp()(0)/(A/cm));
            Jpy[(*nd)->id()]  = static_cast<float>(node_data->Jp()(1)/(A/cm));
            Jpz[(*nd)->id()]  = static_cast<float>(node_data->Jp()(2)/(A/cm));
          }

          psi[(*nd)->id()]  = static_cast<float>(node_data->psi()/V);
          Ec[(*nd)->id()]   = static_cast<float>(node_data->Ec()/eV);
          Ev[(*nd)->id()]   = static_cast<float>(node_data->Ev()/eV);
          qFn[(*nd)->id()]  = static_cast<float>(node_data->qFn()/eV);
          qFp[(*nd)->id()]  = static_cast<float>(node_data->qFp()/eV);
          T  [(*nd)->id()]  = static_cast<float>(node_data->T()/K);

          if(ebm_solution_data)
          {
            Tn [(*nd)->id()]  = static_cast<float>(node_data->Tn()/K);
            Tp [(*nd)->id()]  = static_cast<float>(node_data->Tp()/K);
          }

          if(optical_generation)
            OptG[(*nd)->id()] = static_cast<float>(node_data->OptG()/(concentration_scale/s));

          if(particle_generation)
            PatG[(*nd)->id()] = static_cast<float>(node_data->PatG()/(concentration_scale/s));

          if(ddm_ac_data)
          {
            psi_ac [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->psi_ac()/V);
            n_ac   [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->n_ac()/concentration_scale);
            p_ac   [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->p_ac()/concentration_scale);
            T_ac   [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->T_ac()/K);
            Tn_ac  [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->Tn_ac()/K);
            Tp_ac  [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->Tp_ac()/K);
          }

          if(optical_complex_field)
          {
            OptE_complex[(*nd)->id()] = static_cast<std::complex<float> >(node_data->OptE_complex()/(V/cm));
            OptH_complex[(*nd)->id()] = static_cast<std::complex<float> >(node_data->OptH_complex()/(A/cm));
          }
        }
      }
    }


    write_node_scaler_solution(psi, "psi", grid);
    write_node_scaler_solution(Ec,  "Ec",  grid);
    write_node_scaler_solution(Ev,  "Ev",  grid);
    write_node_scaler_solution(qFn, "elec quasi Fermi level", grid);
    write_node_scaler_solution(qFp, "hole quasi Fermi level", grid);

    if(semiconductor_material)
    {
      write_node_scaler_solution(Na,  "Na", grid);
      write_node_scaler_solution(Nd,  "Nd", grid);
      write_node_scaler_solution(n,   "electron density", grid);
      write_node_scaler_solution(p,   "hole density", grid);
      write_node_scaler_solution(net_doping,  "net_doping", grid);
      write_node_scaler_solution(net_charge,  "net_charge", grid);
      write_node_scaler_solution(mole_x, "mole x", grid);
      write_node_scaler_solution(mole_y, "mole y", grid);
      //write_node_vector_solution(Jnx, Jny, Jnz, "Jn", grid);
      //write_node_vector_solution(Jpx, Jpy, Jpz, "Jp", grid);
    }

    write_node_scaler_solution(T,   "temperature", grid);

    if(ebm_solution_data)
    {
      write_node_scaler_solution(Tn,  "elec temperature", grid);
      write_node_scaler_solution(Tp,  "hole temperature", grid);
    }

    if(ddm_ac_data)
    {
      write_node_complex_solution(psi_ac, "psi AC", grid);
      write_node_complex_solution(n_ac,   "electron AC", grid);
      write_node_complex_solution(p_ac,   "hole AC", grid);
      write_node_complex_solution(T_ac,   "temperature AC", grid);
      write_node_complex_solution(Tn_ac,  "elec temperature AC", grid);
      write_node_complex_solution(Tp_ac,  "hole temperature AC", grid);
    }

    if(optical_complex_field)
    {
      write_node_complex_solution(OptE_complex, "Optical E", grid);
      write_node_complex_solution(OptH_complex, "Optical H", grid);
    }

    if(optical_generation)
      write_node_scaler_solution(OptG,   "Optical Generation", grid);

    if(particle_generation)
      write_node_scaler_solution(PatG,   "Radiation Generation", grid);
  }

  // write cell based data
  {
    std::map<unsigned int, float> Ex,  Ey,  Ez;
    std::map<unsigned int, float> Jnx,  Jny,  Jnz;
    std::map<unsigned int, float> Jpx,  Jpy,  Jpz;

    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      for(unsigned int n=0; n<region->n_cell(); ++n)
      {
        const Elem * elem = region->get_region_elem(n);
        const FVM_CellData * elem_data = region->get_region_elem_data(n);

        Ex[elem->id()]  = static_cast<float>(elem_data->E()(0)/(V/cm));
        Ey[elem->id()]  = static_cast<float>(elem_data->E()(1)/(V/cm));
        Ez[elem->id()]  = static_cast<float>(elem_data->E()(2)/(V/cm));

        Jnx[elem->id()]  = static_cast<float>(elem_data->Jn()(0)/(A/cm));
        Jny[elem->id()]  = static_cast<float>(elem_data->Jn()(1)/(A/cm));
        Jnz[elem->id()]  = static_cast<float>(elem_data->Jn()(2)/(A/cm));

        Jpx[elem->id()]  = static_cast<float>(elem_data->Jp()(0)/(A/cm));
        Jpy[elem->id()]  = static_cast<float>(elem_data->Jp()(1)/(A/cm));
        Jpz[elem->id()]  = static_cast<float>(elem_data->Jp()(2)/(A/cm));
      }
    }

    // write solution for boundary elements, just keep the same as mesh elems
    {
      std::vector<unsigned int>       el;
      std::vector<unsigned short int> sl;
      std::vector<short int>          il;
      unsigned int id = mesh.n_active_elem();
      mesh.boundary_info->build_side_list (el, sl, il);
      for(unsigned int n=0; n<il.size(); ++n, ++id)
      {
        const Elem * elem = mesh.elem(el[n]);

        Ex[id]  = Ex[elem->id()];
        Ey[id]  = Ey[elem->id()];
        Ez[id]  = Ez[elem->id()];

        Jnx[id]  = Jnx[elem->id()];
        Jny[id]  = Jny[elem->id()];
        Jnz[id]  = Jnz[elem->id()];

        Jpx[id]  = Jpx[elem->id()];
        Jpy[id]  = Jpy[elem->id()];
        Jpz[id]  = Jpz[elem->id()];
      }
    }

    write_cell_vector_solution(Ex,  Ey,  Ez,  "electrical field", grid);
    write_cell_vector_solution(Jnx, Jny, Jnz, "electron current", grid);
    write_cell_vector_solution(Jpx, Jpy, Jpz, "hole current", grid);
  }
}


void VTKIO::write_node_scaler_solution(std::map<unsigned int, float> & sol, const std::string & sol_name, vtkUnstructuredGrid*& grid)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if (Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();
    vtk_sol_array->SetName(sol_name.c_str());
    vtk_sol_array->SetNumberOfValues(sol.size());

    // save data into vtk data array
    std::map<unsigned int, float>::iterator it=sol.begin();
    for( ; it!=sol.end(); ++it)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array->InsertValue(it->first, it->second);
    }

    grid->GetPointData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }
}


void VTKIO::write_node_complex_solution(std::map<unsigned int, std::complex<float> > & sol, const std::string & sol_name, vtkUnstructuredGrid*& grid)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if ( Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array_magnitude = vtkFloatArray::New();
    vtkFloatArray *vtk_sol_array_angle     = vtkFloatArray::New();

    vtk_sol_array_magnitude->SetName((sol_name+" abs").c_str());
    vtk_sol_array_magnitude->SetNumberOfValues(sol.size());

    vtk_sol_array_angle->SetName((sol_name+" angle").c_str());
    vtk_sol_array_angle->SetNumberOfValues(sol.size());

    // save data into vtk data array
    std::map<unsigned int, std::complex<float> >::iterator it=sol.begin();
    for( ; it!=sol.end(); ++it)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array_magnitude->InsertValue(it->first, std::abs(it->second));
      vtk_sol_array_angle->InsertValue(it->first, std::arg(it->second));
    }

    grid->GetPointData()->AddArray(vtk_sol_array_magnitude);
    grid->GetPointData()->AddArray(vtk_sol_array_angle);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array_magnitude->Delete();
    vtk_sol_array_angle->Delete();
  }

}


void VTKIO::write_node_vector_solution(std::map<unsigned int, float > & sol_x,
                                       std::map<unsigned int, float > & sol_y,
                                       std::map<unsigned int, float > & sol_z,
                                       const std::string & sol_name, vtkUnstructuredGrid*& grid)
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
    vtk_sol_array->SetNumberOfTuples(sol_x.size());
    vtk_sol_array->SetName(sol_name.c_str());

    // save data into vtk data array
    for( unsigned int i=0; i<sol_x.size(); i++)
    {
      float sol[3]={sol_x[i], sol_y[i], sol_z[i]};
      vtk_sol_array->InsertTuple(i, &sol[0] );
    }

    grid->GetPointData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }
}



void VTKIO::write_cell_scaler_solution(std::map<unsigned int, float> & sol, const std::string & sol_name, vtkUnstructuredGrid*& grid)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if (Genius::processor_id() == 0)
  {
    //create vtk data array
    vtkFloatArray *vtk_sol_array = vtkFloatArray::New();
    vtk_sol_array->SetName(sol_name.c_str());
    vtk_sol_array->SetNumberOfValues(sol.size());

    // save data into vtk data array
    std::map<unsigned int, float>::iterator it=sol.begin();
    for( ; it!=sol.end(); ++it)
    {
      // set scalar value to vtkFloatArray
      vtk_sol_array->InsertValue(it->first, it->second);
    }

    grid->GetCellData()->AddArray(vtk_sol_array);

    //---------------------------------------------------------

    //free the vtk array
    vtk_sol_array->Delete();
  }
}


void VTKIO::write_cell_vector_solution(std::map<unsigned int, float > & sol_x,
                                       std::map<unsigned int, float > & sol_y,
                                       std::map<unsigned int, float > & sol_z,
                                       const std::string & sol_name, vtkUnstructuredGrid*& grid)
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
    vtk_sol_array->SetNumberOfTuples(sol_x.size());
    vtk_sol_array->SetName(sol_name.c_str());

    // save data into vtk data array
    for( unsigned int i=0; i<sol_x.size(); i++)
    {
      float sol[3]={sol_x[i], sol_y[i], sol_z[i]};
      vtk_sol_array->InsertTuple(i, &sol[0] );
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
  // only do it at root node
  if(Genius::processor_id() != 0) return;

  out << "# vtk DataFile Version 3.0"      <<'\n';
  out << "Date Generated by Genius TCAD"   <<'\n';
  out << "ASCII"                           <<'\n';
  out << "DATASET UNSTRUCTURED_GRID"       <<'\n';
  out << "POINTS " << mesh.n_nodes()       << " float" << '\n';

  // write out nodal data
  MeshBase::const_node_iterator nd = mesh.active_nodes_begin();
  MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();
  for (;nd!=nd_end;++nd)
  {
    const Node * node = (*nd);
    out << (*node)(0)/um << " " << (*node)(1)/um << " " << (*node)(2)/um << '\n';
  }

  out << std::endl;
}


void VTKIO::cells_to_vtk(const MeshBase& mesh,  std::ofstream & out)
{

  // only do it at root node
  if(Genius::processor_id() != 0) return;

  int cell_size = mesh.n_active_elem();

  MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();
  for ( ; it != end; ++it)
  {
    const Elem *elem  = (*it);
    cell_size   += elem->n_nodes();
  }

  out<<"CELLS "<<mesh.n_active_elem()<<" "<<cell_size<<'\n';


  for (it = mesh.active_elements_begin(); it != end; ++it)
  {
    const Elem *elem  = (*it);
    out << elem->n_nodes();
    for(unsigned int i=0;i<elem->n_nodes();++i)
    {
      out << " " << elem->node(i);
    }
    out << std::endl;
  }
  out << std::endl;

  out << "CELL_TYPES " << mesh.n_active_elem() << '\n';


  for (it = mesh.active_elements_begin(); it != end; ++it)
  {
    unsigned int celltype;
    Elem *elem  = (*it);
    switch(elem->type())
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
      celltype = 5;  //VTK_TRIANGLE;
      break;// 3
    case TRI6:
      celltype = 22; //VTK_QUADRATIC_TRIANGLE;
      break;// 4
    case QUAD4:
    case QUAD4_FVM:
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
        std::cerr<<"element type "<<elem->type()<<" not implemented"<<std::endl;
        genius_error();
      }
    }

    out << celltype << std::endl;
  }

  out << std::endl;
}


void VTKIO::meshinfo_cell_to_vtk(const MeshBase& mesh, std::ofstream & out)
{

  // only do it at root node
  if(Genius::processor_id() != 0) return;

  //write cell based region and partition info to vtk

  // region information
  out << "SCALARS region float 1"  <<'\n';
  out << "LOOKUP_TABLE default"    <<'\n';
  MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end();
  for ( ; it != end; ++it)
  {
    out << static_cast<float>((*it)->subdomain_id())<<'\n';
  }
  out<<std::endl;

  // material information
  std::map<unsigned int, unsigned int> sub_id_to_material_id;
  for( unsigned int n=0; n<mesh.n_subdomains(); n++)
  {
    const std::string &material = mesh.subdomain_material(n);
    sub_id_to_material_id[n] = Material::get_id_by_material(material);
  }

  out << "SCALARS material float 1"  <<'\n';
  out << "LOOKUP_TABLE default"    <<'\n';
  for (it = mesh.active_elements_begin(); it != end; ++it)
  {
    out << static_cast<float>(sub_id_to_material_id[(*it)->subdomain_id()])<<'\n';
  }
  out<<std::endl;


  // partition information
  out<<"SCALARS partition float 1"  <<'\n';
  out<<"LOOKUP_TABLE default"       <<'\n';
  for (it = mesh.active_elements_begin(); it != end; ++it)
  {
    out<<static_cast<float>((*it)->processor_id())<<'\n';
  }
  out<<std::endl;

}




void VTKIO::meshinfo_node_to_vtk(const MeshBase& mesh, std::ofstream & out)
{
  // only do it at root node
  if(Genius::processor_id() != 0) return;

  out<<"SCALARS node_order float 1" <<'\n';
  out<<"LOOKUP_TABLE default"       <<'\n';
  MeshBase::const_node_iterator nd = mesh.active_nodes_begin();
  MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();
  for (; nd != nd_end; nd++)
  {
    out<<static_cast<float>((*nd)->id())<<'\n' ;
  }
  out<<std::endl;


  out<<"SCALARS node_partition float 1" <<'\n';
  out<<"LOOKUP_TABLE default"       <<'\n';
  nd = mesh.active_nodes_begin();
  for (; nd != nd_end; nd++)
  {
    out<<static_cast<float>((*nd)->processor_id())<<'\n' ;
  }
  out<<std::endl;


  out<<"SCALARS boundary float 1"  <<'\n';
  out<<"LOOKUP_TABLE default"      <<'\n';
  nd = mesh.active_nodes_begin();
  for(;nd!=nd_end;++nd)
  {
    if (mesh.boundary_info->boundary_id(*nd) != BoundaryInfo::invalid_id )
    {
      out<<static_cast<float>(mesh.boundary_info->boundary_id(*nd))<<'\n';
    }
    else
      out<<0.0<<'\n';
  }
  out<<std::endl;
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

  const std::vector<SolverSpecify::SolverType> & solve_histroy = system.solve_histroy();

  for(unsigned int n=0; n<solve_histroy.size(); ++n)
  {
    switch  (solve_histroy[n])
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
    const unsigned int n_elems = mesh.n_active_elem();
    out << "CELL_DATA "<<n_elems     <<'\n';

    meshinfo_cell_to_vtk(mesh, out);

    // solution data
    std::map<unsigned int, float> Ex,  Ey,  Ez;
    std::map<unsigned int, float> Jnx,  Jny,  Jnz;
    std::map<unsigned int, float> Jpx,  Jpy,  Jpz;

    for( unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      for(unsigned int n=0; n<region->n_cell(); ++n)
      {
        const Elem * elem = region->get_region_elem(n);
        const FVM_CellData * elem_data = region->get_region_elem_data(n);

        Ex[elem->id()]  = static_cast<float>(elem_data->E()(0)/(V/cm));
        Ey[elem->id()]  = static_cast<float>(elem_data->E()(1)/(V/cm));
        Ez[elem->id()]  = static_cast<float>(elem_data->E()(2)/(V/cm));

        Jnx[elem->id()]  = static_cast<float>(elem_data->Jn()(0)/(A/cm));
        Jny[elem->id()]  = static_cast<float>(elem_data->Jn()(1)/(A/cm));
        Jnz[elem->id()]  = static_cast<float>(elem_data->Jn()(2)/(A/cm));

        Jpx[elem->id()]  = static_cast<float>(elem_data->Jp()(0)/(A/cm));
        Jpy[elem->id()]  = static_cast<float>(elem_data->Jp()(1)/(A/cm));
        Jpz[elem->id()]  = static_cast<float>(elem_data->Jp()(2)/(A/cm));
      }
    }

    write_cell_vector_solution(Ex,  Ey,  Ez,  "electrical_field", out);
    write_cell_vector_solution(Jnx, Jny, Jnz, "electron_current", out);
    write_cell_vector_solution(Jpx, Jpy, Jpz, "hole_current", out);
  }


  // write node based data
  {
    //write node based boundary info to vtk
    const unsigned int n_nodes = mesh.n_nodes();
    out<<"POINT_DATA "<<n_nodes      <<'\n';

    meshinfo_node_to_vtk(mesh, out);

    // scale back to normal unit
    double concentration_scale = pow(cm, -3);

    // gather data from all the processor
    std::map<unsigned int, float> Na, Nd, net_doping, net_charge, psi, Ec, Ev, qFn, qFp, n, p, T, Tn, Tp; //scalar value
    std::map<unsigned int, float> mole_x, mole_y;
    std::map<unsigned int, float> OptG, PatG;

    std::map<unsigned int, std::complex<float> >psi_ac, n_ac, p_ac, T_ac, Tn_ac, Tp_ac; //complex value
    std::map<unsigned int, std::complex<float> > OptE_complex, OptH_complex;

    std::map<unsigned int, float> Jnx, Jny, Jnz; //vector value
    std::map<unsigned int, float> Jpx, Jpy, Jpz;

    // search all the node belongs to current processor
    MeshBase::const_node_iterator nd = mesh.pid_nodes_begin(Genius::processor_id());
    MeshBase::const_node_iterator nd_end = mesh.pid_nodes_end(Genius::processor_id());
    for (; nd != nd_end; ++nd)
    {
      if( (*nd)->processor_id() != Genius::processor_id() ) continue;

      for( unsigned int r=0; r<system.n_regions(); r++)
      {
        const SimulationRegion * region = system.region(r);
        const FVM_Node * fvm_node = region->region_fvm_node(*nd);
        if(fvm_node==NULL) continue;

        const FVM_NodeData * node_data = region->region_node_data(*nd);

        // if the fvm_node lies on the interface of two material regions,
        // we shall use the node data in the more important region.
        // it just for visualization reason
        if( fvm_node->boundary_id() != BoundaryInfo::invalid_id )
        {
          unsigned int bc_index = system.get_bcs()->get_bc_index_by_bd_id(fvm_node->boundary_id());
          const BoundaryCondition * bc = system.get_bcs()->get_bc(bc_index);
          const FVM_Node * primary_fvm_node = (*bc->region_node_begin(fvm_node->root_node())).second.second;
          genius_assert(primary_fvm_node->root_node()==(*nd));
          node_data = primary_fvm_node->node_data();
        }

        assert ( node_data != NULL );
        {
          if(semiconductor_material)
          {
            Na        [(*nd)->id()]  = static_cast<float>(node_data->Total_Na()/concentration_scale);
            Nd        [(*nd)->id()]  = static_cast<float>(node_data->Total_Nd()/concentration_scale);
            net_doping[(*nd)->id()]  = static_cast<float>(node_data->Net_doping()/concentration_scale);
            net_charge[(*nd)->id()]  = static_cast<float>(node_data->Net_charge()/concentration_scale);
            n         [(*nd)->id()]  = static_cast<float>(node_data->n()/concentration_scale);
            p         [(*nd)->id()]  = static_cast<float>(node_data->p()/concentration_scale);
            mole_x    [(*nd)->id()]  = static_cast<float>(node_data->mole_x());
            mole_y    [(*nd)->id()]  = static_cast<float>(node_data->mole_y());

            Jnx[(*nd)->id()]  = static_cast<float>(node_data->Jn()(0)/(A/cm));
            Jny[(*nd)->id()]  = static_cast<float>(node_data->Jn()(1)/(A/cm));
            Jnz[(*nd)->id()]  = static_cast<float>(node_data->Jn()(2)/(A/cm));

            Jpx[(*nd)->id()]  = static_cast<float>(node_data->Jp()(0)/(A/cm));
            Jpy[(*nd)->id()]  = static_cast<float>(node_data->Jp()(1)/(A/cm));
            Jpz[(*nd)->id()]  = static_cast<float>(node_data->Jp()(2)/(A/cm));
          }

          psi[(*nd)->id()]  = static_cast<float>(node_data->psi()/V);
          Ec[(*nd)->id()]   = static_cast<float>(node_data->Ec()/eV);
          Ev[(*nd)->id()]   = static_cast<float>(node_data->Ev()/eV);
          qFn[(*nd)->id()]  = static_cast<float>(node_data->qFn()/eV);
          qFp[(*nd)->id()]  = static_cast<float>(node_data->qFp()/eV);
          T  [(*nd)->id()]  = static_cast<float>(node_data->T()/K);

          if(ebm_solution_data)
          {
            Tn [(*nd)->id()]  = static_cast<float>(node_data->Tn()/K);
            Tp [(*nd)->id()]  = static_cast<float>(node_data->Tp()/K);
          }

          if(optical_generation)
            OptG[(*nd)->id()] = static_cast<float>(node_data->OptG()/(concentration_scale/s));

          if(particle_generation)
            PatG[(*nd)->id()] = static_cast<float>(node_data->PatG()/(concentration_scale/s));

          if(ddm_ac_data)
          {
            psi_ac [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->psi_ac()/V);
            n_ac   [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->n_ac()/concentration_scale);
            p_ac   [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->p_ac()/concentration_scale);
            T_ac   [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->T_ac()/K);
            Tn_ac  [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->Tn_ac()/K);
            Tp_ac  [(*nd)->id()]   = static_cast<std::complex<float> >(node_data->Tp_ac()/K);
          }

          if(optical_complex_field)
          {
            OptE_complex[(*nd)->id()] = static_cast<std::complex<float> >(node_data->OptE_complex()/(V/cm));
            OptH_complex[(*nd)->id()] = static_cast<std::complex<float> >(node_data->OptH_complex()/(A/cm));
          }
        }
      }
    }

    write_node_scaler_solution(psi, "psi", out);
    write_node_scaler_solution(Ec,  "Ec",  out);
    write_node_scaler_solution(Ev,  "Ev",  out);
    write_node_scaler_solution(qFn, "elec_quasi_Fermi_level", out);
    write_node_scaler_solution(qFp, "hole_quasi_Fermi_level", out);

    if(semiconductor_material)
    {
      write_node_scaler_solution(Na,  "Na", out);
      write_node_scaler_solution(Nd,  "Nd", out);
      write_node_scaler_solution(n,   "electron_density", out);
      write_node_scaler_solution(p,   "hole_density", out);
      write_node_scaler_solution(net_doping,  "net_doping", out);
      write_node_scaler_solution(net_charge,  "net_charge", out);
      write_node_scaler_solution(mole_x, "mole_x", out);
      write_node_scaler_solution(mole_y, "mole_y", out);
      //write_node_vector_solution(Jnx, Jny, Jnz, "Jn", out);
      //write_node_vector_solution(Jpx, Jpy, Jpz, "Jp", out);
    }

    write_node_scaler_solution(T,   "temperature", out);

    if(ebm_solution_data)
    {
      write_node_scaler_solution(Tn,  "elec_temperature", out);
      write_node_scaler_solution(Tp,  "hole_temperature", out);
    }

    if(ddm_ac_data)
    {
      write_node_complex_solution(psi_ac, "psi_AC", out);
      write_node_complex_solution(n_ac,   "electron_AC", out);
      write_node_complex_solution(p_ac,   "hole_AC", out);
      write_node_complex_solution(T_ac,   "temperature_AC", out);
      write_node_complex_solution(Tn_ac,  "elec_temperature_AC", out);
      write_node_complex_solution(Tp_ac,  "hole_temperature_AC", out);
    }

    if(optical_complex_field)
    {
      write_node_complex_solution(OptE_complex, "Optical_E", out);
      write_node_complex_solution(OptH_complex, "Optical_H", out);
    }

    if(optical_generation)
      write_node_scaler_solution(OptG,   "Optical_Generation", out);

    if(particle_generation)
      write_node_scaler_solution(PatG,   "Radiation_Generation", out);
  }

}

void VTKIO::write_node_scaler_solution(std::map<unsigned int, float> & sol, const std::string & sol_name, std::ofstream & out)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if ( Genius::processor_id() == 0)
  {
    out << "SCALARS "<<sol_name<<" float 1" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    std::map<unsigned int, float>::iterator it=sol.begin();
    for( ; it!=sol.end(); ++it)
    {
      out << it->second << std::endl;
    }
    out<<std::endl;
  }
}


void VTKIO::write_node_complex_solution(std::map<unsigned int, std::complex<float> > & sol, const std::string & sol_name, std::ofstream & out)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if ( Genius::processor_id() == 0)
  {
    // out put magnitude
    out << "SCALARS "<<sol_name<<"_abs float 1" << std::endl;
    out << "LOOKUP_TABLE default"   << std::endl;
    std::map<unsigned int, std::complex<float> >::iterator it=sol.begin();
    for( ; it!=sol.end(); ++it)
    {
      out << std::abs(it->second) << std::endl;
    }
    out<<std::endl;


    // out put phase angle
    out << "SCALARS "<<sol_name<<"_angle float 1" << std::endl;
    out << "LOOKUP_TABLE default"   << std::endl;
    for(it=sol.begin(); it!=sol.end(); ++it)
    {
      out << std::arg(it->second) << std::endl;
    }
    out<<std::endl;
  }

}


void VTKIO::write_node_vector_solution(std::map<unsigned int, float > & sol_x,
                                       std::map<unsigned int, float > & sol_y,
                                       std::map<unsigned int, float > & sol_z,
                                       const std::string & sol_name, std::ofstream & out)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol_x);
  Parallel::gather(0, sol_y);
  Parallel::gather(0, sol_z);

  if ( Genius::processor_id() == 0)
  {
    out << "VECTORS "<<sol_name<<" float" << std::endl;
    for( unsigned int i=0; i<sol_x.size(); i++)
    {
      out << sol_x[i] <<" "<< sol_y[i] <<" "<< sol_z[i] <<std::endl;
    }
    out<<std::endl;
  }
}


void VTKIO::write_cell_scaler_solution(std::map<unsigned int, float> & sol, const std::string & sol_name, std::ofstream & out)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol);

  if ( Genius::processor_id() == 0)
  {
    out << "SCALARS "<<sol_name<<" float 1" << std::endl;
    out << "LOOKUP_TABLE default" << std::endl;
    std::map<unsigned int, float>::iterator it=sol.begin();
    for( ; it!=sol.end(); ++it)
    {
      out << it->second << std::endl;
    }
    out<<std::endl;
  }
}

void VTKIO::write_cell_vector_solution(std::map<unsigned int, float > & sol_x,
                                       std::map<unsigned int, float > & sol_y,
                                       std::map<unsigned int, float > & sol_z,
                                       const std::string & sol_name, std::ofstream & out)
{
  // this should run on parallel for all the processor
  Parallel::gather(0, sol_x);
  Parallel::gather(0, sol_y);
  Parallel::gather(0, sol_z);

  if ( Genius::processor_id() == 0)
  {
    out << "VECTORS "<<sol_name<<" float" << std::endl;
    for( unsigned int i=0; i<sol_x.size(); i++)
    {
      out << sol_x[i] <<" "<< sol_y[i] <<" "<< sol_z[i] <<std::endl;
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
  Mesh & mesh = system.mesh();

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

      std::map<const std::string, int> bd_map;
      typedef std::map<const std::string, int>::iterator Bd_It;

      const Mesh::element_iterator el_end = mesh.elements_end();
      for (Mesh::element_iterator el = mesh.elements_begin(); el != el_end; ++el)
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



