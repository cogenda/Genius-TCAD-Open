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

//  $Id: simulation_system.cc,v 1.53 2008/07/09 09:10:08 gdiso Exp $

#include <sstream>
#include <numeric>
#include <queue>
#include <algorithm>

#include "parser.h"
#include "unstructured_mesh.h"
#include "simulation_system.h"
#include "simulation_region.h"
#include "semiconductor_region.h"
#include "insulator_region.h"
#include "conductor_region.h"
#include "resistance_region.h"
#include "vacuum_region.h"
#include "pml_region.h"
#include "parallel.h"
#include "boundary_info.h"
#include "boundary_condition_collector.h"
#include "surface_locator_hub.h"
#include "electrical_source.h"
#include "field_source.h"

#include "vtk_io.h"
#include "vtk2_io.h"
#include "cgns_io.h"
#include "stanford_io.h"
#include "tif3d_io.h"
#include "gdml_io.h"
#include "dfise_io.h"
#include "gmsh_io.h"
#include "unv_io.h"

#include "spice_ckt.h"
#include "location_io.h"

#include "interpolation_2d_csa.h"

#include "perf_log.h"
#include "sync_file.h"


#if defined(HAVE_TR1_UNORDERED_MAP)
  #include <tr1/unordered_map>
#elif defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER) || defined(HAVE_UNORDERED_MAP)
  #include <unordered_map>
#endif


#define __SELF_CHECK__


SimulationSystem::SimulationSystem(MeshBase & mesh)
  : _mesh(mesh), _cylindrical_mesh(false), _distributed_mesh(true), _resistive_metal_mode(false), _block_partition(true),
    _bcs(0), _electrical_source(0),
    _field_source(0), _spice_ckt(0), _global_z_width(false)
{
  // set PhysicalUnit
  PhysicalUnit::set_unit( std::pow(1e18,1.0/3.0) );

  _T_external = 300.0*PhysicalUnit::K;
  _z_width = 1.0*PhysicalUnit::um;
}



SimulationSystem::SimulationSystem(MeshBase & mesh, Parser::InputParser & _decks)
  :  _T_external(300.0), _mesh(mesh), _cylindrical_mesh(false), _distributed_mesh(true), _resistive_metal_mode(false), _block_partition(true),
    _bcs(0), _electrical_source(0),
    _field_source(0), _spice_ckt(0), _global_z_width(false), _z_width(1.0)
{

  MESSAGE<<"Constructing Simulation System...\n"<<std::endl;  RECORD();

  // set physical unit
  double cm = 1e6;
  double s  = 1e12;
  double V  = 1.0;
  double C  = 1.0/1.602176462e-19;
  double K  = 1.0/300;
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();
    if( c.key() == "GLOBAL" )   // find the global card
    {
      if( c.is_parameter_exist("dopingscale") )
        cm = std::pow(c.get_real("dopingscale", 1e18), 1.0/3.0);

      if( c.is_parameter_exist("cm") )
        cm = c.get_real("cm", 1e6);

      if( c.is_parameter_exist("second") )
        s = c.get_real("second", 1e12);

      if( c.is_parameter_exist("volt") )
        V = c.get_real("volt", 1.0);

      if( c.is_parameter_exist("coulomb") )
        C = c.get_real("coulomb", 1.0/1.602176462e-19);

      if( c.is_parameter_exist("kelvin") )
        K = c.get_real("kelvin", 1.0/300);
    }
  }

  // set PhysicalUnit
  PhysicalUnit::set_unit(cm, s, V, C, K);

  _T_external *= PhysicalUnit::K;
  _z_width *= PhysicalUnit::um;

  _coord_to_crystal_coord = TensorValue<double>(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0);

  // parse GLOBAL card again
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();
    if( c.key() == "GLOBAL" )   // find the global card
    {
      if( c.is_parameter_exist("texternal") )
        _T_external   = c.get_real("texternal", 300.0, "latticetemp")*PhysicalUnit::K;

      if( c.is_parameter_exist("z.width") )
      {
        _global_z_width = true;
        _z_width = c.get_real("z.width", 1.0)*PhysicalUnit::um;
      }

      if( c.is_parameter_exist("axis.x") && c.is_parameter_exist("axis.y") && c.is_parameter_exist("axis.z"))
      {
        _axis_x = c.get_array<int>("axis.x");
        _axis_y = c.get_array<int>("axis.y");
        _axis_z = c.get_array<int>("axis.z");
        _build_coord_to_crystal_coord();
      }


      _cylindrical_mesh = c.get_bool("cylindricalmesh", false);
      _distributed_mesh = c.get_bool("distributedmesh", true);
      _resistive_metal_mode = c.get_bool("resistivemetal", false);
      _block_partition = c.get_bool("blockpartition", true);

      double res = c.get_real("leakage.res", 1e100)*PhysicalUnit::V/PhysicalUnit::A;
      double cap = c.get_real("leakage.cap", 0.0)*PhysicalUnit::C/PhysicalUnit::V;
      MetalSimulationRegion::set_aux_parasitic_parameter(std::max(res, 1e-3*PhysicalUnit::V/PhysicalUnit::A), cap);
    }
  }

  if(_cylindrical_mesh) _z_width = 1.0;


  MESSAGE<< "External Temperature = " << _T_external/PhysicalUnit::K << 'K' <<std::endl; RECORD();

  // electrical source
  _electrical_source = new ElectricalSource( _decks );

  // field sources
  _field_source = new FieldSource(*this, _decks);

  // set magnetic field
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();
    if( c.key() == "MAGNETICFIELD" )
    {
      _magnetic_field(0) = c.get_real("bx", 0.0)*PhysicalUnit::V/(PhysicalUnit::m*PhysicalUnit::m)*PhysicalUnit::s;
      _magnetic_field(1) = c.get_real("by", 0.0)*PhysicalUnit::V/(PhysicalUnit::m*PhysicalUnit::m)*PhysicalUnit::s;
      _magnetic_field(2) = c.get_real("bz", 0.0)*PhysicalUnit::V/(PhysicalUnit::m*PhysicalUnit::m)*PhysicalUnit::s;
    }
  }

  // set spice netlist
  for( _decks.begin(); !_decks.end(); _decks.next() )
  {
    Parser::Card c = _decks.get_current_card();
    if( c.key() == "CIRCUIT" && (c.is_parameter_exist("netlist") || c.is_parameter_exist("spice.file")))
    {
      std::string ckt_file = c.get_string("netlist", "", "spice.file");
      std::string local_ckt_file = sync_file(ckt_file.c_str());
      // only the last processor link to ngspice
      if(Genius::processor_id()==Genius::n_processors()-1)
      {
        _spice_ckt = new SPICE_CKT(local_ckt_file);
        _spice_ckt->set_temperature(_T_external/PhysicalUnit::K);
      }
      else // create empty SPICE_CKT
        _spice_ckt = new SPICE_CKT;

      remove(local_ckt_file.c_str());
      // sync information between processors
      _spice_ckt->sync();

      if(c.is_parameter_exist("spice.nodeset"))
      {
        std::string spice_nodeset_file = c.get_string("spice.nodeset", "");
        std::string local_spice_nodeset_file = sync_file(spice_nodeset_file.c_str());

        if(Genius::processor_id()==Genius::n_processors()-1)
          _spice_ckt->import_solution(local_spice_nodeset_file);
        remove(local_spice_nodeset_file.c_str());
      }

      break;
    }
  }

  _bcs = new BoundaryConditionCollector(*this, _decks);
}



SimulationSystem::~SimulationSystem()
{
  for (unsigned int r=0; r<n_regions(); r++)
    delete _simulation_regions[r];
  _simulation_regions.clear();

  delete _bcs;
  delete _electrical_source;
  delete _field_source;
  delete _spice_ckt;
}



void SimulationSystem::_build_coord_to_crystal_coord()
{
  Point crystal_x(1.0, 0.0, 0.0);
  Point crystal_y(0.0, 1.0, 0.0);
  Point crystal_z(0.0, 0.0, 1.0);


  Point x = Point(_axis_x[0], _axis_x[1], _axis_x[2]).unit();
  Point y = Point(_axis_y[0], _axis_y[1], _axis_y[2]).unit();
  Point z = Point(_axis_z[0], _axis_z[1], _axis_z[2]).unit();

  double alpha, beta, gramma;

  if (z==crystal_z)
  {
    alpha=acos(x[0]);
    beta=0;
    gramma=0;
  }
  else if (z == -crystal_z)
  {
    alpha=acos(x[0]);
    beta=3.14159265358979323846;
    gramma=0;
  }
  else
  {
    alpha=acos(-z[1]/sqrt(1-z[2]*z[2]));
    beta=acos(z[2]);
    gramma=acos(Point(cos(alpha), sin(alpha), 0)*x);
  }

 TensorValue<double> D(cos(alpha), sin(alpha), 0,-sin(alpha), cos(alpha), 0,0, 0, 1);
 TensorValue<double> C(1, 0, 0, 0, cos(beta), sin(beta), 0, -sin(beta), cos(beta));
 TensorValue<double> B(cos(gramma), sin(gramma), 0, -sin(gramma), cos(gramma), 0, 0, 0, 1);

 _coord_to_crystal_coord = B*C*D;

 //std::cout<<_coord_to_crystal_coord;
}





double SimulationSystem::z_width() const
{ return (_mesh.mesh_dimension()==2) ? _z_width : 1.0 ;}



unsigned int SimulationSystem::dim() const
{ return _mesh.mesh_dimension(); }


void SimulationSystem::clear(bool clear_mesh)
{
  if(clear_mesh)
    _mesh.clear();

  for (unsigned int r=0; r<n_regions(); r++)
    delete _simulation_regions[r];
  _simulation_regions.clear();

  _bcs->clear();

  _electrical_source->clear_bc_source_map();

  //since we cleared all the solution data, previous solve histroy is meaningless
  _solver_active_history.clear();
}



void SimulationSystem::init_region()
{
  for(unsigned int r=0; r<_simulation_regions.size(); r++)
    _simulation_regions[r]->init(_T_external);
}


void SimulationSystem::reinit_region_after_import()
{
  for(unsigned int r=0; r<_simulation_regions.size(); r++)
    _simulation_regions[r]->reinit_after_import();
}


void SimulationSystem::init_region_post_process()
{
  if( this->empty() ) return;
#if 0
  // check fvm cell quality
  std::vector<Real> region_cell_quality;
  for(unsigned int r=0; r<_simulation_regions.size(); r++)
    region_cell_quality.push_back(_simulation_regions[r]->fvm_cell_quality());
  Parallel::min(region_cell_quality);

  for(unsigned int r=0; r<_simulation_regions.size(); r++)
  {
    Real fvm_cell_quality = region_cell_quality[r];
    std::string region = _simulation_regions[r]->name();
    if(fvm_cell_quality < 1e-4)
    {
      MESSAGE << "Warning: Voronoi cell quality " << fvm_cell_quality << " (the ratio of min/max surface area) of region " << region << " is poor. "
          << "It may cause inaccurate result or even convergence difficults. Try to improve mesh quality if possible!" << std::endl<< std::endl;
      RECORD();
    }
  }
#endif
}


unsigned int SimulationSystem::n_regions() const
{ return  _simulation_regions.size(); }


SimulationRegion * SimulationSystem::region(unsigned int n)
{
  assert(n < _simulation_regions.size());
  return _simulation_regions[n];
}


const SimulationRegion * SimulationSystem::region(unsigned int n) const
{
  assert(n < _simulation_regions.size());
  return _simulation_regions[n];
}


bool SimulationSystem::has_region(const std::string & region_label) const
{
  for(unsigned int r=0; r<n_regions(); r++)
    if( _simulation_regions[r]->name() == region_label )
      return true;
  return false;
}



bool SimulationSystem::has_region_with_material(const std::string & material) const
{
  for(unsigned int r=0; r<n_regions(); r++)
    if( _simulation_regions[r]->material() == material )
      return true;
  return false;
}


SimulationRegion * SimulationSystem::region(const std::string & region_label)
{
  for(unsigned int r=0; r<n_regions(); r++)
    if( _simulation_regions[r]->name() == region_label )
      return _simulation_regions[r];
  return NULL;
}


const SimulationRegion * SimulationSystem::region(const std::string & region_label) const
{
  for(unsigned int r=0; r<n_regions(); r++)
    if( _simulation_regions[r]->name() == region_label )
      return _simulation_regions[r];
  return NULL;
}


bool SimulationSystem::has_single_compound_semiconductor_region() const
{
  for(unsigned int n=0; n<n_regions(); ++n)
  {
    const SimulationRegion * region = this->region(n);
    if(Material::IsSingleCompSemiconductor(region->material()))
      return true;
  }
  return false;
}


bool SimulationSystem::has_complex_compound_semiconductor_region() const
{
  for(unsigned int n=0; n<n_regions(); ++n)
  {
    const SimulationRegion * region = this->region(n);
    if(Material::IsComplexCompSemiconductor(region->material()))
      return true;
  }
  return false;
}


void SimulationSystem::build_simulation_system()
{
  // sync resistive_metal_mode
  Parallel::broadcast(_resistive_metal_mode);

  // each region has its own FVM mesh
  build_region_fvm_mesh();

  // boundary condition
  _bcs->bc_setup();

  // electrical source should konw where is the (electrode) bc
  _electrical_source->link_to_bcs( _bcs );

  // clear point/surface locator to save memory
  _mesh.clear_point_locator();
  _mesh.clear_surface_locator();

  // delete remote elems
  if(_distributed_mesh && !_field_source->request_serial_mesh() && Genius::processor_id() !=0 )
    _mesh.delete_remote_elements(true, true);

  // remove remote object in each region
  for(unsigned int n = 0; n < this->n_regions(); n++)
  {
    SimulationRegion * region = _simulation_regions[n];
    region->remove_remote_object();
  }
}




void SimulationSystem::build_region_fvm_mesh()
{
  START_LOG("build_region_fvm_mesh()", "SimulationSystem");

  MESSAGE<<"Building simulation data structure on all processors..."<<std::endl;  RECORD();

  // we should convert initial mesh elements to FVM element
  // NOTE: for parallel situation, only local elements are converted for saving memory
  // the mesh is prepared after the function call all_fvm_elem ()

  {
    START_LOG("build_region_fvm_mesh(1)", "SimulationSystem");

    UnstructuredMesh & mesh = dynamic_cast<UnstructuredMesh &>(_mesh);

    MESSAGE<<"  Create mesh topological information...";  RECORD();
    // 2d or 3d mesh?
    mesh.count_mesh_dimension();

    mesh.build_mesh_bounding_box();

    // this function will renumber the the node/elem
    mesh.all_first_order();

    // let all the elements find their neighbors
    mesh.find_neighbors();

#if 0
    // reorder the elem/node index by Reverse Cuthill-McKee Algorithm
    std::string err;
    if(!mesh.reorder_elems(err))
    {
      MESSAGE<<err;RECORD();
      genius_error();
    }
#endif
    MESSAGE<<std::endl;  RECORD();


    MESSAGE<<"  Partition mesh...";  RECORD();

    // prepare for partition
    if(_block_partition)
      mesh.subdomain_cluster(this->build_subdomain_cluster());

    // partition the mesh.
    mesh.partition(Genius::n_processors());

    // ok, mesh is prepared
    mesh.set_prepared();

    // remove remote mesh elements when processor_id > 1
    // however, keep boundary elems for later bc setup
    if(_distributed_mesh && !_field_source->request_serial_mesh() && Genius::processor_id() !=0 )
      _mesh.delete_remote_elements(true, false);

    MESSAGE<<std::endl;  RECORD();

    STOP_LOG("build_region_fvm_mesh(1)", "SimulationSystem");



    START_LOG("build_region_fvm_mesh(2)", "SimulationSystem");
    MESSAGE<<"  Create mesh element for finite volume method...";  RECORD();

    if(_cylindrical_mesh)
    {
      std::string error;
      if( mesh.convert_to_cylindrical_fvm_mesh (error) == false )
      {
        MESSAGE<<"  bad mesh."<<std::endl;  RECORD();
        MESSAGE<<"  ERROR:" << error <<std::endl;  RECORD();
        genius_error();
      }
    }
    else
    {
      std::string error;
      if( mesh.convert_to_fvm_mesh (error) == false )
      {
        MESSAGE<<"  bad mesh."<<std::endl;  RECORD();
        MESSAGE<<"  ERROR:" << error <<std::endl;  RECORD();
        genius_error();
      }
    }

    //std::cout<<"FVM MESH " << _mesh.memory_usage()/(1024*1024)<<std::endl;


    MESSAGE<<std::endl;  RECORD();
    STOP_LOG("build_region_fvm_mesh(2)", "SimulationSystem");
  }


  _simulation_regions.resize( _mesh.n_subdomains() );
  std::map<unsigned int,  SimulationRegion *> subdomain_id_to_region_map;

  unsigned int dim = _mesh.mesh_dimension();
  // create simulation region from subdomain information
  for(unsigned int r=0; r<_mesh.n_subdomains(); r++)
  {
    std::string label = _mesh.subdomain_label_by_id(r);
    std::string material = _mesh.subdomain_material(r);
    switch( Material::material_type( material ) )
    {
      case Material::Semiconductor                :
      case Material::SingleCompoundSemiconductor  :
      case Material::ComplexCompoundSemiconductor :
      {
        _simulation_regions[r] = new SemiconductorSimulationRegion(label, material, _T_external, dim, _z_width, _coord_to_crystal_coord);
        break;
      }
      case Material::Insulator     :
      {
        _simulation_regions[r] = new InsulatorSimulationRegion(label, material, _T_external, dim, _z_width);
        break;
      }
      case Material::Conductor     :
      {
        _simulation_regions[r] = new ElectrodeSimulationRegion(label, material, _T_external, dim, _z_width);
        break;
      }
      case Material::Resistance     :
      {
        if(_resistive_metal_mode)
          _simulation_regions[r] = new MetalSimulationRegion(label, material, _T_external, dim, _z_width);
        else
          _simulation_regions[r] = new ElectrodeSimulationRegion(label, material, _T_external, dim, _z_width);
        break;
      }
      case Material::Vacuum       :
      {
        _simulation_regions[r] = new VacuumSimulationRegion(label, material, _T_external, dim, _z_width);
        break;
      }
      case Material::PML       :
      {
        _simulation_regions[r] = new PMLSimulationRegion(label, material, _T_external, dim, _z_width);
        break;
      }
      default:
      {
        MESSAGE<< "ERROR: Region "<< label << " with unsupported material type " << material << std::endl; RECORD();
        genius_error();
      }
    }

    //
    _simulation_regions[r]->set_subdomain_id(r);
    subdomain_id_to_region_map[r] = _simulation_regions[r];

    //set partition weight for each region
    _mesh.set_subdomain_weight(r, Material::material_weight(_mesh.subdomain_material(r)));
  }

  // each region should hold subdomain_id_to_region_map
  SimulationRegion::set_subdomain_id_to_region_map(subdomain_id_to_region_map);


  START_LOG("build_region_fvm_mesh(3)", "SimulationSystem");
  MESSAGE<<"  Building finite volume cells...";  RECORD();

  typedef const Node *                    key_type;
  typedef FVM_Node *                      val_type;
  typedef std::pair<key_type, val_type>   key_val_pair;

#if defined(HAVE_UNORDERED_MAP)
  typedef std::unordered_multimap<key_type, val_type> map_type;
#elif defined(HAVE_TR1_UNORDERED_MAP) || defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER)
  typedef std::tr1::unordered_multimap<key_type, val_type> map_type;
#else
  typedef std::multimap<key_type, val_type>  map_type;
#endif

  map_type _node_to_fvm_node_map;

  // A convenient typedef
  typedef map_type::iterator Iter;

  // search in all the LOCAL element
  MeshBase::element_iterator       el  = _mesh.local_elements_begin();
  const MeshBase::element_iterator end = _mesh.local_elements_end();
  for (; el != end; ++el)
  {
    Elem * elem = *el;
    genius_assert(elem->on_local());

    std::vector<FVM_Node *> elem_fvm_nodes(elem->n_nodes());
    std::vector<FVM_Node *> global_elem_fvm_nodes(elem->n_nodes(), (FVM_Node *)0);

    for (unsigned int n=0; n<elem->n_nodes(); n++)
    {
      Node * node = elem->get_node(n);

      // only set informations for local node
      genius_assert( node->on_local() );

      // create a fvm node from cell's node
      FVM_Node *fvm_node = new FVM_Node( node );
      genius_assert(fvm_node);

      // set the subdomain_id of fvm node the same as the element
      fvm_node->set_subdomain_id( elem->subdomain_id () );

      // the fvm node belongs to which cell, also record the local index of the root node
      fvm_node->add_elem_it_belongs(elem, n) ;

      // set the partial volume associated with this node
      // only on_local FVM elements own non-zero value!
      fvm_node->set_control_volume( elem->partial_volume_truncated(n) );

      // buffer for later use
      elem_fvm_nodes[n] = fvm_node;
      // if this node already has a FVM_Node
      std::pair<Iter, Iter> pos = _node_to_fvm_node_map.equal_range(node);
      while (pos.first != pos.second)
      {
        FVM_Node *  fvm_current = (*pos.first).second;
        if ( fvm_current->subdomain_id () == fvm_node->subdomain_id () )
        {
          global_elem_fvm_nodes[n] = fvm_current;
          break;
        }
        ++pos.first;
      }
    }

    // process FVM_Node neighbor information
    for (unsigned int e=0; e<elem->n_edges(); e++)
    {
      //elem edges
      std::pair<unsigned int, unsigned int> edge_nodes;
      elem->nodes_on_edge(e, edge_nodes);

      // surface area associated with this edge
      Real surface_area = elem->partial_area_with_edge(e);

      FVM_Node * fvm_node1  = elem_fvm_nodes[edge_nodes.first];
      FVM_Node * fvm_node1g = global_elem_fvm_nodes[edge_nodes.first];
      FVM_Node * fvm_node2  = elem_fvm_nodes[edge_nodes.second];
      FVM_Node * fvm_node2g = global_elem_fvm_nodes[edge_nodes.second];

      fvm_node1->add_fvm_node_neighbor(fvm_node2g ? fvm_node2g : fvm_node2, surface_area);
      fvm_node2->add_fvm_node_neighbor(fvm_node1g ? fvm_node1g : fvm_node1, surface_area);
    }

    // at last, we insert this FVM Node into
    // std::multimap<const Node * , FVM_Node *> _node_to_fvm_node_map;
    for (unsigned int n=0; n<elem->n_nodes(); n++)
    {
      FVM_Node * fvm_node   = elem_fvm_nodes[n];
      FVM_Node * fvm_node_g = global_elem_fvm_nodes[n];

      if(fvm_node_g)
      {
        *fvm_node_g +=  *fvm_node;
        delete fvm_node;
        fvm_node = fvm_node_g;
      }
      else
      {
        _node_to_fvm_node_map.insert( std::make_pair(elem->get_node(n), fvm_node) );
      }

      elem->hold_fvm_node( n, fvm_node );
    }

  }






  // prepare ghost node information
  {
    Iter it_fvm_end = _node_to_fvm_node_map.end();
    for(Iter  it_fvm = _node_to_fvm_node_map.begin(); it_fvm != it_fvm_end; ++it_fvm )
    {
      // skip nonlocal fvm_node
      if( !it_fvm->second->on_local() ) continue;

      // the FVM Nodes with same root node. they will be ghost nodes in different region
      std::pair<Iter, Iter> pos = _node_to_fvm_node_map.equal_range( it_fvm->first );

      // insert them into ghost node map. The interface area is set to 0.0 here, will be changed later.
      while (pos.first != pos.second)
      {
        if ( pos.first != it_fvm )
         it_fvm->second->set_ghost_node( (*pos.first).second, (*pos.first).second->subdomain_id (), 0.0 );

        ++pos.first;
      }
    }
  }


  MESSAGE<<std::endl;  RECORD();
  STOP_LOG("build_region_fvm_mesh(3)", "SimulationSystem");


  START_LOG("build_region_fvm_mesh(4)", "SimulationSystem");
  MESSAGE<<"  Building boundary cells...";  RECORD();
  // we scan boundary face to find the area of interface side
  // NOTE here we should use side list of active elements!
  std::vector<unsigned int>       elems;
  std::vector<unsigned short int> sides;
  std::vector<short int>          bds;
  _mesh.boundary_info->build_active_side_list (elems, sides, bds);

  {
    typedef const Node *                    key_type;
    typedef std::pair<unsigned int, Real>   val_type;
    typedef std::pair<key_type, val_type>   key_val_pair;

#if defined(HAVE_UNORDERED_MAP)
    typedef std::unordered_multimap<key_type, val_type> map_type;
#elif defined(HAVE_TR1_UNORDERED_MAP) || defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER)
    typedef std::tr1::unordered_multimap<key_type, val_type> map_type;
#else
    typedef std::multimap<key_type, val_type>  map_type;
#endif

    // map <node , pair<subdomain_id, boundary_area> >
    map_type bd_area_map;
    typedef map_type::iterator Bda_It;

    for (size_t nbd=0; nbd<elems.size(); nbd++ )
    {
      // get the element which has boundary/interface side
      const Elem* elem = _mesh.elem(elems[nbd]);
      if( !elem || !elem->on_local() ) continue;

      // get the side
      AutoPtr<Elem> side (elem->build_side(sides[nbd]));

      // build corresponding FVM elem of the side
      AutoPtr<Elem> fvm_side = Elem::build (Elem::fvm_compatible_type(side->type()), side->parent());

      for (unsigned int v=0; v < side->n_vertices(); v++)
        fvm_side->set_node(v) = side->get_node(v);

      fvm_side->prepare_for_fvm();

      for (unsigned int v=0; v < fvm_side->n_vertices(); v++)
      {
        const Node * node = fvm_side->get_node(v);
        if( !node->on_local() ) continue;

        // if we can find this node exists in bd_area_map
        if ( bd_area_map.find(node) != bd_area_map.end() )
        {

          std::pair<Bda_It, Bda_It> pos = bd_area_map.equal_range(node);

          //if Node has same subdomain_id as the element. assign area to this node
          while (pos.first != pos.second)
          {
            if ( (*pos.first).second.first == elem->subdomain_id() )
            {
              (*pos.first).second.second +=  fvm_side->partial_volume_truncated(v);
              break;
            }
            ++pos.first;
          }
          // not find? insert a new Node
          if (pos.first == pos.second)
            bd_area_map.insert(pos.first, std::make_pair(node, std::make_pair(elem->subdomain_id(), fvm_side->partial_volume_truncated(v))));

        }
        else // not find? insert a new Node
        {
          bd_area_map.insert(std::make_pair(node, std::make_pair(elem->subdomain_id(), fvm_side->partial_volume_truncated(v))));
        }
      }

    }


    // set FVM interface area
    Bda_It it_bd = bd_area_map.begin();
    Bda_It it_bd_end = bd_area_map.end();
    for( ; it_bd != it_bd_end; ++it_bd )
    {
      // the FVM Nodes with same root node. they may be ghost nodes in different region
      std::pair<Iter, Iter> pos = _node_to_fvm_node_map.equal_range( (*it_bd).first );

      // skip nonlocal fvm_node
      if(! (*pos.first).second->on_local() ) continue;

      // insert them into ghost node map.
      while (pos.first != pos.second)
      {
        // set ghost node with the same subdomain id as boundary node
        (*pos.first).second->set_ghost_node_area( (*it_bd).second.first, (*it_bd).second.second );

        ++pos.first;
      }
    }
  }
  MESSAGE<<std::endl;  RECORD();
  STOP_LOG("build_region_fvm_mesh(4)", "SimulationSystem");


  START_LOG("build_region_fvm_mesh(5)", "SimulationSystem");
  MESSAGE<<"  Building norm vector for each interface...";  RECORD();
  // build norm vector of interface
  {
    std::multimap<const Node * , std::pair<unsigned int, VectorValue<Real> > > bd_norm_map;
    typedef std::multimap<const Node * , std::pair<unsigned int, VectorValue<Real> > >::iterator Bdn_It;
    for (size_t nbd=0; nbd<elems.size(); nbd++ )
    {
      // get the element which has boundary/interface side
      const Elem* elem = _mesh.elem(elems[nbd]);

      if(!elem || !elem->on_local()) continue;

      // the vector norm to boundary/interface
      VectorValue<Real> norm = elem->outside_unit_normal(sides[nbd]);

      // get the side
      AutoPtr<Elem> side (elem->build_side(sides[nbd]));
      for (unsigned int v=0; v < side->n_vertices(); v++)
      {
        const Node * node = side->get_node(v);
        bd_norm_map.insert(std::make_pair(node, std::make_pair(elem->subdomain_id(), norm)));
      }
    }

    // classify norm vector to each FVM_Node
    std::map<FVM_Node *, std::vector<VectorValue<Real> > > fvm_norm_vectors;
    Bdn_It it_bd = bd_norm_map.begin();
    for( ; it_bd != bd_norm_map.end(); ++it_bd )
    {
      const Node * node = (*it_bd).first;
      unsigned int sub_id = (*it_bd).second.first;

      // the FVM Nodes with same root node. they may be ghost nodes in different region
      std::pair<Iter, Iter> pos = _node_to_fvm_node_map.equal_range( node );

      // skip fvm_node not on processor
      if(! (*pos.first).second->on_processor() ) continue;

      while (pos.first != pos.second)
      {
        FVM_Node * fvm_node = (*pos.first).second;
        if( sub_id == fvm_node->subdomain_id() )
          fvm_norm_vectors[fvm_node].push_back((*it_bd).second.second);
        ++pos.first;
      }
    }

    // ok, save average norm vector to FVM_Node
    std::map<FVM_Node *, std::vector<VectorValue<Real> > >::iterator it = fvm_norm_vectors.begin();
    for(; it!=fvm_norm_vectors.end(); ++it)
    {
      const std::vector<VectorValue<Real> > & norms = it->second;
      VectorValue<Real> norm;
      for(unsigned int n=0; n<norms.size(); ++n)
      {
        norm += norms[n];
      }
      norm /= norms.size();
      it->first->set_norm(norm.unit(true));
    }
  }

  MESSAGE<<std::endl;  RECORD();
  STOP_LOG("build_region_fvm_mesh(5)", "SimulationSystem");


  START_LOG("build_region_fvm_mesh(6)", "SimulationSystem");
  MESSAGE<<"  Setup simulation regions...";  RECORD();

  // reserve memory for data block
  {
    std::vector<unsigned int> region_cells(this->n_regions(), 0);
    for (el=_mesh.local_elements_begin(); el != end; ++el)
    {
      unsigned int region_index = (*el)-> subdomain_id ();
      region_cells[region_index]++;
    }

    std::vector<unsigned int> region_nodes(this->n_regions(), 0);
    Iter it_fvm_end = _node_to_fvm_node_map.end();
    for(Iter it_fvm = _node_to_fvm_node_map.begin(); it_fvm != it_fvm_end; ++it_fvm )
    {
      //insert the FVM Node into correspoding region
      FVM_Node *fvm_node = (*it_fvm).second;
      unsigned int region_index = fvm_node->subdomain_id();
      if(fvm_node->on_local()) region_nodes[region_index]++;
    }

    for(unsigned int n = 0; n < this->n_regions(); n++)
    {
      _simulation_regions[n]->reserve_data_block(region_cells[n], region_nodes[n]);
    }
  }

  // insert element / FVM_Node pointer into each region
  {
    for (el=_mesh.local_elements_begin(); el != end; ++el)
    {
      unsigned int region_index = (*el)-> subdomain_id ();
      _simulation_regions[region_index]->insert_cell(*el);
    }

    // this is overkill for parallel simulation since
    // we needn't all the nodal information
    Iter it_fvm_end = _node_to_fvm_node_map.end();
    for(Iter it_fvm = _node_to_fvm_node_map.begin(); it_fvm != it_fvm_end; ++it_fvm )
    {
      //insert the FVM Node into correspoding region
      unsigned int region_index = (*it_fvm).second->subdomain_id();
      _simulation_regions[region_index]->insert_fvm_node((*it_fvm).second);
    }
  }


  // ask all the regions to do some pre process
  for(unsigned int n = 0; n < this->n_regions(); n++)
  {
    _simulation_regions[n]->prepare_for_use();
  }
  // pre process should be executed in parallel
  for(unsigned int n = 0; n < this->n_regions(); n++)
  {
    _simulation_regions[n]->prepare_for_use_parallel();
  }

  MESSAGE<<std::endl;  RECORD();
  STOP_LOG("build_region_fvm_mesh(6)", "SimulationSystem");


  START_LOG("build_region_fvm_mesh(7)", "SimulationSystem");
  MESSAGE<<"  Setup hanging node...";  RECORD();
  // set region hanging node. we had make sure that no hanging node on the region interface.
  for(unsigned int n = 0; n < this->n_regions(); n++)
  {
    SimulationRegion * region = _simulation_regions[n];

    const Node * node = NULL;

    // only search local cells to see if they contain hanging nodes.
    SimulationRegion::const_element_iterator it = region->elements_begin();
    SimulationRegion::const_element_iterator it_end = region->elements_end();

    for(; it!=it_end; ++it)
    {
      // set hanging node on element side
      for(unsigned int n=0; n<(*it)->n_sides(); ++n)
        if ( (node = (*it)->is_hanging_node_on_side(n))!=NULL )
          region->add_hanging_node_on_side(node, *it, n);

      // set hanging node on element edge
      for(unsigned int e=0; e<(*it)->n_edges(); ++e)
        if ( (node = (*it)->is_hanging_node_on_edge(e))!=NULL )
          region->add_hanging_node_on_edge(node, *it, e);
    }
  }
  MESSAGE<<std::endl;  RECORD();
  STOP_LOG("build_region_fvm_mesh(7)", "SimulationSystem");


  START_LOG("build_region_fvm_mesh(8)", "SimulationSystem");
  MESSAGE<<"  Setup node distance to nearest surface...";  RECORD();
  SurfaceLocatorHub & surface_locator = _mesh.surface_locator();
  for(unsigned int n = 0; n < this->n_regions(); n++)
  {
    SimulationRegion * region = _simulation_regions[n];
    if(region->type() != SemiconductorRegion) continue;

    SimulationRegion::processor_node_iterator it = region->on_local_nodes_begin();
    SimulationRegion::processor_node_iterator it_end = region->on_local_nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = *it;
      FVM_NodeData * node_data = fvm_node->node_data();

      const Point p = *(fvm_node->root_node());

      Point project_point;
      std::pair<const Elem*, unsigned int> surface_elem_pair = surface_locator(p, n, project_point);
      if( surface_elem_pair.first )
        node_data->dmin() = (p-project_point).size();
    }
  }
  MESSAGE<<std::endl;  RECORD();
  STOP_LOG("build_region_fvm_mesh(8)", "SimulationSystem");


#ifdef __SELF_CHECK__

  double fvm_volume = 0.0;
  for(unsigned int n = 0; n < this->n_regions(); n++)
  {
    const SimulationRegion * region = _simulation_regions[n];
    std::vector<double> vols;
    region->region_node_vol(vols);

    fvm_volume += std::accumulate(vols.begin(), vols.end(), 0.0);
  }


  double mesh_volume = 0.0;
  {
    MeshBase::const_element_iterator       el  = _mesh.this_pid_elements_begin();
    const MeshBase::const_element_iterator end = _mesh.this_pid_elements_end();
    for (; el != end; ++el)
    {
      const Elem * elem = *el;
      mesh_volume += elem->volume();
    }
    Parallel::sum(mesh_volume);
  }

  //std::cout<< fvm_volume << " " << mesh_volume << std::endl;
  //genius_assert( std::abs(fvm_volume - mesh_volume)/mesh_volume < 1e-3 );

#endif


  MESSAGE<<"Simulation data structure build ok.\n"<<std::endl;  RECORD();


  STOP_LOG("build_region_fvm_mesh()", "SimulationSystem");
}



bool SimulationSystem::self_consistent() const
{
  return false;
}



void SimulationSystem::estimate_error (const Parser::Card &c, ErrorVector & error_per_cell) const
{

  // get the refinement depedent variable from input card

  bool v_volume = false, signedlog = false, gradient=true;
  SolutionVariable variable;

  if( c.is_enum_value("variable","volume"))
    v_volume = true;
  else
  {
    variable  = solution_string_to_enum(FormatVariableString(c.get_string("variable", "")));
    if(variable_data_type(variable)!=SCALAR)
    {
      MESSAGE<<"ERROR at " << c.get_fileline() <<" REFINE: Refine variable should be scalar value."<<std::endl; RECORD();
      genius_error();
    }
    signedlog = c.is_enum_value("measure", "signedlog");
  }

  if( c.is_enum_value("evaluation","quantity"))
    gradient = false;

  std::set<std::string> regions;
  for(unsigned int idx=0; idx<c.parameter_size(); idx++)
  {
    Parser::Parameter p = c.get_parameter(idx);
    if ( p.name() == "region" )
      regions.insert( p.get_string() );
  }

  bool check_coordinates=false;
  if( c.is_parameter_exist("x.min") || c.is_parameter_exist("x.max") ||
      c.is_parameter_exist("y.min") || c.is_parameter_exist("y.max") ||
      c.is_parameter_exist("z.min") || c.is_parameter_exist("z.max") )
    check_coordinates=true;

  double xmin = c.get_real("x.min", 0.0, "x.left")*PhysicalUnit::um    - 1e-6;
  double xmax = c.get_real("x.max", 0.0, "x.right")*PhysicalUnit::um   + 1e-6;
  double ymin = c.get_real("y.min", 0.0, "y.top")*PhysicalUnit::um     - 1e-6;
  double ymax = c.get_real("y.max", 0.0, "y.bottom")*PhysicalUnit::um  + 1e-6;
  double zmin = c.get_real("z.min", 0.0, "z.front")*PhysicalUnit::um   - 1e-6;
  double zmax = c.get_real("z.max", 0.0, "z.back")*PhysicalUnit::um    + 1e-6;

  if(xmin>xmax || ymin>ymax || zmin>zmax)
  {
    MESSAGE<<"ERROR at " << c.get_fileline() <<" REFINE: Refine region has incorrect XYZ bound."<<std::endl; RECORD();
    genius_error();
  }


  std::map<unsigned int, ErrorVectorReal> cell_error_map;

  bool refine_flag = true;
  // if regions array is not empty, we only refine region in the regions array!
  if( !regions.empty() ) refine_flag = false;

  // fill gradient of var for all the cells in each region
  for(unsigned int n=0; n<n_regions(); n++)
  {
    const SimulationRegion * region = this->region(n);

    bool region_refine_flag = false;
    // check if this region should be refined.
    if( regions.find(region->name()) != regions.end())
      region_refine_flag = true;

    SimulationRegion::const_element_iterator it = region->elements_begin();
    SimulationRegion::const_element_iterator it_end = region->elements_end();
    for(; it!=it_end; ++it)
    {
      // only process cell belongs to this processor
      if( (*it)->processor_id() != Genius::processor_id() ) continue;

      bool element_refine_flag = true;
      std::vector<PetscScalar> var_vertex;

      double exmin=1e100, exmax=-1e100;
      double eymin=1e100, eymax=-1e100;
      double ezmin=1e100, ezmax=-1e100;
      for(unsigned int nd=0; nd<(*it)->n_nodes(); ++nd)
      {
        const FVM_Node * fvm_node = (*it)->get_fvm_node(nd);
        const FVM_NodeData * fvm_node_data = fvm_node->node_data();
        const Node * node = fvm_node->root_node();

        if ( check_coordinates )
        {
          exmin = (*node)(0) < exmin ? (*node)(0) : exmin;
          exmax = (*node)(0) > exmax ? (*node)(0) : exmax;
          eymin = (*node)(1) < eymin ? (*node)(1) : eymin;
          eymax = (*node)(1) > eymax ? (*node)(1) : eymax;
          ezmin = (*node)(2) < ezmin ? (*node)(2) : ezmin;
          ezmax = (*node)(2) > ezmax ? (*node)(2) : ezmax;
        }

        if (!v_volume)
        {
          PetscScalar var =  fvm_node_data->get_variable_real(variable);
          if(signedlog) var = std::sign(var)*std::log10(1+std::abs(var));
          var_vertex.push_back  ( var );
        }
      }

      if ( check_coordinates )
        if (exmin>xmax || exmax<xmin || eymin>ymax || eymax<ymin || ezmin>zmax || ezmax<zmin)
          element_refine_flag = false;

      if((refine_flag || region_refine_flag) && element_refine_flag)
      {
        if (v_volume)
        {
          cell_error_map[(*it)->id()] = static_cast<ErrorVectorReal>((*it)->volume()/std::pow(PhysicalUnit::um,3.0)); // use cell volume as error
        }
        else
        {
          // compute the gradient
          if(gradient)
          {
            VectorValue<PetscScalar> grad_var   = (*it)->gradient(var_vertex);
            cell_error_map[(*it)->id()] = static_cast<ErrorVectorReal>(grad_var.size()*(*it)->hmax()); // use gradient*hmax as error
          }
          // use cell average value as error
          else
          {
            cell_error_map[(*it)->id()] = std::accumulate(var_vertex.begin(), var_vertex.end(), 0.0)/var_vertex.size();
          }
        }
      }
      else
        cell_error_map[(*it)->id()] = 0.0;
    }
  }

  // gather from all the processors
  Parallel::gather(0, cell_error_map);

  if( Genius::processor_id() == 0)
  {
    // reserve memory for error_per_cell vector
    error_per_cell.resize (_mesh.max_elem_id());
    // fill error_per_cell with 0 as init value
    std::fill(error_per_cell.begin(), error_per_cell.end(), 0.0);

    // fill into error_per_cell
    std::map<unsigned int, ErrorVectorReal>::iterator it = cell_error_map.begin();
    for(; it!=cell_error_map.end(); ++it)
      error_per_cell[(*it).first] = (*it).second;
  }

}



/**
 * fill system data (mesh and node data) into interpolator for later usage
 * type 0 -- linear interpolation
 * type 1 -- signed log interpolation
 * type 2 -- asinh interpolation
 */
void SimulationSystem::fill_interpolator(InterpolationBase *interpolator,
    const std::string & variable_string,
    InterpolationBase::InterpolationType type) const
{
  SolutionVariable variable = solution_string_to_enum(FormatVariableString(variable_string));
  genius_assert(variable!=INVALID_Variable);
  genius_assert(variable_data_type(variable)==SCALAR);

  int group_code = interpolator->set_group_code(variable_string);

  std::map<unsigned int, double> value_map;
  for( unsigned int r=0; r<this->n_regions(); r++)
  {
    const SimulationRegion * region = this->region(r);

    SimulationRegion::const_processor_node_iterator on_processor_nodes_it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator on_processor_nodes_it_end = region->on_processor_nodes_end();
    for(; on_processor_nodes_it!=on_processor_nodes_it_end; ++on_processor_nodes_it)
    {
      const FVM_Node * fvm_node = *on_processor_nodes_it;
      const FVM_NodeData * node_data = fvm_node->node_data();

      // if the fvm_node lies on the interface of two material regions,
      // we shall use the node data in the more important region.
      if( fvm_node->boundary_id() != BoundaryInfo::invalid_id )
      {
        unsigned int bc_index = this->get_bcs()->get_bc_index_by_bd_id(fvm_node->boundary_id());
        const BoundaryCondition * bc = this->get_bcs()->get_bc(bc_index);
        const FVM_Node * primary_fvm_node = (*bc->region_node_begin(fvm_node->root_node())).second.second;
        node_data = primary_fvm_node->node_data();
      }
      if(node_data->is_variable_valid(variable))
        value_map [fvm_node->root_node()->id()] = node_data->get_variable_real(variable);
    }
  }
  Parallel::allgather(value_map);

  interpolator->set_interpolation_type(group_code, type);

  // fill the interpolator
  std::map<unsigned int, double>::const_iterator it = value_map.begin();
  for(; it != value_map.end(); ++it)
    interpolator->add_scatter_data(_mesh.point(it->first), group_code, it->second);

  interpolator->setup(group_code);
}


/**
 * get data from interpolator after mesh refinement
 */
void SimulationSystem::do_interpolation(const InterpolationBase * interpolator , const std::string & variable_string)
{
  SolutionVariable variable = solution_string_to_enum(FormatVariableString(variable_string));
  genius_assert(variable!=INVALID_Variable);
  genius_assert(variable_data_type(variable)==SCALAR);

  int group_code = interpolator->group_code(variable_string);

  // fill gradient of var for all the cells in each region
  for(unsigned int n=0; n<n_regions(); n++)
  {
    SimulationRegion * region = this->region(n);

    SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
    SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      FVM_Node * fvm_node = (*node_it);

      FVM_NodeData * node_data = fvm_node->node_data();
      if(node_data->is_variable_valid(variable))
      {
        double value = interpolator->get_interpolated_value(*(fvm_node->root_node()), group_code);
        node_data->set_variable_real(variable, value);
      }
    }
  }
}



std::vector< std::vector<unsigned int > > SimulationSystem::build_subdomain_cluster()
{
  std::vector<std::vector<unsigned int> > subdomain_adjncy;
  _mesh.subdomain_graph(subdomain_adjncy);

  std::map<Material::MaterialType, std::set<unsigned int> > material_subdomains;
  for(unsigned int s=0; s<_mesh.n_subdomains(); s++)
  {
    Material::MaterialType type = Material::material_type(_mesh.subdomain_material(s));
    material_subdomains[type].insert(s);
  }

  std::map<unsigned int, std::vector<unsigned int> > material_subdomain_cluster;

  std::map<Material::MaterialType, std::set<unsigned int> >::const_iterator material_subdomains_it = material_subdomains.begin();
  for(; material_subdomains_it != material_subdomains.end(); ++material_subdomains_it)
  {
    const std::set<unsigned int> & subdomains = material_subdomains_it->second;
    for(std::set<unsigned int>::const_iterator it=subdomains.begin(); it!=subdomains.end(); ++it)
    {
      unsigned int region_me = *it;
      // find regions with the same material type and connected to me
      std::queue<unsigned int> Q;
      std::set<unsigned int> visit_flag;

      Q.push( region_me );
      while(!Q.empty())
      {
        unsigned int region = Q.front();
        Q.pop();

        visit_flag.insert(region);
        material_subdomain_cluster[region_me].push_back(region);

        const std::vector<unsigned int> & region_neighbors = subdomain_adjncy[region];
        for(unsigned int n=0; n<region_neighbors.size(); ++n)
        {
          unsigned int region_neighbor = region_neighbors[n];
          if( visit_flag.find(region_neighbor) == visit_flag.end() && subdomains.find(region_neighbor) != subdomains.end())
            Q.push(region_neighbor);
        }
      }
    }
  }

  std::set<unsigned int> subdomain_flag;
  std::vector< std::vector<unsigned int > > subdomain_cluster;
  std::map<unsigned int, std::vector<unsigned int> >::const_iterator it = material_subdomain_cluster.begin();
  for( ; it != material_subdomain_cluster.end(); ++it)
  {
    unsigned int region = it->first;
    if(subdomain_flag.find(region) != subdomain_flag.end()) continue;

    const std::vector<unsigned int> & cluster = it->second;
    subdomain_cluster.push_back(cluster);
    subdomain_flag.insert(cluster.begin(), cluster.end());
  }

  /*
  for(unsigned int n=0; n < subdomain_cluster.size(); ++n)
   {
    for(unsigned int m=0; m < subdomain_cluster[n].size(); ++m)
      std::cout<< _mesh.subdomain_label_by_id(subdomain_cluster[n][m] )<<std::endl;

    std::cout<<std::endl;
  }
  */

  return subdomain_cluster;
}



void SimulationSystem::print_info (std::ostream& os) const
{
  os << "Simulation System Information on processor " << Genius::processor_id() << " :" << '\n'
  << " total regions = "         << this->n_regions()           << '\n';

  for(unsigned int n = 0; n < this->n_regions(); n++)
  {
    os << "   region " << _simulation_regions[n]->name()
    << " with material " << _simulation_regions[n]->material()
    << " has "      << _simulation_regions[n]->n_on_processor_cell() << " cells, "
    << _simulation_regions[n]->n_on_processor_node() << " of "
    << _simulation_regions[n]->n_node()<<" total nodes."
    << '\n';
  }

  os << '\n';
}



void SimulationSystem::sync_print_info () const
{
  std::stringstream   ss;
  ss  << "Simulation System Information on processor " << Genius::processor_id() << " :" << '\n';


  std::vector<unsigned int> region_nodes;
  std::vector<unsigned int> region_elems;
  for(unsigned int n = 0; n < this->n_regions(); n++)
  {
    region_elems.push_back( _simulation_regions[n]->n_on_processor_cell() );
    region_nodes.push_back( _simulation_regions[n]->n_on_processor_node() );
  }

  Parallel::sum(region_nodes);
  Parallel::sum(region_elems);

  for(unsigned int n = 0; n < this->n_regions(); n++)
  {
    unsigned int n_on_processor_cell = _simulation_regions[n]->n_on_processor_cell();
    unsigned int n_on_processor_node = _simulation_regions[n]->n_on_processor_node();

    if(n_on_processor_cell==0 && n_on_processor_node==0) continue;

    ss << "   region " << _simulation_regions[n]->name()
    << " with material " << _simulation_regions[n]->material()
    << " has "
    << n_on_processor_cell << " of " << region_elems[n] << " total cells, "
    << n_on_processor_node << " of " << region_nodes[n] << " total nodes."
    << '\n';
  }

  ss << '\n';

  std::string out_info = ss.str();
  std::vector<char> out_info_char(out_info.size());
  std::copy(out_info.begin(), out_info.end(), out_info_char.begin());
  Parallel::allgather(out_info_char);
  out_info_char.push_back(0);
  out_info = std::string(&out_info_char[0]);

  MESSAGE<<out_info<<std::endl;  RECORD();

  //PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", out_info.c_str());
  //PetscSynchronizedFlush(PETSC_COMM_WORLD);
}



size_t SimulationSystem::memory_size() const
{
  size_t counter = sizeof(*this);

  //SimulationRegion
  for(unsigned int n=0; n<_simulation_regions.size(); ++n)
    counter += _simulation_regions[n]->memory_size();


  return counter;

}



void SimulationSystem::export_vtk(const std::string& filename, bool ascii) const
{
  if(!ascii)
  {
#ifdef HAVE_VTK
    std::string file_name = filename;
    // preprocess vtk file extension to make sure it has a ".vtu" format
    if (file_name.rfind(".vtu") > file_name.size())
    {
      // file name has a vtk extension, change it to vtu
      if (file_name.rfind(".vtk") < file_name.size())
      {
        file_name.replace(file_name.rfind("k"), 1, "u");
        MESSAGE<<"Change VTK file extension from .vtk to .vtu to meet XML VTK format." << std::endl; RECORD();
      }

      // if no file extension, add .vtu
      else
      {
        file_name += ".vtu";
        MESSAGE<<"Add VTK file extension .vtu to meet XML VTK format." << std::endl; RECORD();
      }
    }

    MESSAGE<<"Write System to XML VTK file "<< file_name << "...\n" << std::endl; RECORD();
    VTKIO(*this).write (file_name);
#else
    MESSAGE<<"Genius is not compiled with XML VTK support, skip VTK export... "<< std::endl; RECORD();
#endif

  }
  else
  {
    std::string file_name = filename;
    // preprocess vtk file extension to make sure it has a ".vtk" format
    if (file_name.rfind(".vtk") > file_name.size())
    {
      // file name has a vtu extension, change it to vtk
      if (file_name.rfind(".vtu") < file_name.size())
      {
        file_name.replace(file_name.rfind("u"), 1, "k");
        MESSAGE<<"Change VTK file extension .vtu to .vtk since Legacy VTK format is used here." << std::endl; RECORD();
      }
      // if no file extension, add .vtk
      else
      {
        file_name += ".vtk";
        MESSAGE<<"Add VTK file extension .vtk to meet Legacy VTK format." << std::endl; RECORD();
      }
    }

    MESSAGE<<"Write System to Legacy VTK file "<< file_name << "...\n" << std::endl; RECORD();
    VTKIO(*this).write (file_name);
  }


}



void SimulationSystem::export_vtk2(const std::string& filename, const std::vector<std::string> & variables) const
{
#ifdef HAVE_VTK
    std::string file_name = filename;
    // preprocess vtk file extension to make sure it has a ".vtu" format
    if (file_name.rfind(".vtu") > file_name.size())
    {
      // file name has a vtk extension, change it to vtu
      if (file_name.rfind(".vtk") < file_name.size())
      {
        file_name.replace(file_name.rfind("k"), 1, "u");
        MESSAGE<<"Change VTK file extension from .vtk to .vtu to meet XML VTK format." << std::endl; RECORD();
      }
      // if no file extension, add .vtu
      else
      {
        file_name += ".vtu";
        MESSAGE<<"Add VTK file extension .vtu to meet XML VTK format." << std::endl; RECORD();
      }
    }

    MESSAGE<<"Write System to XML VTK file "<< file_name << "...\n" << std::endl; RECORD();
    VTK2IO(*this, variables).write (file_name);
#endif
}


void SimulationSystem::export_node_location(const std::string& filename, const PetscScalar unit, const bool number) const
{
  MESSAGE<<"Write mesh node coordinates to file "<< filename << "...\n" << std::endl; RECORD();

  LocationIO io = LocationIO(*this);
  io.setNumbering (number);
  io.setLengthUnit (unit);
  io.write (filename);
}


void SimulationSystem::export_cgns(const std::string& filename) const
{
  MESSAGE<<"Write System to CGNS file "<< filename << "...\n" << std::endl; RECORD();

  CGNSIO(*this).write (filename);
}


void SimulationSystem::export_ise(const std::string& filename) const
{
  MESSAGE<<"Write System to DF-ISE file "<< filename << "...\n"; RECORD();

  DFISEIO(*this).write (filename);
}


void SimulationSystem::export_tif(const std::string& filename) const
{
  MESSAGE<<"Write System to TIF file "<< filename << "...\n"; RECORD();

  STIFIO(*this, "medici").write (filename);
}


void SimulationSystem::export_gdml_surface(const std::string& filename) const
{
  if(_mesh.mesh_dimension()==2)
  {
    MESSAGE<<"Write geometry information to GDML file requires 3D Mesh." << std::endl;
    return;
  }

  MESSAGE<<"Write geometry information (region) to GDML file "<< filename << "...\n" << std::endl; RECORD();

  GDMLIO(*this).write (filename);
}


void SimulationSystem::import_cgns(const std::string& filename)
{

  MESSAGE<<"Import System from CGNS file "<< filename << "...\n" << std::endl; RECORD();

  CGNSIO(*this).read (filename);
}

void SimulationSystem::import_vtk(const std::string& filename)
{
#ifdef HAVE_VTK
  MESSAGE<<"Import System from VTK file "<< filename << "...\n" << std::endl; RECORD();
  VTKIO(*this).read (filename);
#else
  MESSAGE<<"Genius is not compiled with XML VTK support, skip VTK import... "<< std::endl; RECORD();
#endif
}


void SimulationSystem::import_silvaco(const std::string& filename)
{
  MESSAGE<<"Import System from Silvaco file "<< filename << "...\n" << std::endl; RECORD();

  STIFIO(*this, "silvaco").read (filename);
}



void SimulationSystem::import_silvaco_list(const std::vector<std::string> & filenames)
{
  MESSAGE<<"Import System from Silvaco files...\n" << std::endl; RECORD();

  STIFIO(*this, "silvaco").read (filenames);
}


void SimulationSystem::import_tif(const std::string& filename)
{
  MESSAGE<<"Import System from TIF file "<< filename << "...\n" << std::endl; RECORD();

  STIFIO(*this, "medici").read (filename);
}



void SimulationSystem::import_tif_list(const std::vector<std::string> & filenames)
{
  MESSAGE<<"Import System from TIF files...\n" << std::endl; RECORD();

  STIFIO(*this, "medici").read (filenames);
}


void SimulationSystem::import_tif3d(const std::string& filename)
{
  MESSAGE<<"Import System from TIF3D file "<< filename << "...\n" << std::endl; RECORD();

  TIF3DIO(*this).read (filename);
}


void SimulationSystem::import_suprem(const std::string& filename)
{
  MESSAGE<<"Import System from suprem file "<< filename << "...\n" << std::endl; RECORD();

  STIFIO(*this, "suprem").read (filename);
}


void SimulationSystem::import_ise(const std::string& filename)
{

  MESSAGE<<"Import System from DF-ISE file "<< filename << "..." << std::endl; RECORD();

  DFISEIO(*this).read (filename);

  MESSAGE<< std::endl; RECORD();
}


void SimulationSystem::import_gmsh(const std::string& filename)
{
  MESSAGE<<"Import System from GMSH file "<< filename << "...\n" << std::endl; RECORD();

  GmshIO(*this).read (filename);
}


void SimulationSystem::import_unv(const std::string& filename)
{
  MESSAGE<<"Import System from UNV file "<< filename << "...\n" << std::endl; RECORD();

  UNVIO(*this).read (filename);
}