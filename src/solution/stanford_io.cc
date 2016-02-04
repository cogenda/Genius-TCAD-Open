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
#include<map>
#include <iterator>

// Local includes
#include "silvaco.h"
#include "medici.h"
#include "suprem.h"


#include "stanford_io.h"

#include "mesh_base.h"
#include "boundary_info.h"
#include "mesh_communication.h"
#include "simulation_region.h"
#include "log.h"
#include "parallel.h"
#include "material.h"

#define DEBUG

using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::V;
using PhysicalUnit::K;
using PhysicalUnit::eV;

/**
 * This method implements reading a mesh from 2D TIF format.
 */
void STIFIO::read (const std::string& filename)
{
  StanfordTIF * tif_reader=0;

  if(_format == "silvaco")
    tif_reader = new SilvacoTIF(filename);
  if(_format == "medici")
    tif_reader = new MediciTIF(filename);
  if(_format == "suprem")
    tif_reader = new SupremTIF(filename);

  // read tif file
  int ierr;
  std::string err;
  if( Genius::processor_id() == 0)
  {
    ierr = tif_reader->read(err);
  }
  Parallel::broadcast(ierr);
  if(!ierr)
  {
    delete tif_reader;
    MESSAGE<< err << std::endl; RECORD();
    genius_error();
  }

  //tif_reader->export_tif("debug.tif");
  //tif_reader->export_sup("debug.sup");


  // broadcast
  tif_reader->broadcast();

  if(tif_reader->dim() == 2)
    _import_2d(tif_reader);
  else
    _import_3d_silvaco(tif_reader);


  delete tif_reader;

}



void STIFIO::read (const std::vector<std::string> & files)
{
  std::vector<const StanfordTIF *> readers;

  if(_format == "silvaco")
  {
    for(unsigned int n=0; n<files.size(); ++n)
    {
      StanfordTIF * reader = new SilvacoTIF(files[n]);

      // read tif file
      int ierr;
      std::string err;
      if( Genius::processor_id() == 0)
      {
        ierr = reader->read(err);
      }
      Parallel::broadcast(ierr);
      if(!ierr)
      {
        delete reader;
        MESSAGE<< err << std::endl; RECORD();
        genius_error();
      }

      readers.push_back(reader);
    }
  }

  if(_format == "medici")
  {
    for(unsigned int n=0; n<files.size(); ++n)
    {
      StanfordTIF * reader = new MediciTIF(files[n]);

      // read tif file
      int ierr;
      std::string err;
      if( Genius::processor_id() == 0)
      {
        ierr = reader->read(err);
      }
      Parallel::broadcast(ierr);
      if(!ierr)
      {
        delete reader;
        MESSAGE<< err << std::endl; RECORD();
        genius_error();
      }

      readers.push_back(reader);
    }
  }

  StanfordTIF * tif_reader = StanfordTIF::merge(readers);

  for(unsigned int n=0; n<readers.size(); ++n)
    delete readers[n];


  //tif_reader->export_tif("debug.tif");


  // broadcast
  tif_reader->broadcast();

  if(tif_reader->dim() == 2)
    _import_2d(tif_reader);
  else
    _import_3d_silvaco(tif_reader);



  delete tif_reader;
}




void STIFIO::_import_2d (StanfordTIF * tif_reader)
{
  /*
   * after that, we fill mesh structure with mesh information read from TIF
   */
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  // clear the system
  system.clear();

  // map node * to tif index in TIF::node_array
  std::map<Node *, int>       node_to_tif_index_map;
  // map tif region to mesh region
  std::map<int, int>          tif_region_to_mesh_region;
  std::map<int, int>          mesh_region_to_tif_region;

  if( Genius::processor_id() == 0)
  {
    // fill node location
    std::vector<StanfordTIF::Node_t>::const_iterator tif_node_it = tif_reader->tif_nodes().begin();
    for(int i=0; tif_node_it!=tif_reader->tif_nodes().end(); ++i, ++tif_node_it)
    {
      Node * node = mesh.add_point( Point(tif_node_it->x*um, tif_node_it->y*um, 0) );
      node_to_tif_index_map[node] = i;
    }

    // map bc_index to bc label
    std::map<std::string, std::pair<short int, bool> > bd_map; // <label, <id, user_defind> >
    typedef std::map<std::string, std::pair<short int, bool> >::iterator Bd_It;

    // fill region label and region material
    // at the same time. set all the remaining boundary edges as "region_neumann"
    {
      std::vector<StanfordTIF::Region_t>::const_iterator tif_region_it  = tif_reader->region_array().begin();
      for(int r=0; tif_region_it!=tif_reader->region_array().end(); ++tif_region_it)
      {
        if(tif_region_it->segment) continue;

        std::string material = Material::FormatMaterialString(tif_region_it->material);
        mesh.set_subdomain_label(r, tif_region_it->name );
        mesh.set_subdomain_material(r, material);
        tif_region_to_mesh_region[tif_region_it->index] = r;
        mesh_region_to_tif_region[r] = tif_region_it->index;

        for(unsigned int i=0; i<tif_region_it->boundary.size(); ++i)
        {
          StanfordTIF::Edge_t & edge =  tif_reader->edge(tif_region_it->boundary[i]);
          edge.bc_index = r+1;
        }

        bd_map[tif_region_it->name + "_Neumann"] = std::make_pair(r+1,false);

        ++r;
      }
    }
    mesh.set_n_subdomains() = tif_region_to_mesh_region.size();

    // process explicit defined boundarys. assign bc_index to these boundary edges
    std::vector<StanfordTIF::Region_t>::const_iterator tif_boundary_it =  tif_reader->region_array().begin();
    for(short int bc_index=tif_reader->region_array().size()+1; tif_boundary_it !=  tif_reader->region_array().end(); ++bc_index, ++tif_boundary_it)
    {
      // it is a region, not boundary, skip it
      if(!tif_boundary_it->segment) continue;

      for(unsigned int i=0; i<tif_boundary_it->boundary.size(); ++i)
        tif_reader->edge(tif_boundary_it->boundary[i]).bc_index = bc_index;

      bd_map[tif_boundary_it->name] = std::make_pair(bc_index, true);
    }


    // build edge map
    std::map<StanfordTIF::Edge_t, int, StanfordTIF::lt_edge> edge_table;
    std::vector<StanfordTIF::Edge_t>::const_iterator  tif_edge_it = tif_reader->tif_edges().begin();
    for(; tif_edge_it!=tif_reader->tif_edges().end(); ++tif_edge_it)
    {
      genius_assert( tif_edge_it->bc_index != BoundaryInfo::invalid_id);
      edge_table[*tif_edge_it] = tif_edge_it->bc_index;
    }

    // fill triangles
    std::vector<StanfordTIF::Tri_t>::const_iterator tif_tri_it  = tif_reader->tif_tris().begin();
    for(; tif_tri_it!=tif_reader->tif_tris().end(); ++tif_tri_it)
    {
      Elem* elem = 0;
      bool quad = false;
      if( tif_tri_it->c4 < 0 )
        elem = mesh.add_elem(Elem::build(TRI3).release());
      else
      {
        quad = true;
        elem = mesh.add_elem(Elem::build(QUAD4).release());
      }

      // tri elem node
      elem->set_node(0) = mesh.node_ptr( tif_tri_it->c1 );
      elem->set_node(1) = mesh.node_ptr( tif_tri_it->c2 );
      elem->set_node(2) = mesh.node_ptr( tif_tri_it->c3 );
      if( quad )
        elem->set_node(3) = mesh.node_ptr( tif_tri_it->c4 );

      // which region this tri belongs to
      elem->subdomain_id() = tif_region_to_mesh_region[tif_tri_it->region];

      StanfordTIF::Edge_t edge;

      // process edge0
      edge.point1 = tif_tri_it->c1;
      edge.point2 = tif_tri_it->c2;
      if(edge.point1 > edge.point2) std::swap(edge.point1, edge.point2);
      std::map<StanfordTIF::Edge_t, int, StanfordTIF::lt_edge>::iterator edge_pointer = edge_table.find(edge);
      if( edge_pointer != edge_table.end() )
      {
        // this edge should on external boundary or region interface
        mesh.boundary_info->add_side(elem, 0, edge_pointer->second);
      }

      // process edge1
      edge.point1 = tif_tri_it->c2;
      edge.point2 = tif_tri_it->c3;
      if(edge.point1 > edge.point2) std::swap(edge.point1, edge.point2);
      edge_pointer = edge_table.find(edge);
      if( edge_pointer != edge_table.end() )
      {
        // this edge should on external boundary or region interface
        mesh.boundary_info->add_side(elem, 1, edge_pointer->second);
      }

      // process edge2
      if(!quad)
      {
        edge.point1 = tif_tri_it->c3;
        edge.point2 = tif_tri_it->c1;
        if(edge.point1 > edge.point2) std::swap(edge.point1, edge.point2);
        edge_pointer = edge_table.find(edge);
        if( edge_pointer != edge_table.end() )
        {
          // this edge should on external boundary or region interface
          mesh.boundary_info->add_side(elem, 2, edge_pointer->second);
        }
      }
      else
      {
        edge.point1 = tif_tri_it->c3;
        edge.point2 = tif_tri_it->c4;
        if(edge.point1 > edge.point2) std::swap(edge.point1, edge.point2);
        edge_pointer = edge_table.find(edge);
        if( edge_pointer != edge_table.end() )
        {
          mesh.boundary_info->add_side(elem, 2, edge_pointer->second);
        }

        edge.point1 = tif_tri_it->c4;
        edge.point2 = tif_tri_it->c1;
        if(edge.point1 > edge.point2) std::swap(edge.point1, edge.point2);
        edge_pointer = edge_table.find(edge);
        if( edge_pointer != edge_table.end() )
        {
          mesh.boundary_info->add_side(elem, 3, edge_pointer->second);
        }
      }
    }

    //however, the interface information should be set here

    std::vector<unsigned int>       elems;
    std::vector<unsigned short int> sides;
    std::vector<short int>          bds;

    // get all the boundary element
    mesh.boundary_info->build_side_list (elems, sides, bds);

    //build neighbor information for boundary element. then elem->neighbor() is functional
    mesh.boundary_info->find_neighbors();

    for (size_t nbd=0; nbd<elems.size(); nbd++ )
    {
      // get the element which has boundary/interface side
      const Elem* elem = mesh.elem(elems[nbd]);
      genius_assert(elem->on_boundary() || elem->on_interface());

      //is it an interface side && not a labeled face
      if( elem->neighbor(sides[nbd]) && bds[nbd]<=static_cast<short int>(tif_reader->region_array().size()) )
      {
        // the element and its neighbor should in diffetent subdomain
        unsigned int sbd_id1 = elem->subdomain_id();
        unsigned int sbd_id2 = elem->neighbor(sides[nbd])->subdomain_id();

        // delete the overkilled boundary side
        if (sbd_id1 == sbd_id2)
        {
          mesh.boundary_info->remove(elem, sides[nbd]);
          continue;
        }

        // the side should be an interface side
        genius_assert(elem->on_interface());
        genius_assert(elem->neighbor(sides[nbd])->on_interface());

        //remove the pair-element from boundary
        mesh.boundary_info->remove(elem, sides[nbd]);
        mesh.boundary_info->remove(elem->neighbor(sides[nbd]),
                                   elem->neighbor(sides[nbd])->which_neighbor_am_i(elem));

        // build the label for the interface, which has the form of RegionLabel1_to_RegionLabel2,
        // the two region is alpha ordered.
        std::string bd_label;
        if( tif_reader->region(sbd_id1).name < tif_reader->region(sbd_id2).name)
          bd_label = tif_reader->region(sbd_id1).name + "_to_" + tif_reader->region(sbd_id2).name;
        else
          bd_label = tif_reader->region(sbd_id2).name + "_to_" + tif_reader->region(sbd_id1).name;

        short int bd_index;

        // if the label already exist
        if( bd_map.find(bd_label) != bd_map.end() )
          bd_index = (*bd_map.find(bd_label)).second.first;
        else
        {
          //else, increase bd_index, insert it into bd_map
          bd_index = tif_reader->region_array().size() + bd_map.size() + 1;
          bd_map[bd_label] = std::make_pair(bd_index, false);
        }

        // add pair-element to boundary with new bd_index
        mesh.boundary_info->add_side(elem, sides[nbd], bd_index);
        mesh.boundary_info->add_side(elem->neighbor(sides[nbd]),
                                     elem->neighbor(sides[nbd])->which_neighbor_am_i(elem),
                                     bd_index);
      }
    }

    // after the previous remove and insert operation, the number of boundary id maybe changed. renumber it here.
    mesh.boundary_info->rebuild_ids();


    //write down bd labels
    for(Bd_It bd_it = bd_map.begin(); bd_it != bd_map.end(); ++bd_it)
    {
      mesh.boundary_info->set_label_to_id( bd_it->second.first, bd_it->first, bd_it->second.second);
    }

    // magic number, for 2D mesh, should < 2008
    mesh.magic_num() = 312;

  }


  /*
   * set mesh structure for all processors, and build simulation system
   */

  // broadcast mesh to all the processor
  MeshCommunication mesh_comm;
  mesh_comm.broadcast(mesh);


  // build simulation system
  system.build_simulation_system();
  system.sync_print_info();

  /*
   * after that, set doping infomation here. this should be done for all the processors
   */

  // node id to TIF index map, this should be same for all processors
  std::map<unsigned int, int> node_id_to_tif_index_map;
  // fill node_id_to_tif_index_map by processor 0
  if( Genius::processor_id() == 0)
  {
    std::map<Node *, int>::iterator it = node_to_tif_index_map.begin();
    for(; it != node_to_tif_index_map.end(); ++it)
      node_id_to_tif_index_map[(*it).first->id()] = (*it).second;
  }

  // broadcast node_id_to_tif_index_map and tif_region_to_mesh_region to all the processors
  Parallel::broadcast(node_id_to_tif_index_map , 0);
  Parallel::broadcast(tif_region_to_mesh_region, 0);
  Parallel::broadcast(mesh_region_to_tif_region, 0);


  // ok, we had got enough informations for set up each simulation region

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);

    switch ( region->type() )
    {
    case SemiconductorRegion :
      {
        bool sigle   = Material::IsSingleCompSemiconductor(region->material());
        bool complex = Material::IsComplexCompSemiconductor(region->material());

        const StanfordTIF::SolHead_t & sol_head = tif_reader->sol_head();
        for(int n=0; n<sol_head.sol_num; ++n)
        {
          const std::string  variable = sol_head.sol_name_array[n];
          const std::string  variable_unit = sol_head.solution_unit(variable);
          region->add_variable( SimulationVariable(variable, SCALAR, POINT_CENTER, variable_unit, invalid_uint, true, true)  );
        }

        SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
        SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
        for(; node_it!=node_it_end; ++node_it)
        {
          FVM_Node * fvm_node = (*node_it);
          FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);

          int tif_region_index = mesh_region_to_tif_region[r];
          // tif_node_index is the index in TIF file that this FVM node has
          int tif_node_index = node_id_to_tif_index_map[fvm_node->root_node()->id()];

          // doping
          {
            node_data->Na()   = tif_reader->acceptor(tif_region_index, tif_node_index) * pow(cm, -3);
            node_data->Nd()   = tif_reader->donor(tif_region_index, tif_node_index) * pow(cm, -3);
          }

          // mole fraction
          if(sigle)
          {
            node_data->mole_x()   = tif_reader->mole_x(tif_region_index, tif_node_index);
          }

          if(complex)
          {
            node_data->mole_x()   = tif_reader->mole_x(tif_region_index, tif_node_index);
            node_data->mole_y()   = tif_reader->mole_y(tif_region_index, tif_node_index);
          }

          for(int n=0; n<sol_head.sol_num; ++n)
          {
            const std::string  variable_name = sol_head.sol_name_array[n];
            const SimulationVariable & variable = region->get_variable(variable_name, POINT_CENTER);
            node_data->data<PetscScalar>(variable.variable_index) = tif_reader->solution(n, tif_region_index, tif_node_index)*variable.variable_unit;
          }
        }
        region->init(system.T_external());
        break;
      }
    case InsulatorRegion     :
      {
        region->init(system.T_external());
        break;
      }
    case ElectrodeRegion     :
      {
        region->init(system.T_external());
        break;
      }
    case MetalRegion    :
      {
        region->init(system.T_external());
        break;
      }
    case VacuumRegion        :
      {
        region->init(system.T_external());
        break;
      }
    case PMLRegion           :
      {
        region->init(system.T_external());
        break;
      }
    default: genius_error();
    }
  }
#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  system.init_region_post_process();

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}


void STIFIO::_import_3d_silvaco (StanfordTIF * tif_reader)
{
  /*
   * after that, we fill mesh structure with mesh information read from TIF
   */
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  // clear the system
  system.clear();

  // map node * to tif index in TIF::node_array
  std::map<Node *, std::pair<int, int> >       node_to_tif_index_map;
  std::map<std::pair<int, int>, Node * >       tif_index_to_node_map;

  // map tif region to mesh region
  std::map<int, int>          tif_region_to_mesh_region;


  if( Genius::processor_id() == 0)
  {
    for(int n=0; n<tif_reader->nz(); ++n)
    {
      double z = tif_reader->z(n);
      // fill node location
      std::vector<StanfordTIF::Node_t>::const_iterator tif_node_it = tif_reader->tif_nodes().begin();
      std::vector<StanfordTIF::Node_t>::const_iterator tif_node_it_end = tif_reader->tif_nodes().end();
      for(int i=0; tif_node_it!=tif_node_it_end; ++i, ++tif_node_it)
      {
        Node * node = mesh.add_point( Point(tif_node_it->x*um, tif_node_it->y*um, z*um) );
        node_to_tif_index_map[node] = std::make_pair(i, n);
        tif_index_to_node_map[std::make_pair(i, n)] = node;
      }
    }


    // map bc_index to bc label
    std::map<std::string, std::pair<short int, bool> > bd_map; // <label, <id, user_defind> >
    typedef std::map<std::string, std::pair<short int, bool> >::iterator Bd_It;

    // boundary face map < <edge1, z>, <edge2, z> >
    std::map<std::pair<std::pair<int, int>, std::pair<int, int> >, int> face_map;

    // fill region label and region material
    // at the same time. set all the remaining boundary edges as "region_neumann"
    {
      std::vector<StanfordTIF::Region_t>::const_iterator tif_region_it  = tif_reader->region_array().begin();
      for(int r=0; tif_region_it!=tif_reader->region_array().end(); ++tif_region_it)
      {
        if(tif_region_it->segment) continue;

        std::string material = Material::FormatMaterialString(tif_region_it->material);
        mesh.set_subdomain_label(r, tif_region_it->name );
        mesh.set_subdomain_material(r, material);
        tif_region_to_mesh_region[tif_region_it->index] = r;

        for(unsigned int i=0; i<tif_region_it->boundary.size(); ++i)
        {
          StanfordTIF::Edge_t & edge =  tif_reader->edge(tif_region_it->boundary[i]);
          edge.bc_index = r+1;
        }

        bd_map[tif_region_it->name + "_Neumann"] = std::make_pair(r+1,false);

        r++;
      }


      // process explicit defined boundarys. assign bc_index to these boundary edges
      std::vector<StanfordTIF::Region_t>::const_iterator tif_boundary_it =  tif_reader->region_array().begin();
      for(short int bc_index=tif_reader->region_array().size()+1; tif_boundary_it !=  tif_reader->region_array().end(); ++bc_index, ++tif_boundary_it)
      {
        // it is a region, not boundary, skip it
        if(!tif_boundary_it->segment) continue;

        for(unsigned int i=0; i<tif_boundary_it->boundary.size(); ++i)
          tif_reader->edge(tif_boundary_it->boundary[i]).bc_index = bc_index;

        bd_map[tif_boundary_it->name] = std::make_pair(bc_index, true);
      }
    }







    // fill prisms
    std::vector<StanfordTIF::Prism_t>::const_iterator tif_prism_it  = tif_reader->tif_prisms().begin();
    std::vector<StanfordTIF::Prism_t>::const_iterator tif_prism_it_end  = tif_reader->tif_prisms().end();
    for(; tif_prism_it!=tif_prism_it_end; ++tif_prism_it)
    {
      Elem* elem = mesh.add_elem(Elem::build(PRISM6).release());

      const StanfordTIF::Tri_t & tri = tif_reader->tri(tif_prism_it->tri);
      int z1 = tif_prism_it->z1;
      int z2 = tif_prism_it->z2;

      // prism elem node
      elem->set_node(0) = tif_index_to_node_map[std::make_pair(tri.c1,z1)];
      elem->set_node(1) = tif_index_to_node_map[std::make_pair(tri.c2,z1)];
      elem->set_node(2) = tif_index_to_node_map[std::make_pair(tri.c3,z1)];

      elem->set_node(3) = tif_index_to_node_map[std::make_pair(tri.c1,z2)];
      elem->set_node(4) = tif_index_to_node_map[std::make_pair(tri.c2,z2)];
      elem->set_node(5) = tif_index_to_node_map[std::make_pair(tri.c3,z2)];

      // which region this prism belongs to
      elem->subdomain_id() = tif_prism_it->region;
    }


    mesh.set_n_subdomains() = tif_region_to_mesh_region.size();


    //however, the interface information should be set here

    std::vector<unsigned int>       elems;
    std::vector<unsigned short int> sides;
    std::vector<short int>          bds;

    // get all the boundary element
    mesh.boundary_info->build_side_list (elems, sides, bds);

    //build neighbor information for boundary element. then elem->neighbor() is functional
    mesh.boundary_info->find_neighbors();

    for (size_t nbd=0; nbd<elems.size(); nbd++ )
    {
      // get the element which has boundary/interface side
      const Elem* elem = mesh.elem(elems[nbd]);
      genius_assert(elem->on_boundary() || elem->on_interface());

      //is it an interface side && not a labeled face
      if( elem->neighbor(sides[nbd]) && bds[nbd]<=static_cast<short int>(tif_reader->region_array().size()) )
      {
        // the element and its neighbor should in diffetent subdomain
        unsigned int sbd_id1 = elem->subdomain_id();
        unsigned int sbd_id2 = elem->neighbor(sides[nbd])->subdomain_id();

        // delete the overkilled boundary side
        if (sbd_id1 == sbd_id2)
        {
          mesh.boundary_info->remove(elem, sides[nbd]);
          continue;
        }

        // the side should be an interface side
        genius_assert(elem->on_interface());
        genius_assert(elem->neighbor(sides[nbd])->on_interface());

        //remove the pair-element from boundary
        mesh.boundary_info->remove(elem, sides[nbd]);
        mesh.boundary_info->remove(elem->neighbor(sides[nbd]),
                                   elem->neighbor(sides[nbd])->which_neighbor_am_i(elem));

        // build the label for the interface, which has the form of RegionLabel1_to_RegionLabel2,
        // the two region is alpha ordered.
        std::string bd_label;
        if( tif_reader->region(sbd_id1).name < tif_reader->region(sbd_id2).name)
          bd_label = tif_reader->region(sbd_id1).name + "_to_" + tif_reader->region(sbd_id2).name;
        else
          bd_label = tif_reader->region(sbd_id2).name + "_to_" + tif_reader->region(sbd_id1).name;

        short int bd_index;

        // if the label already exist
        if( bd_map.find(bd_label) != bd_map.end() )
          bd_index = (*bd_map.find(bd_label)).second.first;
        else
        {
          //else, increase bd_index, insert it into bd_map
          bd_index = tif_reader->region_array().size() + bd_map.size() + 1;
          bd_map[bd_label] = std::make_pair(bd_index, false);
        }

        // add pair-element to boundary with new bd_index
        mesh.boundary_info->add_side(elem, sides[nbd], bd_index);
        mesh.boundary_info->add_side(elem->neighbor(sides[nbd]),
                                     elem->neighbor(sides[nbd])->which_neighbor_am_i(elem),
                                     bd_index);
      }
    }

    // after the previous remove and insert operation, the number of boundary id maybe changed. renumber it here.
    mesh.boundary_info->rebuild_ids();


    //write down bd labels
    for(Bd_It bd_it = bd_map.begin(); bd_it != bd_map.end(); ++bd_it)
    {
      mesh.boundary_info->set_label_to_id( bd_it->second.first, bd_it->first, bd_it->second.second);
    }


    // magic number, for 3D mesh, should > 2008
    mesh.magic_num() = 3121;
  }

  Parallel::broadcast(tif_region_to_mesh_region , 0);

  /*
   * set mesh structure for all processors, and build simulation system
   */

  // broadcast mesh to all the processor
  MeshCommunication mesh_comm;
  mesh_comm.broadcast(mesh);


  // build simulation system
  system.build_simulation_system();
  system.sync_print_info();

  /*
   * after that, set doping infomation here. this should be done for all the processors
   */

  // node id to TIF index map, this should be same for all processors
  std::map<unsigned int, std::pair<int, int> > node_id_to_tif_index_map;
  {
    std::vector<unsigned int> node_id;
    std::vector<int> c_array;
    std::vector<int> z_array;
    // fill node_id_to_tif_index_map by processor 0
    if( Genius::processor_id() == 0)
    {
      std::map<Node *, std::pair<int, int> >::iterator it = node_to_tif_index_map.begin();
      for(; it != node_to_tif_index_map.end(); ++it)
      {
        node_id.push_back(it->first->id());
        c_array.push_back(it->second.first);
        z_array.push_back(it->second.second);
      }
    }

    //broadcast node_id_to_tif_index_map and tif_region_to_mesh_region to all the processors
    Parallel::broadcast(node_id, 0);
    Parallel::broadcast(c_array, 0);
    Parallel::broadcast(z_array, 0);

    for(unsigned int n=0; n<node_id.size(); ++n)
      node_id_to_tif_index_map[node_id[n]] = std::make_pair(c_array[n], z_array[n]);
  }



  // ok, we had got enough informations for set up each simulation region
  std::map<std::pair<std::pair<int, int>, int>, unsigned int> solution_map; // map <node_index, z_index, region_index> to data_index
  typedef std::map<std::pair<std::pair<int, int>, int>, unsigned int>::iterator Solution_It;

  for(unsigned int n=0; n<tif_reader->sol_data_array().size(); ++n)
  {
    int tif_region = tif_reader->sol_data(n).region_index;
    int z = tif_reader->sol_data(n).zplane;
    int index = tif_reader->sol_data(n).index;

    int region = tif_region_to_mesh_region.find(tif_region)->second;
    std::pair<std::pair<int, int>, int> key = std::make_pair(std::make_pair(index,z), region);
    solution_map.insert(std::make_pair(key, n));
  }

  const double concentration_scale = pow(cm, -3);

  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);


    switch ( region->type() )
    {
    case SemiconductorRegion :
      {
        bool sigle   = Material::IsSingleCompSemiconductor(region->material());
        bool complex = Material::IsComplexCompSemiconductor(region->material());

        SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
        SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
        for(; node_it!=node_it_end; ++node_it)
        {
          FVM_Node * fvm_node = (*node_it);
          FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);

          // tif_node_index is the index in TIF file that this FVM node has
          std::pair<int, int> tif_node_index = node_id_to_tif_index_map[fvm_node->root_node()->id()];
          int region_index = r;
          std::pair<std::pair<int, int>, int> key = std::make_pair(tif_node_index, region_index);
          if(solution_map.find(key) != solution_map.end())
          {
            unsigned int data_index = solution_map.find(key)->second;

            // doping
            {
              node_data->Na()   = tif_reader->acceptor(data_index) * concentration_scale;
              node_data->Nd()   = tif_reader->donor(data_index) * concentration_scale;
            }

            // mole fraction
            if(sigle)
            {
              node_data->mole_x()   = tif_reader->mole_x(data_index);
            }

            if(complex)
            {
              node_data->mole_x()   = tif_reader->mole_x(data_index);
              node_data->mole_y()   = tif_reader->mole_y(data_index);
            }
          }

        }
        region->init(system.T_external());
        break;
      }
    case InsulatorRegion     :
      {
        region->init(system.T_external());
        break;
      }
    case ElectrodeRegion     :
      {
        region->init(system.T_external());
        break;
      }
    case MetalRegion    :
      {
        region->init(system.T_external());
        break;
      }
    case VacuumRegion        :
      {
        region->init(system.T_external());
        break;
      }
    case PMLRegion           :
      {
        region->init(system.T_external());
        break;
      }
    default: genius_error();
    }
  }




  system.init_region_post_process();

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}




void STIFIO::write (const std::string& filename)
{
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();
  if( system.dim() != 2 ) return;

  const MeshBase & mesh = system.mesh();

  std::map<std::string,  std::vector<unsigned int> >  region_nodes_map;
  std::map<std::string,  std::map<std::string, std::vector<PetscScalar> > > region_solution_map;

  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);
    if(region->type() != SemiconductorRegion) continue;

    region->region_node(region_nodes_map[region->name()]);
    region->get_variable_data("Na", POINT_CENTER, region_solution_map[region->name()]["Na"]) ;
    region->get_variable_data("Nd", POINT_CENTER, region_solution_map[region->name()]["Nd"]) ;
    region->get_variable_data("mole_x", POINT_CENTER, region_solution_map[region->name()]["mole_x"]);
    region->get_variable_data("mole_y", POINT_CENTER, region_solution_map[region->name()]["mole_y"]);
  }


  if(Genius::processor_id() == 0)
  {
    MediciTIF * tif_writer = new MediciTIF();

    for (unsigned int n=0; n<mesh.n_nodes(); n++)
    {
      Point p = mesh.point(n);

      StanfordTIF::Node_t node;
      node.index = n;
      node.x = p[0]/um;
      node.y = p[1]/um;
      node.z = node.h = 0.0;

      tif_writer->tif_nodes().push_back(node);
    }

    //classfy boundary label
    std::vector<unsigned int>       elems;
    std::vector<unsigned short int> sides;
    std::vector<short int>          bds;

    std::set< std::pair<unsigned int, unsigned int> > edges;
    std::map<std::string, std::vector< std::pair<unsigned int, unsigned int> > > bd_edges;
    std::map<unsigned int, std::vector< std::pair<unsigned int, unsigned int> > > region_bd_edges;


    // get all the boundary element
    mesh.boundary_info->build_side_list (elems, sides, bds);
    for(unsigned int b=0; b<bds.size(); ++b)
    {
      short int bd = bds[b];
      std::string label = mesh.boundary_info->get_label_by_id(bd);
      bool user_define = mesh.boundary_info->boundary_id_has_user_defined_label(bd);

      const Elem * elem = mesh.elem(elems[b]);
      unsigned int region = elem->subdomain_id();
      unsigned short int side = sides[b];

      AutoPtr<Elem> side_face = elem->build_side(side);
      unsigned int edge_p1 = side_face->get_node(0)->id();
      unsigned int edge_p2 = side_face->get_node(1)->id();
      if(edge_p1 > edge_p2) std::swap(edge_p1, edge_p2);
      edges.insert(std::make_pair(edge_p1, edge_p2));

      region_bd_edges[region].push_back(std::make_pair(edge_p1, edge_p2));

      if(user_define)
      {
        bd_edges[label].push_back(std::make_pair(edge_p1, edge_p2));
      }
    }


    std::set< std::pair<unsigned int, unsigned int> >::iterator edge_it=edges.begin();
    for(; edge_it!=edges.end(); ++edge_it)
    {
       const std::pair<unsigned int, unsigned int> & edge_points = *edge_it;
       StanfordTIF:: Edge_t edge;
       edge.index = std::distance(edges.begin(), edge_it);
       edge.point1 = edge_points.first;
       edge.point2 = edge_points.second;

       tif_writer->tif_edges().push_back(edge);
    }


    // regions
    for(unsigned int r=0; r<mesh.n_subdomains(); ++r)
    {
       std::string label = mesh.subdomain_label_by_id(r);
       StanfordTIF::Region_t region;
       region.index = r;
       region.name = label;
       region.material = mesh.subdomain_material(r);
       region.segment = false;

       const std::vector< std::pair<unsigned int, unsigned int> > & region_edges = region_bd_edges.find(r)->second;
       for(unsigned int b=0; b<region_edges.size(); ++b)
       {
         const std::pair<unsigned int, unsigned int> & edge_points = region_edges[b];
         unsigned int edge_index = std::distance(edges.begin(),  edges.find(edge_points));
         region.boundary.push_back(edge_index);
       }

       tif_writer->region_array().push_back(region);
    }


    std::map<std::string, std::vector< std::pair<unsigned int, unsigned int> > >::const_iterator bd_edge_it = bd_edges.begin();
    for(int b=0; bd_edge_it != bd_edges.end(); b++, bd_edge_it++)
    {
       StanfordTIF::Region_t region;
       region.index = b;
       region.name = bd_edge_it->first;
       region.material = "surface";
       region.segment = true;

       const std::vector< std::pair<unsigned int, unsigned int> > & seg_edges = bd_edge_it->second;
       for(unsigned int b=0; b<seg_edges.size(); ++b)
       {
         const std::pair<unsigned int, unsigned int> & edge_points = seg_edges[b];
         unsigned int edge_index = std::distance(edges.begin(),  edges.find(edge_points));
         region.boundary.push_back(edge_index);
       }

       tif_writer->region_array().push_back(region);
    }

    MeshBase::const_element_iterator       el  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end();
    for (; el != end; ++el)
    {
      const Elem * elem = *el;
      switch(elem->type())
      {
        case TRI3:
        case TRI3_FVM:
        case TRI3_CY_FVM:
        {
           StanfordTIF::Tri_t tri;
           tri.index = elem->id();
           tri.region = elem->subdomain_id();
           tri.c1 = elem->get_node(0)->id();
           tri.c2 = elem->get_node(1)->id();
           tri.c3 = elem->get_node(2)->id();
           tri.c4 = -1;
           tri.t1 = elem->neighbor(1) ? elem->neighbor(1)->id() : -1024;
           tri.t2 = elem->neighbor(2) ? elem->neighbor(2)->id() : -1024;
           tri.t3 = elem->neighbor(0) ? elem->neighbor(0)->id() : -1024;
           tri.t4 = -1;

           tif_writer->tif_tris().push_back(tri);
           break;
        }
        case QUAD4:
        case QUAD4_FVM:
        case QUAD4_CY_FVM:
        {
           StanfordTIF::Tri_t tri;
           tri.index = elem->id();
           tri.region = elem->subdomain_id();
           tri.c1 = elem->get_node(0)->id();
           tri.c2 = elem->get_node(1)->id();
           tri.c3 = elem->get_node(2)->id();
           tri.c4 = elem->get_node(3)->id();
           tri.t1 = elem->neighbor(0) ? elem->neighbor(0)->id() : -1024;
           tri.t2 = elem->neighbor(1) ? elem->neighbor(1)->id() : -1024;
           tri.t3 = elem->neighbor(2) ? elem->neighbor(2)->id() : -1024;
           tri.t4 = elem->neighbor(3) ? elem->neighbor(3)->id() : -1024;

           tif_writer->tif_tris().push_back(tri);
           break;
        }
        default: genius_error();
      }
    }

    StanfordTIF::SolHead_t & sol_head = tif_writer->sol_head();
    sol_head.sol_num = 4;
    sol_head.sol_name_array.push_back("Net");
    sol_head.sol_name_array.push_back("Total");
    sol_head.sol_name_array.push_back("mole_x");
    sol_head.sol_name_array.push_back("mole_y");
    sol_head.sol_unit_array.push_back("cm^-3");
    sol_head.sol_unit_array.push_back("cm^-3");
    sol_head.sol_unit_array.push_back("1");
    sol_head.sol_unit_array.push_back("1");

    for(unsigned int r=0; r<system.n_regions(); r++)
    {
      const SimulationRegion * region = system.region(r);
      if(region->type() != SemiconductorRegion) continue;

      const std::vector<unsigned int> & nodes = region_nodes_map[region->name()];
      const std::vector<PetscScalar> & na = region_solution_map[region->name()]["Na"];
      const std::vector<PetscScalar> & nd = region_solution_map[region->name()]["Nd"];
      const std::vector<PetscScalar> & mx = region_solution_map[region->name()]["mole_x"];
      const std::vector<PetscScalar> & my = region_solution_map[region->name()]["mole_y"];

      for(unsigned int n=0; n<nodes.size(); ++n)
      {
        StanfordTIF::SolData_t sol;
        sol.index = nodes[n];
        sol.material = region->material();
        sol.data_array.push_back(nd[n]-na[n]);
        sol.data_array.push_back(nd[n]+na[n]);
        sol.data_array.push_back(mx[n]);
        sol.data_array.push_back(my[n]);

        tif_writer->sol_data_array().push_back(sol);
      }
    }

    tif_writer->export_tif(filename);
    delete tif_writer;
  }
}

