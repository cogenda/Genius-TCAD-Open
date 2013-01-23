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

// Local includes
#include "silvaco.h"
#include "medici.h"

#include "stanford_io.h"

#include "mesh_base.h"
#include "boundary_info.h"
#include "mesh_communication.h"
#include "simulation_region.h"
#include "log.h"
#include "parallel.h"
#include "material.h"

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


  // broadcast
  tif_reader->broadcast();

  if(tif_reader->dim() == 2)
    read_2d(tif_reader);
  else
    read_3d_silvaco(tif_reader);

}



void STIFIO::read_2d (StanfordTIF * tif_reader)
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

        std::string material = tif_region_it->material;
        mesh.set_subdomain_label(r, tif_region_it->name );
        mesh.set_subdomain_material(r, material);
        tif_region_to_mesh_region[tif_region_it->index] = r;

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
      Elem* elem = mesh.add_elem(Elem::build(TRI3).release());

      // tri elem node
      elem->set_node(0) = mesh.node_ptr( tif_tri_it->c1 );
      elem->set_node(1) = mesh.node_ptr( tif_tri_it->c2 );
      elem->set_node(2) = mesh.node_ptr( tif_tri_it->c3 );

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
        genius_assert( tif_tri_it->t3<0 || tif_tri_it->region != tif_reader->tri(tif_tri_it->t3).region );
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
        genius_assert( tif_tri_it->t1<0 || tif_tri_it->region != tif_reader->tri(tif_tri_it->t1).region );
        mesh.boundary_info->add_side(elem, 1, edge_pointer->second);
      }

      // process edge2
      edge.point1 = tif_tri_it->c3;
      edge.point2 = tif_tri_it->c1;
      if(edge.point1 > edge.point2) std::swap(edge.point1, edge.point2);
      edge_pointer = edge_table.find(edge);
      if( edge_pointer != edge_table.end() )
      {
        // this edge should on external boundary or region interface
        genius_assert( tif_tri_it->t2<0 || tif_tri_it->region != tif_reader->tri(tif_tri_it->t2).region );
        mesh.boundary_info->add_side(elem, 2, edge_pointer->second);
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
  Parallel::broadcast(tif_region_to_mesh_region , 0);



  // ok, we had got enough informations for set up each simulation region
  std::map<std::pair<int, int>, unsigned int> solution_map; // map <node_index, region_index> to data_index
  typedef std::map<std::pair<int, int>, unsigned int>::iterator Solution_It;

  for(unsigned int n=0; n<tif_reader->sol_data_array().size(); ++n)
  {
    int tif_region = tif_reader->sol_data(n).region_index;
    int region = tif_region_to_mesh_region.find(tif_region)->second;
    std::pair<int, int> key = std::make_pair(tif_reader->sol_data(n).index, region);
    solution_map.insert(std::make_pair(key, n));
  }


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
          int tif_node_index = node_id_to_tif_index_map[fvm_node->root_node()->id()];
          int region_index = r;
          std::pair<int, int> key = std::make_pair(tif_node_index, region_index);
          assert(solution_map.find(key) != solution_map.end());
          unsigned int data_index = solution_map.find(key)->second;

          // doping
          {
            node_data->Na()   = tif_reader->acceptor(data_index) * pow(cm, -3);
            node_data->Nd()   = tif_reader->donor(data_index) * pow(cm, -3);
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


  delete tif_reader;


  system.init_region_post_process();

}


void STIFIO::read_3d_silvaco (StanfordTIF * tif_reader)
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

        std::string material = tif_region_it->material;
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


  delete tif_reader;


  system.init_region_post_process();


}







