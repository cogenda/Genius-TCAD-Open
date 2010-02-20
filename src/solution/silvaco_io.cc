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
#include "silvaco_io.h"

#include "mesh.h"
#include "boundary_info.h"
#include "mesh_communication.h"
#include "parallel.h"
#include "simulation_system.h"
#include "show_mesh_2d.h"


using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::V;
using PhysicalUnit::K;
using PhysicalUnit::eV;

/**
 * This method implements reading a mesh from a specified file
 * in Silvaco TIF format.
 */
void STIFIO::read (const std::string& filename)
{
  SilvacoTIF tif_reader(filename);

  // read silvaco tif file
  int ierr;
  if( Genius::processor_id() == 0)
    ierr = tif_reader.read();

  Parallel::broadcast(ierr);
  if(ierr) return;

  /*
   * after that, we fill mesh structure with mesh information read from TIF
   */
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  Mesh & mesh = system.mesh();

  // clear the system
  system.clear();

  // map node * to tif index in TIF::node_array
  std::map<Node *, int>       node_to_tif_index_map;

  if( Genius::processor_id() == 0)
  {
    // fill node location
    std::vector<SilvacoTIF::Node_t>::const_iterator tif_node_it = tif_reader.tif_nodes().begin();
    for(int i=0; tif_node_it!=tif_reader.tif_nodes().end(); ++i, ++tif_node_it)
    {
      Node * node = mesh.add_point( Point(tif_node_it->x*um, tif_node_it->y*um, 0) );
      node_to_tif_index_map[node] = i;
    }

    // map bc_index to bc label
    std::map<const std::string, short int> bd_map;
    typedef std::map<const std::string, short int>::iterator Bd_It;


    // fill region label and region material
    // at the same time. set all the remaining boundary edges as "region_neumann"
    std::map<int, int> tif_region_to_mesh_region;
    std::vector<SilvacoTIF::Region_t>::const_iterator tif_region_it  = tif_reader.region_array().begin();
    for(int r=0; tif_region_it!=tif_reader.region_array().end(); ++tif_region_it)
    {
      if(tif_region_it->segment) continue;

      std::string material = tif_region_it->material;
      mesh.set_subdomain_label(r, tif_region_it->name );
      mesh.set_subdomain_material(r, material);
      tif_region_to_mesh_region[tif_region_it->index] = r;

      for(unsigned int i=0; i<tif_region_it->boundary.size(); ++i)
      {
        SilvacoTIF::Edge_t & edge =  tif_reader.edge(tif_region_it->boundary[i]);
        edge.bc_index = r+1;
      }

      bd_map[tif_region_it->name + "_Neumann"] = r+1;

      ++r;
    }
    mesh.set_n_subdomains() = tif_region_to_mesh_region.size();

    // process explicit defined boundarys. assign bc_index to these boundary edges
    std::vector<SilvacoTIF::Region_t>::iterator tif_boundary_it =  tif_reader.region_array().begin();
    for(short int bc_index=tif_reader.region_array().size()+1; tif_boundary_it !=  tif_reader.region_array().end(); ++bc_index, ++tif_boundary_it)
    {
      // it is a region, not boundary, skip it
      if(!tif_boundary_it->segment) continue;

      for(unsigned int i=0; i<tif_boundary_it->boundary.size(); ++i)
        tif_reader.edge(tif_boundary_it->boundary[i]).bc_index = bc_index;

      bd_map[tif_boundary_it->name] = bc_index;
    }


    // build edge map
    std::map<SilvacoTIF::Edge_t, int, SilvacoTIF::lt_edge> edge_table;
    std::vector<SilvacoTIF::Edge_t>::iterator  tif_edge_it = tif_reader.tif_edges().begin();
    for(; tif_edge_it!=tif_reader.tif_edges().end(); ++tif_edge_it)
    {
      genius_assert( tif_edge_it->bc_index != BoundaryInfo::invalid_id);
      edge_table[*tif_edge_it] = tif_edge_it->bc_index;
    }

    // fill triangles
    std::vector<SilvacoTIF::Tri_t>::const_iterator tif_tri_it  = tif_reader.tif_tris().begin();
    for(; tif_tri_it!=tif_reader.tif_tris().end(); ++tif_tri_it)
    {
      Elem* elem = mesh.add_elem(Elem::build(TRI3).release());

      // tri elem node
      elem->set_node(0) = mesh.node_ptr( tif_tri_it->c1 );
      elem->set_node(1) = mesh.node_ptr( tif_tri_it->c2 );
      elem->set_node(2) = mesh.node_ptr( tif_tri_it->c3 );

      // which region this tri belongs to
      elem->subdomain_id() = tif_region_to_mesh_region[tif_tri_it->region];

      SilvacoTIF::Edge_t edge;

      // process edge0
      edge.point1 = tif_tri_it->c1;
      edge.point2 = tif_tri_it->c2;
      if(edge.point1 > edge.point2) std::swap(edge.point1, edge.point2);
      std::map<SilvacoTIF::Edge_t, int, SilvacoTIF::lt_edge>::iterator edge_pointer = edge_table.find(edge);
      if( edge_pointer != edge_table.end() )
      {
        // this edge should on external boundary or region interface
        genius_assert( (*tif_tri_it).t3<0 || (*tif_tri_it).region != tif_reader.tri((*tif_tri_it).t3).region );
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
        genius_assert( (*tif_tri_it).t1<0 || (*tif_tri_it).region != tif_reader.tri((*tif_tri_it).t1).region );
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
        genius_assert( (*tif_tri_it).t2<0 || (*tif_tri_it).region != tif_reader.tri((*tif_tri_it).t2).region );
        mesh.boundary_info->add_side(elem, 2, edge_pointer->second);
      }
    }

    //however, the interface information should be set here

    std::vector<unsigned int>       elems;
    std::vector<unsigned short int> sides;
    std::vector<short int>          bds;

    // get all the boundary element
    mesh.boundary_info->build_side_list (elems, sides, bds);

    //build neighbor information for mesh. then elem->neighbor() is functional
    mesh.find_neighbors();

    for (size_t nbd=0; nbd<elems.size(); nbd++ )
    {
      // get the element which has boundary/interface side
      const Elem* elem = mesh.elem(elems[nbd]);
      genius_assert(elem->on_boundary() || elem->on_interface());

      //is it an interface side && not a label face
      if( elem->neighbor(sides[nbd])!=NULL && bds[nbd]<=static_cast<short int>(tif_reader.region_array().size()) )
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
        if( tif_reader.region(sbd_id1).name < tif_reader.region(sbd_id2).name)
          bd_label = tif_reader.region(sbd_id1).name + "_to_" + tif_reader.region(sbd_id2).name;
        else
          bd_label = tif_reader.region(sbd_id2).name + "_to_" + tif_reader.region(sbd_id1).name;

        short int bd_index;

        // if the label already exist
        if( bd_map.find(bd_label) != bd_map.end() )
          bd_index = (*bd_map.find(bd_label)).second;
        else
        {
          //else, increase bd_index, insert it into bd_map
          bd_index = tif_reader.region_array().size() + bd_map.size() + 1;
          bd_map.insert(std::pair<const std::string, short int>(bd_label,bd_index));
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
    Bd_It bd_it = bd_map.begin();
    for(; bd_it != bd_map.end(); ++bd_it)
    {
      mesh.boundary_info->set_label_to_id( (*bd_it).second, (*bd_it).first );
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

  // broadcast node_id_to_tif_index_map to all the processors
  Parallel::broadcast(node_id_to_tif_index_map , 0);

  // broadcast SolHead_t to all processors
  Parallel::broadcast(tif_reader.sol_head().sol_num, 0);
  if(Genius::processor_id() != 0)
    tif_reader.sol_head().sol_name_array.resize(tif_reader.sol_head().sol_num);
  for(int n=0; n<tif_reader.sol_head().sol_num; ++n)
    Parallel::broadcast(tif_reader.sol_head().sol_name_array[n]);

  //broadcast SolData to all processors
  unsigned int n_solution = tif_reader.sol_data_array().size();
  Parallel::broadcast(n_solution, 0);
  if(Genius::processor_id() != 0)
    tif_reader.sol_data_array().resize(n_solution);
  for(unsigned int n=0; n < n_solution; ++n)
  {
    Parallel::broadcast(tif_reader.sol_data(n).index);
    Parallel::broadcast(tif_reader.sol_data(n).material);
    Parallel::broadcast(tif_reader.sol_data(n).data_array);
  }

  // ok, we had got enough informations for set up each simulation region
  std::multimap<int,  SilvacoTIF::SolData_t> solution_map;
  typedef std::multimap<int,  SilvacoTIF::SolData_t>::iterator Solution_It;
  for(unsigned int n=0; n<n_solution; ++n)
    solution_map.insert(std::pair<int,  SilvacoTIF::SolData_t>(tif_reader.sol_data(n).index, tif_reader.sol_data(n)));

  unsigned int donor        = tif_reader.sol_head().solution_index("Donor");
  unsigned int acceptor     = tif_reader.sol_head().solution_index("Acceptor");

  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);

    switch ( region->type() )
    {
    case SemiconductorRegion :
      {
        SimulationRegion::node_iterator node_it = region->nodes_begin();
        SimulationRegion::node_iterator node_it_end = region->nodes_end();
        for(; node_it!=node_it_end; ++node_it)
        {
          FVM_Node * fvm_node = (*node_it).second;
          if( !fvm_node->root_node()->on_local() ) continue;

          FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);

          // tif_node_index is the index in TIF file that this FVM node has
          int tif_node_index = node_id_to_tif_index_map[(*node_it).first];
          // however, one TIF node may has several solution data when it lies on region interface, we should find corrent
          // solution data for this region which has the same material string as this region.
          std::pair<Solution_It, Solution_It> sol_it_pair = solution_map.equal_range(tif_node_index);
          Solution_It sol_it = sol_it_pair.first;
          for(; sol_it!=sol_it_pair.second; ++sol_it)
            if((*sol_it).second.material == region->material()) break;

          if( sol_it != sol_it_pair.second && sol_it != solution_map.end())
          {
            // doping
            if(donor!=invalid_uint && acceptor!=invalid_uint)
            {
              node_data->Na()   = (*sol_it).second.data_array[acceptor] * pow(cm, -3);
              node_data->Nd()   = (*sol_it).second.data_array[donor   ] * pow(cm, -3);
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
    case ConductorRegion     :
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

}












