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
#include "tif_data.h"
#include "tif_io.h"
#include "mesh_base.h"
#include "boundary_info.h"
#include "mesh_communication.h"
#include "simulation_region.h"
#include "parallel.h"
#include "material.h"

using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::V;
using PhysicalUnit::K;
using PhysicalUnit::eV;

/**
 * This method implements reading a mesh from a specified file
 * in TIF format.
 */
void TIFIO::read (const std::string& filename)
{
  /*
   * first, we call yyparse to read tif file
   */
  if( Genius::processor_id() == 0)
  {
    TIF::yyin = fopen(filename.c_str(), "r");
    genius_assert( TIF::yyin != NULL );

    genius_assert( !TIF::yyparse() );

    fclose(TIF::yyin);
  }

  /*
   * after that, we fill mesh structure with mesh information read from TIF
   */
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  // clear the system
  system.clear();

  // map node * to tif index in TIF::node_array
  std::map<Node *, int>       node_to_tif_index_map;

  if( Genius::processor_id() == 0)
  {
    // fill node location
    std::vector<TIF::Node_t>::iterator tif_node_it = TIF::node_array.begin();
    for(int i=0; tif_node_it!=TIF::node_array.end(); ++i, ++tif_node_it)
    {
      Node * node = mesh.add_point( Point((*tif_node_it).x*um, (*tif_node_it).y*um, 0) );
      node_to_tif_index_map[node] = i;
    }

    // map bc_index to bc label
    std::map<std::string, std::pair<short int, bool> > bd_map;
    typedef std::map<std::string, std::pair<short int, bool> >::iterator Bd_It;


    // fill region label and region material
    // at the same time. set all the remaining boundary edges as "region_neumann"
    mesh.set_n_subdomains() = TIF::region_array.size();
    std::vector<TIF::Region_t>::iterator tif_region_it  = TIF::region_array.begin();
    for(int r=1; tif_region_it!=TIF::region_array.end(); ++r, ++tif_region_it)
    {
      std::string material = (*tif_region_it).type;
      mesh.set_subdomain_label( (*tif_region_it).index, (*tif_region_it).name );
      mesh.set_subdomain_material( (*tif_region_it).index, material);

      for(unsigned int i=0; i<(*tif_region_it).boundary.size(); ++i)
      {
        TIF::Edge_t & edge =  TIF::edge_array[(*tif_region_it).boundary[i]];
        edge.bc_index = r;
      }

      bd_map[(*tif_region_it).name + "_Neumann"] = std::make_pair(r, false);
    }

    // process explicit defined boundarys. assign bc_index to these boundary edges
    std::vector<TIF::Interface_t>::iterator tif_interface_it =  TIF::interface_array.begin();
    for(short int bc_index=TIF::region_array.size()+1; tif_interface_it !=  TIF::interface_array.end(); ++bc_index, ++tif_interface_it)
    {
      // it is a region, not boundary, skip it
      if( (*tif_interface_it).region >= 0)
      {
        genius_assert((*tif_interface_it).boundary.size()==0);
        continue;
      }

      for(unsigned int i=0; i<(*tif_interface_it).boundary.size(); ++i)
        TIF::edge_array[ (*tif_interface_it).boundary[i] ].bc_index = bc_index;

      bd_map[(*tif_interface_it).name] = std::make_pair(bc_index, true);
    }


    // build edge map
    std::map<TIF::Edge_t, int, TIF::lt_edge> edge_table;
    std::vector<TIF::Edge_t>::iterator  tif_edge_it = TIF::edge_array.begin();
    for(; tif_edge_it!=TIF::edge_array.end(); ++tif_edge_it)
    {
      genius_assert( (*tif_edge_it).bc_index != BoundaryInfo::invalid_id);
      edge_table[*tif_edge_it] = (*tif_edge_it).bc_index;
    }

    // fill triangles
    std::vector<TIF::Tri_t>::iterator tif_tri_it  = TIF::tri_array.begin();
    for(; tif_tri_it!=TIF::tri_array.end(); ++tif_tri_it)
    {
      Elem* elem = mesh.add_elem(Elem::build(TRI3).release());

      // tri elem node
      elem->set_node(0) = mesh.node_ptr( tif_tri_it->c1 );
      elem->set_node(1) = mesh.node_ptr( tif_tri_it->c2 );
      elem->set_node(2) = mesh.node_ptr( tif_tri_it->c3 );

      // which region this tri belongs to
      elem->subdomain_id() = (*tif_tri_it).region;

      std::vector< std::pair<unsigned int, unsigned int> > edges;
      edges.push_back( std::make_pair(tif_tri_it->c1, tif_tri_it->c2) );
      edges.push_back( std::make_pair(tif_tri_it->c2, tif_tri_it->c3) );
      edges.push_back( std::make_pair(tif_tri_it->c3, tif_tri_it->c1) );

      for(unsigned int n=0; n<edges.size(); ++n)
      {
        TIF::Edge_t edge;
        edge.point1 = edges[n].first;
        edge.point2 = edges[n].second;
        if( edge.point1 > edge.point2 )
          std::swap(edge.point1, edge.point2);

        std::map<TIF::Edge_t, int, TIF::lt_edge>::iterator edge_pointer = edge_table.find(edge);
        if( edge_pointer != edge_table.end() )
          mesh.boundary_info->add_side(elem, n, edge_pointer->second);
      }
    }


    // fill quadrangles
    std::vector<TIF::Quad_t>::iterator tif_quad_it  = TIF::quad_array.begin();
    for(; tif_quad_it!=TIF::quad_array.end(); ++tif_quad_it)
    {
      Elem* elem = mesh.add_elem(Elem::build(QUAD4).release());

      // quadrangle elem node
      elem->set_node(0) = mesh.node_ptr( tif_quad_it->c1 );
      elem->set_node(1) = mesh.node_ptr( tif_quad_it->c2 );
      elem->set_node(2) = mesh.node_ptr( tif_quad_it->c3 );
      elem->set_node(3) = mesh.node_ptr( tif_quad_it->c4 );

      // which region this tri belongs to
      elem->subdomain_id() = (*tif_quad_it).region;

      std::vector< std::pair<unsigned int, unsigned int> > edges;
      edges.push_back( std::make_pair(tif_quad_it->c1, tif_quad_it->c2) );
      edges.push_back( std::make_pair(tif_quad_it->c2, tif_quad_it->c3) );
      edges.push_back( std::make_pair(tif_quad_it->c3, tif_quad_it->c4) );
      edges.push_back( std::make_pair(tif_quad_it->c4, tif_quad_it->c1) );

      for(unsigned int n=0; n<edges.size(); ++n)
      {
        TIF::Edge_t edge;
        edge.point1 = edges[n].first;
        edge.point2 = edges[n].second;
        if( edge.point1 > edge.point2 )
          std::swap(edge.point1, edge.point2);

        std::map<TIF::Edge_t, int, TIF::lt_edge>::iterator edge_pointer = edge_table.find(edge);
        if( edge_pointer != edge_table.end() )
          mesh.boundary_info->add_side(elem, n, edge_pointer->second);
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

      // is it a neumann boundary?
      if( elem->neighbor(sides[nbd]) == NULL && bds[nbd]<=static_cast<short int>(TIF::region_array.size()) )
      {
        unsigned int sbd_id = elem->subdomain_id();
        assert(elem->on_boundary());
        //remove the element from boundary
        mesh.boundary_info->remove(elem, sides[nbd]);

        std::string bd_label = TIF::region_array[sbd_id].name + "_Neumann";

        short int bd_index;

        // if the label already exist
        if( bd_map.find(bd_label) != bd_map.end() )
          bd_index = (*bd_map.find(bd_label)).second.first;
        else
        {
          //else, increase bd_index, insert it into bd_map
          bd_index = sbd_id + 1;
          bd_map[bd_label] = std::make_pair(bd_index, false);
        }

        // add pair-element to boundary with new bd_index
        mesh.boundary_info->add_side(elem, sides[nbd], bd_index);
      }

      //is it an interface side && not a label face
      if( elem->neighbor(sides[nbd])!=NULL && bds[nbd]<=static_cast<short int>(TIF::region_array.size()) )
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
        if( TIF::region_array[sbd_id1].name < TIF::region_array[sbd_id2].name)
          bd_label = TIF::region_array[sbd_id1].name + "_to_" + TIF::region_array[sbd_id2].name;
        else
          bd_label = TIF::region_array[sbd_id2].name + "_to_" + TIF::region_array[sbd_id1].name;

        short int bd_index;

        // if the label already exist
        if( bd_map.find(bd_label) != bd_map.end() )
          bd_index = (*bd_map.find(bd_label)).second.first;
        else
        {
          //else, increase bd_index, insert it into bd_map
          bd_index = TIF::interface_array.size() + TIF::region_array.size() + bd_map.size() + 1;
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
    Bd_It bd_it = bd_map.begin();
    for(; bd_it != bd_map.end(); ++bd_it)
    {
      mesh.boundary_info->set_label_to_id( (*bd_it).second.first, (*bd_it).first, (*bd_it).second.second );
    }

    // magic number, for 2D mesh, should < 2008
    mesh.magic_num() = 192;
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
   * after that, set doping infomation, this should be done for all the processors
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
  Parallel::broadcast(TIF::sol_head.sol_num, 0);
  if(Genius::processor_id() != 0)
    TIF::sol_head.sol_name_array.resize(TIF::sol_head.sol_num);
  for(int n=0; n<TIF::sol_head.sol_num; ++n)
    Parallel::broadcast(TIF::sol_head.sol_name_array[n]);

  //broadcast SolData to all processors
  unsigned int n_solution = TIF::sol_data.size();
  Parallel::broadcast(n_solution, 0);
  if(Genius::processor_id() != 0)
    TIF::sol_data.resize(n_solution);
  for(unsigned int n=0; n < n_solution; ++n)
  {
    Parallel::broadcast(TIF::sol_data[n].index);
    Parallel::broadcast(TIF::sol_data[n].material);
    Parallel::broadcast(TIF::sol_data[n].data_array);
  }


  // ok, we had got enough informations for set up each simulation region
  std::multimap<int,  TIF::SolData_t> solution_map;
  typedef std::multimap<int,  TIF::SolData_t>::iterator Solution_It;
  for(unsigned int n=0; n<TIF::sol_data.size(); ++n)
    solution_map.insert(std::pair<int,  TIF::SolData_t>(TIF::sol_data[n].index, TIF::sol_data[n]));

  unsigned int net          = TIF::sol_head.solution_index("Net");
  unsigned int total        = TIF::sol_head.solution_index("Total");

  unsigned int psi          = TIF::sol_head.solution_index("v");
  unsigned int elec_density = TIF::sol_head.solution_index("n");
  unsigned int hole_density = TIF::sol_head.solution_index("p");

  unsigned int latt_temp    = TIF::sol_head.solution_index("tl");
  unsigned int elec_temp    = TIF::sol_head.solution_index("tn");
  unsigned int hole_temp    = TIF::sol_head.solution_index("tp");

  unsigned int mole_x       = TIF::sol_head.solution_index("mole_x");
  unsigned int mole_y       = TIF::sol_head.solution_index("mole_y");

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
          // however, one TIF node may has several solution data when it lies on region interface, we should find corrent
          // solution data for this region which has the same material string as this region.
          std::pair<Solution_It, Solution_It> sol_it_pair = solution_map.equal_range(tif_node_index);
          Solution_It sol_it = sol_it_pair.first;
          for(; sol_it!=sol_it_pair.second; ++sol_it)
            if((*sol_it).second.material == region->material()) break;

          if( sol_it != sol_it_pair.second && sol_it != solution_map.end())
          {
            // field solution
            if(psi!=invalid_uint && elec_density!=invalid_uint && hole_density!=invalid_uint )
            {
              node_data->n()   = (*sol_it).second.data_array[elec_density] * pow(cm, -3);
              node_data->p()   = (*sol_it).second.data_array[hole_density] * pow(cm, -3);
              node_data->psi() = (*sol_it).second.data_array[psi         ] * V;
            }

            if(latt_temp != invalid_uint)
              node_data->T()   = (*sol_it).second.data_array[latt_temp   ] * K;
            else
              node_data->T()   = system.T_external();

            if(elec_temp != invalid_uint)
              node_data->Tn()  = (*sol_it).second.data_array[elec_temp   ] * K;
            else
              node_data->Tn()  = system.T_external();

            if(hole_temp != invalid_uint)
              node_data->Tp()  = (*sol_it).second.data_array[hole_temp   ] * K;
            else
              node_data->Tp()  = system.T_external();

            // doping
            if(net!=invalid_uint && total!=invalid_uint)
            {
              double net_doping   = (*sol_it).second.data_array[net  ] * pow(cm, -3);
              double total_doping = (*sol_it).second.data_array[total] * pow(cm, -3);

              node_data->Na()   = 0.5*std::abs(total_doping - net_doping);
              node_data->Nd()   = 0.5*std::abs(total_doping + net_doping);
            }

            // mole fraction
            if(sigle && mole_x!=invalid_uint)
            {
              node_data->mole_x()   = (*sol_it).second.data_array[mole_x];
            }

            if(complex && mole_x!=invalid_uint && mole_y!=invalid_uint)
            {
              node_data->mole_x()   = (*sol_it).second.data_array[mole_x];
              node_data->mole_y()   = (*sol_it).second.data_array[mole_y];
            }
          }
        }
#if 0
        // after import previous solutions, we re-init region here
        if(psi!=invalid_uint && elec_density!=invalid_uint && hole_density!=invalid_uint )
          region->reinit_after_import();
        else // or we have to set initial value to semiconductor node
          region->init(system.T_external());
#endif
        region->init(system.T_external());
        break;
      }
    case InsulatorRegion     :
      {
        SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
        SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
        for(; node_it!=node_it_end; ++node_it)
        {
          FVM_Node * fvm_node = (*node_it);

          FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);

          // tif_node_index is the index in TIF file that this FVM node has
          int tif_node_index = node_id_to_tif_index_map[fvm_node->root_node()->id()];
          // however, one TIF node may has several solution data when it lies on region interface, we should find corrent
          // solution data for this region which has the same material string as this region.
          std::pair<Solution_It, Solution_It> sol_it_pair = solution_map.equal_range(tif_node_index);
          Solution_It sol_it = sol_it_pair.first;
          for(; sol_it!=sol_it_pair.second; ++sol_it)
            if((*sol_it).second.material == region->material()) break;

          if( sol_it != sol_it_pair.second && sol_it != solution_map.end())
          {
            // field solution
            if(psi != invalid_uint)
              node_data->psi() = (*sol_it).second.data_array[psi         ] * V;

            if(latt_temp != invalid_uint)
              node_data->T()   = (*sol_it).second.data_array[latt_temp   ] * K;
            else
              node_data->T()   = system.T_external();
          }
        }
#if 0
        if(psi!=invalid_uint)
          region->reinit_after_import();
        else
          region->init(system.T_external());
#endif
        region->init(system.T_external());
        break;
      }
    case ElectrodeRegion     :
      {
        SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
        SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
        for(; node_it!=node_it_end; ++node_it)
        {
          FVM_Node * fvm_node = (*node_it);

          FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);

          // tif_node_index is the index in TIF file that this FVM node has
          int tif_node_index = node_id_to_tif_index_map[fvm_node->root_node()->id()];
          // however, one TIF node may has several solution data when it lies on region interface, we should find corrent
          // solution data for this region which has the same material string as this region.
          std::pair<Solution_It, Solution_It> sol_it_pair = solution_map.equal_range(tif_node_index);
          Solution_It sol_it = sol_it_pair.first;
          for(; sol_it!=sol_it_pair.second; ++sol_it)
            if((*sol_it).second.material == region->material()) break;

          if( sol_it != sol_it_pair.second && sol_it != solution_map.end())
          {
            // field solution
            if(psi != invalid_uint)
              node_data->psi() = (*sol_it).second.data_array[psi         ] * V;

            if(latt_temp != invalid_uint)
              node_data->T()   = (*sol_it).second.data_array[latt_temp   ] * K;
            else
              node_data->T()   = system.T_external();
          }
        }
#if 0
        if(psi!=invalid_uint)
          region->reinit_after_import();
        else
          region->init(system.T_external());
#endif
        region->init(system.T_external());
        break;
      }
    case MetalRegion     :
      {
        SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
        SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
        for(; node_it!=node_it_end; ++node_it)
        {
          FVM_Node * fvm_node = (*node_it);

          FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);

          // tif_node_index is the index in TIF file that this FVM node has
          int tif_node_index = node_id_to_tif_index_map[fvm_node->root_node()->id()];
          // however, one TIF node may has several solution data when it lies on region interface, we should find corrent
          // solution data for this region which has the same material string as this region.
          std::pair<Solution_It, Solution_It> sol_it_pair = solution_map.equal_range(tif_node_index);
          Solution_It sol_it = sol_it_pair.first;
          for(; sol_it!=sol_it_pair.second; ++sol_it)
            if((*sol_it).second.material == region->material()) break;

          if( sol_it != sol_it_pair.second && sol_it != solution_map.end())
          {
            // field solution
            if(psi != invalid_uint)
              node_data->psi() = (*sol_it).second.data_array[psi         ] * V;

            if(latt_temp != invalid_uint)
              node_data->T()   = (*sol_it).second.data_array[latt_temp   ] * K;
            else
              node_data->T()   = system.T_external();
          }
        }
#if 0
        if(psi!=invalid_uint)
          region->reinit_after_import();
        else
          region->init(system.T_external());
#endif
        region->init(system.T_external());
        break;
      }
    case VacuumRegion     :
      {
        region->init(system.T_external());
        break;
      }
    case PMLRegion        :
      {
        region->init(system.T_external());
        break;
      }
    default: genius_error();
    }
  }

  system.init_region_post_process();

  // we can free TIF data now
  TIF::clear();
}












