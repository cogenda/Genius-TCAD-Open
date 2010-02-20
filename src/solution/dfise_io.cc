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
#include "dfise_io.h"
#include "dfise.h"
#include "mesh.h"
#include "boundary_info.h"
#include "elem.h"
#include "parallel.h"
#include "mesh_communication.h"
#include "simulation_system.h"


#if   defined(HAVE_HASH_MAP)
# include <hash_map>
#elif defined(HAVE_EXT_HASH_MAP)
# include <ext/hash_map>
#endif

using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::V;
using PhysicalUnit::K;
using PhysicalUnit::eV;

/**
 * This method implements reading a mesh from a specified file
 * in DF-ISE format.
 */
void DFISEIO::read (const std::string& filename)
{

  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  Mesh & mesh = system.mesh();

  // clear the system
  system.clear();

  DFISE::DFISE_MESH ise_reader;

  // map node * to dfise node index
  std::map<Node *, unsigned int>       node_to_dfise_node_index_map;

  if( Genius::processor_id() == 0)
  {
    ise_reader.parse_dfise(filename);
    //ise_reader.export_vtk("debug.vtk");

    const DFISE::INFO & grid_info = ise_reader.get_grid_info();
    const DFISE::GRID & grid      = ise_reader.get_grid();

    /*
     * fill node location
     */
    for(unsigned int n=0; n<grid.Vertices.size(); ++n)
    {
      Node * node = mesh.add_point(grid.Vertices[n]*um);
      node_to_dfise_node_index_map[node] = n;
    }

    /*
     * set elements
     */

    // hold the boundary elems here
    std::vector<std::pair<Elem *, int> >boundary_elems;

    // search for all the dfise elements
    for(unsigned int n=0; n<grid.Elements.size(); ++n)
    {
      const DFISE::Element & element = grid.Elements[n];
      switch(element.elem_type)
      {
      case EDGE2    :
        {
          //edge2 should always be boundary elems
          genius_assert(grid_info.dimension == 2);

          Elem* elem = Elem::build(EDGE2).release();
          // fill elem node
          elem->set_node(0) = mesh.node_ptr( element.vertices[0] );
          elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
          boundary_elems.push_back(std::make_pair(elem , element.region_index));
          break;
        }
      case TRI3     :
        {
          Elem* elem = Elem::build(TRI3).release();
          // fill elem node
          elem->set_node(0) = mesh.node_ptr( element.vertices[0] );
          elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
          elem->set_node(2) = mesh.node_ptr( element.vertices[2] );

          // which region this elem belongs to
          if(grid_info.dimension == 2)
          {
            mesh.add_elem(elem);
            elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);
          }
          else
          {
            genius_assert(grid_info.dimension == 3);
            boundary_elems.push_back(std::make_pair(elem , element.region_index));
          }
          break;
        }
      case QUAD4    :
        {
          Elem* elem = Elem::build(QUAD4).release();

          // fill elem node
          elem->set_node(0) = mesh.node_ptr( element.vertices[0] );
          elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
          elem->set_node(2) = mesh.node_ptr( element.vertices[2] );
          elem->set_node(3) = mesh.node_ptr( element.vertices[3] );

          // which region this elem belongs to
          if(grid_info.dimension == 2)
          {
            mesh.add_elem(elem);
            elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);
          }
          else
          {
            genius_assert(grid_info.dimension == 3);
            boundary_elems.push_back(std::make_pair(elem , element.region_index));
          }
          break;
        }
      case TET4     :
        {
          genius_assert(grid_info.dimension == 3);
          Elem* elem = mesh.add_elem(Elem::build(TET4).release());

          // fill elem node
          elem->set_node(0) = mesh.node_ptr( element.vertices[2] );
          elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
          elem->set_node(2) = mesh.node_ptr( element.vertices[0] );
          elem->set_node(3) = mesh.node_ptr( element.vertices[3] );

          // which region this elem belongs to
          elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);
          break;
        }
      case PYRAMID5 :
        {
          genius_assert(grid_info.dimension == 3);
          Elem* elem = mesh.add_elem(Elem::build(PYRAMID5).release());

          // fill elem node
          elem->set_node(0) = mesh.node_ptr( element.vertices[3] );
          elem->set_node(1) = mesh.node_ptr( element.vertices[2] );
          elem->set_node(2) = mesh.node_ptr( element.vertices[1] );
          elem->set_node(3) = mesh.node_ptr( element.vertices[0] );
          elem->set_node(4) = mesh.node_ptr( element.vertices[4] );

          // which region this elem belongs to
          elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);

          break;
        }
      case PRISM6   :
        {
          genius_assert(grid_info.dimension == 3);
          Elem* elem = mesh.add_elem(Elem::build(PRISM6).release());

          // fill elem node
          elem->set_node(0) = mesh.node_ptr( element.vertices[2] );
          elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
          elem->set_node(2) = mesh.node_ptr( element.vertices[0] );
          elem->set_node(3) = mesh.node_ptr( element.vertices[3] );
          elem->set_node(4) = mesh.node_ptr( element.vertices[4] );
          elem->set_node(5) = mesh.node_ptr( element.vertices[5] );

          // which region this elem belongs to
          elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);
          genius_assert(elem->fvm_compatible_test());
          break;
        }
      case HEX8     :
        {
          genius_assert(grid_info.dimension == 3);
          Elem* elem = mesh.add_elem(Elem::build(HEX8).release());

          // fill elem node, the df-ise node order is different to libmesh order!
          // HEX8:  7        6                  HEX8: 6        7
          //        o--------o                        o--------o
          //       /:       /|                       /:       /|
          //      / :      / |                      / :      / |
          //   4 /  :   5 /  |                   3 /  :   5 /  |
          //    o--------o   |                    o--------o   |
          //    |   o....|...o 2                  |   o....|...o 4
          //    |  .3    |  /                     |  .2    |  /
          //    | .      | /                      | .      | /
          //    |.       |/                       |.       |/
          //    o--------o                        o--------o
          //    0        1                        0        1
          //      libmesh                           df-ise
          // and vertices order is 2 4 1 0, 3 5 7 6
          elem->set_node(0) = mesh.node_ptr( element.vertices[3] );
          elem->set_node(1) = mesh.node_ptr( element.vertices[2] );
          elem->set_node(2) = mesh.node_ptr( element.vertices[1] );
          elem->set_node(3) = mesh.node_ptr( element.vertices[0] );
          elem->set_node(4) = mesh.node_ptr( element.vertices[4] );
          elem->set_node(5) = mesh.node_ptr( element.vertices[5] );
          elem->set_node(6) = mesh.node_ptr( element.vertices[6] );
          elem->set_node(7) = mesh.node_ptr( element.vertices[7] );

          // which region this elem belongs to
          elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);
          genius_assert(elem->fvm_compatible_test());
          break;
        }
      default       : genius_error();
      }
    }

    //process boundary elems
    {
      typedef unsigned int                    key_type;
      typedef std::pair<Elem*, unsigned char> val_type;
      typedef std::pair<key_type, val_type>   key_val_pair;

#if   defined(HAVE_HASH_MAP)
      typedef std::hash_multimap<key_type, val_type> map_type;
#elif defined(HAVE_EXT_HASH_MAP)
# if    (__GNUC__ == 3) && (__GNUC_MINOR__ == 0) // gcc 3.0
      typedef std::hash_multimap<key_type, val_type> map_type;
# elif (__GNUC__ >= 3)                          // gcc 3.1 & newer
      typedef __gnu_cxx::hash_multimap<key_type, val_type> map_type;
# else
      // XLC and who knows what other compilers get here.
      // Try the most standard thing we can:
      typedef std::multimap<key_type, val_type>  map_type;
# endif
#else
      typedef std::multimap<key_type, val_type>  map_type;
#endif

      // A map from side keys to corresponding elements & side numbers
      map_type side_to_elem_map;

      const Mesh::element_iterator el_end = mesh.elements_end();
      for (Mesh::element_iterator el = mesh.elements_begin(); el != el_end; ++el)
      {
        Elem* elem = *el;
        for (unsigned int ms=0; ms<elem->n_neighbors(); ms++)
        {
          // Get the key for the side of this element
          const unsigned int key = elem->key(ms);

          key_val_pair kvp;
          kvp.first         = key;
          kvp.second.first  = elem;
          kvp.second.second = ms;

          side_to_elem_map.insert (kvp);
        }
      }

      // find the boundary_elems belongs to which elem/side pair
      for(unsigned int n=0; n<boundary_elems.size(); ++n)
      {
        Elem *boundary_elem   = boundary_elems[n].first;
        int   boundary_region = boundary_elems[n].second; //which dfise region this boundary_elem belongs to

        const unsigned int key = boundary_elem->key();

        // Look for elements that have an identical side key
        std::pair<map_type::iterator, map_type::iterator>  bounds = side_to_elem_map.equal_range(key);
        assert (bounds.first != bounds.second);

        // May be multiple keys, check all the possible elements which _might_ be neighbors.
        while (bounds.first != bounds.second)
        {
          // Get the potential element
          Elem* elem = bounds.first->second.first;

          // Get the side for the neighboring element
          const unsigned int ns = bounds.first->second.second;
          const AutoPtr<DofObject> elem_side(elem->side(ns));
          assert (elem_side.get() != NULL);

          if(*elem_side == *boundary_elem)
          {
            //skip interface elem, which will be processed later
            if(grid_info.region_to_boundaryregion(boundary_region)>=0)
              mesh.boundary_info->add_side(elem, ns, grid_info.region_to_boundaryregion(boundary_region));
          }
          ++bounds.first;
        }
      }

      //free the boundary_elems
      for(unsigned int n=0; n<boundary_elems.size(); ++n)
        delete boundary_elems[n].first;

      // map bc_index to bc label
      std::map<const std::string, int> bd_map;
      typedef std::map<const std::string, int>::iterator Bd_It;
      for(unsigned int n=0; n<grid_info.n_boundary_regions(); ++n)
      {
        bd_map[grid_info.boundary_label(n)]=n;
      }

      //NOTE: we havn't process all the interface and region external boundaries
      //build neighbor information for mesh. then elem->neighbor() is functional
      mesh.find_neighbors();

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
            std::string bd_label = grid_info.fieldregion_label(elem->subdomain_id()) + "_Neumann" ;
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
            std::string region1 = grid_info.fieldregion_label(sbd_id1);
            std::string region2 = grid_info.fieldregion_label(sbd_id2);
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



    // fill region label and region material
    mesh.set_n_subdomains() = grid_info.n_field_regions();
    for(unsigned int r=0; r<grid_info.nb_regions; ++r)
    {
      if(!grid_info.is_boundary_region(r) && !grid_info.is_interface_region(r))
      {
        mesh.set_subdomain_label( grid_info.region_to_fieldregion(r), grid_info.regions[r] );
        mesh.set_subdomain_material( grid_info.region_to_fieldregion(r), grid_info.materials[r]);
      }
    }

    // magic number, for 2D mesh, should < 2008
    if(grid_info.dimension == 2)
      mesh.magic_num() = 535;
    if(grid_info.dimension == 3)
      mesh.magic_num() = 40011;


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

  // broadcast dataset to all the processors
  {
    unsigned int n_datasets = ise_reader.n_datasets();
    Parallel::broadcast(n_datasets);

    //for each dataset
    for(unsigned int n=0; n<n_datasets; ++n)
    {
      DFISE::DATASET * dataset = NULL;
      if(Genius::processor_id() == 0)
        dataset = ise_reader.get_dataset(n);
      else
        dataset = new DFISE::DATASET;

      // broadcast critical data in this dataset
      Parallel::broadcast(dataset->name);
      Parallel::broadcast(dataset->Regions);
      Parallel::broadcast(dataset->Scalar_Values);
      Parallel::broadcast(dataset->node_to_value_index_map);

      if(Genius::processor_id() != 0)
        ise_reader.add_dataset(dataset);
    }

  }


  // set node id to dfise node index map, this should be same for all the processors
  std::map<unsigned int, unsigned int> node_id_to_dfise_node_index_map;
  // fill node_id_to_dfise_node_id_map by processor 0
  if( Genius::processor_id() == 0)
  {
    std::map<Node *, unsigned int>::iterator it = node_to_dfise_node_index_map.begin();
    for(; it != node_to_dfise_node_index_map.end(); ++it)
      node_id_to_dfise_node_index_map[(*it).first->id()] = (*it).second;
  }
  // broadcast node_id_to_dfise_node_index_map to all the processors
  Parallel::broadcast(node_id_to_dfise_node_index_map);


  // ok, we had got enough informations for set up each simulation region
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);

    switch ( region->type() )
    {
    case SemiconductorRegion :
      {
        SimulationRegion::node_iterator node_it = region->nodes_begin();
        SimulationRegion::node_iterator node_it_end = region->nodes_end();

        bool set_doping_by_concentration = ise_reader.is_value_exist("DopingConcentration", r);

        for(; node_it!=node_it_end; ++node_it)
        {
          FVM_Node * fvm_node = (*node_it).second;
          if( !fvm_node->root_node()->on_local() ) continue;

          FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);
          unsigned int dfise_node_index = node_id_to_dfise_node_index_map[fvm_node->root_node()->id()];

          if(set_doping_by_concentration)
          {
            double DopingConcentration = ise_reader.get_scaler_value("DopingConcentration", r, dfise_node_index)* pow(cm, -3);
            node_data->Na()   = DopingConcentration < 0 ? std::abs(DopingConcentration) : 0;
            node_data->Nd()   = DopingConcentration > 0 ? std::abs(DopingConcentration) : 0;
          }
          else
          {
            //Phosphorus PhosphorusActiveConcentration and etc.
            node_data->P()    = ise_reader.get_scaler_value_by_fuzzy("Phosphorus", r, dfise_node_index)* pow(cm, -3);;
            //Arsenic ArsenicActiveConcentration and etc.
            node_data->As()   = ise_reader.get_scaler_value_by_fuzzy("Arsenic", r, dfise_node_index)* pow(cm, -3);;
            //Antimony and etc
            node_data->Sb()   = ise_reader.get_scaler_value_by_fuzzy("Antimony", r, dfise_node_index)* pow(cm, -3);;
            //Boron BoronConcentration BoronActiveConcentration
            node_data->B()    = ise_reader.get_scaler_value_by_fuzzy("Boron", r, dfise_node_index)* pow(cm, -3);;
          }
        }

        // after import previous solutions, we re-init region here
        region->init(system.T_external());

        break;
      }
    case InsulatorRegion     :
      {
        SimulationRegion::node_iterator node_it = region->nodes_begin();
        SimulationRegion::node_iterator node_it_end = region->nodes_end();
        for(; node_it!=node_it_end; ++node_it)
        {
          FVM_Node * fvm_node = (*node_it).second;
          if( !fvm_node->root_node()->on_local() ) continue;

          FVM_NodeData * node_data = fvm_node->node_data();
          genius_assert(node_data);
        }

        region->init(system.T_external());

        break;
      }
    case ConductorRegion     :
      {
        SimulationRegion::node_iterator node_it = region->nodes_begin();
        SimulationRegion::node_iterator node_it_end = region->nodes_end();
        for(; node_it!=node_it_end; ++node_it)
        {
          FVM_Node * fvm_node = (*node_it).second;
          if( !fvm_node->root_node()->on_local() ) continue;

          FVM_NodeData * node_data = fvm_node->node_data();
          genius_assert(node_data);
        }
        region->init(system.T_external());

        break;
      }

    case VacuumRegion     :
      {
        region->init(system.T_external());
        break;
      }

    default: genius_error();
    }
  }


}












