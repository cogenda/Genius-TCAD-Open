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
#include<set>

// Local includes
#include "dfise_io.h"
#include "dfise.h"
#include "mesh_base.h"
#include "boundary_info.h"
#include "elem.h"
#include "parallel.h"
#include "mesh_communication.h"
#include "simulation_system.h"
#include "simulation_region.h"
#include "boundary_condition_collector.h"




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
  MeshBase & mesh = system.mesh();

  // clear the system
  system.clear();

  ise_reader = new DFISE::DFISE_MESH;

  if( Genius::processor_id() == 0)
  {
    ise_reader->parse_dfise(filename);

    const DFISE::INFO & grid_info = ise_reader->get_grid_info();
    const DFISE::GRID & grid      = ise_reader->get_grid();

    // fill node location
    VectorValue<double> translate(grid.translate);
    TensorValue<double> transform(grid.transform);
    for(unsigned int n=0; n<grid.Vertices.size(); ++n)
    {
      VectorValue<double> p(grid.Vertices[n][0], grid.Vertices[n][1], grid.Vertices[n][2]);
      p = transform*p + translate;
      Node * node = mesh.add_point(p*um);
      _node_to_dfise_node_index_map[node] = n;
    }


    //set elements
    for(unsigned int n=0; n<grid.Elements.size(); ++n)
    {
      _set_mesh_element(grid_info, grid.Elements[n]);
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

    //build neighbor information for mesh. then elem->neighbor() is functional
    mesh.find_neighbors();

    _set_boundary(grid_info);

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
  unsigned int n_datasets = ise_reader->n_datasets();
  Parallel::broadcast(n_datasets);

  //for each dataset
  for(unsigned int n=0; n<n_datasets; ++n)
  {
    DFISE::DATASET * dataset = NULL;
    if(Genius::processor_id() == 0)
      dataset = ise_reader->get_dataset(n);
    else
      dataset = new DFISE::DATASET;

    // broadcast critical data in this dataset
    Parallel::broadcast(dataset->name);
    Parallel::broadcast(dataset->Regions);
    Parallel::broadcast(dataset->Scalar_Values);
    Parallel::broadcast(dataset->node_to_value_index_map);

    if(Genius::processor_id() != 0)
      ise_reader->add_dataset(dataset);
  }


  // set node id to dfise node index map, this should be same for all the processors
  if( Genius::processor_id() == 0)
  {
    std::map<Node *, unsigned int>::iterator it = _node_to_dfise_node_index_map.begin();
    for(; it != _node_to_dfise_node_index_map.end(); ++it)
      _node_id_to_dfise_node_index_map[(*it).first->id()] = (*it).second;
  }

  // broadcast _node_id_to_dfise_node_index_map to all the processors
  Parallel::broadcast(_node_id_to_dfise_node_index_map);

  // setup each region
  _set_region();

  delete ise_reader;
}



void DFISEIO::_set_mesh_element( const DFISE::INFO & grid_info, const DFISE::Element &element)
{
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  unsigned int dim = grid_info.dimension;
  switch(element.elem_code)
  {
      case 1    :
      {
        //edge2 should always be boundary elems
        genius_assert(dim == 2);

        Elem* elem = Elem::build(EDGE2).release();
        // fill elem node
        elem->set_node(0) = mesh.node_ptr( element.vertices[0] );
        elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
        _boundary_elem_regions.push_back(std::make_pair(elem , element.region_index));
        break;
      }
      case 2     :
      {
        Elem* elem = Elem::build(TRI3).release();
        // fill elem node
        elem->set_node(0) = mesh.node_ptr( element.vertices[0] );
        elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
        elem->set_node(2) = mesh.node_ptr( element.vertices[2] );

        // which region this elem belongs to
        if(dim == 2)
        {
          mesh.add_elem(elem);
          elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);
        }
        else
        {
          genius_assert(dim == 3);
          _boundary_elem_regions.push_back(std::make_pair(elem , element.region_index));
        }
        break;
      }
      case 3    :
      {
        Elem* elem = Elem::build(QUAD4).release();

        // fill elem node
        elem->set_node(0) = mesh.node_ptr( element.vertices[0] );
        elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
        elem->set_node(2) = mesh.node_ptr( element.vertices[2] );
        elem->set_node(3) = mesh.node_ptr( element.vertices[3] );

        // which region this elem belongs to
        if(dim == 2)
        {
          mesh.add_elem(elem);
          elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);
          genius_assert(elem->fvm_compatible_test());
        }
        else
        {
          genius_assert(dim == 3);
          _boundary_elem_regions.push_back(std::make_pair(elem , element.region_index));
        }
        break;
      }
      case 5     :
      {
        genius_assert(dim == 3);
        Elem* elem = mesh.add_elem(Elem::build(TET4).release());

        // fill elem node
        elem->set_node(0) = mesh.node_ptr( element.vertices[0] );
        elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
        elem->set_node(2) = mesh.node_ptr( element.vertices[2] );
        elem->set_node(3) = mesh.node_ptr( element.vertices[3] );

        // which region this elem belongs to
        elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);
        break;
      }
      case 6 :
      {
        genius_assert(dim == 3);
        Elem* elem = mesh.add_elem(Elem::build(PYRAMID5).release());

        // fill elem node
        elem->set_node(0) = mesh.node_ptr( element.vertices[0] );
        elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
        elem->set_node(2) = mesh.node_ptr( element.vertices[2] );
        elem->set_node(3) = mesh.node_ptr( element.vertices[3] );
        elem->set_node(4) = mesh.node_ptr( element.vertices[4] );

        // which region this elem belongs to
        elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);
        genius_assert(elem->fvm_compatible_test());
        break;
      }
      case 7   :
      {
        genius_assert(dim == 3);
        Elem* elem = mesh.add_elem(Elem::build(PRISM6).release());
        // fill elem node
        elem->set_node(0) = mesh.node_ptr( element.vertices[0] );
        elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
        elem->set_node(2) = mesh.node_ptr( element.vertices[2] );
        elem->set_node(3) = mesh.node_ptr( element.vertices[3] );
        elem->set_node(4) = mesh.node_ptr( element.vertices[4] );
        elem->set_node(5) = mesh.node_ptr( element.vertices[5] );
        // which region this elem belongs to
        elem->subdomain_id() = grid_info.region_to_fieldregion(element.region_index);
        genius_assert(elem->fvm_compatible_test());
        break;
      }
      case 8     :
      {
        genius_assert(dim == 3);
        Elem* elem = mesh.add_elem(Elem::build(HEX8).release());

        elem->set_node(0) = mesh.node_ptr( element.vertices[0] );
        elem->set_node(1) = mesh.node_ptr( element.vertices[1] );
        elem->set_node(2) = mesh.node_ptr( element.vertices[2] );
        elem->set_node(3) = mesh.node_ptr( element.vertices[3] );
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



void DFISEIO::_set_boundary(const DFISE::INFO & grid_info)
{
  SimulationSystem & system = FieldInput<SimulationSystem>::system();
  MeshBase & mesh = system.mesh();

  //process boundary elems

  typedef unsigned int                    key_type;
  typedef std::pair<Elem*, unsigned char> val_type;
  typedef std::pair<key_type, val_type>   key_val_pair;

#if defined(HAVE_UNORDERED_MAP)
  typedef std::unordered_multimap<key_type, val_type> map_type;
#elif defined(HAVE_TR1_UNORDERED_MAP) || defined(HAVE_TR1_UNORDERED_MAP_WITH_STD_HEADER)
  typedef std::tr1::unordered_multimap<key_type, val_type> map_type;
#else
  typedef std::multimap<key_type, val_type>  map_type;
#endif

  // A map from side keys to corresponding elements & side numbers
  map_type side_to_elem_map;

  const MeshBase::element_iterator el_end = mesh.elements_end();
  for (MeshBase::element_iterator el = mesh.elements_begin(); el != el_end; ++el)
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
  for(unsigned int n=0; n<_boundary_elem_regions.size(); ++n)
  {
    Elem *boundary_elem   = _boundary_elem_regions[n].first;
    int   boundary_region = _boundary_elem_regions[n].second; //which dfise region this boundary_elem belongs to

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
        //genius_assert(elem->on_boundary(ns) || elem->on_interface() );
        //skip interface elem, which will be processed later
        if(grid_info.region_to_boundaryregion(boundary_region)>=0)
          mesh.boundary_info->add_side(elem, ns, grid_info.region_to_boundaryregion(boundary_region));
      }

      ++bounds.first;
    }
  }

  //free the boundary_elems
  for(unsigned int n=0; n<_boundary_elem_regions.size(); ++n)
    delete _boundary_elem_regions[n].first;

  // map bc_index to bc label
  std::map<std::string, int> bd_map;
  std::map<int, std::string> inverse_bd_map;
  typedef std::map<std::string, int>::iterator Bd_It;
  for(unsigned int n=0; n<grid_info.n_boundary_regions(); ++n)
  {
    bd_map[grid_info.boundary_label(n)]=n;
    inverse_bd_map[n] = grid_info.boundary_label(n);
  }


  short int new_bd_id_begin = bd_map.size();

  // a dfise boundary may cover several regions, which breaks genius boundary rule, fix it here
  {
    std::vector<unsigned int>       boundary_el;
    std::vector<unsigned short int> boundary_sl;
    std::vector<short int>          boundary_il;
    mesh.boundary_info->build_side_list (boundary_el, boundary_sl, boundary_il);

    // <bd-index, <region_index, region_index>, elem-side> >
    typedef std::map<short int, std::map<std::pair<unsigned int, unsigned int>, std::vector<std::pair<unsigned int, unsigned int> > > > BoundaryRegions;
    BoundaryRegions boundary_regions;
    for(unsigned int n=0; n<boundary_el.size(); n++)
    {
      const Elem * elem = mesh.elem(boundary_el[n]);
      unsigned int region_1 = elem->subdomain_id();
      unsigned int region_2 = invalid_uint;
      if(elem->neighbor(boundary_sl[n]))
      {
        region_2 = elem->neighbor(boundary_sl[n])->subdomain_id();
        if(region_1 > region_2) std::swap(region_1, region_2);
      }
      std::pair<unsigned int, unsigned int> region_pair(region_1, region_2);
      boundary_regions[boundary_il[n]][region_pair].push_back(std::make_pair(boundary_el[n], boundary_sl[n]));
    }

    BoundaryRegions::const_iterator it = boundary_regions.begin();
    for( ; it != boundary_regions.end(); ++it)
    {
      // this boundary only conver one region (pair), that is right
      if(it->second.size() == 1 ) continue;

      // we should split this boundary by regions
      std::string boundary_label = inverse_bd_map[it->first];
      bd_map.erase(boundary_label);

      const std::map<std::pair<unsigned int, unsigned int>, std::vector<std::pair<unsigned int, unsigned int> > > & region_boundary_elems = it->second;
      std::map<std::pair<unsigned int, unsigned int>, std::vector<std::pair<unsigned int, unsigned int> > >::const_iterator region_boundary_elems_it = region_boundary_elems.begin();
      for(; region_boundary_elems_it != region_boundary_elems.end(); ++region_boundary_elems_it)
      {
        unsigned int region_1 = region_boundary_elems_it->first.first;
        unsigned int region_2 = region_boundary_elems_it->first.second;

        std::string new_boundary_label;
        if(region_2 == invalid_uint)
          new_boundary_label = boundary_label + '_' + mesh.subdomain_label_by_id(region_1);
        else
        {
          if(mesh.subdomain_label_by_id(region_1) < mesh.subdomain_label_by_id(region_2))
            new_boundary_label = boundary_label + '_' + mesh.subdomain_label_by_id(region_1) + "_to_" + mesh.subdomain_label_by_id(region_2);
          else
            new_boundary_label = boundary_label + '_' + mesh.subdomain_label_by_id(region_2) + "_to_" + mesh.subdomain_label_by_id(region_1);
        }

        short int new_bd_id = new_bd_id_begin++;
        bd_map[new_boundary_label] = new_bd_id;

        const std::vector<std::pair<unsigned int, unsigned int> > & boundary_elems = region_boundary_elems_it->second;
        for(unsigned int n=0; n<boundary_elems.size(); ++n)
        {
          const Elem * elem = mesh.elem(boundary_elems[n].first);
          unsigned int side = boundary_elems[n].second;
          mesh.boundary_info->remove( elem, side );
          mesh.boundary_info->add_side(elem, side, new_bd_id);
        }
      }
    }

  }



  //NOTE: we havn't process all the interface and region external boundaries, do it here

  std::set<short int> dummy_bds;
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
          bd_id = new_bd_id_begin++;
          bd_map[bd_label] = bd_id;
          dummy_bds.insert(bd_id);
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
          bd_id = new_bd_id_begin++;
          bd_map[bd_label] = bd_id;
          dummy_bds.insert(bd_id);
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
    short int bd_id = bd_it->second;
    bool user_defined = (dummy_bds.find(bd_id) == dummy_bds.end() );
    mesh.boundary_info->set_label_to_id( bd_it->second, bd_it->first, user_defined );
  }

}


void DFISEIO::_set_region()
{
  SimulationSystem & system = FieldInput<SimulationSystem>::system();

  const double doping_scale = pow(cm, -3);

  // ok, we had got enough informations for set up each simulation region
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = system.region(r);

    switch ( region->type() )
    {
        case SemiconductorRegion :
        {
          bool has_net_doping_concentration = ise_reader->is_value_exist("DopingConcentration", r);
          bool has_total_doping_concentration = ise_reader->is_value_exist("TotalConcentration", r);

          bool species_doping_concentration = ( ise_reader->is_value_exist_fuzzy("Boron", r) ||
                                                ise_reader->is_value_exist_fuzzy("Phosphorus", r) ||
                                                ise_reader->is_value_exist_fuzzy("Arsenic", r) ||
                                                ise_reader->is_value_exist_fuzzy("Antimony", r) );
          unsigned int Boron_index=invalid_uint, Phosphorus_index=invalid_uint, Arsenic_index=invalid_uint, Antimony_index=invalid_uint;
          // add variables to this region
          if( species_doping_concentration )
          {
            Boron_index      = region->add_variable( SimulationVariable("Boron", SCALAR, POINT_CENTER, "cm^-3", invalid_uint, true, true)  );
            Phosphorus_index = region->add_variable( SimulationVariable("Phosphorus", SCALAR, POINT_CENTER, "cm^-3", invalid_uint, true, true)  );
            Arsenic_index    = region->add_variable( SimulationVariable("Arsenic", SCALAR, POINT_CENTER, "cm^-3", invalid_uint, true, true)  );
            Antimony_index   = region->add_variable( SimulationVariable("Antimony", SCALAR, POINT_CENTER, "cm^-3", invalid_uint, true, true)  );
          }

          SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
          SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
          for(; node_it!=node_it_end; ++node_it)
          {
            FVM_Node * fvm_node = (*node_it);

            FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);
            unsigned int dfise_node_index = _node_id_to_dfise_node_index_map[fvm_node->root_node()->id()];

            double NetConcentration = ise_reader->get_scaler_value("DopingConcentration", r, dfise_node_index)* doping_scale;
            double TotalConcentration = ise_reader->get_scaler_value("TotalConcentration", r, dfise_node_index)* doping_scale;

            double Phosphorus=0.0;
            if( ise_reader->is_value_exist("PhosphorusActiveConcentration", r) )
              Phosphorus = ise_reader->get_scaler_value("PhosphorusActiveConcentration", r, dfise_node_index)* doping_scale;
            else if( ise_reader->is_value_exist("PhosphorusConcentration", r) )
              Phosphorus = ise_reader->get_scaler_value("PhosphorusConcentration", r, dfise_node_index)* doping_scale;

            double Arsenic=0.0;
            if( ise_reader->is_value_exist("ArsenicActiveConcentration", r) )
              Arsenic = ise_reader->get_scaler_value("ArsenicActiveConcentration", r, dfise_node_index)* doping_scale;
            else if( ise_reader->is_value_exist("ArsenicConcentration", r) )
              Arsenic = ise_reader->get_scaler_value("ArsenicConcentration", r, dfise_node_index)* doping_scale;

            double Antimony=0.0;
            if( ise_reader->is_value_exist("AntimonyActiveConcentration", r) )
              Antimony = ise_reader->get_scaler_value("AntimonyActiveConcentration", r, dfise_node_index)* doping_scale;
            else if( ise_reader->is_value_exist("AntimonyConcentration", r) )
              Antimony = ise_reader->get_scaler_value("AntimonyConcentration", r, dfise_node_index)* doping_scale;

            double Boron=0.0;
            if( ise_reader->is_value_exist("BoronActiveConcentration", r) )
              Boron = ise_reader->get_scaler_value("BoronActiveConcentration", r, dfise_node_index)* doping_scale;
            else if( ise_reader->is_value_exist("BoronConcentration", r) )
              Boron = ise_reader->get_scaler_value("BoronConcentration", r, dfise_node_index)* doping_scale;


            // Have species doping concentration and don't have DopingConcentration
            if(species_doping_concentration && !has_net_doping_concentration)
            {
              node_data->data<Real>(Phosphorus_index) = Phosphorus;
              node_data->data<Real>(Arsenic_index)    = Arsenic;
              node_data->data<Real>(Antimony_index)   = Antimony;
              node_data->data<Real>(Boron_index)      = Boron;
              node_data->Na()   = Boron;
              node_data->Nd()   = Phosphorus + Arsenic + Antimony;
              continue;
            }

            // both have DopingConcentration and TotalConcentration
            if( has_total_doping_concentration )
            {
              node_data->Na()   = 0.5*(TotalConcentration - NetConcentration);
              node_data->Nd()   = 0.5*(TotalConcentration + NetConcentration);
              continue;
            }

            // check if species_doping_concentration match net_doping_concentration

            if(species_doping_concentration && std::abs((Phosphorus + Arsenic + Antimony) - Boron - NetConcentration) < 1e-2*std::abs(NetConcentration) )
            {
              node_data->data<Real>(Phosphorus_index) = Phosphorus;
              node_data->data<Real>(Arsenic_index)    = Arsenic;
              node_data->data<Real>(Antimony_index)   = Antimony;
              node_data->data<Real>(Boron_index)      = Boron;
              node_data->Na()   = Boron;
              node_data->Nd()   = Phosphorus + Arsenic + Antimony;
              continue;
            }

            // at last, we have to set doping by only NetConcentration
            node_data->Na()   = NetConcentration < 0 ? -NetConcentration : 0;
            node_data->Nd()   = NetConcentration > 0 ?  NetConcentration : 0;
          }

          // after import previous solutions, we re-init region here
          region->init(system.T_external());

          break;
        }
        case InsulatorRegion     :
        case ElectrodeRegion     :
        case MetalRegion         :
        {
          SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
          SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
          for(; node_it!=node_it_end; ++node_it)
          {
            FVM_Node * fvm_node = (*node_it);
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

  system.init_region_post_process();
}






//------------------------------------------------------------------------------

//------------------------------------------------------------------------------







void DFISEIO::write (const std::string& filename)
{
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();
  const MeshBase & mesh = system.mesh();

  // prepare data, we have to create DopingConcentration and NetConcentration
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    SimulationRegion * region = const_cast<SimulationRegion *>(system.region(r));
    if(region->type() == SemiconductorRegion)
    {
      unsigned int  doping_index   = region->add_variable( SimulationVariable("DopingConcentration", SCALAR, POINT_CENTER, "cm^-3", invalid_uint, true, true)  );
      unsigned int  net_index      = region->add_variable( SimulationVariable("NetConcentration", SCALAR, POINT_CENTER, "cm^-3", invalid_uint, true, true)  );

      SimulationRegion::local_node_iterator node_it = region->on_local_nodes_begin();
      SimulationRegion::local_node_iterator node_it_end = region->on_local_nodes_end();
      for(; node_it!=node_it_end; ++node_it)
      {
        FVM_Node * fvm_node = (*node_it);

        FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data);
        node_data->data<Real>(doping_index) = node_data->Total_doping();
        node_data->data<Real>(net_index) = node_data->Net_doping();
      }
    }
  }

  std::vector<DFISE::DATASET *> datasets;
  datasets.push_back(_build_dataset("DopingConcentration", "DopingConcentration"));
  datasets.push_back(_build_dataset("NetConcentration", "NetConcentration"));

  if( Genius::processor_id() == 0)
  {
    ise_writer = new DFISE::DFISE_MESH;

    _build_region();

    _build_element();


    // fill grid info
    DFISE::INFO & grid_info = ise_writer->get_grid_info();
    {
      grid_info.version = 1.0;
      grid_info.type = DFISE::INFO::grid;
      grid_info.dimension = mesh.mesh_dimension();
      grid_info.nb_vertices = mesh.n_nodes();
      grid_info.nb_edges = _edges.size();
      grid_info.nb_faces = _faces.size();
      grid_info.nb_elements = _elems.size();
      grid_info.nb_regions = _regions.size();
      grid_info.regions = _regions;
      grid_info.materials = _materials;
    }


    // fill grid data
    DFISE::GRID & grid = ise_writer->get_grid();
    {
      grid.dimension = mesh.mesh_dimension();

      // vertex
      MeshBase::const_node_iterator nd = mesh.active_nodes_begin();
      MeshBase::const_node_iterator nd_end = mesh.active_nodes_end();
      for ( ; nd!=nd_end; ++nd )
      {
        const Node * node = (*nd);
        DFISE::Point p;
        {
          p[0] =  (*node)(0)/um; //scale to um
          p[1] =  (*node)(1)/um; //scale to um
          p[2] =  (*node)(2)/um; //scale to um
        }
        grid.Vertices.push_back(p);
      }

      // edge
      for(unsigned int n=0; n<_edges.size(); ++n)
      {
        grid.add_edge( _edges[n] );
      }

      // face

      for(unsigned int n=0; n<_faces.size(); ++n)
      {
        const Elem * face = _faces[n];
        grid.Faces.push_back(_edge_index(face));
      }

      // location
      if(mesh.mesh_dimension() == 2)
      {
        for(unsigned int n=0; n<_edges.size(); ++n)
        {
          grid.Locations.push_back(_location[n]);
        }
      }
      if(mesh.mesh_dimension() == 3)
      {
        for(unsigned int n=0; n<_faces.size(); ++n)
        {
          grid.Locations.push_back(_location[n]);
        }
      }

      // element
      for(unsigned int n=0; n<_elems.size(); ++n)
      {
        const Elem * elem = _elems[n];
        DFISE::Element DFISE_Elem;
        DFISE_Elem.elem_code = _elem_code(elem->type());

        switch( Elem::dim(elem->type()) )
        {
            case 1 :
            {
              DFISE_Elem.faces.push_back(elem->get_node(0)->id());
              DFISE_Elem.faces.push_back(elem->get_node(1)->id());
              break;
            }
            case 2 :  DFISE_Elem.faces = _edge_index(elem); break;
            case 3 :  DFISE_Elem.faces = _face_index(elem); break;
        }
        grid.Elements.push_back(DFISE_Elem);
      }

      // regions
      grid.regions = _regions;
      grid.materials = _materials;
      grid.region_elements = _region_elements;
    }


    // fill data info
    DFISE::INFO & data_info = ise_writer->get_dataset_info();
    {
      data_info.version = 1.0;
      data_info.type = DFISE::INFO::dataset;
      data_info.dimension = mesh.mesh_dimension();
      data_info.nb_vertices = mesh.n_nodes();
      data_info.nb_edges = _edges.size();
      data_info.nb_faces = _faces.size();
      data_info.nb_elements = _elems.size();
      data_info.nb_regions = _regions.size();
    }

    // set dataset, dataset will be freed by ise_writer
    for(unsigned int n=0; n<datasets.size(); ++n)
    {
      data_info.datasets.push_back( datasets[n]->name );
      data_info.functions.push_back( datasets[n]->function );
      ise_writer->add_dataset(datasets[n], true);
    }

    ise_writer->write_dfise(filename);
    delete ise_writer;
  }

  for(unsigned int n=0; n<datasets.size(); ++n)
    delete datasets[n];

  _write_on_exit();

}





void DFISEIO::_build_region()
{
  genius_assert( Genius::is_first_processor() );

  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();
  const MeshBase & mesh = system.mesh();

  std::map<short int, std::vector< std::pair<const Elem *, unsigned int> > > boundary_elem_side_map;
  mesh.boundary_info->active_elem_with_boundary_id(boundary_elem_side_map);

  // region elements
  MeshBase::const_element_iterator el_end = mesh.elements_end();
  for (MeshBase::const_element_iterator el = mesh.elements_begin(); el != el_end; ++el)
  {
    const Elem* elem = *el;
    _elems.push_back(elem);
  }

  // region name/materials
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);
    _regions.push_back(region->name());
    _materials.push_back(region->material());
  }

  // boundary cells, labels and material "Contact"
  for(unsigned int b=0; b<system.get_bcs()->n_bcs(); b++)
  {
    const BoundaryCondition * bc = system.get_bcs()->get_bc(b);
    if(bc->boundary_type()==BOUNDARY && bc->is_electrode())
    {
      _regions.push_back(bc->label());
      _materials.push_back("Contact");

      const std::vector< std::pair<const Elem *, unsigned int> > & boundary_elems
          = boundary_elem_side_map.find(bc->boundary_id())->second;

      for( unsigned int n=0; n< boundary_elems.size(); n++)
      {
        const std::pair<const Elem *, unsigned int> & elem_pair = boundary_elems[n];
        Elem * bc_elem = elem_pair.first->build_side(elem_pair.second, false).release();
        bc_elem->subdomain_id() = _regions.size()-1;
        _elems.push_back(bc_elem);
        _boundary_elems.push_back(bc_elem);
      }
    }
  }

  // statistic region elements
  _region_elements.resize(_regions.size());
  for(unsigned int n=0; n<_elems.size(); ++n)
  {
    const Elem * elem = _elems[n];
    _region_elements[elem->subdomain_id()].push_back(n);
  }

}



void DFISEIO::_build_element()
{
  genius_assert( Genius::is_first_processor() );

  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();
  const MeshBase & mesh = system.mesh();

  // build element
  MeshBase::const_element_iterator el_end = mesh.elements_end();
  for (MeshBase::const_element_iterator el = mesh.elements_begin(); el != el_end; ++el)
  {
    const Elem* elem = *el;

    for (unsigned int ms=0; ms<elem->n_sides(); ms++)
    {
      char location = 'i';
      if( elem->on_boundary(ms) ) location='e';
      if( elem->on_interface(ms) ) location='f';

      // faces
      if( elem->dim() == 3)
      {
        AutoPtr<Elem> side_elem = elem->build_side(ms, false);
        bool face_not_exist = _add_face(side_elem.get()); // return true when face not exist
        int face_index = side_elem->id(); // id is the face index
        if(face_not_exist)
        {
          // ok, a new face is added to _faces, set the location of it
          _location.push_back(location);
          // also, we know the face index is positive
          _elem_faces[elem].push_back(face_index);
          // release the autoprt, so the face will not be delete here
          side_elem.release();
        }
        else
        {
          // there is already face exist, which comes from my neighbor elem,
          // as a result, the existing face order should be inversed
          _elem_faces[elem].push_back(-face_index-1);
        }
      }

      // edges
      if( elem->dim() == 2)
      {
        std::pair<unsigned int, unsigned int> nodes;
        elem->nodes_on_edge(ms, nodes);
        if(_add_edge(elem->get_node(nodes.first)->id(), elem->get_node(nodes.second)->id()))
        {
          // ok, a new edge is added to _edges, set the location of it
          _location.push_back(location);
        }
      }
    }

    // extra edge for coords
    if( elem->dim() == 3)
    {
      for (unsigned int ns=0; ns<elem->n_edges(); ns++)
      {
        std::pair<unsigned int, unsigned int> nodes;
        elem->nodes_on_edge(ns, nodes);
        _add_edge(elem->get_node(nodes.first)->id(), elem->get_node(nodes.second)->id());
      }
    }
  }
}



bool DFISEIO::_add_face( Elem *face )
{
  const unsigned int key = face->key();

  typedef _multimap_type::iterator key_elem_it;
  std::pair<key_elem_it, key_elem_it>  bounds = _faces_map.equal_range(key);
  while (bounds.first != bounds.second)
  {
    const Elem * elem = bounds.first->second;
    if( *elem == *face ) // this face already exist
    {
      face->set_id() = elem->id();
      return false;
    }
    ++bounds.first;
  }

  // face not exist, add it to faces
  if( bounds.first == bounds.second )
  {
    face->set_id() = _faces.size();
    _faces.push_back(face);
    _faces_map.insert( std::make_pair(key, face) );
    return true;
  }

  return false;
}


bool DFISEIO::_add_edge( int n1, int n2 )
{
  std::pair<int, int> edge(n1, n2);
  if( _edges_map.find(edge) == _edges_map.end() ) // edge not exist
  {
    _edges_map.insert( std::make_pair(edge, static_cast<int>(_edges_map.size()) ) );
    _edges.push_back(edge);
    return true;
  }
  return false;
}



std::vector<int> DFISEIO::_edge_index(const Elem *face) const
{
  std::vector<int> edge_index;
  for (unsigned int ms=0; ms<face->n_sides(); ms++)
  {
    std::pair<unsigned int, unsigned int> nodes;
    face->nodes_on_edge(ms, nodes);

    std::pair<int, int> face_edge(face->get_node(nodes.first)->id(), face->get_node(nodes.second)->id());

    std::map< std::pair<int, int>, int, _lt_pair_int>::const_iterator it  = _edges_map.find(face_edge);
    assert( it != _edges_map.end() );

    int index = it->second;
    const std::pair<int, int> & edge = it->first;

    if( edge.first == face_edge.first )
    {
      assert(edge.second == face_edge.second);
      edge_index.push_back(index);
    }
    if( edge.first == face_edge.second )
    {
      assert(edge.second == face_edge.first);
      edge_index.push_back(-index-1);
    }
  }

  return edge_index;
}



std::vector<int> DFISEIO::_face_index(const Elem *cell) const
{
  // we had precomputed face_index to mesh cells
  if( _elem_faces.find(cell) != _elem_faces.end() )
  {
    const std::vector<int> & index = _elem_faces.find(cell)->second;
    std::vector<unsigned int> order;
    cell->side_order( ISE,  order);
    std::vector<int> face_index(index.size());
    for(unsigned int n=0; n<order.size(); ++n)
      face_index[n] = index[order[n]];
    return face_index;
  }

  // however, for boundary cells, we have to work it out
  typedef _multimap_type::const_iterator key_elem_it;
  std::vector<int> face_index;

  std::vector<unsigned int> order;
  cell->side_order( ISE,  order);
  for (unsigned int ms=0; ms<order.size(); ms++)
  {
    const Point norm = cell->outside_unit_normal(order[ms]);

    AutoPtr<Elem> face_elem = cell->build_side(order[ms]);
    const unsigned int key = face_elem->key();

    std::pair<key_elem_it, key_elem_it>  bounds = _faces_map.equal_range(key);
    while (bounds.first != bounds.second)
    {
      const Elem * elem = bounds.first->second;
      if( *elem == *face_elem ) break;
      ++bounds.first;
    }

    assert (bounds.first != bounds.second);
    const Elem * face = bounds.first->second;
    int index = face->id();

    Point face_norm;
    {
      Point p1 = face->point(0);
      Point p2 = face->point(1);
      Point p3 = face->point(2);
      face_norm = ((p2 - p1).cross(p3 - p2)).unit();
    }

    if( norm*face_norm > 0 )
      face_index.push_back(index);
    else
      face_index.push_back(-index-1);
  }
  return face_index;
}


int DFISEIO::_elem_code(ElemType t) const
{
  switch(t)
  {
      case EDGE2:
      case EDGE2_FVM:    return 1;
      case TRI3:
      case TRI3_FVM:     return 2;
      case QUAD4:
      case QUAD4_FVM:    return 3;
      case TET4:
      case TET4_FVM:     return 5;
      case PYRAMID5:
      case PYRAMID5_FVM: return 6;
      case PRISM6:
      case PRISM6_FVM:   return 7;
      case HEX8:
      case HEX8_FVM:     return 8;
      default: return 0;
  }
  return 0;
}



DFISE::DATASET * DFISEIO::_build_dataset( const std::string & variable,  const std::string & function)
{
  std::vector<std::string> validity;
  std::map<unsigned int, double> scalar;
  const SimulationSystem & system = FieldOutput<SimulationSystem>::system();
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);
    if(region->has_variable(variable, POINT_CENTER))
    {
      validity.push_back(region->name());

      std::vector<Real> data;
      std::vector<unsigned int> index;
      region->get_variable_data<Real>(variable, POINT_CENTER, data);
      region->region_node(index);

      for(unsigned int n=0; n<index.size(); ++n)
        scalar.insert( std::make_pair(index[n], data[n]) );
    }
  }

  if( Genius::processor_id() == 0)
  {
    DFISE::DATASET * dataset = new DFISE::DATASET;
    dataset->name = variable;
    dataset->function = function;
    dataset->type = DFISE::DATASET::scalar;
    dataset->dimension = 1;
    dataset->location = DFISE::DATASET::vertex;
    dataset->validity = validity;
    dataset->n_data = scalar.size();
    std::map<unsigned int, double>::const_iterator it = scalar.begin();
    for(; it != scalar.end(); ++it)
      dataset->Scalar_Values.push_back(it->second);
    return dataset;
  }

  return 0;
}


void DFISEIO::_write_on_exit()
{
  for(unsigned int n=0; n<_faces.size(); ++n)
  {
    delete _faces[n];
  }

  for(unsigned int n=0; n<_boundary_elems.size(); ++n)
    delete _boundary_elems[n];
}

