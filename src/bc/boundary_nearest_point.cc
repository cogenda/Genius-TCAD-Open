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
#include "mesh_base.h"
#include "mesh_tools.h"
#include "boundary_info.h"
#include "simulation_system.h"
#include "simulation_region.h"
#include "boundary_condition_is.h"
#include "surface_locator_hub.h"
#include "boundary_condition_collector.h"


using PhysicalUnit::nm;

// find the nearset point on gate/charge boundary
void InsulatorSemiconductorInterfaceBC::_find_nearest_points_in_gate_region()
{

  SimulationSystem  & system  = this->system();
  const MeshBase & mesh = system.mesh();

  const std::vector<const Node *> & nodes = this->nodes();

  // build the bounding box for this boundary
  std::pair<Point, Point> boundary_bounding_box;
  {
    Point min(1.e30,   1.e30,  1.e30);
    Point max(-1.e30, -1.e30, -1.e30);

    for (unsigned int n=0; n<nodes.size(); n++)
      for (unsigned int i=0; i<3; i++)
      {
        min(i) = std::min(min(i), (*nodes[n])(i));
        max(i) = std::max(max(i), (*nodes[n])(i));
      }
    boundary_bounding_box = std::make_pair(min, max);
  }

  // the nearst nodes <Node *, region>, some of them may not on processor, we must set them as on local later
  std::multimap<Node *, unsigned int> target_nodes;

  // build fast surface locator
  SurfaceLocatorHub & surface_locator = mesh.surface_locator();


  for(unsigned int r=0; r<system.n_regions(); ++r)
  {
    const SimulationRegion * region = system.region(r);
    // only consider nearset point on conductor or resistance region
    if( region->type() != MetalRegion &&  region->type() != ElectrodeRegion ) continue;

    const std::pair<Point, Point> & region_bounding_box = region->boundingbox();

    Real distance = MeshTools::minimal_distance(region_bounding_box, boundary_bounding_box);
    if( distance < 1e-6*nm || distance > 30*nm ) continue; // skip adjacent and far away regions


    for(unsigned int n=0; n<nodes.size(); ++n)
    {
      // skip node not on processor
      if( !nodes[n]->on_processor() ) continue;

      const Point p = *nodes[n];
      // point is far from boundingbox of region, skip it
      if( MeshTools::minimal_distance(region_bounding_box, p) >  30*nm ) continue;

      Point project_point;
      std::pair<const Elem*, unsigned int> surface_elem_pair = surface_locator(p, r, project_point, 30*nm);
      if( surface_elem_pair.first == NULL ) continue;

      // ok, which bc the nearset point on?
      short int bd_id = mesh.boundary_info->boundary_id(surface_elem_pair.first, surface_elem_pair.second);
      const BoundaryCondition * bc = system.get_bcs()->get_bc_by_bd_id(bd_id);

      if( bc->bc_type() == ChargedContact ||
          bc->bc_type() == GateContact    ||
          bc->bc_type() == IF_Insulator_Metal )
      {
        // save it
        NearestPoint loc;
        loc.region = const_cast<SimulationRegion *>(region);
        loc.bc = const_cast<BoundaryCondition *>(bc);
        loc.elem = surface_elem_pair.first;
        loc.side = surface_elem_pair.second;
        loc.p = project_point;
        _node_nearest_point_map.insert(std::make_pair(nodes[n], loc));

        // save the target nodes
        AutoPtr<Elem> side = loc.elem->build_side(loc.side);
        for(unsigned int i=0; i<side->n_nodes(); ++i)
        {
          Node * n = side->get_node(i);
          target_nodes.insert(std::make_pair(n, loc.elem->subdomain_id()));
        }
      }

    }
  }


  // set all target node not on this processor as local, create FVM_Node .
  std::set<SimulationRegion *> modified_regions;
  std::multimap<Node *, unsigned int>::iterator node_it = target_nodes.begin();
  for( ; node_it != target_nodes.end(); ++node_it)
  {
    Node * node = node_it->first;
    if( node->on_local() ) continue;

    node->on_local() = true;

    // update fvm_node in the region
    SimulationRegion * region = system.region(node_it->second);
    region->insert_fvm_node(new FVM_Node( node ));

    // mark regions need to be updated
    modified_regions.insert(region);
  }

  // rebuild region fvm node array
  std::set<SimulationRegion *>::iterator region_it = modified_regions.begin();
  for( ; region_it != modified_regions.end(); ++region_it)
  {
    (*region_it)->rebuild_region_fvm_node_list();
  }

}
