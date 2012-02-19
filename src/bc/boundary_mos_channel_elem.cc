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
#include <queue>
#include <algorithm>

#include "mesh_base.h"
#include "mesh_tools.h"
#include "boundary_info.h"
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "boundary_condition_gate.h"
#include "boundary_condition_ir.h"
#include "boundary_condition_simplegate.h"
#include "boundary_condition_is.h"
#include "surface_locator_hub.h"
#include "boundary_condition_collector.h"
#include "parallel.h"

using PhysicalUnit::nm;


void ResistanceInsulatorBC::_find_mos_channel_elem()
{
#if 0
  SimulationSystem  & system  = this->system();
  const MeshBase & mesh = system.mesh();

  // build fast surface locator
  SurfaceLocatorHub & surface_locator = mesh.surface_locator();

  // gate region
  const SimulationRegion * gate_region = bc_regions().second;
  assert( gate_region->type() == MetalRegion );
  const std::pair<Point, Point> & gate_region_bounding_box = gate_region->boundingbox();

  for(unsigned int r=0; r<system.n_regions(); ++r)
  {
    SimulationRegion * region = system.region(r);
    // only consider nearset point on semiconductor region
    if( region->type() != SemiconductorRegion ) continue;
    if( region->is_neighbor(gate_region) ) continue; // not gate, but ohmic/schottky contact

    const std::pair<Point, Point> & region_bounding_box = region->boundingbox();
    Real distance = MeshTools::minimal_distance(region_bounding_box, gate_region_bounding_box);
    if( distance > 100*nm ) continue; // skip far away regions

    SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<SemiconductorSimulationRegion *>(region);
    SimulationRegion::element_iterator elem_it = region->elements_begin();
    SimulationRegion::element_iterator elem_it_end = region->elements_end();
    for(; elem_it != elem_it_end; elem_it++)
    {
      const Elem * e = *elem_it;
      if( !e->on_processor() ) continue;

      const Point centroid = e->centroid();
      // point is far from boundingbox of gate region, skip it
      if( std::abs(MeshTools::minimal_distance(gate_region_bounding_box, centroid)) >  100*nm ) continue;

      Point project_point;
      std::pair<const Elem*, unsigned int> surface_elem_pair = surface_locator(centroid, gate_region->subdomain_id(), project_point, 100*nm);
      if( surface_elem_pair.first != NULL )
        semiconductor_region->_elem_in_mos_channel.insert(e);

    }

  }
#endif
}



void GateContactBC::_find_mos_channel_elem()
{
#if 0
  SimulationSystem  & system  = this->system();
  const MeshBase & mesh = system.mesh();

  // build fast surface locator
  SurfaceLocatorHub & surface_locator = mesh.surface_locator();


  // build the bounding box for this boundary
  std::pair<Point, Point> boundary_bounding_box;
  {
    const std::vector<const Node *> & nodes = this->nodes();

    std::vector<Real> min(3, 1.e30);
    std::vector<Real> max(3, -1.e30);

    for (unsigned int n=0; n<nodes.size(); n++)
      for (unsigned int i=0; i<3; i++)
      {
        min[i] = std::min(min[i], (*nodes[n])(i));
        max[i] = std::max(max[i], (*nodes[n])(i));
      }
    Parallel::min(min);
    Parallel::max(max);
    boundary_bounding_box = std::make_pair(Point(&min[0]), Point(&max[0]));
  }


  for(unsigned int r=0; r<system.n_regions(); ++r)
  {
    SimulationRegion * region = system.region(r);
    // only consider nearset point on semiconductor region
    if( region->type() != SemiconductorRegion ) continue;

    const std::pair<Point, Point> & region_bounding_box = region->boundingbox();
    Real distance = MeshTools::minimal_distance(region_bounding_box, boundary_bounding_box);
    if( distance > 100*nm ) continue; // skip far away regions

    SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<SemiconductorSimulationRegion *>(region);

    SimulationRegion::element_iterator elem_it = region->elements_begin();
    SimulationRegion::element_iterator elem_it_end = region->elements_end();
    for(; elem_it != elem_it_end; elem_it++)
    {
      const Elem * e = *elem_it;
      if( !e->on_processor() ) continue;

      const Point centroid = e->centroid();
      // point is far from boundingbox of gate region, skip it
      if( std::abs(MeshTools::minimal_distance(boundary_bounding_box, centroid)) >  100*nm ) continue;

      Point project_point;
      std::pair<const Elem*, unsigned int> surface_elem_pair = surface_locator(centroid, this->boundary_id(), project_point, 100*nm);
      if( surface_elem_pair.first != NULL )
        semiconductor_region->_elem_in_mos_channel.insert(e);
    }
  }
#endif
}


void SimpleGateContactBC::_find_mos_channel_elem()
{
  // also very slow
  SimulationSystem  & system  = this->system();
  const MeshBase & mesh = system.mesh();

  SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<SemiconductorSimulationRegion *>(bc_regions().first);

  std::vector<const Elem *> el;
  std::vector<unsigned int> sl;
  mesh.boundary_info->active_elem_with_boundary_id (el, sl, this->boundary_id());

  for(unsigned int n=0; n<el.size(); ++n)
  {
    const Elem * elem = el[n];
    if( !elem->on_processor() ) continue;

    Point origin_point = elem->centroid();

    // do breadth-first search
    std::set<const Elem *> visited_elem;
    std::queue<const Elem *> queue;

    queue.push(elem);
    visited_elem.insert(elem);

    while(!queue.empty())
    {
      const Elem * current_elem = queue.front();
      queue.pop();

      if( current_elem->on_processor() && (current_elem->centroid() - origin_point).size() < 20*nm )
        semiconductor_region->_elem_in_mos_channel.insert(current_elem);
      for(unsigned int i=0; i<current_elem->n_sides(); ++i)
      {
        const Elem * neighbor_elem = current_elem->neighbor(i);
        if( !neighbor_elem ) continue;

        if( visited_elem.find(neighbor_elem) == visited_elem.end() )
        {
          if( (neighbor_elem->centroid() - origin_point).size() < 20*nm )
            queue.push(neighbor_elem);
        }

        visited_elem.insert(neighbor_elem);
      }
    }
  }
}


void InsulatorSemiconductorInterfaceBC::_find_mos_channel_elem()
{
#if 0
  // too slow
  SimulationSystem  & system  = this->system();
  const MeshBase & mesh = system.mesh();

  // build the bounding box for this boundary
  std::pair<Point, Point> boundary_bounding_box;
  {
    const std::vector<const Node *> & nodes = this->nodes();

    std::vector<Real> min(3, 1.e30);
    std::vector<Real> max(3, -1.e30);

    for (unsigned int n=0; n<nodes.size(); n++)
      for (unsigned int i=0; i<3; i++)
    {
      min[i] = std::min(min[i], (*nodes[n])(i));
      max[i] = std::max(max[i], (*nodes[n])(i));
    }
    Parallel::min(min);
    Parallel::max(max);
    boundary_bounding_box = std::make_pair(Point(&min[0]), Point(&max[0]));
  }

  // build fast surface locator
  SurfaceLocatorHub & surface_locator = mesh.surface_locator();

  SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<SemiconductorSimulationRegion *>(bc_regions().first);
  SimulationRegion::element_iterator elem_it = semiconductor_region->elements_begin();
  SimulationRegion::element_iterator elem_it_end = semiconductor_region->elements_end();
  for(; elem_it != elem_it_end; elem_it++)
  {
    const Elem * e = *elem_it;
    if( !e->on_processor() ) continue;

    const Point centroid = e->centroid();
    // point is far from boundingbox of this region, skip it
    if( std::abs(MeshTools::minimal_distance(boundary_bounding_box, centroid)) >  40*nm ) continue;

    Point project_point;
    std::pair<const Elem*, unsigned int> surface_elem_pair = surface_locator(centroid, this->boundary_id(), project_point, 30*nm);
    if( surface_elem_pair.first != NULL )
      semiconductor_region->_elem_in_mos_channel.insert(e);
  }
#endif

#if 0
  // fast, but may lost some elements
  SimulationSystem  & system  = this->system();
  const MeshBase & mesh = system.mesh();

  SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<SemiconductorSimulationRegion *>(bc_regions().first);

  std::vector<const Elem *> el;
  std::vector<unsigned int> sl;
  mesh.boundary_info->active_elem_with_boundary_id (el, sl, this->boundary_id());

  std::set<const Elem *> visited_elem;
  for(unsigned int n=0; n<el.size(); ++n)
  {
    const Elem * elem = el[n];
    if( !elem->on_processor() ) continue;

    if(visited_elem.find(elem) != visited_elem.end()) continue;
    visited_elem.clear();
    visited_elem.insert(elem);

    Point origin_point = elem->centroid();

    // do breadth-first search
    std::queue<const Elem *> queue;
    queue.push(elem);
    while(!queue.empty())
    {
      const Elem * current_elem = queue.front();
      queue.pop();

      if( current_elem->on_processor() && (current_elem->centroid() - origin_point).size() < 30*nm )
        semiconductor_region->_elem_in_mos_channel.insert(current_elem);
      for(unsigned int i=0; i<current_elem->n_sides(); ++i)
      {
        const Elem * neighbor_elem = current_elem->neighbor(i);
        if( !neighbor_elem ) continue;

        if( visited_elem.find(neighbor_elem) == visited_elem.end() )
        {
          //if( (neighbor_elem->centroid() - origin_point).size() < 30*nm )
            queue.push(neighbor_elem);
        }

        visited_elem.insert(neighbor_elem);
      }
    }
  }
#endif


  SimulationSystem  & system  = this->system();
  const MeshBase & mesh = system.mesh();

  SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<SemiconductorSimulationRegion *>(bc_regions().first);

  std::vector<const Elem *> el;
  std::vector<unsigned int> sl;
  mesh.boundary_info->active_elem_with_boundary_id (el, sl, this->boundary_id());

  for(unsigned int n=0; n<el.size(); ++n)
  {
    const Elem * elem = el[n];
    if( !elem->on_processor() ) continue;
    if( elem->subdomain_id() != semiconductor_region->subdomain_id() ) continue;

    std::set<const Elem *> visited_elem;
    visited_elem.insert(elem);

    Point origin_point = elem->centroid();
    Point dir = -elem->outside_unit_normal(sl[n]);
    Real radii = 1.5*elem->hmax();


    // do breadth-first search
    std::queue<const Elem *> queue;
    queue.push(elem);
    while(!queue.empty())
    {
      const Elem * current_elem = queue.front();
      queue.pop();

      semiconductor_region->_elem_in_mos_channel.insert(current_elem);

      for(unsigned int i=0; i<current_elem->n_sides(); ++i)
      {
        const Elem * neighbor_elem = current_elem->neighbor(i);
        if( !neighbor_elem ) continue;
        if( !neighbor_elem->on_processor() ) continue;
        if( neighbor_elem->subdomain_id() != semiconductor_region->subdomain_id() ) continue;

        if( visited_elem.find(neighbor_elem) == visited_elem.end() )
        {
          Point vec_neighbor = neighbor_elem->centroid() - origin_point;
          if( vec_neighbor*dir > 0 && vec_neighbor*dir < 20*nm && (vec_neighbor - vec_neighbor*dir*dir).size() < radii )
            queue.push(neighbor_elem);
        }

        visited_elem.insert(neighbor_elem);
      }
    }
  }

}

