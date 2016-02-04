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
#include <numeric>


#include "dose_rate.h"
#include "simulation_system.h"
#include "simulation_region.h"
#include "boundary_info.h"
#include "mesh_tools.h"
#include "parallel.h"
#include "physical_unit.h"


using PhysicalUnit::um;
using PhysicalUnit::mm;
using PhysicalUnit::kg;


std::vector<double> DoseRate::density;
std::vector<double> DoseRate::constrain;
double DoseRate::min_leaf=0.0;


DoseRate::DoseRate(const SimulationSystem & system )
  : _system(system), _mesh(system.mesh())
{

  DoseRate::build_region_density(system);

  DoseRate::build_mesh_constrain(system.mesh());

  min_leaf = 1.0*mm;

  // mesh bounding box
  std::pair<Point, Point> bbox = MeshTools::global_bounding_box(_mesh);

  const Point cent = 0.5*(bbox.second + bbox.first);
  const Real  rad  = 0.5*(bbox.second - bbox.first).size();

  Point low  = cent - Point(rad, rad, rad);
  Point high = cent + Point(rad, rad, rad);

  OcTreeDataDoseRate * root_data = new OcTreeDataDoseRate();
  this->mesh_elem(root_data->elems);
  root_data->calculate_density();

  _octree = new OcTree(low, high, root_data, 3);
  _octree->refine();

}


DoseRate::~DoseRate()
{
  delete _octree;
}



void DoseRate::refine()
{
  _octree->refine();
}



void DoseRate::energy_deposite(const Point &p1, const Point &p2, double e)
{
  const double length = (p2-p1).size();
  if(length == 0.0) return;

  std::vector<std::pair<OcTreeNode, double> > result;
  _octree->intersect(p1, p2, result);

  for(unsigned int n=0; n<result.size(); ++n)
  {
    const OcTreeNode & node = result[n].first;
    const double l = result[n].second;

    OcTreeDataDoseRate * dose_rate = dynamic_cast<OcTreeDataDoseRate *>(node.data());
    dose_rate->electron_energy += e*l/length;
  }
}



void DoseRate::sync_energy_deposite()
{
  std::vector<double> energy;
  OcTree::tree_leaf_iterator leaf_it = _octree->begin_leaf();
  for ( ; leaf_it != _octree->end_leaf(); ++leaf_it )
  {
    OcTreeDataDoseRate * dose_rate = dynamic_cast<OcTreeDataDoseRate *>(leaf_it->data());
    energy.push_back(dose_rate->electron_energy);
  }
  Parallel::sum(energy);

  leaf_it = _octree->begin_leaf();
  for (unsigned int n=0 ; leaf_it != _octree->end_leaf(); ++leaf_it, ++n )
  {
    OcTreeDataDoseRate * dose_rate = dynamic_cast<OcTreeDataDoseRate *>(leaf_it->data());
    dose_rate->electron_energy = energy[n];
  }
}



void DoseRate::clear_energy_deposite()
{
  OcTree::tree_leaf_iterator leaf_it = _octree->begin_leaf();
  for ( ; leaf_it != _octree->end_leaf(); ++leaf_it )
  {
    OcTreeDataDoseRate * dose_rate = dynamic_cast<OcTreeDataDoseRate *>(leaf_it->data());
    dose_rate->electron_energy = 0.0;
  }
}


double DoseRate::total_energy() const
{
  double energy = 0.0;
  OcTree::tree_leaf_iterator leaf_it = _octree->begin_leaf();
  for ( ; leaf_it != _octree->end_leaf(); ++leaf_it )
  {
    OcTreeDataDoseRate * dose_rate = dynamic_cast<OcTreeDataDoseRate *>(leaf_it->data());
    energy += dose_rate->electron_energy;
  }
  return energy;
}



double DoseRate::energy_deposite_density(const Point & p) const
{
  OcTree::tree_iterator_base it = _octree->find_leaf_has_point(p);
  double volume = it->volume();

  OcTreeDataDoseRate * dose_rate = dynamic_cast<OcTreeDataDoseRate *>(it->data());
  if(dose_rate->weighted_density > 0.0)
    return dose_rate->electron_energy/volume/dose_rate->weighted_density;
  return 0.0;
}



void DoseRate::export_vtk(const std::string &file)
{
  /*
  double energy = 0.0;
  for(double x=0.0*um; x<=3.0*um; x+=0.001*um)
    for(double z=0.0*um; z<=0.2*um; z+=0.001*um)
  {
      energy_deposite(Point(x, 0.1*um, z), Point(x, 2*um, z), 100.0);
      energy += 100.0;
  }


  sync_energy_deposite();
  */

  if(Genius::processor_id() == 0)
    _octree->export_vtk(file);
}




void DoseRate::particle_endpoint(const Point &p)
{
  OcTree::tree_iterator_base it = _octree->find_leaf_has_point(p);
  if( it.node==0 ) return;
  OcTreeDataDoseRate * dose_rate = dynamic_cast<OcTreeDataDoseRate *>(it->data());
  dose_rate->electron_endpoint.push_back(p);
}


void DoseRate::clear_particle_endpoint()
{
  OcTree::tree_leaf_iterator leaf_it = _octree->begin_leaf();
  for ( ; leaf_it != _octree->end_leaf(); ++leaf_it )
  {
    OcTreeDataDoseRate * dose_rate = dynamic_cast<OcTreeDataDoseRate *>(leaf_it->data());
    dose_rate->electron_endpoint.clear();
  }
}



void DoseRate::build_mesh_constrain(const MeshBase & mesh)
{
  constrain.clear();
  constrain.resize(mesh.n_subdomains(), 1e30);

  std::vector<unsigned int>       el;
  std::vector<unsigned short int> sl;
  std::vector<short int>          il;
  mesh.boundary_info->build_on_processor_side_list(el, sl, il);
  for(unsigned int n=0; n<el.size(); n++)
  {
    const Elem * elem = mesh.elem(el[n]);
    constrain[elem->subdomain_id()] = std::min(constrain[elem->subdomain_id()] ,elem->hmin());
  }

  Parallel::min(constrain);
}


void DoseRate::build_region_density(const SimulationSystem & system)
{
  density.clear();
  // density of each region
  for(unsigned int r=0; r<system.n_regions(); r++)
  {
    const SimulationRegion * region = system.region(r);
    density.push_back(region->get_density());
  }
}



void DoseRate::mesh_elem(std::vector<CElem> &elems) const
{
  std::vector<unsigned int> id;
  std::vector<double> location;
  std::vector<unsigned int> subdomain;

  MeshBase::const_element_iterator       el  = _mesh.this_pid_elements_begin();
  const MeshBase::const_element_iterator end = _mesh.this_pid_elements_end();
  for (; el != end; ++el)
  {
    const Elem * elem = *el;
    Point center = elem->centroid();

    id.push_back(elem->id());
    location.push_back(center[0]);
    location.push_back(center[1]);
    location.push_back(center[2]);
    subdomain.push_back(elem->subdomain_id());
  }

  Parallel::allgather(id);
  Parallel::allgather(location);
  Parallel::allgather(subdomain);

  for(unsigned int n=0; n<id.size(); ++n)
  {
    CElem celem;
    celem.id = id[n];
    celem.center = Point(location[3*n+0], location[3*n+1], location[3*n+2]);
    celem.subdomain = subdomain[n];
    elems.push_back(celem);
  }
}


//------------------------------------------------------------------------------------------------

std::vector<OcTreeDataBase *> DoseRate::OcTreeDataDoseRate::subdivide(const OcTreeNode & node)
{
  std::vector<OcTreeDataDoseRate *> child;
  for(unsigned int n=0; n<8; ++n)
    child.push_back(new OcTreeDataDoseRate());

  for( unsigned int e=0; e<elems.size(); ++e)
  {
    const Point & center = elems[e].center;
    OcTreeLocation loc = node.get_sub_location(center);
    if(loc.key() != -1)
      child[loc.key()]->elems.push_back(elems[e]);
  }


  for( unsigned int e=0; e<electron_endpoint.size(); ++e)
  {
    const Point & p = electron_endpoint[e];
    OcTreeLocation loc = node.get_sub_location(p);
    if(loc.key() != -1)
      child[loc.key()]->electron_endpoint.push_back(p);
  }


  for(unsigned int n=0; n<child.size(); ++n)
  {
    child[n]->calculate_density();
    if(child[n]->weighted_density == 0.0)
      child[n]->weighted_density = this->weighted_density;
  }

  for(unsigned int n=0; n<child.size(); ++n)
    child[n]->electron_energy = electron_energy/8.0;

  elems.clear();
  electron_endpoint.clear();
  electron_energy = 0.0;

  std::vector<OcTreeDataBase *> ret;
  for(unsigned int n=0; n<child.size(); ++n)
    ret.push_back((OcTreeDataBase *)child[n]);
  return ret;
}


bool DoseRate::OcTreeDataDoseRate::refine(const OcTreeNode & node, const std::vector<OcTreeNode> &neighbor) const
{
  // too many elems
  if(elems.size() > 30) return true;


  if(elems.size()>1 && node.volume() > std::pow(DoseRate::min_leaf, 3.0)) return true;


  // size constrain
  std::set<unsigned int> subdomain;
  for(unsigned int n=0; n<elems.size(); ++n)
    subdomain.insert(elems[n].subdomain);
  if(subdomain.size() >= 2)
  {
    double length = 1e30;
    std::set<unsigned int>::const_iterator  it = subdomain.begin();
    for(; it != subdomain.end(); ++it)
      length = std::min(length, DoseRate::constrain[*it]);

    if( node.volume() > length*length*length ) return true;
  }


  // stop electron


  // no refine needed
  return false;
}


double DoseRate::OcTreeDataDoseRate::value(const std::string &) const
{
  if(weighted_density > 0.0)
    return electron_energy/weighted_density;
  return 0.0;
}


void DoseRate::OcTreeDataDoseRate::calculate_density()
{
  std::vector<unsigned int> n_elems(DoseRate::density.size(), 0);
  for(unsigned int n=0; n<elems.size(); n++)
    n_elems[elems[n].subdomain]++;

  unsigned int total_elems = std::accumulate(n_elems.begin(), n_elems.end(), 0 );

  weighted_density = 0.0;
  for(unsigned int n=0; n<DoseRate::density.size(); ++n)
    weighted_density += DoseRate::density[n]*n_elems[n]/(1e-10+total_elems);

}


