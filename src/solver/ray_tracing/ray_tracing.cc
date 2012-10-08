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

#include <stack>
#include <iomanip>
#include <numeric>


#include "sphere.h"
#include "mesh_base.h"
#include "boundary_info.h"
#include "semiconductor_region.h"
#include "mesh_tools.h"
#include "field_source.h"
#include "light_lenses.h"
#include "ray_tracing/light_thread.h"
#include "ray_tracing/object_tree.h"
#include "ray_tracing/ray_tracing.h"
#include "ray_tracing/anti_reflection_coating.h"
#include "parallel.h"


using PhysicalUnit::um;
using PhysicalUnit::J;
using PhysicalUnit::s;
using PhysicalUnit::cm;
using PhysicalUnit::m;
using PhysicalUnit::eps0;
using PhysicalUnit::mu0;
using PhysicalUnit::h;



RayTraceSolver::RayTraceSolver(SimulationSystem & system, const Parser::Card & c)
    : SolverBase(system), _card(c), surface_elem_tree(0)
{
  system.record_active_solver(this->solver_type());
}


int RayTraceSolver::create_solver()
{
  MESSAGE<< '\n' << "Ray tracing Solver init... ";
  RECORD();

  // since we need serial mesh here
  MeshBase & mesh = _system.mesh();
  genius_assert(mesh.is_serial());

  // do necessary precomputation for fast ray tracing
  surface_elem_tree = new ObjectTree(mesh);
  build_elems_node_map();
  build_elems_edge_map();
  build_boundary_elems_map();
  build_anti_reflection_coating_surface_map();
  build_elem_carrier_density();
  // parse input deck
  define_lenses();
  create_rays();

  MESSAGE<< _total_rays <<" rays for each wave length."<<std::endl;
  RECORD();

  return 0;
}


int RayTraceSolver::solve()
{
  START_LOG("solve()", "RayTraceSolver");

  // for each wavelentgh
  for(unsigned int n=0; n<_optical_sources.size(); ++n)
  {
    double lamda     = _optical_sources[n].wave_length;
    double intensity = _optical_sources[n].power;
    double power     = _dim==2 ? intensity*_wave_plane.min_dist : intensity*_wave_plane.ray_area();

    build_region_refractive_index(lamda);

    // clear and re-create the array to record energy deposition
    _band_absorption_energy_in_elem.clear();
    _band_absorption_energy_in_elem.resize(_system.mesh().n_elem(), 0.0);

    _total_absorption_energy_in_elem.clear();
    _total_absorption_energy_in_elem.resize(_system.mesh().n_elem(), 0.0);

    MESSAGE<< "  process light of " /*<< std::setiosflags(std::ios::fixed)*/  << lamda/um << " um";
    RECORD();

    //process all the rays
    unsigned int n_on_processor_rays = _wave_plane.n_on_processor_rays();
    for(unsigned int k=0; k<n_on_processor_rays; ++k)
    {
      // create ray
      LightThread * light = new  LightThread(_wave_plane.ray_start_point(k),
                                             _wave_plane.norm,
                                             _wave_plane.E_dir,
                                             lamda,
                                             power,
                                             power
                                            );

      if(!_lenses->empty())  light = (*_lenses) << light;

      // call function ray_tracing to process a single ray
      ray_tracing(light);

      //indicator
      if(k%(1+(n_on_processor_rays)/20)==0) // +1 for prevent divide by zero error
      {
        MESSAGE<< ".";
        RECORD();
      }
#if defined(HAVE_FENV_H) && defined(DEBUG)
      genius_assert( !fetestexcept(FE_INVALID) );
#endif

    }

    // gather energy deposit from all the processors
    Parallel::sum(_band_absorption_energy_in_elem);
    Parallel::sum(_total_absorption_energy_in_elem);

    // convert energy deposit to carrier optical generation
    optical_generation(n);

    MESSAGE<< "ok" <<std::endl;
    RECORD();
  }
#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  STOP_LOG("solve()", "RayTraceSolver");

  return 0;
}


int RayTraceSolver::destroy_solver()
{
  delete surface_elem_tree;

  {
    std::map<const Elem *, std::vector<const Elem *>,  lt_edge>::iterator it = _elems_shared_this_edge.begin();
    for(; it!=_elems_shared_this_edge.end(); ++it)
      delete it->first;
    _elems_shared_this_edge.clear();

  }

  {
    std::map<const Elem *, std::vector<std::pair<const Elem*, unsigned int> >, lt_edge>::iterator it = _boundary_edge_to_elem_side_map.begin();
    for(; it!=_boundary_edge_to_elem_side_map.end(); ++it)
      delete it->first;
    _boundary_edge_to_elem_side_map.clear();

  }

  {
    std::map<short int, ARCoatings *>::const_iterator it = _arc_surface.begin();
    for( ; it != _arc_surface.end(); ++it)
      delete it->second;

  }


  MESSAGE<< "Ray tracing Solver finished.\n" <<std::endl;
  RECORD();

  return 0;
}


//---------------------------------------------------------------


void RayTraceSolver::define_lenses()
{
  std::vector<std::string> lenses;
  if(_card.is_parameter_exist("lenses"))
    lenses = _card.get_array<std::string>("lenses");

  for(unsigned int idx=0; idx<_card.parameter_size(); idx++)
  {
    Parser::Parameter p = _card.get_parameter(idx);
    if ( p.name() == "lens" )
    {
      lenses.push_back(p.get_string());
    }
  }

  _lenses = _system.get_field_source()->light_lenses();
  _lenses->set_lenses(lenses);
}


void RayTraceSolver::create_rays()
{

  // optical wave is defined by command line
  if(_card.is_parameter_exist("lambda")||_card.is_parameter_exist("wavelength"))
  {
    OpticalSource source;
    source.wave_length = _card.get_real("lambda", 0.532, "wavelength")*um;// wave length of light source
    source.power       = _card.get_real("intensity", 0.0)*J/s/cm/cm;      // incident power
    source.eta_auto    = !_card.is_parameter_exist("quan.eff");
    source.eta         = _card.get_real("quan.eff", 1.0);
    _optical_sources.push_back(source);
  }

  // read optical wave from spectrum file
  if(_card.is_parameter_exist("spectrumfile"))
  {
    parse_spectrum_file(_card.get_string("spectrumfile", ""));
  }


  const MeshBase &mesh = _system.mesh();

  // ray direction
  Point dir;
  if(_card.is_parameter_exist("k"))
  {
    std::vector<double> k = _card.get_array<double>("k");
    k.resize(3, 0.0);
    dir = Point(&k[0]).unit();
  }
  else
  {
    double phi   = _card.get_real("k.phi",   0.0)/180.0*3.14159265358979323846;
    double theta = _card.get_real("k.theta", 0.0)/180.0*3.14159265358979323846;
    dir = Point(sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta));
  }

  // E vector direction
  Point E_dir;
  if(_card.is_parameter_exist("e"))
  {
    std::vector<double> e = _card.get_array<double>("e");
    e.resize(3, 0.0);
    E_dir = Point(&e[0]).unit();
  }
  else
  {
    double phi   = _card.get_real("e.phi",  90.0)/180.0*3.14159265358979323846;
    double theta = _card.get_real("e.theta", 0.0)/180.0*3.14159265358979323846;
    E_dir = Point(sin(phi)*cos(theta), cos(phi), sin(phi)*sin(theta));
  }

  // check if E//k
  if( std::abs(dir.dot(E_dir)) > 1.0 - 1e-3 )
  {
    MESSAGE<<"\nERROR at " <<_card.get_fileline()<< " Ray Tracing: E vector not perpendicular to wave vector." << std::endl; RECORD();
    genius_error();
  }
  // compute the E
  E_dir = (E_dir - dir*E_dir*dir).unit();

  // determine the ray density -- the distance between rays
  double min_dist = 1e+10;
  MeshBase::const_element_iterator       el  = mesh.elements_begin();
  const MeshBase::const_element_iterator el_end = mesh.elements_end();
  for (; el != el_end; ++el)
  {
    double hmin = (*el)->hmin();
    if( hmin<min_dist ) min_dist = hmin;
  }

  // the mesh bounding sphere (circle in 2D)
  Sphere bounding_sphere = MeshTools::bounding_sphere(mesh);
  // consider lenses
  if( !_lenses->empty() )
  {
    Sphere lenses_bounding_sphere = _lenses->bounding_sphere();
    bounding_sphere = MeshTools::bounding_sphere(bounding_sphere, lenses_bounding_sphere);
  }

  // the ray start plane
  if(_card.is_parameter_exist("ray.center"))
  {
    std::vector<double> c = _card.get_array<double>("ray.center");
    c.resize(3, 0.0); // if use only gives one cood, pad it
    _wave_plane.center    = Point(&c[0])*um;
    if( bounding_sphere.below_surface(_wave_plane.center) )
      _wave_plane.center  = _wave_plane.center - 2*bounding_sphere.radius()*dir;
  }
  else
    _wave_plane.center    = bounding_sphere.center() - 2*bounding_sphere.radius()*dir;

  _wave_plane.norm     = dir;
  _wave_plane.E_dir    = E_dir;
  if(_card.is_parameter_exist("ray.radius"))
    _wave_plane.R      = _card.get_real("ray.radius",  100.0)*um;
  else
    _wave_plane.R      = bounding_sphere.radius();
  if(_card.is_parameter_exist("ray.distance"))
    _wave_plane.min_dist = _card.get_real("ray.distance", 0.1)*um;
  else
    _wave_plane.min_dist = min_dist/_card.get_real("ray.density", 10.0);

  // this vector stores all the ray start points
  std::vector<Point> ray_start_points;

  //3D mesh
  if(surface_elem_tree->is_octree())
  {
    // compute all the start point of the ray
    const Point y_axis(0,1,0);

    // two direction vector in ray start plane
    Point d1, d2;

    // if ray dir is parallel to y axis
    if((y_axis.cross(dir)).size()<=1e-9)
    {
      d1 = Point(1, 0, 0); // x-axis
      d2 = dir.cross(d1);
    }
    else // ray is not parallel to y axis
    {
      d1 = ((y_axis.cross(dir)).cross(dir)).unit();
      d2 = dir.cross(d1);
    }

    // how many rays in one of the direction
    int n_rays_half = int(ceil(_wave_plane.R/_wave_plane.min_dist));
    for(int i=-n_rays_half; i<=n_rays_half; ++i)
      for(int j=-n_rays_half; j<=n_rays_half; ++j)
      {
        Point s = _wave_plane.center + i*_wave_plane.min_dist*d1 + j*_wave_plane.min_dist*d2;
        if(surface_elem_tree->hit_boundbox(s, dir))
          ray_start_points.push_back(s);
      }
    _dim = 3;
  }

  //2D mesh
  if(surface_elem_tree->is_quadtree())
  {
    const Point z_axis(0,0,1);
    Point d = dir.cross(z_axis);

    int n_rays_half = int(ceil(_wave_plane.R/_wave_plane.min_dist));
    for(int i=-n_rays_half; i<=n_rays_half; ++i)
    {
      Point s = _wave_plane.center + i*_wave_plane.min_dist*d;
      if(surface_elem_tree->hit_boundbox(s, dir))
        ray_start_points.push_back(s);
    }

    _dim = 2;
  }

  // compute start point on all the processors
  _total_rays = ray_start_points.size();
  unsigned int on_process_rays = ray_start_points.size()/Genius::n_processors();
  unsigned int begin = on_process_rays*Genius::processor_id();
  unsigned int end   = std::min(begin+on_process_rays, ray_start_points.size());
  for(unsigned int n=begin; n<end; ++n)
    _wave_plane.ray_start_points.push_back(ray_start_points[n]);

}



void RayTraceSolver::parse_spectrum_file(const std::string & filename)
{
  // only processor 0 read the spectrum file
  std::vector<OpticalSource> _opt_srcs;
  if(Genius::processor_id() == 0)
  {
    std::ifstream in(filename.c_str());
    assert(in.good());

    std::string line;
    while(!in.eof())
    {
      std::getline(in, line);

      // skip empty lines
      if(line.size()==0) continue;
      // skip the line begin with '#'
      if(line.find('#')==0) continue;

      std::stringstream ss;
      ss << line;

      OpticalSource source;
      // read wave length and power information from line
      ss >> source.wave_length;
      ss >> source.power;
      // test if it contains extra parameter for quantum efficiency
      ss >> source.eta;
      // if no quantum efficiency parameter find, set eta_auto flag to true
      source.eta_auto = ss.fail();

      _opt_srcs.push_back(source);
    }
    in.close();
  }

  // broadcast to all the processors
  unsigned int n_sources = _opt_srcs.size();
  if (Genius::n_processors() > 1)
  {
    Parallel::broadcast(n_sources);
    if (Genius::processor_id() != 0)
      _opt_srcs.resize(n_sources);
    for(unsigned int n=0; n<n_sources; ++n)
    {
      std::vector<double> package(4);
      if (Genius::processor_id() == 0)
      {
        package[0] = _opt_srcs[n].wave_length;
        package[1] = _opt_srcs[n].power;
        package[2] = double(_opt_srcs[n].eta_auto);
        package[3] = _opt_srcs[n].eta;
      }

      Parallel::broadcast(package);

      if (Genius::processor_id() != 0)
      {
        _opt_srcs[n].wave_length = package[0];
        _opt_srcs[n].power       = package[1];
        _opt_srcs[n].eta_auto    = package[2]>0.5;
        _opt_srcs[n].eta         = package[3];
      }
    }
  }

  if(n_sources<2)
  {
    MESSAGE<< "ERROR: Spectrum file "<< filename << "should contain at least two spectrums.\n" <<std::endl;
    RECORD();
    genius_error();
  }

  //post process


  double total_power = 0.0;
  for(unsigned int n=0; n<n_sources; ++n)
  {
    double lambda_distance;
    if(n==0)
      lambda_distance = 0.5*(_opt_srcs[1].wave_length -_opt_srcs[0].wave_length);
    else if(n==n_sources-1)
      lambda_distance = 0.5*(_opt_srcs[n_sources-1].wave_length -_opt_srcs[n_sources-2].wave_length);
    else
      lambda_distance = 0.5*(_opt_srcs[n+1].wave_length -_opt_srcs[n-1].wave_length);

    OpticalSource source;
    source.wave_length = _opt_srcs[n].wave_length*um;
    source.power = (_opt_srcs[n].power*lambda_distance)*J/s/(cm*cm);
    source.eta_auto    = _opt_srcs[n].eta_auto;
    source.eta         = _opt_srcs[n].eta;
    _optical_sources.push_back(source);

    total_power += _opt_srcs[n].power*lambda_distance;
  }

  MESSAGE<<"Read "<<n_sources<<" spectrums from file "<<filename<<".\n"
  <<"Total power intensity is approximate "<<std::setprecision(3)<<std::setw(10)<<total_power<<" W/(cm^2)"
  <<std::endl;

  RECORD();

}

void RayTraceSolver::build_region_refractive_index(double lamda)
{
  _region_refractive_index.clear();

  for(unsigned int n=0; n<_system.n_regions(); ++n)
  {
    SimulationRegion* region = _system.region(n);
    Complex r_index = region->get_optical_refraction(lamda);

    std::pair<double, double> refract_index = std::make_pair(r_index.real(), r_index.imag());

    _region_refractive_index[n] = refract_index;
  }

  // set env refractive index
  _region_refractive_index[invalid_uint] = std::make_pair(1.0, 0.0);
}


void RayTraceSolver::build_elem_carrier_density()
{
  std::vector<unsigned int> elem_id;
  std::map<unsigned int, double> elem_n;
  std::map<unsigned int, double> elem_p;

  for(unsigned int r=0; r<_system.n_regions(); ++r)
  {
    const SimulationRegion* region = _system.region(r);
    if( region->type() != SemiconductorRegion ) continue;

    SimulationRegion::const_element_iterator it = region->elements_begin();
    SimulationRegion::const_element_iterator it_end = region->elements_end();
    for(; it != it_end; ++it)
    {
      const Elem * elem = *it;
      if( !elem->on_processor() ) continue;

      double n=0, p=0;
      for(unsigned int i=0; i<elem->n_nodes(); ++i)
      {
        const FVM_Node * fvm_node = elem->get_fvm_node(i);
        const FVM_NodeData * node_data = fvm_node->node_data();
        n += node_data->n()/elem->n_nodes();
        p += node_data->p()/elem->n_nodes();
      }
      elem_id.push_back(elem->id());
      elem_n.insert(std::make_pair(elem->id(), n));
      elem_p.insert(std::make_pair(elem->id(), p));
    }
  }

  Parallel::allgather(elem_id);
  Parallel::allgather(elem_n);
  Parallel::allgather(elem_p);

  genius_assert(elem_id.size() == elem_n.size());
  genius_assert(elem_id.size() == elem_p.size());

  for(unsigned int i=0; i<elem_id.size(); ++i)
  {
    double n = elem_n.find(elem_id[i])->second;
    double p = elem_p.find(elem_id[i])->second;
    _elem_carrier_density[elem_id[i]] = std::make_pair(n, p);
  }
}


double RayTraceSolver::get_free_carrier_absorption(const Elem* elem, double wavelength) const
{
  const SimulationRegion * region =  _system.region(elem->subdomain_id());
  if( region->type() == SemiconductorRegion )
  {
    const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
    std::pair<double, double> carrier = _elem_carrier_density.find(elem->id())->second;
    return semiconductor_region->material()->optical->FreeCarrierAbsorption(wavelength, carrier.first, carrier.second, _system.T_external() );
  }
  return 0.0;
}



void RayTraceSolver::build_elems_node_map()
{
  MeshTools::build_nodes_to_elem_map(_system.mesh(), _elems_shared_this_node);
}


void RayTraceSolver::build_elems_edge_map()
{
  const MeshBase &mesh = _system.mesh();
  MeshBase::const_element_iterator       el  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();

  for (; el != end; ++el)
    for (unsigned int e=0; e<(*el)->n_edges(); e++)
    {
      bool edge_exist = false;
      const Elem * edge_elem = (*el)->build_edge(e).release();
      if(_elems_shared_this_edge.find(edge_elem) != _elems_shared_this_edge.end())
        edge_exist = true;

      _elems_shared_this_edge[edge_elem].push_back(*el);

      if(edge_exist) delete edge_elem;
    }
}



void RayTraceSolver::build_anti_reflection_coating_surface_map()
{
  // build anti-reflection coatings
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    const BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    if( bc->has_scalar_array("coatings"))
    {
      ARCoatings * arc = new ARCoatings;
      arc->inner_region = bc->bc_subdomains().first;
      arc->layers = 0;
      std::vector<double>  raw = bc->scalar_array("coatings");
      for(unsigned int i=0; i <raw.size(); )
      {
        arc->layer_refractive_index.push_back( std::complex<double>(raw[i++],  -raw[i++]) );
        arc->layer_thickness.push_back(raw[i++]*um);
        arc->layers++;
      }

      _arc_surface.insert( std::make_pair(bc->boundary_id(), arc) );
    }
  }

}


void RayTraceSolver::build_boundary_elems_map()
{

  const MeshBase &mesh = _system.mesh();
  std::vector<unsigned int>        el;
  std::vector<unsigned short int>  sl;
  std::vector<short int>           il;

  mesh.boundary_info->build_side_list(el, sl, il);
  for(unsigned int n=0; n<el.size(); ++n)
  {
    const Elem * boundary_elem = mesh.elem(el[n]);
    AutoPtr<Elem> side = boundary_elem->build_side(sl[n], false);

    for(unsigned int nd=0; nd<side->n_nodes(); ++nd)
    {
      const Node * node = side->get_node(nd);
      _boundary_node_to_elem_side_map[node].push_back(std::make_pair(boundary_elem, sl[n]));
    }

    for(unsigned int ne=0; ne<side->n_edges(); ++ne)
    {
      bool edge_exist = false;
      const Elem * edge = side->build_side(ne, false).release();
      if(_boundary_edge_to_elem_side_map.find(edge)!=_boundary_edge_to_elem_side_map.end())
        edge_exist = true;

      _boundary_edge_to_elem_side_map[edge].push_back(std::make_pair(boundary_elem, sl[n]));

      if(edge_exist) delete edge;
    }

    const BoundaryCondition * bc =  _system.get_bcs()->get_bc_by_bd_id(il[n]);
    genius_assert(bc);
    _elem_to_surface_map[std::make_pair(boundary_elem, sl[n])] = bc;
  }

}


bool RayTraceSolver::is_surface(const Elem *elem, unsigned int side) const
{
  return ( _elem_to_surface_map.find(std::make_pair(elem, side)) != _elem_to_surface_map.end() );
}


bool RayTraceSolver::is_full_reflect_surface(const Elem *elem, unsigned int side) const
{
  const BoundaryCondition * bc = _elem_to_surface_map.find(std::make_pair(elem, side))->second;
  if(bc->has_flag("reflection") && bc->flag("reflection"))
    return true;
  return false;
}


const ARCoatings * RayTraceSolver::is_anti_reflection_coating_surface(const Elem *elem, unsigned int side) const
{
  const BoundaryCondition * bc = _elem_to_surface_map.find(std::make_pair(elem, side))->second;
  std::map<short int, ARCoatings *>::const_iterator pair = _arc_surface.find(bc->boundary_id());
  if( pair == _arc_surface.end() )
    return 0;
  return pair->second;
}




const Elem * RayTraceSolver::ray_hit(const Point &p, const Point &dir, const std::vector<const Elem *> &elems, IntersectionResult & result) const
{
  const Elem * hit_elem = NULL;
  double dist = 1e10;
  for(unsigned int n=0; n<elems.size(); ++n)
  {
    const Elem * elem = elems[n];
    IntersectionResult ret;
    elem->ray_hit(p, dir, ret, _dim);
    if(ret.hit_points.size()<2) continue; //missed or only one intersection point
    if(ret.hit_points[0].t<-1e-10) continue; //when t<0, it is not ray hit..
    if(ret.hit_points[0].t<dist)
    {
      dist     = ret.hit_points[0].t;
      hit_elem = elem;
      result   = ret;
    }
  }
  return hit_elem;
}


const Elem * RayTraceSolver::ray_hit(const Point &p, const Point &dir, const Elem *elem, IntersectionResult & result) const
{
  elem->ray_hit(p, dir, result, _dim);

  //missed
  if(result.state==Missed) return NULL;

  // normal result
  if(result.hit_points.size()==2)  return elem;

  // the ray is tangent with this elem, we shoud find more detailed information
  if(result.hit_points.size()==1)
  {
    Hit_Point hit_point = result.hit_points[0];
    result.hit_points.clear();
    switch(hit_point.point_location)
    {
    case  on_face  :
    case  on_side  : return elem; break; //impossible!
    case  on_edge  :
      {
        unsigned int edge_index = hit_point.mark;
        AutoPtr<Elem> edge = elem->build_edge(edge_index);
        const std::vector<const Elem *> & elems = _elems_shared_this_edge.find(edge.get())->second;
        const Elem * next_elem = this->ray_hit(p, dir, elems, result);
        assert(next_elem!=elem);
        return next_elem;
      }
    case on_vertex :
      {
        unsigned int vertex_index = hit_point.mark;
        const Node * node = elem->get_node(vertex_index);
        const std::vector<const Elem *> & elems = _elems_shared_this_node[node->id()];
        const Elem * next_elem = this->ray_hit(p, dir, elems, result);
        assert(next_elem!=elem);
        return next_elem;
      }
    }
  }

  //we should never reach here
  genius_error();
  return NULL;
}



void RayTraceSolver::ray_tracing(LightThread *ray)
{

  // use stack to save all the rays (origin and secondary)
  std::stack<LightThread *> ray_stack;
  ray_stack.push(ray);

  while(!ray_stack.empty())
  {
    LightThread * current_ray = ray_stack.top();
    ray_stack.pop();

    if(current_ray==NULL) continue;

    // the ray doesn't hit any elem yet?
    if(current_ray->hit_elem==NULL)
    {
      // find the first element this ray hit
      const Elem * elem = surface_elem_tree->hit(current_ray);
      if(elem==NULL) {delete current_ray; continue;}

      current_ray->hit_elem = elem;

      assert(elem->on_boundary());
      // get the intersection result
      elem->ray_hit(current_ray->start_point(), current_ray->dir(), current_ray->result, _dim);

      // the first intersection point
      assert(current_ray->result.hit_points.size());

      std::vector<std::pair<const Elem*, unsigned int> > hit_elems;
      Hit_Point & hit_point = current_ray->result.hit_points[0];
      switch(hit_point.point_location)
      {
      case on_face   :  //the ray-mesh intersection point is on elem face
      case on_side   :  //2d case
        {
          if(elem->on_boundary(hit_point.mark)) //safe guard
            hit_elems.push_back( std::make_pair(elem, hit_point.mark) );
          break;
        }
      case on_edge   : //the ray-mesh intersection point is on elem edge, we must split the ray
        {
          assert(_dim==3); //only valid for 3d mesh
          unsigned int edge_index = hit_point.mark;
          AutoPtr<Elem> edge = elem->build_edge(edge_index);
          assert(_boundary_edge_to_elem_side_map.find(edge.get())!=_boundary_edge_to_elem_side_map.end());
          hit_elems = _boundary_edge_to_elem_side_map.find(edge.get())->second;
          break;
        }
      case on_vertex :  //the ray-mesh intersection point is on elem vertex, we must split the ray
        {
          unsigned int vertex_index = hit_point.mark;
          const Node * current_node = elem->get_node(vertex_index);
          assert(_boundary_node_to_elem_side_map.find(current_node)!=_boundary_node_to_elem_side_map.end());
          hit_elems = _boundary_node_to_elem_side_map.find(current_node)->second;
          break;
        }
      }

      current_ray->power() = current_ray->power()/(hit_elems.size()+1e-10);
      for(unsigned int n=0; n<hit_elems.size(); ++n)
      {
        const Elem * boundary_elem = hit_elems[n].first;
        unsigned int side = hit_elems[n].second;
        assert(is_surface(boundary_elem, side));

        Point p = hit_point.p;
        Point norm = boundary_elem->outside_unit_normal(side);
        //the surface norm should has a angle >90 degree to ray dir
        if(norm.dot(current_ray->dir()) > -1e-10) continue;

        // if reflect surface
        if(is_full_reflect_surface(boundary_elem, side))
        {
          LightThread * reflect_ray = current_ray->reflection(p, norm);
          // some stupid skill: shift the reflect ray to prevent it hit this elem again
          reflect_ray->start_point() = reflect_ray->start_point() + 1e-6*reflect_ray->dir();
          // the reflect ray hit the mesh again?
          const Elem *surface_elem = surface_elem_tree->hit(reflect_ray);
          if(surface_elem && surface_elem!=elem)
          {
            reflect_ray->hit_elem = this->ray_hit(reflect_ray->start_point(), reflect_ray->dir(), surface_elem, reflect_ray->result);
            if(reflect_ray->hit_elem)
              ray_stack.push(reflect_ray);
            else
              delete reflect_ray;
          }
          else
            delete reflect_ray;
          continue;
        }


        double n1 = get_refractive_index_re(boundary_elem->subdomain_id(side));
        double n2 = get_refractive_index_re(boundary_elem->subdomain_id());

        // if we have a anti reflection coating
        const ARCoatings * arc = is_anti_reflection_coating_surface(boundary_elem, side);
        bool inv = (arc && arc->inner_region != boundary_elem->subdomain_id());

        //generate reflect/refract rays
        std::pair<LightThread *, LightThread *> ray_pair = current_ray->interface_light_gen_linear_polarized(p, norm, n1, n2, arc, inv);

        // for refract ray
        LightThread *refract_ray = ray_pair.second;
        if(refract_ray)
        {
          refract_ray->hit_elem = this->ray_hit(refract_ray->start_point(), refract_ray->dir(), boundary_elem, refract_ray->result);
          if(refract_ray->hit_elem)
            ray_stack.push(refract_ray);
          else
            delete refract_ray;
        }

        // for reflect ray
        LightThread *reflect_ray = ray_pair.first;
        if(reflect_ray)
        {
          // some stupid skill: shift the reflect ray to prevent it hit this elem again
          reflect_ray->start_point() = reflect_ray->start_point() + 1e-8*reflect_ray->dir();
          // the refract ray hit the mesh again?
          const Elem *surface_elem = surface_elem_tree->hit(reflect_ray);
          if(surface_elem && surface_elem!=elem)
          {
            reflect_ray->hit_elem = this->ray_hit(reflect_ray->start_point(), reflect_ray->dir(), surface_elem, reflect_ray->result);
            if(reflect_ray->hit_elem)
              ray_stack.push(reflect_ray);
            else
              delete reflect_ray;
          }
          else
            delete reflect_ray;
        }
      }

      delete current_ray;
      continue;
    }


    //assert(current_ray->result.hit_points.size()==2);
    if( current_ray->result.hit_points.size() != 2 )
    {
      // FIXME, should not happen...
      delete current_ray;
      continue;
    }

    // calculate energy deposit
    const Elem * elem = current_ray->hit_elem;
    Hit_Point  end_point = current_ray->result.hit_points[1];


    double a_band = 4*3.14159265358979*this->get_refractive_index_im(elem->subdomain_id())/current_ray->wavelength();
    double a_tail = 0.0;
    double a_fc   = this->get_free_carrier_absorption(elem, current_ray->wavelength());

    std::vector<double> energy_deposit = current_ray->advance_to(end_point.p, a_band, a_tail, a_fc);
    double total_energy_deposit = std::accumulate(energy_deposit.begin(), energy_deposit.end(), 0.0);

    switch(current_ray->result.state)
    {
      // all the energy deposited in this elem
    case Intersect_Body :
      _band_absorption_energy_in_elem[elem->id()] += energy_deposit[0];
      _total_absorption_energy_in_elem[elem->id()] += total_energy_deposit;
      break;
      // two elem shares the energy deposite
    case On_Face        :
      {
        _band_absorption_energy_in_elem[elem->id()] += 0.5*energy_deposit[0];
        _total_absorption_energy_in_elem[elem->id()] += 0.5*total_energy_deposit;
        unsigned int side = current_ray->result.mark;
        const Elem * neighbor = elem->neighbor(side);
        if(neighbor)
        {
          _band_absorption_energy_in_elem[neighbor->id()] += 0.5*energy_deposit[0];
          _total_absorption_energy_in_elem[neighbor->id()] += 0.5*total_energy_deposit;
        }
        break;
      }
      // all the elems have this edge shares the deposited energy
    case Overlap_Edge   :
      {
        unsigned int e = current_ray->result.mark;
        AutoPtr<Elem> edge = elem->build_edge(e);
        const std::vector<const Elem *> & elems = _elems_shared_this_edge.find(edge.get())->second;
        assert(elems.size());
        for(unsigned int n=0; n<elems.size(); ++n)
        {
          _band_absorption_energy_in_elem[elems[n]->id()] += energy_deposit[0]/elems.size();
          _total_absorption_energy_in_elem[elems[n]->id()] += total_energy_deposit/elems.size();
        }
        break;
      }
      //we should never reach here
    default: genius_error();
    }


    if(current_ray->is_dead())
    { delete current_ray; continue; }

    // safe guard: when the number of rays in stack exceed 1000, we may fall into endless loop
    // force to exit
    if(ray_stack.size()>1000)
    {
      while(!ray_stack.empty())
      {
        LightThread * current_ray = ray_stack.top();
        ray_stack.pop();
        delete current_ray;
      }
      return;
    }

    // find next ray elem intersection

    current_ray->result.hit_points.clear();
    switch(end_point.point_location)
    {
    case  on_face   : //3d
    case  on_side   : //2d
      {
        unsigned int side = end_point.mark;
        const Elem * next_elem = elem->neighbor(side);
        if(next_elem && next_elem->subdomain_id() == elem->subdomain_id())
        {
          current_ray->hit_elem = next_elem;
          next_elem->ray_hit(current_ray->start_point(), current_ray->dir(), current_ray->result, _dim);
          ray_stack.push(current_ray);
        }
        else //we are on material interface
        {
          // if reflect surface
          if(is_surface(elem, side) && is_full_reflect_surface(elem, side))
          {  delete current_ray; continue; }

          Point p = end_point.p;
          Point norm = - elem->outside_unit_normal(side);
          double n1 = get_refractive_index_re(elem->subdomain_id());
          double n2 = get_refractive_index_re(elem->subdomain_id(side));
          //generate reflect/refract rays
          std::pair<LightThread *, LightThread *> ray_pair = current_ray->interface_light_gen_linear_polarized(p, norm, n1, n2);

          // for refract ray
          LightThread *refract_ray = ray_pair.second;
          if(refract_ray)
          {
            refract_ray->hit_elem = next_elem;
            if(next_elem)
              next_elem->ray_hit(refract_ray->start_point(), refract_ray->dir(), refract_ray->result, _dim);
            else
              refract_ray->start_point() = refract_ray->start_point() + 1e-8*refract_ray->dir();
            ray_stack.push(refract_ray);
          }

          // for reflect ray
          LightThread *reflect_ray = ray_pair.first;
          if(reflect_ray)
          {
            reflect_ray->hit_elem = elem;
            elem->ray_hit(reflect_ray->start_point(), reflect_ray->dir(), reflect_ray->result, _dim);
            assert(reflect_ray->result.state!=Missed);
            ray_stack.push(reflect_ray);
          }
          delete current_ray;
        }
      }
      break;
    case  on_edge   :
      {
        unsigned int edge_index = end_point.mark;
        AutoPtr<Elem> edge = elem->build_edge(edge_index);
        // the edge is not on boundary
        if(_boundary_edge_to_elem_side_map.find(edge.get())==_boundary_edge_to_elem_side_map.end())
        {
          const std::vector<const Elem *> & edge_elems = _elems_shared_this_edge.find(edge.get())->second;
          const std::vector<const Elem *> & node1_elems = _elems_shared_this_node[edge->node(0)];
          const std::vector<const Elem *> & node2_elems = _elems_shared_this_node[edge->node(1)];
          std::set<const Elem *> elems_set;
          elems_set.insert(edge_elems.begin(), edge_elems.end());
          elems_set.insert(node1_elems.begin(), node1_elems.end());
          elems_set.insert(node2_elems.begin(), node2_elems.end());
          std::vector<const Elem *> elems;
          elems.insert(elems.end(), elems_set.begin(), elems_set.end());

          const Elem * next_elem = this->ray_hit(current_ray->start_point(), current_ray->dir(), elems, current_ray->result);
          if(next_elem && next_elem!=elem)
          {
            current_ray->hit_elem = next_elem;
            ray_stack.push(current_ray);
            continue;
          }
        }
        else // the edge is on boundary
        {
          Point p = end_point.p;
          unsigned int edge_index = end_point.mark;
          AutoPtr<Elem> edge = elem->build_edge(edge_index);
          const std::vector<std::pair<const Elem*, unsigned int> > & elems = _boundary_edge_to_elem_side_map.find(edge.get())->second;

          //split current_ray, each edge on boundary shoud be shared by 2 elem
          unsigned int effective_faces = 2;
          current_ray->power() = current_ray->power()/effective_faces;

          for(unsigned int n=0; n<elems.size(); ++n)
          {
            const Elem * boundary_elem = elems[n].first;
            unsigned int side = elems[n].second;

            // if reflect surface
            if(is_surface(boundary_elem, side) && is_full_reflect_surface(boundary_elem, side)) continue;

            Point norm = boundary_elem->outside_unit_normal(side);
            //the surface norm should has a angle >90 degree to ray dir
            if(norm.dot(current_ray->dir()) > -1e-10) continue;

            double n1 = get_refractive_index_re(boundary_elem->subdomain_id(side));
            double n2 = get_refractive_index_re(boundary_elem->subdomain_id());
            //generate reflect/refract rays
            std::pair<LightThread *, LightThread *> ray_pair = current_ray->interface_light_gen_linear_polarized(p, norm, n1, n2);

            // for refract ray
            LightThread *refract_ray = ray_pair.second;
            if(refract_ray)
            {
              refract_ray->hit_elem = this->ray_hit(refract_ray->start_point(), refract_ray->dir(), boundary_elem, refract_ray->result);
              if(refract_ray->hit_elem)
                ray_stack.push(refract_ray);
              else
                delete refract_ray;
            }

            // for reflect ray
            LightThread *reflect_ray = ray_pair.first;
            if(reflect_ray)
            {
              reflect_ray->hit_elem = this->ray_hit(reflect_ray->start_point(), reflect_ray->dir(), elem, reflect_ray->result);
              if(reflect_ray->hit_elem)
                ray_stack.push(reflect_ray);
              else
                delete reflect_ray;
            }
          }

          delete current_ray;
        }
      }
      break;
    case  on_vertex :
      {
        unsigned int vertex_index = end_point.mark;
        const Node * node = elem->get_node(vertex_index);
        std::vector<const Elem *> & elems = _elems_shared_this_node[node->id()];
        // the node is not on boundary
        if( _boundary_node_to_elem_side_map.find(node)==_boundary_node_to_elem_side_map.end())
        {
          const Elem * next_elem = this->ray_hit(current_ray->start_point(), current_ray->dir(), elems, current_ray->result);
          assert(next_elem && next_elem!=elem);
          current_ray->hit_elem = next_elem;
          ray_stack.push(current_ray);
          continue;
        }
        else
        {
          Point p = end_point.p;
          unsigned int vertex_index = end_point.mark;
          const Node * current_node = elem->get_node(vertex_index);
          const std::vector<std::pair<const Elem*, unsigned int> > & elems = _boundary_node_to_elem_side_map.find(current_node)->second;

          //split current_ray
          unsigned int effective_faces = 0;
          for(unsigned int n=0; n<elems.size(); ++n)
          {
            Point norm = elems[n].first->outside_unit_normal(elems[n].second);
            //the surface norm should has a angle >90 degree to ray dir
            if(norm.dot(current_ray->dir()) > -1e-10) continue;
            effective_faces++;
          }
          current_ray->power() = current_ray->power()/effective_faces;

          for(unsigned int n=0; n<elems.size(); ++n)
          {
            const Elem * boundary_elem = elems[n].first;
            unsigned int side = elems[n].second;

            // if reflect surface
            if(is_surface(boundary_elem, side) && is_full_reflect_surface(boundary_elem, side)) continue;

            Point norm = boundary_elem->outside_unit_normal(side);
            //the surface norm should has a angle >90 degree to ray dir
            if(norm.dot(current_ray->dir()) > -1e-10) continue;

            double n1 = get_refractive_index_re(boundary_elem->subdomain_id(side));
            double n2 = get_refractive_index_re(boundary_elem->subdomain_id());
            //generate reflect/refract rays
            std::pair<LightThread *, LightThread *> ray_pair = current_ray->interface_light_gen_linear_polarized(p, norm, n1, n2);

            // for refract ray
            LightThread *refract_ray = ray_pair.second;
            if(refract_ray)
            {
              refract_ray->hit_elem = this->ray_hit(refract_ray->start_point(), refract_ray->dir(), boundary_elem, refract_ray->result);
              if(refract_ray->hit_elem)
                ray_stack.push(refract_ray);
              else
                delete refract_ray;
            }

            // for reflect ray
            LightThread *reflect_ray = ray_pair.first;
            if(reflect_ray)
            {
              reflect_ray->hit_elem = this->ray_hit(reflect_ray->start_point(), reflect_ray->dir(), elem, reflect_ray->result);
              if(reflect_ray->hit_elem)
                ray_stack.push(reflect_ray);
              else
                delete reflect_ray;
            }
          }
          delete current_ray;
        }
      }
      break;
      //we should never reach here
    default: genius_error();
    }

  }
}



void RayTraceSolver::optical_generation(unsigned int n)
{
  const MeshBase &mesh = _system.mesh();

  double c = 1.0/sqrt(eps0*mu0);
  double lamda    = _optical_sources[n].wave_length;
  double quan_eff = _optical_sources[n].eta;

  for (unsigned int i=0; i<_band_absorption_energy_in_elem.size(); i++)
  {
    const Elem * elem = mesh.elem(i);
    if(!elem->on_local()) continue; //skip nonlocal elements

    SimulationRegion* elem_region = _system.region(elem->subdomain_id());

    // only semiconductor region generating carriers
    if(elem_region->type() == SemiconductorRegion)
    {
      double Eg = elem_region->get_optical_Eg(elem_region->T_external());
      double E_photon = h*c/lamda;
      if(_optical_sources[n].eta_auto)
      {
        // calculate optical gen quantum efficiency
        quan_eff = floor(E_photon/Eg);
      }

      double gen = _band_absorption_energy_in_elem[i]/E_photon*quan_eff;
      double heat = _band_absorption_energy_in_elem[i] - Eg*gen;

      double volumn = 0;
      for(unsigned int nd=0; nd<elem->n_nodes(); nd++)
      {
        volumn += elem->partial_volume_truncated(nd);
      }

      for(unsigned int nd=0; nd<elem->n_nodes(); nd++)
      {
        const Node* node = elem->get_node(nd);
        double vol_ratio = elem->partial_volume_truncated(nd)/(volumn+1e-10);

        FVM_Node* fvm_node = elem_region->region_fvm_node(node);
        assert(fvm_node->node_data());

        fvm_node->node_data()->OptG() += gen*vol_ratio/fvm_node->volume();
        fvm_node->node_data()->OptQ() += heat*vol_ratio/fvm_node->volume();
      }
    }
  }


  for (unsigned int i=0; i<_total_absorption_energy_in_elem.size(); i++)
  {
    const Elem * elem = mesh.elem(i);
    if(!elem->on_local()) continue; //skip nonlocal elements

    SimulationRegion* elem_region = _system.region(elem->subdomain_id());

    double energy = _total_absorption_energy_in_elem[i];
    if(energy==0.0) continue;

    double volumn = 0;
    for(unsigned int nd=0; nd<elem->n_nodes(); nd++)
    {
      volumn += elem->partial_volume_truncated(nd);
    }

    for(unsigned int nd=0; nd<elem->n_nodes(); nd++)
    {
      const Node* node = elem->get_node(nd);
      double vol_ratio = elem->partial_volume_truncated(nd)/(volumn+1e-10);

      FVM_Node* fvm_node = elem_region->region_fvm_node(node);
      assert(fvm_node->node_data());

      fvm_node->node_data()->OptE() += energy*vol_ratio/fvm_node->volume();
    }
  }

}


