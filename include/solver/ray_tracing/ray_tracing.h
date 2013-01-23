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

#ifndef __ray_tracing_h__
#define __ray_tracing_h__

//C++ include
#include <string>
#include <vector>
#include <map>
#include <utility>

//local include
#include "elem.h"
#include "edge_edge2.h"
#include "parser.h"
#include "solver_base.h"


class ObjectTree;
class LightThread;
class LightLenses;
class ARCoatings;

/**
 * Ray tracing program  to calculate photogeneration carriers, which is
 * needed when simulating optical device.
 */

class RayTraceSolver : public SolverBase
{
public:
  /**
   * the constructor of RayTraceSolver, take system as parameter
   */
  RayTraceSolver(SimulationSystem & system, const Parser::Card & c);

  /**
   * destructor, do nothing
   */
  ~RayTraceSolver() {}

  /**
   * @return the solver type
   */
  virtual SolverSpecify::SolverType solver_type() const
  {return SolverSpecify::RAY_TRACE;}

  /**
   * virtual function, create the solver
   */
  virtual int create_solver();

  /**
   * virtual function, do the solve process
   */
  virtual int solve();

  /**
   * virtual function, destroy the solver, release internal data
   */
  virtual int destroy_solver();


private:

  /**
   * parameters for ray tracing solver
   */
  const Parser::Card & _card;

  /**
   * read wave spectrum from file
   */
  void parse_spectrum_file(const std::string &);

  /**
   * if we do 2d ray tracing
   */
  unsigned int _dim;

  /**
   * contains refractive_index of each region for a special lamda
   */
  std::map<unsigned int, std::pair<double, double> > _region_refractive_index;

  /**
   * build _region_refractive_index
   */
  void build_region_refractive_index(double lamda);

  /**
   * after build the _region_refractive_index map, we can access refractive index
   */
  double get_refractive_index_re(unsigned int sub_id) const
  { return _region_refractive_index.find(sub_id)->second.first; }

  /**
   * after build the _region_refractive_index map, we can access refractive index
   */
  double get_refractive_index_im(unsigned int sub_id) const
  {return _region_refractive_index.find(sub_id)->second.second;}


  /**
   * carrier density in each semiconductor elem
   */
  std::map<unsigned int, std::pair<double, double> > _elem_carrier_density;

  /**
   * build _elem_carrier_density
   */
  void build_elem_carrier_density();

  /**
   * @return free carrier absorption of given elem
   */
  double get_free_carrier_absorption(const Elem*, double) const;

  /**
   * record all the elements which contains this Node as its vertex
   */
  std::vector<std::vector<const Elem*> > _elems_shared_this_node;

  /**
   * build _elems_shared_this_node map
   */
  void build_elems_node_map();

  /**
   * the '<' operator used in std::map _elems_shared_this_edge
   */
  struct lt_edge
  {
    bool operator()(const Elem * e1, const Elem * e2) const
    {
      unsigned int e1_p1 = e1->node(0);
      unsigned int e1_p2 = e1->node(1);
      if(e1_p1>e1_p2) std::swap(e1_p1, e1_p2);

      unsigned int e2_p1 = e2->node(0);
      unsigned int e2_p2 = e2->node(1);
      if(e2_p1>e2_p2) std::swap(e2_p1, e2_p2);

      if(e1_p2!=e2_p2) return e1_p2<e2_p2;
      return e1_p1<e2_p1;
    }
  };

  /**
   * record all the elements which contains this edge
   */
  std::map<const Elem *, std::vector<const Elem *>,  lt_edge> _elems_shared_this_edge;

  /**
   * build _elems_shared_this_edge map
   */
  void build_elems_edge_map();

  /**
   * record all the boundary/interface node to boundary/interface elem-side pair
   */
  std::map<const Node *, std::vector<std::pair<const Elem*, unsigned int> > >        _boundary_node_to_elem_side_map;

  /**
   * record all the boundary/interface edge to boundary/interface elem-side pair
   */
  std::map<const Elem *, std::vector<std::pair<const Elem*, unsigned int> >, lt_edge> _boundary_edge_to_elem_side_map;

  /**
   * map elem-side to boundary condition
   */
  std::map< std::pair<const Elem *, unsigned int>, const BoundaryCondition * >  _elem_to_surface_map;

  /**
   * build boundary/interface node/edge map and _full_reflect_surface map
   */
  void build_boundary_elems_map();

  /**
   * anti-reflection coatings
   */
  std::map<short int, ARCoatings *>  _arc_surface;

  /**
   * build anti-reflection coating surface
   */
  void build_anti_reflection_coating_surface_map();


  /**
   * @return true if the elem-side is a boundary surface
   */
  bool is_surface(const Elem *, unsigned int) const;

  /**
   * @return true if the elem-side is full reflection
   */
  bool is_full_reflect_surface(const Elem *, unsigned int) const;

  /**
   * @return anti-reflection coating surface, NULL for normal surface
   */
  const ARCoatings * is_anti_reflection_coating_surface(const Elem *, unsigned int) const;

  /**
   * find the first elem ray hit among elements in vector
   */
  const Elem * ray_hit(const Point &p, const Point &dir, const std::vector<const Elem *> &, IntersectionResult&) const;

  /**
   * find the first elem ray hit among elements in vector
   */
  const Elem * ray_hit(const Point &p, const Point &dir, const Elem *, IntersectionResult&) const;

  /**
   * an element tree for fast ray and surface elem intersection determination
   */
  ObjectTree *surface_elem_tree;

  /**
   * when the light source can be considered as plane wave, this struct stores the plane norm to wave direction.
   * we will build a bounding sphere(C,R) of the mesh, then we build the plane with plane_norm = light_direction
   * and plane_center = sphere_center - 2*R*light_direction
   * then we can generate all the start point of the ray and distribute them to all the processors
   */
  struct WavePlane
  {
    /**
     * plane center
     */
    Point center;

    /**
     * plane norm , the same as light direction
     */
    Point norm;

    /**
     * electrical field direction of incident wave
     */
    Point E_dir;

    /**
     * radius of bounding sphere
     */
    double R;

    /**
     * minial distance between 2 rays
     */
    double min_dist;

    /**
     * the effetive area of the ray. ray energy = intensity*ray_area
     */
    double ray_area() const
    { return min_dist*min_dist; }

    /**
     * the start point of each ray, each takes three double value as (X,Y,Z)
     * only record the rays own by local processor
     */
    std::vector<Point> ray_start_points;

    unsigned int n_on_processor_rays() const
    { return static_cast<unsigned int>(ray_start_points.size()); }

    const Point & ray_start_point(unsigned int n) const
    { return ray_start_points[n]; }
  };

  WavePlane _wave_plane;

  /**
   * total ray number for each wave length
   */
  unsigned int _total_rays;

  /**
   * the parameters of optical source
   */
  struct OpticalSource
  {
     double  wave_length;// wave length
     double  power;      // wave power
     double  eta;        // quantum efficiency at this wave length
     bool    eta_auto;   // if it is true, determine quantum efficiency by photon energy and material bandgap
  };

  /**
   * all the optical sources are listed here
   */
  std::vector<OpticalSource> _optical_sources;

  /**
   * create the light source, ray direction and how many rays should be traced
   * if parallel simulation trigged, distribute the task to all the processors
   */
  void create_rays();

  /**
   * lenses system
   */
  LightLenses * _lenses;

  /**
   * create the llenses
   */
  void define_lenses();

  /**
   * do ray tracing of a single ray
   */
  void ray_tracing(LightThread *);

  /**
   * save the energy deposit. for parallel simulation, we must gather this vector
   * from all the processors (call Parallel::sum(_band_absorption_energy_in_elem))
   */
  std::vector<double> _band_absorption_energy_in_elem;

  /**
   * save the energy deposit. for parallel simulation, we must gather this vector
   * from all the processors (call Parallel::sum(_total_absorption_energy_in_elem))
   */
  std::vector<double> _total_absorption_energy_in_elem;

  /**
   * convert energy to optical Generation
   */
  void optical_generation(unsigned int );

  /**
   * calculate total optical carrier generation
   * incident optical energy and total deposite energy
   */
  void statistic() const;
};

#endif

