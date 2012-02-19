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

#include <cstdlib>
#include <dlfcn.h>
#include "boundary_condition.h"
#include "solver_base.h"
#include "injection_hook.h"


InjectionHook::InjectionHook(SolverBase & solver, const std::string & name, void *param)
    : Hook ( solver, name ), electron_inject_density(0), hole_inject_density(0)
{
  std::string injection_boundary;  // find the boundary label given in HOOK command
  std::string injection_current;   // load the injection current expression from external so file

  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for(std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
      parm_it != parm_list.end(); parm_it++)
  {
    if(parm_it->name() == "boundary" && parm_it->type() == Parser::STRING)
      injection_boundary = parm_it->get_string();
    if(parm_it->name() == "current" && parm_it->type() == Parser::STRING)
      injection_current = parm_it->get_string();
  }
  if( injection_boundary.empty() )
  {
    MESSAGE<< "Injection Error: A boundary should be given for carrier inject" <<std::endl;RECORD();
    genius_error();
  }

  // ok, we now load the injection current
  if( injection_current.empty() )
  {
    MESSAGE<< "Injection Error: A dll file should be given for calculating inject current" <<std::endl;RECORD();
    genius_error();
  }

  dll_file = dlopen(injection_current.c_str(), RTLD_LAZY );
  if(dll_file==NULL)
  {
    MESSAGE<< "Injection Error: Can not open dll file "<<injection_current<< " to load inject current" <<std::endl;RECORD();
    genius_error();
  }

  electron_inject_density = (external_funxtion)dlsym (dll_file, "electron_inject_density");
  if(electron_inject_density==NULL)
  {
    MESSAGE<< "Injection Error: Can not load electron inject current from dll file "<<injection_current<<"."<<std::endl;RECORD();
    genius_error();
  }

  hole_inject_density = (external_funxtion)dlsym (dll_file, "hole_inject_density");
  if(hole_inject_density==NULL)
  {
    MESSAGE<< "Injection Error: Can not load hole inject current from dll file "<<injection_current<<"."<<std::endl;RECORD();
    genius_error();
  }


  // find the injection boundary by boundary label
  SimulationSystem &system = _solver.get_system();
  BoundaryConditionCollector * bcs = system.get_bcs();
  _injection_bc = bcs->get_bc(injection_boundary);

  if( _injection_bc == NULL )
  {
    MESSAGE<< "Injection Error: Can not find boundary " << injection_boundary << " for carrier inject" <<std::endl;RECORD();
    genius_error();
  }

  // only NeumannBoundary and IF_Insulator_Semiconductor can be used for injection
  if( _injection_bc->bc_type() != NeumannBoundary && _injection_bc->bc_type() != IF_Insulator_Semiconductor )
  {
    MESSAGE<< "Injection Error: Inject boundary " << injection_boundary
        << " should be neumann or semiconductor-insulator interface" <<std::endl;
    RECORD();
    genius_error();
  }

  // loop all the FVM nodes belogs to this boundary
  BoundaryCondition::node_iterator node_it = _injection_bc->nodes_begin();
  BoundaryCondition::node_iterator end_it = _injection_bc->nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = _injection_bc->region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = _injection_bc->region_node_end(*node_it);
    for(; rnode_it!=end_rnode_it; ++rnode_it  )
    {
      SimulationRegion * region = (*rnode_it).second.first;
      if( region->type() != SemiconductorRegion  ) continue;

      InjectionSurface injection_surface;
      FVM_Node * fvm_node = (*rnode_it).second.second;
      injection_surface.fvm_node = fvm_node;

      // find the surface norm and area
      FVM_Node::fvm_element_iterator elem_it = fvm_node->elem_begin();
      FVM_Node::fvm_element_iterator elem_it_end = fvm_node->elem_end();
      for(; elem_it!=elem_it_end; ++ elem_it)
      {
        const Elem * elem = (*elem_it).first;
        unsigned int side_index = (*elem_it).second;
        Point surface_norm = elem->outside_unit_normal(side_index);
        injection_surface.surface_norm.push_back(surface_norm);

        // build corresponding FVM elem of the side
        AutoPtr<Elem> side = elem->build_side(side_index);
        AutoPtr<Elem> fvm_side = Elem::build (Elem::fvm_compatible_type(side->type()), side->parent());
        for (unsigned int v=0; v < side->n_vertices(); v++)
          fvm_side->set_node(v) = side->get_node(v);
        fvm_side->prepare_for_fvm();

        for (unsigned int v=0; v < fvm_side->n_vertices(); v++)
        {
          const Node * node = fvm_side->get_node(v);
          if(node == fvm_node->root_node())
            injection_surface.surface_area.push_back(fvm_side->partial_volume_truncated(v));
        }
      }
      injection_surface.surface_patches = fvm_node->partly_has_n_elem();

      assert(injection_surface.surface_norm.size() == injection_surface.surface_area.size());
      _injection_surfaces.push_back(injection_surface);
    }
  }
}


InjectionHook::~InjectionHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void InjectionHook::on_init()
{}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void InjectionHook::pre_solve()
{
  // current time (in the unit of second)
  const Real t = SolverSpecify::clock/PhysicalUnit::s;
  // unit of injection current A/cm^2
  const Real scale = PhysicalUnit::A/PhysicalUnit::cm/PhysicalUnit::cm;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  const Real zwidth = _injection_bc->z_width();

  std::vector<InjectionSurface>::iterator it = _injection_surfaces.begin();
  std::vector<InjectionSurface>::iterator it_end = _injection_surfaces.end();
  for(; it!=it_end; ++it)
  {
    FVM_Node * fvm_node = (*it).fvm_node;
    FVM_NodeData * fvm_node_data =  fvm_node->node_data();
    Point p = (*(fvm_node->root_node()))/PhysicalUnit::um;

    Real EIn=0, HIn=0;
    for(unsigned int n=0; n<(*it).surface_patches; ++n)
    {
      // the norm vector of surface
      Point norm = (*it).surface_norm[n];

      // here we consider both 2D/3D case
      Real  area = (*it).surface_area[n]*zwidth;

      // ok, compute injection current here
      EIn += area*electron_inject_density(p.x(), p.y(), p.z(), norm.x(), norm.y(),norm.z(), t)*scale;
      HIn += area*hole_inject_density(p.x(), p.y(), p.z(), norm.x(), norm.y(),norm.z(), t)*scale;
    }

    fvm_node_data->EIn() = EIn;
    fvm_node_data->HIn() = HIn;
  }
}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void InjectionHook::post_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void InjectionHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void InjectionHook::on_close()
{
  std::vector<InjectionSurface>::iterator it = _injection_surfaces.begin();
  std::vector<InjectionSurface>::iterator it_end = _injection_surfaces.end();
  for(; it!=it_end; ++it)
  {
    FVM_Node * fvm_node = (*it).fvm_node;
    FVM_NodeData * fvm_node_data =  fvm_node->node_data();
    fvm_node_data->EIn() = 0.0;
    fvm_node_data->HIn() = 0.0;
  }

  dlclose( dll_file );
}



#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new InjectionHook ( solver, name, fun_data );
  }

}

#endif
