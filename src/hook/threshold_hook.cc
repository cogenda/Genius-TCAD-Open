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
#include <iomanip>

#include "mesh_base.h"
#include "solver_base.h"
#include "threshold_hook.h"
#include "parallel.h"

/*
 * usage: HOOK Load=threshold string<region>=(region_name) real<e.field>=(threshold_value) bool<interrupt>=(true|false)
 * <region> specifies which region will be used for threshold evaluation. If ommited, genius will calculate threshold in all regions
 * <e.field> electrical magnitude, in V/cm
 * <interrupt> indicate if genius will exit when exceeding the given threshold
 */

/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
ThresholdHook::ThresholdHook ( SolverBase & solver, const std::string & name, void * param)
    : Hook ( solver, name ), _violate_threshold(false), _stop_when_violate_threshold(false)
{
  const SimulationSystem & system = get_solver().get_system();

  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for ( std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
        parm_it != parm_list.end(); parm_it++ )
  {
    if( parm_it->name() == "region" )
    {
      _region = parm_it->get_string();
      if( system.region(_region) == NULL )
      {
        if( Genius::is_first_processor() )
          std::cerr<<"ThresholdHook: Invalid given region "<< _region  <<  " to be monitor." << std::endl;
        return;
      }
      continue;
    }
    if( parm_it->name() == "x.min" )
    {
      _lower_bound.x() = parm_it->get_real()*PhysicalUnit::um;
      continue;
    }
    if( parm_it->name() == "y.min" )
    {
      _lower_bound.y() = parm_it->get_real()*PhysicalUnit::um;
      continue;
    }
    if( parm_it->name() == "z.min" )
    {
      _lower_bound.z() = parm_it->get_real()*PhysicalUnit::um;
      continue;
    }
    if( parm_it->name() == "x.max" )
    {
      _upper_bound.x() = parm_it->get_real()*PhysicalUnit::um;
      continue;
    }
    if( parm_it->name() == "y.max" )
    {
      _upper_bound.y() = parm_it->get_real()*PhysicalUnit::um;
      continue;
    }
    if( parm_it->name() == "z.max" )
    {
      _upper_bound.z() = parm_it->get_real()*PhysicalUnit::um;
      continue;
    }

    if( parm_it->name() == "prefix"  )
    {
      _threshold_prefix = parm_it->get_string();
      continue;
    }

    if( parm_it->name() == "interrupt" )
    {
      _stop_when_violate_threshold = parm_it->get_bool();
      continue;
    }

    // the variable to be monitor
    {
      SolutionVariable var = solution_string_to_enum (parm_it->name());
      if( var == INVALID_Variable )
      {
        if( Genius::is_first_processor() )
          std::cerr<<"ThresholdHook: Invalid given variable "<< parm_it->name()  <<  " to be monitor." << std::endl;
        return;
      }

      genius_assert( parm_it->type() == Parser::REAL);

      switch(var)
      {
          case E_FIELD     : _vector_variable_threshold_map[var] = parm_it->get_real() * (PhysicalUnit::V/PhysicalUnit::cm); break;
          case TEMPERATURE : _scalar_variable_threshold_map[var] = parm_it->get_real() * (PhysicalUnit::K); break;
          default      : break;
      }
    }
  }

}


/*----------------------------------------------------------------------
 * destructor, close file
 */
ThresholdHook::~ThresholdHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void ThresholdHook::on_init()
{}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void ThresholdHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void ThresholdHook::post_solve()
{
  if(_scalar_variable_threshold_map.find(TEMPERATURE) !=  _scalar_variable_threshold_map.end())
    _check_T_threshold();

  if(_vector_variable_threshold_map.find(E_FIELD) !=  _vector_variable_threshold_map.end())
    _check_E_threshold();
}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void ThresholdHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void ThresholdHook::on_close()
{}


void ThresholdHook::_check_T_threshold()
{}


void ThresholdHook::_check_E_threshold()
{
  const SimulationSystem & system = get_solver().get_system();
  const MeshBase & mesh = system.mesh();

  const Real E_threshold = _vector_variable_threshold_map[E_FIELD];

  bool box = _is_bound_box_valid();

  Real E_magnitude = 0;
  unsigned int cell = invalid_uint;

  for( unsigned int n=0; n<system.n_regions(); ++n )
  {
    const SimulationRegion * region = system.region(n);
    if( !_region.empty() && region->name() != _region ) continue;

    for(unsigned int e=0; e<region->n_cell(); ++e)
    {
      const Elem * elem = region->get_region_elem(e);

      if(box && !_in_bound_box(elem->centroid())) continue;

      const FVM_CellData * elem_data = region->get_region_elem_data(e);
      const Real E = elem_data->E().size();
      if( E > E_magnitude )
      {
        E_magnitude = E;
        cell = elem->id();
      }
    }
  }

  std::map<Real, unsigned int> order;
  if(cell != invalid_uint)
    order.insert(std::make_pair(E_magnitude, cell));

  Parallel::allgather(order);

  if(!order.empty())
  {
    E_magnitude = order.rbegin()->first;
    cell = order.rbegin()->second;

    AutoPtr<Elem> elem = mesh.elem_clone(cell);

    if( Genius::is_first_processor() )
    {
      Point location = elem->centroid()/(PhysicalUnit::um);
      std::cout << "Threshold "<< _threshold_prefix  << ": Max E magnitude " << E_magnitude/(PhysicalUnit::V/PhysicalUnit::cm) << " V/cm "
      <<"at (" << location[0] <<", " <<location[1] <<", "<<location[2] <<")" << std::endl;
    }

    if( E_magnitude > E_threshold )
    {
      if( Genius::is_first_processor() )
      {
        std::cout << "           which exceed threshold " << E_threshold/(PhysicalUnit::V/PhysicalUnit::cm) << " V/cm !" << std::endl;
      }
      if( _violate_threshold == false )
      {
        system.export_vtk ( _threshold_prefix + "device_violate_E_threshold.vtu", false );
        system.export_cgns ( _threshold_prefix + "device_violate_E_threshold.cgns" );
        _violate_threshold = true;
      }

      if( _violate_threshold && _stop_when_violate_threshold )
      {
        Parallel::verify(_violate_threshold);
        std::abort();
      }
    }
  }
  else
  {
    if( Genius::is_first_processor() )
    {
      std::cerr<<"ThresholdHook: no solution exist in given region/bound box." << std::endl;
    }
  }


}

bool ThresholdHook::_is_bound_box_valid()
{
  return _lower_bound.x()<_upper_bound.x() || _lower_bound.y()<_upper_bound.y() || _lower_bound.z()<_upper_bound.z();
}

bool ThresholdHook::_in_bound_box(const Point &p)
{
  return p.x() >= _lower_bound.x() && p.x() <= _upper_bound.x() &&
         p.y() >= _lower_bound.y() && p.y() <= _upper_bound.y() &&
         p.z() >= _lower_bound.z() && p.z() <= _upper_bound.z();
}


#ifndef CYGWIN

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new ThresholdHook ( solver, name, fun_data );
  }

}

#endif

