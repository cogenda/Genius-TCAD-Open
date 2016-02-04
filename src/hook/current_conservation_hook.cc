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

#include <string>
#include <cstdlib>
#include <ctime>
#include <iomanip>

#include "solver_base.h"
#include "current_conservation_hook.h"
#include "parallel.h"


/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
CurrentConservationHook::CurrentConservationHook(SolverBase & solver, const std::string & name, void * param)
    : Hook(solver, name), _probe_file(SolverSpecify::out_prefix + ".current.conservation.dat")
{
  _p_solver = & solver;

  if ( !Genius::processor_id() )
    _out.open(_probe_file.c_str());

}


/*----------------------------------------------------------------------
 * destructor, close file
 */
CurrentConservationHook::~CurrentConservationHook()
{
  if ( !Genius::processor_id() )
    _out.close();
}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void CurrentConservationHook::on_init()
{
  if ( !Genius::processor_id() )
  {
    time_t          _time;
    time(&_time);

    _out << "# Title: Gnuplot File Created by Genius TCAD Simulation" << std::endl;
    _out << "# Date: " << ctime(&_time) << std::endl;
    _out << "# Plotname: current conservation monitor"  << std::endl;
    _out << "# Variables: " << std::endl;
    int cCnt=0;


    // Transient
    if(SolverSpecify::Type==SolverSpecify::TRANSIENT)
    {
      _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << "Time [s]"  << std::endl;
    }


    // for each region
    for( unsigned int r=0; r<_p_solver->get_system().n_regions(); r++)
    {
      const SimulationRegion * region = _p_solver->get_system().region(r);

      if( region->type() == SemiconductorRegion )
      {
        _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << region->name() << " [A] " << std::endl;
        _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << region->name() << " relative "   << std::endl;
      }

      if( region->type() == InsulatorRegion )
      {
        _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << region->name() << " [A] " << std::endl;
        _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << region->name() << " relative "   << std::endl;
      }

      if( region->type() == MetalRegion )
      {
        _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << region->name() << " [A] " << std::endl;
        _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << region->name() << " relative "   << std::endl;
      }

    }

    _out << std::endl;
  }

}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void CurrentConservationHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void CurrentConservationHook::post_solve()
{
  // for each region

  std::vector<double> region_current_sum;
  std::vector<double> region_current_abs_sum;

  for( unsigned int r=0; r<_p_solver->get_system().n_regions(); r++)
  {
    const SimulationRegion * region = _p_solver->get_system().region(r);

    if( region->type() == SemiconductorRegion )
    {
      std::pair<double, double> current = _current_conservation_semiconductor(region);
      region_current_sum.push_back(current.first);
      region_current_abs_sum.push_back(current.second);
    }

    if( region->type() == InsulatorRegion )
    {
      std::pair<double, double> current = _current_conservation_insulator(region);
      region_current_sum.push_back(current.first);
      region_current_abs_sum.push_back(current.second);
    }

    if( region->type() == MetalRegion )
    {
      std::pair<double, double> current = _current_conservation_metal(region);
      region_current_sum.push_back(current.first);
      region_current_abs_sum.push_back(current.second);
    }
  }


  if ( !Genius::processor_id() )
  {
    // set the float number precision
    _out.precision(6);

    // set output width and format
    _out<< std::scientific << std::right;

    _out<<' ';

    // Transient
    if(SolverSpecify::Type==SolverSpecify::TRANSIENT)
    {
      _out << std::setw(15) << SolverSpecify::clock/PhysicalUnit::s;
    }

    for( unsigned int r=0; r<_p_solver->get_system().n_regions(); r++)
    {
      _out << std::setw(20) << region_current_sum[r]/PhysicalUnit::A;
      _out << std::setw(20) << std::abs(region_current_sum[r])/(region_current_abs_sum[r]+1e-30*PhysicalUnit::A);
    }

    _out << std::endl;
  }


}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void CurrentConservationHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void CurrentConservationHook::on_close()
{}



std::pair<double,double> CurrentConservationHook::_current_conservation_semiconductor(const SimulationRegion *region)
{
  double I = 0.0;
  double I_abs = 0.0;

  const std::map<short int, BoundaryCondition *> & boundaries = region->region_boundaries();

  std::map<short int, BoundaryCondition *>::const_iterator it = boundaries.begin();
  for( ; it!= boundaries.end(); it++)
  {
    double current = 0.0;
    BoundaryCondition * bc = it->second;
    if( bc->bc_type() == OhmicContact )    current = bc->ext_circuit()->current();
    if( bc->bc_type() == SchottkyContact ) current = bc->ext_circuit()->current();
    if( bc->bc_type() == IF_Metal_Ohmic )  current = bc->current();
    if( bc->bc_type() == IF_Metal_Schottky )  current = bc->current();
    if( bc->bc_type() == HomoInterface )
    {
      current += bc->scalar("electron_current");
      current += bc->scalar("hole_current");
      current += bc->scalar("displacement_current");
    }
    if( bc->bc_type() == IF_Insulator_Semiconductor )
    {
      current = bc->scalar("interface_current");
    }

    I += current;
    I_abs += std::abs(current);
  }

  return std::make_pair(I, I_abs);

}


std::pair<double,double> CurrentConservationHook::_current_conservation_insulator(const SimulationRegion *region)
{
  double I_displacement = 0.0;
  double I_displacement_abs = 0.0;

  SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
  SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
  for(; it!=it_end; ++it)
  {
    const FVM_Node * fvm_node = *it;
    const FVM_NodeData * fvm_node_data = fvm_node->node_data();
    PetscScalar V_insulator = fvm_node_data->psi();

    if( fvm_node->on_boundary() )
    {
      if(SolverSpecify::TimeDependent == true)
      {
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
        {
          const FVM_Node *nb_node = (*nb_it).first;
          const FVM_NodeData * nb_node_data = nb_node->node_data();
          // the psi of neighbor node
          PetscScalar V_nb = nb_node_data->psi();
          // distance from nb node to this node
          PetscScalar distance = fvm_node->distance(nb_node);
          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = fvm_node->cv_surface_area(nb_node);
          PetscScalar dEdt;
          if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false) //second order
          {
            PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
            dEdt = ( (2-r)/(1-r)*(V_insulator-V_nb)
                     - 1.0/(r*(1-r))*(fvm_node_data->psi_last()-nb_node_data->psi_last())
                     + (1-r)/r*(fvm_node_data->psi_old()-nb_node_data->psi_old()))/distance/(SolverSpecify::dt_last+SolverSpecify::dt);
          }
          else//first order
          {
            dEdt = ((V_insulator-V_nb)-(fvm_node_data->psi_last()-nb_node_data->psi_last()))/distance/SolverSpecify::dt;
          }

          I_displacement += cv_boundary*fvm_node_data->eps()*dEdt;
          I_displacement_abs += std::abs(cv_boundary*fvm_node_data->eps()*dEdt);
        }
      }
    }

  }

  Parallel::sum(I_displacement);
  Parallel::sum(I_displacement_abs);

  return std::make_pair(I_displacement, I_displacement_abs);

}



std::pair<double,double> CurrentConservationHook::_current_conservation_metal(const SimulationRegion *region)
{
  double I_conductance = 0.0;
  double I_conductance_abs = 0.0;

  const double sigma = region->get_conductance();

  SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
  SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
  for(; it!=it_end; ++it)
  {
    const FVM_Node * fvm_node = *it;
    const FVM_NodeData * fvm_node_data = fvm_node->node_data();
    PetscScalar V_metal = fvm_node_data->psi();

    if( fvm_node->on_boundary() )
    {
      {
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
        {
          const FVM_Node *nb_node = (*nb_it).first;
          const FVM_NodeData * nb_node_data = nb_node->node_data();
          // the psi of neighbor node
          PetscScalar V_nb = nb_node_data->psi();
          // distance from nb node to this node
          PetscScalar distance = fvm_node->distance(nb_node);
          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = std::abs(fvm_node->cv_surface_area(nb_node));

          I_conductance += sigma*cv_boundary*(V_metal-V_nb)/distance;
          I_conductance_abs += std::abs(sigma*cv_boundary*(V_metal-V_nb)/distance);
        }
      }
    }
  }

  Parallel::sum(I_conductance);
  Parallel::sum(I_conductance_abs);

  return std::make_pair(I_conductance, I_conductance_abs);
}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook (SolverBase & solver, const std::string & name, void * fun_data)
  {
    return new CurrentConservationHook(solver, name, fun_data );
  }

}

#endif

