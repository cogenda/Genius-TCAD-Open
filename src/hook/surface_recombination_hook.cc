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
#include <ctime>
#include <cstdlib>
#include <iomanip>

#include "mesh_base.h"
#include "solver_base.h"
#include "semiconductor_region.h"
#include "surface_recombination_hook.h"
#include "parallel.h"


#ifdef WINDOWS
  #include <io.h>      // for windows _access function
#else
  #include <unistd.h>  // for POSIX access function
#endif


/*
 * usage: HOOK Load=surface_recombination string<output>=file_name
 */

/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
SurfaceRecombHook::SurfaceRecombHook ( SolverBase & solver, const std::string & name, void * param)
  : Hook ( solver, name ),_file("surface_recombination.dat")
{
  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for ( std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
        parm_it != parm_list.end(); parm_it++ )
  {
    if( parm_it->name() == "output" )
    {
      _file = parm_it->get_string();
    }
  }

  if ( !Genius::processor_id() )
  {
#ifdef WINDOWS
    bool file_exist = ( _access( (char *) _file.c_str(),  04 ) == 0 );
#else
    bool file_exist = ( access( _file.c_str(),  R_OK ) == 0 );
#endif

    if(file_exist && SolverSpecify::out_append)
      _out.open(_file.c_str(), std::ios::app);
    else
    {
      _out.open(_file.c_str(), std::ios::trunc);
      _write_gnuplot_head();
    }
  }
}


/*----------------------------------------------------------------------
 * destructor, close file
 */
SurfaceRecombHook::~SurfaceRecombHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void SurfaceRecombHook::on_init()
{}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void SurfaceRecombHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void SurfaceRecombHook::post_solve()
{
  const SimulationSystem &system = get_solver().get_system();
  const BoundaryConditionCollector * bcs = system.get_bcs();

  for(unsigned int n=0; n<bcs->n_bcs(); n++)
  {
    const BoundaryCondition * bc = bcs->get_bc(n);

    if(bc->bc_type()==IF_Insulator_Semiconductor)
    {
      double recomb_current;
      double recomb_current_trap_n;
      double recomb_current_trap_p;
      _surface_recombination(bc, recomb_current, recomb_current_trap_n, recomb_current_trap_p);
      _out << std::setw(25) << recomb_current/PhysicalUnit::A;
      _out << std::setw(25) << recomb_current_trap_n/PhysicalUnit::A;
      _out << std::setw(25) << recomb_current_trap_p/PhysicalUnit::A;
    }

    if(bc->bc_type()==NeumannBoundary && bc->flag("surface.recombination") == true)
    {
      double recomb_current;
      double recomb_current_trap_n;
      double recomb_current_trap_p;
      _surface_recombination(bc, recomb_current, recomb_current_trap_n, recomb_current_trap_p);
      _out << std::setw(25) << recomb_current/PhysicalUnit::A;
      _out << std::setw(25) << recomb_current_trap_n/PhysicalUnit::A;
      _out << std::setw(25) << recomb_current_trap_p/PhysicalUnit::A;
    }
  }


  _out << std::endl;
}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void SurfaceRecombHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void SurfaceRecombHook::on_close()
{
  if ( !Genius::processor_id() )
    _out.close();
}


//----------------------------------------------------------------------

void  SurfaceRecombHook::_surface_recombination(const BoundaryCondition * bc, double & s, double &s_trap_e, double &s_trap_h)
{
  PetscScalar recomb_current = 0.0;
  PetscScalar recomb_current_trap_n = 0.0;
  PetscScalar recomb_current_trap_p = 0.0;

  BoundaryCondition::const_node_iterator node_it = bc->nodes_begin();
  BoundaryCondition::const_node_iterator end_it = bc->nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::const_region_node_iterator  rnode_it     = bc->region_node_begin(*node_it);
    BoundaryCondition::const_region_node_iterator  end_rnode_it = bc->region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node =  (*rnode_it).second.second;

      if(region->type() == SemiconductorRegion)
      {
        const SemiconductorSimulationRegion * sregion = (const SemiconductorSimulationRegion *) region;
        const FVM_NodeData * node_data = fvm_node->node_data();

        double n   =  node_data->n();                         // electron density
        double p   =  node_data->p();                         // hole density
        double T   = node_data->T();
        // process interface fixed charge density
        double boundary_area = fvm_node->outside_boundary_surface_area();

        {
          // surface recombination
          Material::MaterialSemiconductor *mt =  sregion->material();
          mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);
          double GSurf = - mt->band->R_Surf(p, n, T) * boundary_area; //generation due to SRH

          recomb_current += GSurf*bc->z_width();
        }

        if( sregion->advanced_model().Trap )
        {
          // consider charge trapping in semiconductor bulk (bulk_flag=true)

          // call the Trap MPI to calculate trap occupancy using the local carrier densities and lattice temperature
          PetscScalar ni = sregion->material()->band->nie(p, n, T);
          sregion->material()->trap->Calculate(false,p,n,ni,T);

          // calculate the rates of electron and hole capture
          PetscScalar TrapElec = sregion->material()->trap->ElectronTrapRate(false,n,ni,T) * boundary_area;
          PetscScalar TrapHole = sregion->material()->trap->HoleTrapRate    (false,p,ni,T) * boundary_area;

          recomb_current_trap_n += TrapElec*bc->z_width();
          recomb_current_trap_p += TrapHole*bc->z_width();
        }

      }
    }
  }

  Parallel::sum(recomb_current);
  Parallel::sum(recomb_current_trap_n);
  Parallel::sum(recomb_current_trap_p);

  s = recomb_current;
  s_trap_e = recomb_current_trap_n;
  s_trap_h = recomb_current_trap_p;
}


void  SurfaceRecombHook::_write_gnuplot_head()
{

  // prepare the file head
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    time_t    _time;
    // get simulation time
    time(&_time);

    // write file head
    _out << "# Title: Gnuplot File for surface recombination current" << std::endl;
    _out << "# Created by Genius TCAD Simulation" << std::endl;
    _out << "# Date: " << ctime(&_time) << std::endl;

    _out << "# Variables: " << std::endl;

    if( SolverSpecify::Type==SolverSpecify::DCSWEEP       ||
        SolverSpecify::Type == SolverSpecify::STEADYSTATE ||
        SolverSpecify::Type == SolverSpecify::OP          ||
        SolverSpecify::Type==SolverSpecify::TRACE         ||
        SolverSpecify::Type==SolverSpecify::TRANSIENT )
    {
      unsigned int n_var = 0;
      // if transient simulation, we need to record time
      if ( SolverSpecify::Type == SolverSpecify::TRANSIENT )
      {
        _out << '#' <<'\t' << ++n_var <<'\t' << "time" << " [s]"<< std::endl;
        _out << '#' <<'\t' << ++n_var <<'\t' << "time_step" << " [s]"<< std::endl;
      }

      // record electrode IV information
      const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
        // search for all the bcs
      for(unsigned int n=0; n<bcs->n_bcs(); n++)
      {
        const BoundaryCondition * bc = bcs->get_bc(n);

        if(bc->bc_type()==IF_Insulator_Semiconductor)
        {
          std::string bc_label = bc->label();

          _out << '#' <<'\t' << ++n_var <<'\t' << bc_label << "  surface recombination current [A]"<< std::endl;
          _out << '#' <<'\t' << ++n_var <<'\t' << bc_label << "  electron trap current [A]"<< std::endl;
          _out << '#' <<'\t' << ++n_var <<'\t' << bc_label << "  hole trap current [A]"<< std::endl;
        }

        if(bc->bc_type()==NeumannBoundary && bc->flag("surface.recombination") == true)
        {
          std::string bc_label = bc->label();

          _out << '#' <<'\t' << ++n_var <<'\t' << bc_label << "  surface recombination current [A]"<< std::endl;
          _out << '#' <<'\t' << ++n_var <<'\t' << bc_label << "  electron trap current [A]"<< std::endl;
          _out << '#' <<'\t' << ++n_var <<'\t' << bc_label << "  hole trap current [A]"<< std::endl;
        }
      }
    }


  }
}



#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new SurfaceRecombHook ( solver, name, fun_data );
  }

}

#endif

