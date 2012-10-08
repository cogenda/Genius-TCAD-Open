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
#include <time.h>
#include <string>
#include <cstdlib>
#include <iomanip>

#include "solver_base.h"
#include "tunneling_hook.h"

using PhysicalUnit::A;
using PhysicalUnit::V;
using PhysicalUnit::cm;

/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
TunnelingHook::TunnelingHook ( SolverBase & solver, const std::string & name, void * param )
    : Hook ( solver, name ), _output_prefix ( SolverSpecify::out_prefix )
{
  if ( !Genius::processor_id() )
  {
    _out.open((_output_prefix+".tunneling_current.dat").c_str());

    // write file head
    _out << "# Title: Gnuplot File Created by Genius TCAD Simulation" << std::endl << std::endl;

    _out << "# Plotname: tunneling current" << std::endl;
    // write variables
    _out << "# Variables: " << std::endl;

    unsigned int n_var = 0;
    const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
    for(unsigned int n=0; n<bcs->n_bcs(); n++)
    {
      const BoundaryCondition * bc = bcs->get_bc(n);
      if( bc->is_electrode() )
      {
        std::string bc_label = bc->label();
        if(!bc->electrode_label().empty())
          bc_label = bc->electrode_label();
        _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_Vapp"      << " [V]"<< std::endl;
      }

      if ( bc->bc_type() == IF_Insulator_Semiconductor )
      {
        _out << '#' <<'\t' << ++n_var <<'\t' << bc->label() + "_J_CBET"      << " [A*cm^-2]"<< std::endl;
        _out << '#' <<'\t' << ++n_var <<'\t' << bc->label() + "_J_VBHT"      << " [A*cm^-2]"<< std::endl;
        _out << '#' <<'\t' << ++n_var <<'\t' << bc->label() + "_J_VBET"      << " [A*cm^-2]"<< std::endl;
        _out << '#' <<'\t' << ++n_var <<'\t' << bc->label() + "_J"           << " [A*cm^-2]"<< std::endl;
      }
    }

    _out << std::endl;
  }

}


/*----------------------------------------------------------------------
 * destructor, close file
 */
TunnelingHook::~TunnelingHook()
{
}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void TunnelingHook::on_init()
{
}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void TunnelingHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void TunnelingHook::post_solve()
{
  if ( !Genius::processor_id() )
  {
    // set the float number precision
    _out.precision(6);

    // set output width and format
    _out<< std::scientific << std::right;

    const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
    for(unsigned int n=0; n<bcs->n_bcs(); n++)
    {
      const BoundaryCondition * bc = bcs->get_bc(n);

      if( bc->is_electrode() )
      {
        _out << std::setw(15) << bc->ext_circuit()->potential()/V;
      }

      if ( bc->bc_type() == IF_Insulator_Semiconductor )
      {
        _out << std::setw(15) << bc->scalar("J_CBET")/(A/cm/cm)
             << std::setw(15) << bc->scalar("J_VBHT")/(A/cm/cm)
             << std::setw(15) << bc->scalar("J_VBET")/(A/cm/cm)
             << std::setw(15) << (bc->scalar("J_VBHT") - bc->scalar("J_CBET") - bc->scalar("J_VBET"))/(A/cm/cm);
      }
    }

    _out << std::endl;
  }
}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void TunnelingHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void TunnelingHook::on_close()
{}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new TunnelingHook ( solver, name, fun_data );
  }

}

#endif

