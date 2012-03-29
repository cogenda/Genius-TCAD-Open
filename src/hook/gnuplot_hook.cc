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

//  $Id: gnuplot_hook.cc,v 1.7 2008/07/09 12:25:19 gdiso Exp $

#include <string>
#include <cstdlib>
#include <iomanip>

#include "solver_base.h"
#include "gnuplot_hook.h"
#include "spice_ckt.h"
#include "mxml.h"
#include "MXMLUtil.h"

/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
GnuplotHook::GnuplotHook(SolverBase & solver, const std::string & name, void * file)
    : Hook(solver, name), _input_file((const char *)file),
    _gnuplot_file(SolverSpecify::out_prefix + ".dat"), _mixA(false)
{

  SolverSpecify::SolverType solver_type = this->get_solver().solver_type();

  // if we are called by mixA solver?
  if(solver_type == SolverSpecify::DDML1MIXA  ||
     solver_type == SolverSpecify::DDML2MIXA ||
     solver_type == SolverSpecify::EBML3MIXA
    )
    _mixA = true;


  if ( !Genius::processor_id() )
  {
#ifdef WINDOWS
    bool file_exist = ( _access( (char *) _gnuplot_file.c_str(),  04 ) == 0 );
#else
    bool file_exist = ( access( _gnuplot_file.c_str(),  R_OK ) == 0 );
#endif

    if(file_exist && SolverSpecify::out_append)
      _out.open(_gnuplot_file.c_str(), std::ios::app);
    else
    {
      _out.open(_gnuplot_file.c_str(), std::ios::trunc);
      _write_gnuplot_head();
    }
  }

}


/*----------------------------------------------------------------------
 * destructor, close file
 */
GnuplotHook::~GnuplotHook()
{

}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void GnuplotHook::on_init()
{}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void GnuplotHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void GnuplotHook::post_solve()
{

  // save electrode IV
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    // set the float number precision
    _out.precision(6);

    // set output width and format
    _out<< std::scientific << std::right;

    if( SolverSpecify::Type==SolverSpecify::DCSWEEP ||
        SolverSpecify::Type==SolverSpecify::TRACE   ||
        SolverSpecify::Type==SolverSpecify::TRANSIENT )
    {
      // if transient simulation, we need to record time
      if (SolverSpecify::Type == SolverSpecify::TRANSIENT)
      {
        _out << SolverSpecify::clock/PhysicalUnit::s << '\t';
        _out << std::setw(15) << SolverSpecify::dt/PhysicalUnit::s;
      }


      if( _mixA )
      {
        const SPICE_CKT * spice_ckt = this->get_solver().get_system().get_circuit();
        for(unsigned int n=0; n<spice_ckt->n_ckt_nodes(); n++)
        {
          _out << std::setw(15) << spice_ckt->get_solution(n);
        }
      }
      else
      {
        // search for all the bc
        const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
        for(unsigned int n=0; n<bcs->n_bcs(); n++)
        {
          const BoundaryCondition * bc = bcs->get_bc(n);

          // electrode
          if( bc->is_electrode() )
          {
            //record vapp, electrode potential and electrode current
            _out << std::setw(15) << bc->ext_circuit()->Vapp()/PhysicalUnit::V;
            _out << std::setw(15) << bc->ext_circuit()->potential()/PhysicalUnit::V;
            _out << std::setw(15) << bc->ext_circuit()->current()/PhysicalUnit::A;

            if( bc->bc_type() == OhmicContact )
            {
              //_out << std::setw(15) << bc->ext_circuit()->current_displacement()/PhysicalUnit::A;
              _out << std::setw(15) << bc->ext_circuit()->current_electron()/PhysicalUnit::A;
              _out << std::setw(15) << bc->ext_circuit()->current_hole()/PhysicalUnit::A;
            }
            continue;
          }

          if( bc->has_current_flow() )
          {
            _out << std::setw(15) << bc->current()/PhysicalUnit::A;
          }

          if( bc->bc_type() == IF_Metal_Ohmic || bc->bc_type() == IF_Metal_Schottky)
          {
            _out << std::setw(15) << bc->psi()/PhysicalUnit::V;
          }

          // charge integral interface
          if( bc->bc_type() == ChargeIntegral )
          {
            _out << std::setw(15) << bc->scalar("qf")/PhysicalUnit::C;
            _out << std::setw(15) << bc->psi()/PhysicalUnit::V;
          }
        }
      }

    }


    if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
    {
      PetscScalar omega = 2*3.14159265358979323846*SolverSpecify::Freq*PhysicalUnit::s;
      _out << SolverSpecify::Freq*PhysicalUnit::s << '\t';

      // search for all the bc
      const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
      const BoundaryCondition * ac_bc = bcs->get_bc(SolverSpecify::Electrode_ACScan[0]);

      for(unsigned int n=0; n<bcs->n_bcs(); n++)
      {
        const BoundaryCondition * bc = bcs->get_bc(n);
        // skip bc which is not electrode
        if( !bc->is_electrode() ) continue;

        // DC potential and current
        _out << std::setw(15) << bc->ext_circuit()->potential()/PhysicalUnit::V;
        _out << std::setw(15) << bc->ext_circuit()->current()/PhysicalUnit::A;

        //record electrode potential and electrode current for AC simulation
        _out << std::setw(15) << std::abs(bc->ext_circuit()->potential_ac())/PhysicalUnit::V;
        _out << std::setw(15) << std::arg(bc->ext_circuit()->potential_ac());

        _out << std::setw(15) << std::abs(bc->ext_circuit()->current_ac())/PhysicalUnit::A;
        _out << std::setw(15) << std::arg(bc->ext_circuit()->current_ac());

        std::complex<PetscScalar> Y;
        if(ac_bc->ext_circuit()->potential_ac() != 0.0)
          Y = (bc->ext_circuit()->current_ac()/PhysicalUnit::A)/(ac_bc->ext_circuit()->potential_ac()/PhysicalUnit::V);
        else
          Y = 0.0;
        _out << std::setw(15) << Y.real();
        _out << std::setw(15) << Y.imag()/omega;
      }
    }

    _out << std::endl;

    ////----
    {
      mxml_node_t *eSolution = get_solver().current_dom_solution_elem();
      if (eSolution)
      {
        mxml_node_t *eOutput  = mxmlFindElement(eSolution, eSolution, "output", NULL, NULL, MXML_DESCEND_FIRST);
        mxml_node_t *eGnuplot = mxmlNewElement(eOutput, "gnuplot");
        mxml_node_t *eFile    = mxmlNewElement(eGnuplot, "file");
        mxmlAdd(eFile, MXML_ADD_AFTER, NULL, MXMLQVariant::makeQVString(_gnuplot_file));
      }
    }
  }

}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void GnuplotHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void GnuplotHook::on_close()
{
  if ( !Genius::processor_id() )
    _out.close();
}



void  GnuplotHook::_write_gnuplot_head()
{

  // prepare the file head
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    // get simulation time
    time(&_time);

    // write file head
    _out << "# Title: Gnuplot File Created by Genius TCAD Simulation" << std::endl;
    _out << "# Date: " << ctime(&_time) << std::endl;

    switch (SolverSpecify::Type)
    {
      case SolverSpecify::DCSWEEP   :
        _out << "# Plotname: DC transfer characteristic" << std::endl; break;
      case SolverSpecify::TRACE     :
        _out << "# Plotname: DC curve trace" << std::endl; break;
      case SolverSpecify::TRANSIENT :
        _out << "# Plotname: Transient Analysis" << std::endl; break;
      case SolverSpecify::ACSWEEP   :
        _out << "# Plotname: AC small signal Analysis" << std::endl; break;
        default: break;
    }

    // write variables
    _out << "# Variables: " << std::endl;

    if( SolverSpecify::Type==SolverSpecify::DCSWEEP ||
        SolverSpecify::Type==SolverSpecify::TRACE   ||
        SolverSpecify::Type==SolverSpecify::TRANSIENT )
    {
      unsigned int n_var = 0;
      // if transient simulation, we need to record time
      if ( SolverSpecify::Type == SolverSpecify::TRANSIENT )
      {
        _out << '#' <<'\t' << ++n_var <<'\t' << "time" << " [s]"<< std::endl;
        _out << '#' <<'\t' << ++n_var <<'\t' << "time_step" << " [s]"<< std::endl;
      }

      if( _mixA ) // mix mode
      {
        const SPICE_CKT * spice_ckt = this->get_solver().get_system().get_circuit();
        for(unsigned int n=0; n<spice_ckt->n_ckt_nodes(); n++)
        {
          if(spice_ckt->is_voltage_node(n))
            _out << '#' <<'\t' << ++n_var <<'\t' << spice_ckt->ckt_node_name(n)   << " [V]"<< std::endl;
          else
            _out << '#' <<'\t' << ++n_var <<'\t' << spice_ckt->ckt_node_name(n)   << " [A]"<< std::endl;
        }
      }
      else
      {
        // record electrode IV information
        const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
        // search for all the bcs
        for(unsigned int n=0; n<bcs->n_bcs(); n++)
        {
          const BoundaryCondition * bc = bcs->get_bc(n);
          // for electrode
          if( bc->is_electrode() )
          {
            std::string bc_label = bc->label();
            if(!bc->electrode_label().empty())
              bc_label = bc->electrode_label();

            _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_Vapp"      << " [V]"<< std::endl;
            _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_potential" << " [V]"<< std::endl;
            _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_current"   << " [A]"<< std::endl;

            if( bc->bc_type() == OhmicContact )
            {
              //_out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_displacement_current"   << " [A]"<< std::endl;
              _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_electron_current"   << " [A]"<< std::endl;
              _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_hole_current"   << " [A]"<< std::endl;
            }

            continue;
          }

          if( bc->has_current_flow() )
          {
            std::string bc_label = bc->label();
            _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_current"   << " [A]"<< std::endl;
          }

          if( bc->bc_type() == IF_Metal_Ohmic || bc->bc_type() == IF_Metal_Schottky)
          {
            std::string bc_label = bc->label();
            _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_average_potential"   << " [V]"<< std::endl;
          }

          // charge integral interface
          if( bc->bc_type() == ChargeIntegral )
          {
            std::string bc_label = bc->label();
            _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_Q"          << " [C]"<< std::endl;
            _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_potential"  << " [V]"<< std::endl;
          }
        }
      }

    }

    if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
    {
      unsigned int n_var = 0;

      _out << '#' <<'\t' << ++n_var <<'\t' << "frequency" << " [Hz]"<< std::endl;

      // record electrode IV information
      const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();

      const BoundaryCondition * ac_bc = bcs->get_bc(SolverSpecify::Electrode_ACScan[0]);
      std::string ac_bc_label = ac_bc->label();

      // search for all the bcs
      for(unsigned int n=0; n<bcs->n_bcs(); n++)
      {
        const BoundaryCondition * bc = bcs->get_bc(n);
        // skip bc which is not electrode
        if( !bc->is_electrode() ) continue;
        //
        std::string bc_label = bc->label();
        if(!bc->electrode_label().empty())
          bc_label = bc->electrode_label();
        // DC
        _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_dc_potential" << " [V]"<< std::endl;
        _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_dc_current" << " [A]"<< std::endl;
        // AC
        _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_ac_potential_magnitude" << " [V]"<< std::endl;
        _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_ac_potential_angle"     << "    "<< std::endl;

        _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_ac_current_magnitude"   << " [A]"<< std::endl;
        _out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_ac_current_angle"       << "    "<< std::endl;

        _out << '#' <<'\t' << ++n_var <<'\t' << "G_" + ac_bc_label + '_' +bc_label   << " [S]"<< std::endl;
        _out << '#' <<'\t' << ++n_var <<'\t' << "C_" + ac_bc_label + '_' +bc_label   << " [F]"<< std::endl;
      }
    }

    _out << std::endl;

  }

}



#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook (SolverBase & solver, const std::string & name, void * fun_data)
  {
    return new GnuplotHook(solver, name, fun_data );
  }

}

#endif

