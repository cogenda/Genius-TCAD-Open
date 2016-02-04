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
#include <ctime>
#include <string>
#include <cstdlib>
#include <iomanip>

#include "config.h"
#ifdef WINDOWS
#else
#include <unistd.h>
#endif

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
    _gnuplot_file(SolverSpecify::out_prefix + ".dat"), _ddm(false), _mixA(false)
{

  SolverSpecify::SolverType solver_type = this->get_solver().solver_type();

  // if we are called by ddm solver?
  if(solver_type == SolverSpecify::DDML1 || solver_type == SolverSpecify::DENSITY_GRADIENT  ||
     solver_type == SolverSpecify::DDML2 || solver_type == SolverSpecify::EBML3
    )
    _ddm = true;

  // if we are called by mixA solver?
  if(solver_type == SolverSpecify::DDML1MIXA || solver_type == SolverSpecify::DDML1MIX  ||
     solver_type == SolverSpecify::DDML2MIXA || solver_type == SolverSpecify::DDML2MIX  ||
     solver_type == SolverSpecify::EBML3MIXA || solver_type == SolverSpecify::EBML3MIX
    )
    _mixA = true;

  // FIXME temporary fix
  if( SolverSpecify::Type==SolverSpecify::OP )
    _gnuplot_file = SolverSpecify::out_prefix + ".op.dat";

  if( SolverSpecify::Type==SolverSpecify::STEADYSTATE )
    _gnuplot_file = SolverSpecify::out_prefix + ".steady.dat";


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
    _out.precision(14);

    // set output width and format
    _out<< std::scientific << std::right;

    if( SolverSpecify::Type==SolverSpecify::DCSWEEP       ||
        SolverSpecify::Type == SolverSpecify::STEADYSTATE ||
        SolverSpecify::Type == SolverSpecify::OP          ||
        SolverSpecify::Type==SolverSpecify::TRACE         ||
        SolverSpecify::Type==SolverSpecify::TRANSIENT )
    {
      // if transient simulation, we need to record time
      if (SolverSpecify::Type == SolverSpecify::TRANSIENT)
      {
        _out << SolverSpecify::clock/PhysicalUnit::s << '\t';
        _out << std::setw(25) << SolverSpecify::dt/PhysicalUnit::s;
      }


      if( _mixA )
      {
        const SPICE_CKT * spice_ckt = this->get_solver().get_system().get_circuit();
        for(unsigned int n=0; n<spice_ckt->n_ckt_nodes(); n++)
        {
          _out << std::setw(25) << spice_ckt->get_solution(n);
        }

        const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
        for(unsigned int n=0; n<bcs->n_bcs(); n++)
        {
          const BoundaryCondition * bc = bcs->get_bc(n);

          // electrode
          if( bc->has_current_flow() )
          {
            _out << std::setw(25) << bc->current()/PhysicalUnit::A;
          }
        }

      }
      else
      {
        // search for all the bc
        double power = 0;
        
        const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
        for(unsigned int n=0; n<bcs->n_bcs(); n++)
        {
          const BoundaryCondition * bc = bcs->get_bc(n);

          // electrode
          if( bc->is_electrode() )
          {
            //record vapp, electrode potential and electrode current
            _out << std::setw(25) << bc->ext_circuit()->Vapp()/PhysicalUnit::V;
            _out << std::setw(25) << bc->ext_circuit()->potential()/PhysicalUnit::V;
            _out << std::setw(25) << bc->ext_circuit()->current()/PhysicalUnit::A;

            if( bc->bc_type() == OhmicContact )
            {
              //_out << std::setw(25) << bc->ext_circuit()->current_displacement()/PhysicalUnit::A;
              _out << std::setw(25) << bc->ext_circuit()->current_electron()/PhysicalUnit::A;
              _out << std::setw(25) << bc->ext_circuit()->current_hole()/PhysicalUnit::A;
            }
            
            power += bc->ext_circuit()->Vapp()*bc->ext_circuit()->current();
            
            continue;
          }

          if( bc->has_current_flow() )
          {
            _out << std::setw(25) << bc->current()/PhysicalUnit::A;
          }

          if( bc->bc_type() == IF_Metal_Ohmic || bc->bc_type() == IF_Metal_Schottky)
          {
            _out << std::setw(25) << bc->psi()/PhysicalUnit::V;
          }

          // charge integral interface
          if( bc->bc_type() == ChargeIntegral )
          {
            _out << std::setw(25) << bc->scalar("qf")/PhysicalUnit::C;
            _out << std::setw(25) << bc->psi()/PhysicalUnit::V;
          }

          // current pass though homo interface
          if( _ddm && bc->bc_type() == HomoInterface)
          {
            _out << std::setw(25) << bc->scalar("electron_current")/PhysicalUnit::A;
            _out << std::setw(25) << bc->scalar("hole_current")/PhysicalUnit::A;
            _out << std::setw(25) << bc->scalar("displacement_current")/PhysicalUnit::A;
          }
        }
        
        _out<< std::setw(25) << power/(PhysicalUnit::V*PhysicalUnit::A);
      }

    }


    if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
    {
      PetscScalar omega = 2*3.14159265358979323846*SolverSpecify::Freq*PhysicalUnit::s;
      _out << SolverSpecify::Freq*PhysicalUnit::s << '\t';

      // search for all the bc
      const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();

      for(unsigned int n=0; n<bcs->n_bcs(); n++)
      {
        const BoundaryCondition * bc = bcs->get_bc(n);
        // skip bc which is not electrode
        if( bc->is_electrode() )
        {

          // DC potential and current
          _out << std::setw(25) << bc->ext_circuit()->potential()/PhysicalUnit::V;
          _out << std::setw(25) << bc->ext_circuit()->current()/PhysicalUnit::A;

          //record electrode potential and electrode current for AC simulation
          _out << std::setw(25) << bc->ext_circuit()->Vac()/PhysicalUnit::V;

          _out << std::setw(25) << bc->ext_circuit()->potential_ac().real()/PhysicalUnit::V;
          _out << std::setw(25) << bc->ext_circuit()->potential_ac().imag()/PhysicalUnit::V;

          _out << std::setw(25) << bc->ext_circuit()->current_ac().real()/PhysicalUnit::A;
          _out << std::setw(25) << bc->ext_circuit()->current_ac().imag()/PhysicalUnit::A;

          std::complex<PetscScalar> Y;
          Y = (bc->ext_circuit()->current_ac()/PhysicalUnit::A)/(SolverSpecify::VAC/PhysicalUnit::V);
          _out << std::setw(25) << Y.real();
          _out << std::setw(25) << Y.imag()/omega;

          continue;
        }

        if( bc->has_current_flow() )
        {
          _out << std::setw(25) << bc->current()/PhysicalUnit::A;
        }
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

//----------------------------------------------------------------------

void  GnuplotHook::_write_gnuplot_head()
{

  // prepare the file head
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    time_t          _time;

    // get simulation time
    time(&_time);

    // write file head
    _out << "# Title: Gnuplot File Created by Genius TCAD Simulation" << std::endl;
    _out << "# Date: " << ctime(&_time) << std::endl;

    switch (SolverSpecify::Type)
    {
      case SolverSpecify::DCSWEEP   :
        _out << "# Plotname: DC transfer characteristic" << std::endl; break;
      case SolverSpecify::OP   :
        _out << "# Plotname: OP" << std::endl; break;
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
        _out << '#' <<'\t' << ++n_var <<'\t' << "Time" << " [s]"<< std::endl;
        _out << '#' <<'\t' << ++n_var <<'\t' << "TimeStep" << " [s]"<< std::endl;
      }

      if( _mixA ) // mix mode
      {
        const SPICE_CKT * spice_ckt = this->get_solver().get_system().get_circuit();


        for(unsigned int n=0; n<spice_ckt->n_ckt_nodes(); n++)
        {
          std::string electrode = spice_ckt->ckt_node_name(n);

          if(spice_ckt->is_voltage_node(n))
            _out << '#' <<'\t' << ++n_var <<'\t' << "V(" + electrode + ")"  << " [V]"<< std::endl;
          else
            _out << '#' <<'\t' << ++n_var <<'\t' << "I(" + electrode + ")"  << " [A]"<< std::endl;
        }

        // record electrode IV information
        const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
        // search for all the bcs
        for(unsigned int n=0; n<bcs->n_bcs(); n++)
        {
          const BoundaryCondition * bc = bcs->get_bc(n);

          if( bc->has_current_flow() )
          {
            std::string bc_label = bc->label();
            if(!bc->electrode_label().empty())
              bc_label = bc->electrode_label() + '[' + bc_label + ']';
            _out << '#' <<'\t' << ++n_var <<'\t' << "I(" + bc_label + ")"   << " [A]"<< std::endl;
          }
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
              bc_label = bc->electrode_label() + '[' + bc_label + ']';

            _out << '#' <<'\t' << ++n_var <<'\t' << "Vapp(" + bc_label + ")"  << " [V]"<< std::endl;
            _out << '#' <<'\t' << ++n_var <<'\t' << "P(" + bc_label + ")"     << " [V]"<< std::endl;
            _out << '#' <<'\t' << ++n_var <<'\t' << "I(" + bc_label + ")"     << " [A]"<< std::endl;

            if( bc->bc_type() == OhmicContact )
            {
              //_out << '#' <<'\t' << ++n_var <<'\t' << bc_label + "_displacement_current"   << " [A]"<< std::endl;
              _out << '#' <<'\t' << ++n_var <<'\t' << "Ie(" + bc_label + ")"   << " [A]"<< std::endl;
              _out << '#' <<'\t' << ++n_var <<'\t' << "Ih(" + bc_label + ")"   << " [A]"<< std::endl;
            }

            continue;
          }

          if( bc->has_current_flow() )
          {
            std::string bc_label = bc->label();
            _out << '#' <<'\t' << ++n_var <<'\t' << "I(" + bc_label + ")"   << " [A]"<< std::endl;
          }

          if( bc->bc_type() == IF_Metal_Ohmic || bc->bc_type() == IF_Metal_Schottky)
          {
            std::string bc_label = bc->label();
            _out << '#' <<'\t' << ++n_var <<'\t' << "V(" + bc_label + ")"   << " [V]"<< std::endl;
          }

          // charge integral interface
          if( bc->bc_type() == ChargeIntegral )
          {
            std::string bc_label = bc->label();
            _out << '#' <<'\t' << ++n_var <<'\t' << "Q(" + bc_label + ")"          << " [C]"<< std::endl;
            _out << '#' <<'\t' << ++n_var <<'\t' << "V(" + bc_label + ")"  << " [V]"<< std::endl;
          }

          if( _ddm && bc->bc_type() == HomoInterface)
          {
            std::string bc_label = bc->label();
            _out << '#' <<'\t' << ++n_var <<'\t' << "Ie(" + bc_label + ")"   << " [A]"<< std::endl;
            _out << '#' <<'\t' << ++n_var <<'\t' << "Ih(" + bc_label + ")"   << " [A]"<< std::endl;
            _out << '#' <<'\t' << ++n_var <<'\t' << "Idisp(" + bc_label + ")"   << " [A]"<< std::endl;
          }
        }
        
        _out << '#' <<'\t' << ++n_var <<'\t' << "Power" << " [W]"<< std::endl;
      }

    }

    if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
    {
      unsigned int n_var = 0;

      _out << '#' <<'\t' << ++n_var <<'\t' << "frequency" << " [Hz]"<< std::endl;

      // record electrode IV information
      const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();

      std::string ac_bc_label = SolverSpecify::Electrode_ACScan[0];

      // search for all the bcs
      for(unsigned int n=0; n<bcs->n_bcs(); n++)
      {
        const BoundaryCondition * bc = bcs->get_bc(n);

        if( bc->is_electrode() )
        {
          //
          std::string bc_label = bc->label();
          if(!bc->electrode_label().empty())
            bc_label = bc->electrode_label() + '[' + bc_label + ']';
          // DC
          _out << '#' <<'\t' << ++n_var <<'\t' << "Vdc(" + bc_label + ")" << " [V]"<< std::endl;
          _out << '#' <<'\t' << ++n_var <<'\t' << "Idc(" + bc_label + ")" << " [A]"<< std::endl;

          // AC
          _out << '#' <<'\t' << ++n_var <<'\t' << "Vac(" + bc_label + ")" << " [V]"<< std::endl;

          _out << '#' <<'\t' << ++n_var <<'\t' << "Pac.real(" + bc_label + ")" << " [V]"<< std::endl;
          _out << '#' <<'\t' << ++n_var <<'\t' << "Pac.imag(" + bc_label + ")" << " [V]"<< std::endl;

          _out << '#' <<'\t' << ++n_var <<'\t' << "Iac.real(" + bc_label + ")" << " [A]"<< std::endl;
          _out << '#' <<'\t' << ++n_var <<'\t' << "Iac.imag(" + bc_label + ")" << " [A]"<< std::endl;

          _out << '#' <<'\t' << ++n_var <<'\t' << "G(" + ac_bc_label + ',' +bc_label   << ") [S]"<< std::endl;
          _out << '#' <<'\t' << ++n_var <<'\t' << "C(" + ac_bc_label + ',' +bc_label   << ") [F]"<< std::endl;

          continue;
        }

        if( bc->has_current_flow() )
        {
          std::string bc_label = bc->label();
          _out << '#' <<'\t' << ++n_var <<'\t' << "Idc(" + bc_label + ")"   << " [A]"<< std::endl;
        }
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

