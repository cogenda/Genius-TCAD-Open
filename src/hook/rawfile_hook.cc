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

//  $Id: rawfile_hook.cc,v 1.9 2008/07/09 12:25:19 gdiso Exp $

#include <string>
#include <cstdlib>
#include <iomanip>

#include "solver_base.h"
#include "rawfile_hook.h"
#include "spice_ckt.h"

/*----------------------------------------------------------------------
 * constructor, open the rawfile for writing
 */
RawFileHook::RawFileHook(SolverBase & solver, const std::string & name, void * file)
    : Hook(solver, name), _input_file((const char *)file), _raw_file(SolverSpecify::out_prefix + ".raw"), _mixA(false), _n_values(0)
{
  if ( !Genius::processor_id() )
    _out.open(_raw_file.c_str());

  SolverSpecify::SolverType solver_type = this->get_solver().solver_type();

  // if we are called by mixA solver?
  if(solver_type == SolverSpecify::DDML1MIXA ||
      solver_type == SolverSpecify::DDML2MIXA ||
      solver_type == SolverSpecify::EBML3MIXA
    )
    _mixA = true;

}


/*----------------------------------------------------------------------
 * destructor, close the raw file
 */
RawFileHook::~RawFileHook()
{

}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void RawFileHook::on_init()
{
  // prepare the buffer info
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    // get simulation time
    time(&_time);

    if( SolverSpecify::Type==SolverSpecify::DCSWEEP ||
        SolverSpecify::Type==SolverSpecify::TRACE   ||
        SolverSpecify::Type==SolverSpecify::TRANSIENT )
    {
      // if transient simulation, we need to record time
      if ( SolverSpecify::Type == SolverSpecify::TRANSIENT )
      {
        _variables.push_back( std::pair<std::string, std::string>("time", "time") );
        _variables.push_back( std::pair<std::string, std::string>("time_step", "time") );
      }

      if( !_mixA )
      {
        // record electrode IV information
        const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
        // search for all the bcs
        for(unsigned int n=0; n<bcs->n_bcs(); n++)
        {
          const BoundaryCondition * bc = bcs->get_bc(n);
          // electrode
          if( bc->is_electrode() )
          {
            std::string bc_label = bc->label();
            if(!bc->electrode_label().empty())
              bc_label = bc->electrode_label();

            _variables.push_back( std::pair<std::string, std::string>(bc_label + "_Vapp", "voltage") );
            _variables.push_back( std::pair<std::string, std::string>(bc_label + "_potential", "voltage")  );
            _variables.push_back( std::pair<std::string, std::string>(bc_label + "_current", "current") );
            continue;
          }

          if( bc->has_current_flow() )
          {
            std::string bc_label = bc->label();
            _variables.push_back( std::pair<std::string, std::string>(bc_label + "_current", "current") );
          }


          if( bc->bc_type() == IF_Metal_Ohmic || bc->bc_type() == IF_Metal_Schottky)
          {
            std::string bc_label = bc->label();
            _variables.push_back( std::pair<std::string, std::string>(bc_label + "_average_potential", "voltage") );
          }

          // charge integral interface
          if( bc->bc_type() == ChargeIntegral )
          {
            std::string bc_label = bc->label();
            _variables.push_back( std::pair<std::string, std::string>(bc_label + "_Q", "charge") );
            _variables.push_back( std::pair<std::string, std::string>(bc_label + "_potential", "voltage")  );
            continue;
          }
        }
      }
      else
      {
        const SPICE_CKT * spice_ckt = this->get_solver().get_system().get_circuit();
        for(unsigned int n=0; n<spice_ckt->n_ckt_nodes(); n++)
        {
          if(spice_ckt->is_voltage_node(n))
            _variables.push_back( std::pair<std::string, std::string>(spice_ckt->ckt_node_name(n), "voltage") );
          else
            _variables.push_back( std::pair<std::string, std::string>(spice_ckt->ckt_node_name(n), "current") );
        }
      }

    }

    if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
    {

      _variables.push_back( std::pair<std::string, std::string>("frequency", "Hz") );

      // record electrode IV information
      const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
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
        //
        _variables.push_back( std::pair<std::string, std::string>(bc_label + "_potential_magnitude", "voltage") );
        _variables.push_back( std::pair<std::string, std::string>(bc_label + "_potential_angle",     "")  );
        _variables.push_back( std::pair<std::string, std::string>(bc_label + "_current_magnitude",   "current") );
        _variables.push_back( std::pair<std::string, std::string>(bc_label + "_current_angle",       "") );
      }
    }

    _values.resize( _variables.size() );

  }

}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void RawFileHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void RawFileHook::post_solve()
{

  // save electrode IV
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    // variable counter
    unsigned int i=0;

    if( SolverSpecify::Type==SolverSpecify::DCSWEEP ||
        SolverSpecify::Type==SolverSpecify::TRACE   ||
        SolverSpecify::Type==SolverSpecify::TRANSIENT )
    {
      // if transient simulation, we need to record time
      if (SolverSpecify::Type == SolverSpecify::TRANSIENT)
      {
        _values[i++].push_back( SolverSpecify::clock/PhysicalUnit::s );
        _values[i++].push_back( SolverSpecify::dt/PhysicalUnit::s );
      }

      if( !_mixA )
      {
        // search for all the bc
        const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
        for(unsigned int n=0; n<bcs->n_bcs(); n++)
        {
          const BoundaryCondition * bc = bcs->get_bc(n);
          // electrode
          if( bc->is_electrode() )
          {
            _values[i++].push_back( bc->ext_circuit()->Vapp()/PhysicalUnit::V );
            _values[i++].push_back( bc->ext_circuit()->potential()/PhysicalUnit::V );
            _values[i++].push_back( bc->ext_circuit()->current()/PhysicalUnit::A );
            continue;
          }

          if( bc->has_current_flow() )
          {
            _values[i++].push_back( bc->current()/PhysicalUnit::A );
          }

          if( bc->bc_type() == IF_Metal_Ohmic || bc->bc_type() == IF_Metal_Schottky)
          {
            _values[i++].push_back( bc->psi()/PhysicalUnit::V );
          }

          // charge integral interface
          if( bc->bc_type() == ChargeIntegral )
          {
            _values[i++].push_back( bc->scalar("qf")/PhysicalUnit::C );
            _values[i++].push_back( bc->psi()/PhysicalUnit::V );
          }
        }
      }
      else
      {
        const SPICE_CKT * spice_ckt = this->get_solver().get_system().get_circuit();
        for(unsigned int n=0; n<spice_ckt->n_ckt_nodes(); n++)
        {
          if(spice_ckt->is_voltage_node(n))
            _values[i++].push_back(spice_ckt->get_solution(n));
          else
            _values[i++].push_back(spice_ckt->get_solution(n));
        }
      }
    }

    if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
    {
      //record frequency
      _values[i++].push_back( SolverSpecify::Freq*PhysicalUnit::s );

      // record electrode IV information
      const BoundaryConditionCollector * bcs = this->get_solver().get_system().get_bcs();
      // search for all the bcs
      for(unsigned int n=0; n<bcs->n_bcs(); n++)
      {
        const BoundaryCondition * bc = bcs->get_bc(n);
        // skip bc which is not electrode
        if( !bc->is_electrode() ) continue;
        //
        _values[i++].push_back( std::abs(bc->ext_circuit()->potential_ac())/PhysicalUnit::V );
        _values[i++].push_back( std::arg(bc->ext_circuit()->potential_ac()) );
        _values[i++].push_back( std::abs(bc->ext_circuit()->current_ac())/PhysicalUnit::A );
        _values[i++].push_back( std::arg(bc->ext_circuit()->current_ac()) );
      }
    }

    if(i)  _n_values++;
  }

}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void RawFileHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void RawFileHook::on_close()
{
  // write spice raw file
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    // write raw file head
    _out << "Title: SPICE Raw File Created by Genius TCAD Simulation" << std::endl;
    _out << "Date: " << ctime(&_time) << std::endl;

    switch (SolverSpecify::Type)
    {
        case SolverSpecify::DCSWEEP :
          _out << "Plotname: DC transfer characteristic" << std::endl; break;
        case SolverSpecify::TRACE     :
          _out << "Plotname: DC curve trace" << std::endl; break;
        case SolverSpecify::TRANSIENT :
          _out << "Plotname: Transient Analysis" << std::endl; break;
        case SolverSpecify::ACSWEEP   :
          _out << "Plotname: AC small signal Analysis" << std::endl; break;
        default: break;
    }

    _out <<  "Flags: real" << std::endl;

    _out <<  "No. Variables: " << _variables.size()    << std::endl;
    _out <<  "No. Points: "    << _n_values  << '\n'   << std::endl;

    // write variables
    _out << "Variables:"<<std::endl;
    for(unsigned int n=0; n<_variables.size(); n++)
    {
      _out << '\t' << n << '\t' << _variables[n].first << '\t' << _variables[n].second << std::endl;
    }

    _out << std::endl;

    // write values
    _out << "Values:"<<std::endl;

    _out << std::setprecision(15) << std::scientific << std::right;

    for(unsigned int i=0; i<_n_values; i++)
    {
      _out << " " << i;
      for(unsigned int n=0; n<_variables.size(); n++)
        _out  << '\t' << std::setw(25) << _values[n][i] << std::endl;
    }

    _out.close();
  }

}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook (SolverBase & solver, const std::string & name, void * fun_data)
  {
    return new RawFileHook(solver, name, fun_data );
  }

}

#endif

