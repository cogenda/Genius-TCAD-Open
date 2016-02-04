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
#include "spice_ckt.h"
#include "probe_hook.h"
#include "parallel.h"


/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
ProbeHook::ProbeHook(SolverBase & solver, const std::string & name, void * param)
    : Hook(solver, name), _probe_file(SolverSpecify::out_prefix + ".probe")
{
  _p_solver = & solver;
  _p_fvm_node   = 0;
  _min_loc = invalid_uint;

  const std::vector<Parser::Parameter> & parm_list = *((std::vector<Parser::Parameter> *)param);
  for(std::vector<Parser::Parameter>::const_iterator parm_it = parm_list.begin();
      parm_it != parm_list.end(); parm_it++)
  {
    if(parm_it->name() == "x" && parm_it->type() == Parser::REAL)
      _pp(0)=parm_it->get_real() * PhysicalUnit::um;
    if(parm_it->name() == "y" && parm_it->type() == Parser::REAL)
      _pp(1)=parm_it->get_real() * PhysicalUnit::um;
    if(parm_it->name() == "z" && parm_it->type() == Parser::REAL)
      _pp(2)=parm_it->get_real() * PhysicalUnit::um;
    if(parm_it->name() == "file" && parm_it->type() == Parser::STRING)
      _probe_file = parm_it->get_string();
    if(parm_it->name() == "region" && parm_it->type() == Parser::STRING)
      _region = parm_it->get_string();
    if(parm_it->name() == "material" && parm_it->type() == Parser::STRING)
      _material = parm_it->get_string();
  }

  const SimulationSystem & system = _p_solver->get_system();
  if(!_region.empty() && !system.has_region(_region)) _region.clear();
  if(!_material.empty() && !system.has_region_with_material(_material)) _material.clear();

  if ( !Genius::processor_id() )
    _out.open(_probe_file.c_str());

}


/*----------------------------------------------------------------------
 * destructor, close file
 */
ProbeHook::~ProbeHook()
{
  if ( !Genius::processor_id() )
    _out.close();
}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void ProbeHook::on_init()
{

  double min_dis = 1e100;
  for( unsigned int r=0; r<_p_solver->get_system().n_regions(); r++)
  {
    const SimulationRegion * region = _p_solver->get_system().region(r);
    if(!_region.empty() && region->name() != _region) continue;
    if(!_material.empty() && region->material() != _material) continue;

    SimulationRegion::const_processor_node_iterator node_it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator node_it_end = region->on_processor_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      const FVM_Node * fvm_node = *node_it;
      const Node * node = fvm_node->root_node();

      double dis = ((*node)-_pp).size();
      if(dis<min_dis)
      {
        min_dis = dis;
        _p_fvm_node = fvm_node;
      }
    }
  }

  // after this call, the _min_loc contains processor_id with minimal min_dis
  Parallel::min_loc(min_dis, _min_loc);

  double x=0,y=0,z=0;
  std::string region;
  std::vector<std::string> var_name;
  int n_var;

  if(_p_fvm_node)
  {
    x = (*_p_fvm_node->root_node())(0);
    y = (*_p_fvm_node->root_node())(1);
    z = (*_p_fvm_node->root_node())(2);
    region = _p_solver->get_system().region(_p_fvm_node->subdomain_id())->label();
  }

  Parallel::broadcast(x, _min_loc);
  Parallel::broadcast(y, _min_loc);
  Parallel::broadcast(z, _min_loc);
  Parallel::broadcast(region, _min_loc);

  if (Genius::processor_id() == _min_loc)
  {
    const FVM_NodeData * node_data = _p_fvm_node->node_data();
    switch (node_data->type())
    {
      case FVM_NodeData::SemiconductorData:
        if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
        {
          n_var = 6;
          var_name.push_back("psi_real [V]");
          var_name.push_back("psi_imag [V]");

          var_name.push_back("n_real [cm^-3]");
          var_name.push_back("n_imag [cm^-3]");

          var_name.push_back("p_real [cm^-3]");
          var_name.push_back("p_imag [cm^-3]");
        }
        else
        {
          n_var = 5;
          var_name.push_back("psi [V]");
          var_name.push_back("n [cm^-3]");
          var_name.push_back("p [cm^-3]");
          var_name.push_back("Optical G [cm^-3/s]");
          var_name.push_back("Partical G [cm^-3/s]");
        }
        break;

      case FVM_NodeData::InsulatorData:
      case FVM_NodeData::ConductorData:
      case FVM_NodeData::ResistanceData:
        if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
        {
          n_var = 2;
          var_name.push_back("psi_real [V]");
          var_name.push_back("psi_imag [V]");
        }
        else
        {
          n_var = 1;
          var_name.push_back("psi [V]");
        }
        break;
      default:
        n_var = 0;
    }

  }

  Parallel::broadcast(n_var, _min_loc);
  var_name.resize(n_var);
  for(int i=0; i<n_var; i++)
    Parallel::broadcast(var_name[i], _min_loc);

  if ( !Genius::processor_id() )
  {
    time_t          _time;
    time(&_time);

    _out << "# Title: Gnuplot File Created by Genius TCAD Simulation" << std::endl;
    _out << "# Date: " << ctime(&_time) << std::endl;
    _out << "# Plotname: Probe"  << std::endl;
    _out << "# Node location (um): ";
    _out << "x=" << x/PhysicalUnit::um << "\ty=" << y/PhysicalUnit::um << "\tz=" << z/PhysicalUnit::um << std::endl;
    _out << "# Node in region: " << region << std::endl;
    _out << "# Variables: " << std::endl;
    int cCnt=0;

    // DC Sweep
    if( SolverSpecify::Type==SolverSpecify::DCSWEEP || SolverSpecify::Type==SolverSpecify::TRACE)
    {
      if ( SolverSpecify::Electrode_VScan.size() )
        for(unsigned int n=0; n<SolverSpecify::Electrode_VScan.size(); ++n)
          _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << SolverSpecify::Electrode_VScan[n]  << " [V]"<< std::endl;

      if ( SolverSpecify::Electrode_IScan.size() )
        for(unsigned int n=0; n<SolverSpecify::Electrode_IScan.size(); ++n)
          _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << SolverSpecify::Electrode_IScan[n]  << " [A]"<< std::endl;
    }

    // Transient
    if(SolverSpecify::Type==SolverSpecify::TRANSIENT)
    {
      _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << "Time [s]"  << std::endl;
    }

    // AC
    if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
    {
      _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << "Frequency [Hz]"  << std::endl;
    }

    for(int i=0; i<n_var; i++)
    {
      _out << '#' << std::setw(10) << ++cCnt << std::setw(30) << var_name[i]  << std::endl;
    }
    _out << std::endl;
  }

}

/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void ProbeHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void ProbeHook::post_solve()
{
  int n_var;
  std::vector<double> var;
  if (Genius::processor_id() == _min_loc)
  {
    const FVM_NodeData * node_data = _p_fvm_node->node_data();

    switch (node_data->type())
    {
      case FVM_NodeData::SemiconductorData:
        if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
        {
          n_var = 6;
          var.push_back(node_data->psi_ac().real()/PhysicalUnit::V);
          var.push_back(node_data->psi_ac().imag()/PhysicalUnit::V);

          var.push_back(node_data->n_ac().real()/std::pow(PhysicalUnit::cm, -3));
          var.push_back(node_data->n_ac().imag()/std::pow(PhysicalUnit::cm, -3));

          var.push_back(node_data->p_ac().real()/std::pow(PhysicalUnit::cm, -3));
          var.push_back(node_data->p_ac().imag()/std::pow(PhysicalUnit::cm, -3));
        }
        else
        {
          n_var = 5;
          var.push_back(node_data->psi()/PhysicalUnit::V);
          var.push_back(node_data->n()/std::pow(PhysicalUnit::cm, -3));
          var.push_back(node_data->p()/std::pow(PhysicalUnit::cm, -3));

          var.push_back(node_data->PatG()/(std::pow(PhysicalUnit::cm, -3)/PhysicalUnit::s));
          var.push_back(node_data->OptG()/(std::pow(PhysicalUnit::cm, -3)/PhysicalUnit::s));
        }
        break;

      case FVM_NodeData::InsulatorData:
      case FVM_NodeData::ConductorData:
      case FVM_NodeData::ResistanceData:
        if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
        {
          n_var = 2;
          var.push_back(node_data->psi_ac().real()/PhysicalUnit::V);
          var.push_back(node_data->psi_ac().imag()/PhysicalUnit::V);
        }
        else
        {
          n_var = 1;
          var.push_back(node_data->psi()/PhysicalUnit::V);
        }
        break;
      default:
        n_var = 0;
    }
  }


  Parallel::broadcast(n_var,_min_loc);
  if(n_var)
    Parallel::broadcast(var,_min_loc);

  SPICE_CKT * ckt = _p_solver->get_system().get_circuit();

  if ( !Genius::processor_id() )
  {
    // set the float number precision
    _out.precision(6);

    // set output width and format
    _out<< std::scientific << std::right;

    _out<<' ';

    // DC Sweep
    if( SolverSpecify::Type==SolverSpecify::DCSWEEP || SolverSpecify::Type==SolverSpecify::TRACE )
    {
      if ( SolverSpecify::Electrode_VScan.size() )
      {
        for(unsigned int n=0; n<SolverSpecify::Electrode_VScan.size(); ++n)
        {
          _out << std::setw(20) << SolverSpecify::Electrode_VScan_Voltage/PhysicalUnit::V;
        }
      }
      
      if ( SolverSpecify::Electrode_IScan.size() )
      {
        for(unsigned int n=0; n<SolverSpecify::Electrode_IScan.size(); ++n)
        {
          _out << std::setw(20) << SolverSpecify::Electrode_IScan_Current/PhysicalUnit::A;
        }
      }
        
    }

    // Transient
    if(SolverSpecify::Type==SolverSpecify::TRANSIENT)
    {
      _out << std::setw(15) << SolverSpecify::clock/PhysicalUnit::s;
    }

    // AC
    if( SolverSpecify::Type==SolverSpecify::ACSWEEP)
    {
      PetscScalar omega = 2*3.14159265358979323846*SolverSpecify::Freq*PhysicalUnit::s;
      _out << std::setw(15) << SolverSpecify::Freq*PhysicalUnit::s;
    }

    for(unsigned int i=0; i<var.size(); i++)
      _out << std::setw(15) << var[i];

    _out << std::endl;
  }

}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void ProbeHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void ProbeHook::on_close()
{}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook (SolverBase & solver, const std::string & name, void * fun_data)
  {
    return new ProbeHook(solver, name, fun_data );
  }

}

#endif

