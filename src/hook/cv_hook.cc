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
#include <numeric>

#include "solver_base.h"
#include "cv_hook.h"
#include "parallel.h"


/*----------------------------------------------------------------------
 * constructor, open the rawfile for writing
 */
CVHook::CVHook(SolverBase & solver, const std::string & name, void * file)
    : Hook(solver, name), _input_file((const char *)file), _raw_file(SolverSpecify::out_prefix + ".cv"), _out(_raw_file.c_str()), _n_values(0)
{}


/*----------------------------------------------------------------------
 * destructor, close the raw file
 */
CVHook::~CVHook()
{ _out.close(); }


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void CVHook::on_init()
{
  // prepare the buffer info

  // get simulation time
  time(&_time);

  if (SolverSpecify::Type == SolverSpecify::DCSWEEP)
  {
    BoundaryConditionCollector * bcs = _solver.get_system().get_bcs();
    for(unsigned int n=0; n<bcs->n_bcs(); n++)
    {
      BoundaryCondition * bc = bcs->get_bc(n);

      if(bc->label() == SolverSpecify::Electrode_VScan[0])
      {
        _sweep_electrode = bc->label();
      }
      if( bc->bc_type() == GateContact )
      {
        _gate_electrodes.push_back (bc->label());
      }
    }
    _gate_charge.resize( _gate_electrodes.size() );
  }

}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void CVHook::pre_solve()
{
}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void CVHook::post_solve()
{

  // available in DCSWEEP mode only
  if (SolverSpecify::Type == SolverSpecify::DCSWEEP)
  {
    BoundaryConditionCollector * bcs = _solver.get_system().get_bcs();
    unsigned int elec_count=0;
    for(unsigned int n=0; n<bcs->n_bcs(); n++)
    {
      BoundaryCondition * bc = bcs->get_bc(n);

      // Save sweep voltage
      if(bc->label() == SolverSpecify::Electrode_VScan[0])
      {
        _vsweep.push_back(bc->ext_circuit()->Vapp()/PhysicalUnit::V);
      }

      // skip bc which is not gate
      if( bc->bc_type() != GateContact ) continue;

      std::vector<double> flux_buffer;

      // device in Z dimension. for 3D mesh, z_width() should return 1.0.
      PetscScalar flux_scale = bc->z_width();

      BoundaryCondition::const_node_iterator node_it = bc->nodes_begin();
      BoundaryCondition::const_node_iterator end_it  = bc->nodes_end();

      for(; node_it != end_it; ++node_it)
      {
        // skip node not belonging to this processor
        if ( (*node_it)->processor_id() != Genius::processor_id() ) continue;

        // iterate over all fvm_node associated to *node_it
        BoundaryCondition::region_node_iterator rnode_it = bc->region_node_begin(*node_it);
        BoundaryCondition::region_node_iterator end_rnode_it = bc->region_node_end(*node_it);

        std::vector<SimulationRegion *> regions;
        std::vector<FVM_Node *> fvm_nodes;

        //
        for (unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it)
        {
          regions.push_back( (*rnode_it).second.first );
          fvm_nodes.push_back( (*rnode_it).second.second );

          switch ( regions[i]->type() )
          {
            case InsulatorRegion:
              {
                FVM_NodeData * node_data = fvm_nodes[i]->node_data();

                FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
                for( ; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it )
                {
                  FVM_Node *nb_node = (*nb_it).second;
                  FVM_NodeData * nb_node_data = nb_node->node_data();

                  PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();
                  PetscScalar cv_area = fvm_nodes[i]->cv_surface_area(nb_node->root_node());
                  PetscScalar EFlux = (node_data->psi()-nb_node_data->psi())/distance;

                  flux_buffer.push_back(cv_area*node_data->eps()*EFlux*flux_scale);
                }
                break;
              }
            case ConductorRegion:
              {
                break;
              }
            default: genius_error(); //we should never reach here
          }
        }
      }
      Parallel::allgather(flux_buffer);
      PetscScalar charge = std::accumulate(flux_buffer.begin(), flux_buffer.end(), 0.0);
      _gate_charge[elec_count++].push_back(charge/PhysicalUnit::C);
    }
    if (elec_count) _n_values++;
  }
}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void CVHook::post_iteration()
{
}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void CVHook::on_close()
{
  // write spice raw file
  // only root processor do this command
  if ( !Genius::processor_id() )
  {
    if (SolverSpecify::Type == SolverSpecify::DCSWEEP)
    {
      // write raw file head
      _out << "# Title: CV Curve Created by Genius TCAD Simulation" << std::endl;
      _out << "# Date: " << ctime(&_time) << std::endl;

      //_out << std::setprecision(15) << std::scientific << std::right;
      _out << "#\t1\t"<< _sweep_electrode << " [V]" << std::endl;
      for(unsigned int n=0; n<_gate_charge.size(); n++)
        _out << "#\t" << n+2 << "\t"<< _gate_electrodes[n] << " [F]" << std::endl;

      _out << std::endl;

      {
        unsigned int i=0;
        _out << _vsweep[i];
        for (unsigned int n=0; n<_gate_charge.size(); n++)
          _out << "\t" << (_gate_charge[n][i+1]-_gate_charge[n][i])/(_vsweep[i+1]-_vsweep[i]);
        _out << std::endl;
      }

      for(unsigned int i=1; i<_n_values-1; i++)
      {
        _out << _vsweep[i];
        for(unsigned int n=0; n<_gate_charge.size(); n++)
        {
          double hl = _vsweep[i-1]-_vsweep[i];
          double hr = _vsweep[i+1]-_vsweep[i];
          double c1 = hr/hl/(hr-hl);
          double c2 = -(hr+hl)/hl/hr;
          double c3 = -hl/hr/(hr-hl);

          _out  << '\t' <<  c1*_gate_charge[n][i-1]
                          + c2*_gate_charge[n][i]
                          + c3*_gate_charge[n][i+1];
        }
        _out << std::endl;
      }

      {
        unsigned int i= _n_values-1;
        _out << _vsweep[i];
        for (unsigned int n=0; n<_gate_charge.size(); n++)
          _out << "\t" << (_gate_charge[n][i]-_gate_charge[n][i-1])/(_vsweep[i]-_vsweep[i-1]);
        _out << std::endl;
      }

    }
  }
}


#ifndef CYGWIN

// dll interface
extern "C"
{
  Hook* get_hook (SolverBase & solver, const std::string & name, void * fun_data)
  {
    return new CVHook(solver, name, fun_data );
  }

}

#endif


