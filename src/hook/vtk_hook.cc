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

#include "solver_base.h"
#include "vtk_hook.h"
#include "spice_ckt.h"
#include "MXMLUtil.h"


/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
VTKHook::VTKHook ( SolverBase & solver, const std::string & name, void * )
    : Hook ( solver, name ), _vtk_prefix ( SolverSpecify::out_prefix ), _mixA ( false )
{
  this->count  =0;
  this->_t_step=0;
  this->_v_step=0;
  this->_i_step=0;
  this->_t_last=0;
  this->_v_last=0;
  this->_i_last=0;
  const SimulationSystem &system = get_solver().get_system();

  std::ostringstream vtk_filename;
  vtk_filename << _vtk_prefix << ( this->count++ ) << ".vtu";
  system.export_vtk ( vtk_filename.str(), false );

  SolverSpecify::SolverType solver_type = this->get_solver().solver_type();
  // if we are called by mixA solver?
  if ( solver_type == SolverSpecify::DDML1MIXA ||
       solver_type == SolverSpecify::DDML2MIXA ||
       solver_type == SolverSpecify::EBML3MIXA
     )
    _mixA = true;
}


/*----------------------------------------------------------------------
 * destructor, close file
 */
VTKHook::~VTKHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void VTKHook::on_init()
{
  std::vector<Parser::Parameter> parm_list = SolverSpecify::Hook_Parameters["vtk"];
  for ( std::vector<Parser::Parameter>::iterator parm_it = parm_list.begin();
        parm_it != parm_list.end(); parm_it++ )
  {
    if ( parm_it->name() == "tstep" && parm_it->type() == Parser::REAL )
      _t_step=parm_it->get_real() * PhysicalUnit::s;
    if ( parm_it->name() == "vstep" && parm_it->type() == Parser::REAL )
      _v_step=parm_it->get_real() * PhysicalUnit::V;
    if ( parm_it->name() == "istep" && parm_it->type() == Parser::REAL )
      _i_step=parm_it->get_real() * PhysicalUnit::A;
  }

  time_sequence.clear();

}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void VTKHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void VTKHook::post_solve()
{

  std::ostringstream vtk_filename;
  if ( SolverSpecify::Type==SolverSpecify::DCSWEEP && SolverSpecify::Electrode_VScan.size() )
  {
    double Vscan = 0;

    // DDM solver
    if ( !_mixA )
    {
      const BoundaryConditionCollector * bcs = _solver.get_system().get_bcs();
      const BoundaryCondition * bc = bcs->get_bc ( SolverSpecify::Electrode_VScan[0] );
      Vscan = bc->ext_circuit()->Vapp();
    }
    // MIXA solver
    else
    {
      SPICE_CKT * spice_ckt = _solver.get_system().get_circuit();
      Vscan = spice_ckt->get_voltage_from ( SolverSpecify::Electrode_VScan[0] );
    }


    if ( std::fabs ( Vscan - this->_v_last ) >= this->_v_step )
    {
      const SimulationSystem &system = get_solver().get_system();

      vtk_filename << _vtk_prefix << ( this->count++ ) << ".vtu";
      system.export_vtk ( vtk_filename.str(),false );

      time_sequence.push_back ( std::make_pair ( Vscan/PhysicalUnit::V, vtk_filename.str() ) );
      _v_last = Vscan;
    }
  }

  if ( SolverSpecify::Type==SolverSpecify::DCSWEEP && SolverSpecify::Electrode_IScan.size() )
  {
    // DDM solver only
    assert ( !_mixA );

    const BoundaryConditionCollector * bcs = _solver.get_system().get_bcs();
    const BoundaryCondition * bc = bcs->get_bc ( SolverSpecify::Electrode_IScan[0] );
    double Iscan = bc->ext_circuit()->Iapp();

    if ( std::fabs ( Iscan - this->_i_last ) >= this->_i_step )
    {
      const SimulationSystem &system = get_solver().get_system();

      vtk_filename << _vtk_prefix << ( this->count++ ) << ".vtu";
      system.export_vtk ( vtk_filename.str(),false );

      time_sequence.push_back ( std::make_pair ( Iscan/PhysicalUnit::A, vtk_filename.str() ) );
      _i_last = Iscan;
    }
  }


  if ( SolverSpecify::Type==SolverSpecify::TRANSIENT )
  {
    if ( SolverSpecify::clock - this->_t_last >= this->_t_step )
    {
      const SimulationSystem &system = get_solver().get_system();

      vtk_filename << _vtk_prefix << ( this->count++ ) << ".vtu";
      system.export_vtk ( vtk_filename.str(), false );

      time_sequence.push_back ( std::make_pair ( SolverSpecify::clock/PhysicalUnit::ps, vtk_filename.str() ) );
      _t_last = SolverSpecify::clock;

    }
  }

  ////----
  if ( !vtk_filename.str().empty() )
  {
    mxml_node_t *eSolution = get_solver().current_dom_solution_elem();
    if ( eSolution )
    {
      mxml_node_t *eOutput  = mxmlFindElement ( eSolution, eSolution, "output", NULL, NULL, MXML_DESCEND_FIRST );
      mxml_node_t *eVtk = mxmlNewElement ( eOutput, "vtk" );
      mxml_node_t *eFile    = mxmlNewElement ( eVtk, "file" );
      mxmlAdd ( eFile, MXML_ADD_AFTER, NULL, MXMLQVariant::makeQVString ( vtk_filename.str() ) );
    }
  }
}



/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void VTKHook::post_iteration()
{}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void VTKHook::on_close()
{
  if ( time_sequence.size() ==0 ) return;

  if ( !Genius::processor_id() )
  {
    std::string pvd_filename;
    pvd_filename = _vtk_prefix + ".pvd";

    std::ofstream   out ( pvd_filename.c_str() );

    out<< "<?xml version=\"1.0\"?>" <<std::endl;
    out<< "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"<<std::endl;
    out<< "<Collection>" <<std::endl;

    for ( unsigned int n=0; n<time_sequence.size(); ++n )
      out<<"  <DataSet timestep=\""<<time_sequence[n].first<<"\" group=\"\" part=\"0\" file=\""<<time_sequence[n].second<<"\"/>"<<std::endl;

    out<<"  </Collection>"<<std::endl;
    out<<"</VTKFile>"<<std::endl;

    out.close();
  }
}


#ifndef CYGWIN

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new VTKHook ( solver, name, fun_data );
  }

}

#endif

