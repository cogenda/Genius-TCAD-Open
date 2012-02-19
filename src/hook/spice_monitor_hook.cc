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

#include <iomanip>

#include "solver_base.h"
#include "spice_monitor_hook.h"
#include "spice_ckt.h"
#include "mixA_solver.h"

/*----------------------------------------------------------------------
 * constructor, open the file for writing
 */
SpiceMonitorHook::SpiceMonitorHook ( SolverBase & solver, const std::string & name, void * param )
    : Hook ( solver, name )
{
  this->_mix_solver = false;
  this->solution_count=0;
  this->iteration_count=0;
}


/*----------------------------------------------------------------------
 * destructor, close file
 */
SpiceMonitorHook::~SpiceMonitorHook()
{}


/*----------------------------------------------------------------------
 *   This is executed before the initialization of the solver
 */
void SpiceMonitorHook::on_init()
{
  this->_mix_solver = false;

  if( _solver.solver_type() == SolverSpecify::DDML1MIXA ||
      _solver.solver_type() == SolverSpecify::DDML2MIXA ||
      _solver.solver_type() == SolverSpecify::EBML3MIXA
    )
  {
    _mix_solver = true;
  }


}



/*----------------------------------------------------------------------
 *   This is executed previously to each solution step.
 */
void SpiceMonitorHook::pre_solve()
{}



/*----------------------------------------------------------------------
 *  This is executed after each solution step.
 */
void SpiceMonitorHook::post_solve()
{
  this->solution_count++;
  this->iteration_count=0;
}


/*----------------------------------------------------------------------
 *  This is executed before each (nonlinear) iteration
 */
void SpiceMonitorHook::pre_iteration()
{
#if 0
  const SimulationSystem &system = _solver.get_system();
  const SPICE_CKT * circuit = system.get_circuit();

  if(Genius::is_last_processor())
  {
    for(unsigned int n=0; n<circuit->n_ckt_nodes(); ++n)
    {
      int global_row;
      std::vector<int> global_col;
      std::vector<double> values;
      circuit->ckt_matrix_row(n, global_row, global_col, values);

      for(unsigned int c=0; c<global_col.size(); ++c)
        std::cout<< '(' << global_row <<','<< global_col[c] << ')' << '=' << values[c] << ' ';
      std::cout<<std::endl;
    }
  }
#endif
}


/*----------------------------------------------------------------------
 *  This is executed after each (nonlinear) iteration
 */
void SpiceMonitorHook::post_iteration(void * _f, void * _x, void * _dx, void * _w, bool & , bool &)
{
  const SimulationSystem &system = _solver.get_system();
  const SPICE_CKT * circuit = system.get_circuit();

  Vec f  = Vec(_f);  // previous function residual
  Vec x  = Vec(_x);  // previous iterate value
  Vec dx = Vec(_dx); // new search direction and length
  Vec w  = Vec(_w);  // current candidate iterate

  PetscScalar    *ff;
  PetscScalar    *xx;
  PetscScalar    *dxx;
  PetscScalar    *ww;

  VecGetArray(f, &ff);
  VecGetArray(x, &xx);
  VecGetArray(dx, &dxx);
  VecGetArray(w, &ww);


  if(Genius::is_last_processor())
  {

    std::cout << std::left << std::setw(30) << "  Node" << std::setw(15) << "Variable" << std::setw(15) << "Update"  << std::endl;
    std::cout << std::left << std::setw(30) << "  ----" << std::setw(15) << "--------" << std::setw(15) << "------"  << std::endl;

    for(unsigned int n=0; n<circuit->n_ckt_nodes(); ++n)
    {
      unsigned int array_offset_x = circuit->array_offset_x(n);

      std::cout << "  " << std::left << std::setw(30) << circuit->ckt_node_name(n)
          << std::setw(15) << xx[array_offset_x]
          << std::setw(15) <<-dxx[array_offset_x]
      << std::endl;
    }
    std::cout<<std::endl;
  }

  VecRestoreArray(f, &ff);
  VecRestoreArray(x, &xx);
  VecRestoreArray(dx, &dxx);
  VecRestoreArray(w, &ww);


  const MixASolverBase & mix_solver = dynamic_cast<MixASolverBase &>(_solver);

  std::string matrix_prefix = SolverSpecify::out_prefix+".monitor";
  std::ostringstream matrix_filename;
  matrix_filename << matrix_prefix << '.' << this->solution_count<< '.' << this->iteration_count << ".mat";

  mix_solver.dump_spice_matrix_petsc( matrix_filename.str());


  this->iteration_count++;

}



/*----------------------------------------------------------------------
 * This is executed after the finalization of the solver
 */
void SpiceMonitorHook::on_close()
{}


#ifdef DLLHOOK

// dll interface
extern "C"
{
  Hook* get_hook ( SolverBase & solver, const std::string & name, void * fun_data )
  {
    return new SpiceMonitorHook ( solver, name, fun_data );
  }
}

#endif

