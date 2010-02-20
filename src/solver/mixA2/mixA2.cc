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



#include "mixA2/mixA2.h"
#include "spice_ckt.h"
#include "parallel.h"



using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::cm;
using PhysicalUnit::K;


/*------------------------------------------------------------------
 * create nonlinear solver contex and adjust some parameters
 */
int MixA2Solver::create_solver()
{

  MESSAGE<< '\n' << "Advanced Mixed Simulation with DDM Level 2 init..." << std::endl;
  RECORD();

  return MixASolverBase::create_solver();
}




/*------------------------------------------------------------------
 * set initial value to solution vector and scaling vector
 */
int MixA2Solver::pre_solve_process(bool load_solution)
{
  if(load_solution)
  {
    // for all the regions
    // NOTE we use DDM2_Fill_Value here!
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DDM2_Fill_Value(x, L);
    }

    spice_fill_value(x, L);

    VecAssemblyBegin(x);
    VecAssemblyBegin(L);

    VecAssemblyEnd(x);
    VecAssemblyEnd(L);
  }

  return MixASolverBase::pre_solve_process(load_solution);
}





/*------------------------------------------------------------------
 * the main solve routine. which is under the control of ngspice
 */
int MixA2Solver::solve()
{

  START_LOG("MixA2Solver_SNES()", "MixA2Solver");

  switch( SolverSpecify::Type )
  {
  case SolverSpecify::OP :
    solve_dcop(); break;

  case SolverSpecify::DCSWEEP :
    solve_dcsweep(); break;

  case SolverSpecify::TRANSIENT :
    solve_transient(); break;

  default: genius_error();
  }

  STOP_LOG("MixA2Solver_SNES()", "MixA2Solver");

  return 0;
}




/*------------------------------------------------------------------
 * restore the solution to each region
 */
int MixA2Solver::post_solve_process()
{

  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  VecGetArray(lx, &lxx);

  //search for all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Update_Solution(lxx);
  }

  VecRestoreArray(lx, &lxx);

  _circuit->save_solution();

  return MixASolverBase::post_solve_process();
}



/*------------------------------------------------------------------
 * return the extra dofs of spice circuit
 */
unsigned int MixA2Solver::extra_dofs() const
{
  return _circuit->n_ckt_nodes();
}


/*------------------------------------------------------------------
 * load previous state into solution vector
 */
int MixA2Solver::diverged_recovery()
{
  // for all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Fill_Value(x, L);
  }

  //load spice previous solution
  if(Genius::is_last_processor())
    _circuit->restore_solution();

  spice_fill_value(x, L);

  VecAssemblyBegin(x);
  VecAssemblyBegin(L);

  VecAssemblyEnd(x);
  VecAssemblyEnd(L);

  return 0;
}


/*------------------------------------------------------------------
 * Potential Newton Damping
 */
void MixA2Solver::potential_damping(Vec x, Vec y, Vec w, PetscTruth *changed_y, PetscTruth *changed_w)
{

  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;

  VecGetArray(x, &xx);  // previous iterate value
  VecGetArray(y, &yy);  // new search direction and length
  VecGetArray(w, &ww);  // current candidate iterate

  PetscScalar dV_max = 0.0; // the max changes of psi
  PetscScalar onePerCMC = 1.0*std::pow(cm,-3);
  PetscScalar T_external = this->get_system().T_external();

  // the lattice temperature limit
  // I think 50K under T_external is ok even for semiconductor cooling.
  PetscScalar TemperatureLimit = T_external - 50*K;

  // we should find dV_max;
  // first, we find in local
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = _system.region(n);
    if( region->type() != SemiconductorRegion ) continue;

    SimulationRegion::const_node_iterator it = region->nodes_begin();
    SimulationRegion::const_node_iterator it_end = region->nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = (*it).second;

      //if this node NOT belongs to this processor or not valid, continue
      if( !fvm_node->on_processor() || !fvm_node->is_valid()) continue;

      // we konw the fvm_node->local_offset() is psi in semiconductor region
      unsigned int local_offset = fvm_node->local_offset();
      if(  std::abs(yy[local_offset]) > dV_max)
        dV_max = std::abs(yy[local_offset]);

      //prevent negative carrier density
      if ( ww[local_offset+1] < onePerCMC )
        ww[local_offset+1]= onePerCMC;
      if ( ww[local_offset+2] < onePerCMC )
        ww[local_offset+2]= onePerCMC;

      // the lattice temperature limit
      if ( ww[local_offset+3] < TemperatureLimit )
        ww[local_offset+3] = TemperatureLimit;

    }
  }

  // for parallel situation, we should find the dv_max in global.
  Parallel::max( dV_max );

  if( dV_max > 1e-6 )
  {
    // compute logarithmic potential damping factor f;
    PetscScalar Vt = kb*T_external/e;
    PetscScalar f = log(1+dV_max/Vt)/(dV_max/Vt);

    // do newton damping here
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      // only consider semiconductor region
      const SimulationRegion * region = _system.region(n);
      //      if( region->type() != SemiconductorRegion ) continue;

      SimulationRegion::const_node_iterator it = region->nodes_begin();
      SimulationRegion::const_node_iterator it_end = region->nodes_end();
      for(; it!=it_end; ++it)
      {
        const FVM_Node * fvm_node = (*it).second;
        //if this node NOT belongs to this processor or not valid, continue
        if( !fvm_node->on_processor() || !fvm_node->is_valid()) continue;

        unsigned int local_offset = fvm_node->local_offset();

        ww[local_offset] = xx[local_offset] - f*yy[local_offset];

      }
    }
  }

  //damping spice nodal voltage update
  if(Genius::is_last_processor())
  {
    PetscScalar max_update=0;
    for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
      if(_circuit->is_voltage_node(n))
      {
        PetscScalar update = yy[_circuit->array_offset(n)];
        if(update > max_update)
          max_update = update;
      }

    if (max_update>1)
    {
      PetscScalar damp_factor = 1/max_update;
      if (damp_factor<0.1) damp_factor=0.1;
      for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
        if(_circuit->is_voltage_node(n))
        {
          unsigned int array_offset = _circuit->array_offset(n);
          ww[array_offset] = xx[array_offset] - damp_factor*yy[array_offset];
        }
    }
  }

  VecRestoreArray(x, &xx);
  VecRestoreArray(y, &yy);
  VecRestoreArray(w, &ww);

  *changed_y = PETSC_FALSE;
  *changed_w = PETSC_TRUE;

  return;
}



/*------------------------------------------------------------------
 * Bank-Rose Newton Damping
 */
void MixA2Solver::bank_rose_damping(Vec , Vec , Vec , PetscTruth *changed_y, PetscTruth *changed_w)
{
  *changed_y = PETSC_FALSE;
  *changed_w = PETSC_FALSE;

  return;
}



/*------------------------------------------------------------------
 * positive density Newton Damping
 */
void MixA2Solver::positive_density_damping(Vec x, Vec y, Vec w, PetscTruth *changed_y, PetscTruth *changed_w)
{

  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;

  VecGetArray(x, &xx);  // previous iterate value
  VecGetArray(y, &yy);  // new search direction and length
  VecGetArray(w, &ww);  // current candidate iterate

  int changed_flag=0;
  PetscScalar onePerCMC = 1.0*std::pow(cm,-3);
  PetscScalar T_external = this->get_system().T_external();

  // the lattice temperature limit
  // I think 50K under T_external is ok even for semiconductor cooling.
  PetscScalar TemperatureLimit = T_external - 50*K;

  // do newton damping here
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = _system.region(n);
    if( region->type() != SemiconductorRegion ) continue;

    SimulationRegion::const_node_iterator it = region->nodes_begin();
    SimulationRegion::const_node_iterator it_end = region->nodes_end();

    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = (*it).second;
      //if this node NOT belongs to this processor or not valid, continue
      if( !fvm_node->on_processor() || !fvm_node->is_valid()) continue;

      unsigned int local_offset = fvm_node->local_offset();


      // the maximum potential update is limited to 1V
      if ( fabs(yy[local_offset]) > 1.0 )
      { ww[local_offset] = xx[local_offset] - std::sign(yy[local_offset])*1.0; changed_flag = 1; }

      //prevent negative carrier density
      if ( ww[local_offset+1] < onePerCMC )
      {
        ww[local_offset+1]= onePerCMC;
        changed_flag=1;
      }
      if ( ww[local_offset+2] < onePerCMC )
      {
        ww[local_offset+2]= onePerCMC;
        changed_flag=1;
      }

      // the lattice temperature limit
      if ( ww[local_offset+3] < TemperatureLimit )
      {
        ww[local_offset+3] = TemperatureLimit;
        changed_flag=1;
      }

    }
  }

  VecRestoreArray(x, &xx);
  VecRestoreArray(y, &yy);
  VecRestoreArray(w, &ww);

  //synch changed_flag, if it is not zero, the vector is changed
  Parallel::sum( changed_flag );

  if(changed_flag)
  {
    *changed_y = PETSC_FALSE;
    *changed_w = PETSC_TRUE;
  }

  return;
}



void MixA2Solver::projection_positive_density_check(Vec x, Vec xo)
{
  PetscScalar    *xx;
  PetscScalar    *oo;

  VecGetArray(x, &xx);
  VecGetArray(xo, &oo);

  PetscScalar onePerCMC = 1.0*std::pow(cm,-3);
  PetscScalar T_external = this->get_system().T_external();

  // the lattice temperature limit
  // I think 50K under T_external is ok even for semiconductor cooling.
  PetscScalar TemperatureLimit = T_external - 50*K;

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = _system.region(n);
    if( region->type() != SemiconductorRegion ) continue;

    SimulationRegion::const_node_iterator it = region->nodes_begin();
    SimulationRegion::const_node_iterator it_end = region->nodes_end();
    for(; it!=it_end; it++)
    {
      const FVM_Node * fvm_node = (*it).second;
      //if this node NOT belongs to this processor or not valid, continue
      if( !fvm_node->on_processor() || !fvm_node->is_valid()) continue;

      unsigned int local_offset = fvm_node->local_offset();

      // prevent negative carrier density
      if ( xx[local_offset+1] < onePerCMC ) xx[local_offset+1]= onePerCMC;
      if ( xx[local_offset+2] < onePerCMC ) xx[local_offset+2]= onePerCMC;

      // lattice temperature limit
      if ( xx[local_offset+3] < TemperatureLimit ) xx[local_offset+3] = TemperatureLimit;
    }
  }

  VecRestoreArray(x,&xx);
  VecRestoreArray(xo,&oo);
}


/*------------------------------------------------------------------
 * evaluate local truncation error
 */
PetscScalar MixA2Solver::LTE_norm()
{

  // time steps
  PetscScalar hn  = SolverSpecify::dt;
  PetscScalar hn1 = SolverSpecify::dt_last;
  PetscScalar hn2 = SolverSpecify::dt_last_last;

  // relative error
  PetscScalar eps_r = SolverSpecify::TS_rtol;
  // abs error
  PetscScalar eps_a = SolverSpecify::TS_atol;

  VecZeroEntries(xp);
  VecZeroEntries(LTE);

  // get the predict solution vector and LTE vector
  if(SolverSpecify::TS_type == SolverSpecify::BDF1)
  {
    VecAXPY(xp, 1+hn/hn1, x_n);
    VecAXPY(xp, -hn/hn1, x_n1);
    VecAXPY(LTE, hn/(hn+hn1), x);
    VecAXPY(LTE, -hn/(hn+hn1), xp);
  }
  else if(SolverSpecify::TS_type == SolverSpecify::BDF2)
  {
    PetscScalar cn  = 1+hn*(hn+2*hn1+hn2)/(hn1*(hn1+hn2));
    PetscScalar cn1 = -hn*(hn+hn1+hn2)/(hn1*hn2);
    PetscScalar cn2 = hn*(hn+hn1)/(hn2*(hn1+hn2));

    VecAXPY(xp, cn,  x_n);
    VecAXPY(xp, cn1, x_n1);
    VecAXPY(xp, cn2, x_n2);
    VecAXPY(LTE, hn/(hn+hn1+hn2),  x);
    VecAXPY(LTE, -hn/(hn+hn1+hn2), xp);
  }

  int N=0;
  PetscScalar r;


  // with LTE vector and relative & abs error, we get the error estimate here
  PetscScalar    *xx, *ll;
  VecGetArray(x, &xx);
  VecGetArray(LTE, &ll);

  // error estimate for each region
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    const SimulationRegion * region = _system.region(n);

    switch ( region->type() )
    {
    case SemiconductorRegion :
      {
        SimulationRegion::const_node_iterator it = region->nodes_begin();
        SimulationRegion::const_node_iterator it_end = region->nodes_end();
        for(; it!=it_end; ++it)
        {
          const FVM_Node * fvm_node = (*it).second;
          //if this node NOT belongs to this processor, continue
          if( !fvm_node->on_processor() ) continue;

          unsigned int local_offset = fvm_node->local_offset();

          ll[local_offset+0] = 0;
          ll[local_offset+1] = ll[local_offset+1]/(eps_r*xx[local_offset+1]+eps_a); // current jn
          ll[local_offset+2] = ll[local_offset+2]/(eps_r*xx[local_offset+2]+eps_a); // current jp
          ll[local_offset+3] = ll[local_offset+3]/(eps_r*xx[local_offset+3]+eps_a); // temperature
        }
        N += 3*region->n_on_processor_node();
        break;
      }
    case InsulatorRegion :
    case ConductorRegion :
      {
        SimulationRegion::const_node_iterator it = region->nodes_begin();
        SimulationRegion::const_node_iterator it_end = region->nodes_end();
        for(; it!=it_end; ++it)
        {
          const FVM_Node * fvm_node = (*it).second;
          //if this node NOT belongs to this processor, continue
          if( !fvm_node->on_processor() ) continue;

          unsigned int local_offset = fvm_node->local_offset();

          ll[local_offset+0] = 0;
          ll[local_offset+1] = ll[local_offset+1]/(eps_r*xx[local_offset+1]+eps_a); // temperature
        }
        N += region->n_on_processor_node();
        break;
      }
    case VacuumRegion:
      break;
    default: genius_error();
    }
  }


  // error estimate for spice ckt
  if(Genius::is_last_processor())
  {
    for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
    {
      unsigned int array_offset = _circuit->array_offset(n);
      ll[array_offset] = ll[array_offset]/(eps_r*xx[array_offset]+eps_a);
    }
    N += _circuit->n_ckt_nodes();
  }

  VecRestoreArray(x, &xx);
  VecRestoreArray(LTE, &ll);

  VecNorm(LTE, NORM_2, &r);

  //for parallel situation, we should sum N of all the processor
  Parallel::sum( N );

  if( N>0 )
  {
    r /= sqrt(static_cast<double>(N));
    return r;
  }

  return 1.0;

}


#if PETSC_VERSION_LE(2,3,3)
#include "private/snesimpl.h"
#endif
void MixA2Solver::error_norm()
{
  PetscScalar    *xx;
  PetscScalar    *ff;

  // scatte global solution vector x to local vector lx
  // this is not necessary since it had already done in function evaluation
  //VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  //VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  // scatte global function vector f to local vector lf
#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
  VecScatterBegin(scatter, f, lf, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, f, lf, INSERT_VALUES, SCATTER_FORWARD);
#endif

#if PETSC_VERSION_LE(2,3,3)
  VecScatterBegin(scatter, snes->vec_func_always, lf, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, snes->vec_func_always, lf, INSERT_VALUES, SCATTER_FORWARD);
#endif

  VecGetArray(lx, &xx);  // solution value
  VecGetArray(lf, &ff);  // function value

  // do clear
  potential_norm        = 0;
  electron_norm         = 0;
  hole_norm             = 0;
  temperature_norm      = 0;
  elec_temperature_norm = 0;
  hole_temperature_norm = 0;

  poisson_norm              = 0;
  elec_continuity_norm      = 0;
  hole_continuity_norm      = 0;
  heat_equation_norm        = 0;
  elec_energy_equation_norm = 0;
  hole_energy_equation_norm = 0;
  electrode_norm            = 0;

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = _system.region(n);

    SimulationRegion::const_node_iterator it = region->nodes_begin();
    SimulationRegion::const_node_iterator it_end = region->nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = (*it).second;
      //if this node NOT belongs to this processor, continue
      if( !fvm_node->on_processor() ) continue;

      unsigned int offset = fvm_node->local_offset();

      switch ( region->type() )
      {
      case SemiconductorRegion :
        {
          potential_norm   += xx[offset+0]*xx[offset+0];
          electron_norm    += xx[offset+1]*xx[offset+1];
          hole_norm        += xx[offset+2]*xx[offset+2];
          temperature_norm += xx[offset+3]*xx[offset+3];

          poisson_norm         += ff[offset+0]*ff[offset+0];
          elec_continuity_norm += ff[offset+1]*ff[offset+1];
          hole_continuity_norm += ff[offset+2]*ff[offset+2];
          heat_equation_norm   += ff[offset+3]*ff[offset+3];
          break;
        }
      case InsulatorRegion :
        {
          potential_norm   += xx[offset]*xx[offset];
          temperature_norm += xx[offset+1]*xx[offset+1];

          poisson_norm        += ff[offset]*ff[offset];
          heat_equation_norm  += ff[offset+1]*ff[offset+1];
          break;
        }
      case ConductorRegion :
        {
          potential_norm   += xx[offset]*xx[offset];
          temperature_norm += xx[offset+1]*xx[offset+1];

          poisson_norm        += ff[offset]*ff[offset];
          heat_equation_norm  += ff[offset+1]*ff[offset+1];
          break;
        }
      case VacuumRegion:
        break;
      default: genius_error();
      }
    }
  }

  if(Genius::is_last_processor())
    electrode_norm = _circuit->ckt_residual_norm2();

  // sum of variable value on all processors
  parallel_only();
  Parallel::sum(potential_norm);
  Parallel::sum(electron_norm);
  Parallel::sum(hole_norm);
  Parallel::sum(temperature_norm);

  Parallel::sum(poisson_norm);
  Parallel::sum(elec_continuity_norm);
  Parallel::sum(hole_continuity_norm);
  Parallel::sum(heat_equation_norm);

  // sqrt to get L2 norm
  potential_norm   = sqrt(potential_norm);
  electron_norm    = sqrt(electron_norm);
  hole_norm        = sqrt(hole_norm);
  temperature_norm = sqrt(temperature_norm);

  poisson_norm        = sqrt(poisson_norm);
  elec_continuity_norm = sqrt(elec_continuity_norm);
  hole_continuity_norm = sqrt(hole_continuity_norm);
  heat_equation_norm  = sqrt(heat_equation_norm);

  Parallel::broadcast(electrode_norm, Genius::last_processor_id());

  VecRestoreArray(lx, &xx);
  VecRestoreArray(lf, &ff);

}




///////////////////////////////////////////////////////////////
// provide function and jacobian evaluation for DDML2 solver //
///////////////////////////////////////////////////////////////





/*------------------------------------------------------------------
 * evaluate the residual of function f at x
 */
void MixA2Solver::build_petsc_sens_residual(Vec x, Vec r)
{

  START_LOG("MixA2Solver_Residual()", "MixA2Solver");

  // scatte global solution vector x to local vector lx
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray(lx, &lxx);

  // clear old data
  VecZeroEntries (r);

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

  // evaluate governing equations of DDML2 in all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Function(lxx, r, add_value_flag);
  }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate time derivative if necessary
  if(SolverSpecify::TimeDependent == true)
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DDM2_Time_Dependent_Function(lxx, r, add_value_flag);
    }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // process hanging node here
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Function_Hanging_Node(lxx, r, add_value_flag);
  }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  build_spice_function(lxx, r, add_value_flag);

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate governing equations of MixA2 for all the boundaries
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->MixA_DDM2_Function(lxx, r, add_value_flag);
  }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);

  this->ground_spice_0_node(r, add_value_flag);

  // assembly the function Vec
  VecAssemblyBegin(r);
  VecAssemblyEnd(r);

  // scale the function vec
  PetscScalar *ff,*scale;
  // get function array and scale array.
  VecGetArray(r, &ff);
  // L is the scaling vector, the Jacobian evaluate function may dynamically update it.
  VecGetArray(L, &scale);

  // scale it!
  for(unsigned int n=0; n<n_local_dofs; n++)
    ff[n] *= scale[n];

  // restore back
  VecRestoreArray(r, &ff);
  VecRestoreArray(L, &scale);

  STOP_LOG("MixA2Solver_Residual()", "MixA2Solver");

}



/*------------------------------------------------------------------
 * evaluate the Jacobian J of function f at x
 */
void MixA2Solver::build_petsc_sens_jacobian(Vec x, Mat *, Mat *)
{

  START_LOG("MixA2Solver_Jacobian()", "MixA2Solver");

  // scatte global solution vector x to local vector lx
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray(lx, &lxx);

  MatZeroEntries(J);

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

  // evaluate Jacobian matrix of governing equations of DDML2 in all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Jacobian(lxx, &J, add_value_flag);
  }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate Jacobian matrix of time derivative if necessary
  if(SolverSpecify::TimeDependent == true)
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DDM2_Time_Dependent_Jacobian(lxx, &J, add_value_flag);
    }


#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  build_spice_jacobian(lxx, &J, add_value_flag);

  // before first assemble, resereve none zero pattern for each boundary

  if( !jacobian_matrix_first_assemble )
  {
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      bc->MixA_DDM2_Jacobian_Reserve(&J, add_value_flag);
    }

    jacobian_matrix_first_assemble = true;

    // after that, we do not allow zero insert/add to matrix
#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
    genius_assert(!MatSetOption(J, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
#endif

#if PETSC_VERSION_LE(2,3,3)
    genius_assert(!MatSetOption(J, MAT_IGNORE_ZERO_ENTRIES));
#endif

  }

  // process hanging node here
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Jacobian_Hanging_Node(lxx, &J, add_value_flag);
  }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate Jacobian matrix of governing equations of DDML2 for all the boundaries
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->MixA_DDM2_Jacobian(lxx, &J, add_value_flag);
  }


#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);

  // assembly the matrix
  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (J, MAT_FINAL_ASSEMBLY);

  this->ground_spice_0_node(&J);

  //scaling the matrix
  MatDiagonalScale(J, L, PETSC_NULL);

  // we use the reciprocal of matrix diagonal as scaling value
  // this will take place at next call
  //MatGetDiagonal(J, L);
  //VecReciprocal(L);

  //MatView(J, PETSC_VIEWER_DRAW_WORLD);
  //getchar();

  STOP_LOG("MixA2Solver_Jacobian()", "MixA2Solver");


}


