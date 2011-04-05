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


#include "ddm2/ddm2.h"
#include "parallel.h"
#include "petsc_utils.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::cm;
using PhysicalUnit::K;
using PhysicalUnit::V;



/*------------------------------------------------------------------
 * create nonlinear solver contex and adjust some parameters
 */
int DDM2Solver::create_solver()
{

  MESSAGE<< '\n' << "DDM Solver Level 2 init..."  << std::endl;
  RECORD();

  return DDMSolverBase::create_solver();
}




/*------------------------------------------------------------------
 * set initial value to solution vector and scaling vector
 */
int DDM2Solver::pre_solve_process(bool load_solution)
{
  if(load_solution)
  {
    // for all the regions
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DDM2_Fill_Value(x, L);
    }


    // for all the bcs
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      bc->DDM2_Fill_Value(x, L);
    }

    VecAssemblyBegin(x);
    VecAssemblyBegin(L);

    VecAssemblyEnd(x);
    VecAssemblyEnd(L);
  }

  DDMSolverBase::pre_solve_process(load_solution);

  return 0;
}





/*------------------------------------------------------------------
 * wrap to each solve implimention
 */
int DDM2Solver::solve()
{

  START_LOG("DDM2Solver_SNES()", "DDM2Solver");

  switch( SolverSpecify::Type )
  {
      case SolverSpecify::EQUILIBRIUM :
      solve_equ(); break;

      case SolverSpecify::STEADYSTATE:
      solve_steadystate(); break;

      case SolverSpecify::DCSWEEP:
      solve_dcsweep(); break;

      case SolverSpecify::TRANSIENT:
      solve_transient();break;

      case SolverSpecify::TRACE:
      solve_iv_trace();break;

      default:
      MESSAGE<< '\n' << "DDM2Solver: Unsupported solve type."; RECORD();
      genius_error();
      break;
  }


  STOP_LOG("DDM2Solver_SNES()", "DDM2Solver");

  return 0;
}




/*------------------------------------------------------------------
 * restore the solution to each region
 */
int DDM2Solver::post_solve_process()
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

  // update bcs
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DDM2_Update_Solution(lxx);
  }

  VecRestoreArray(lx, &lxx);


  DDMSolverBase::post_solve_process();

  return 0;
}



/*------------------------------------------------------------------
 * write the (intermediate) solution to each region
 */
void DDM2Solver::flush_system()
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
}



/*------------------------------------------------------------------
 * load previous state into solution vector
 */
int DDM2Solver::diverged_recovery()
{
  // for all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Fill_Value(x, L);
  }

  // for all the bcs
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DDM2_Fill_Value(x, L);
  }

  VecAssemblyBegin(x);
  VecAssemblyBegin(L);

  VecAssemblyEnd(x);
  VecAssemblyEnd(L);

  return 0;
}


#if 0
/**
 * function for line search pre check.
 */
void DDM2Solver::sens_line_search_pre_check(Vec , Vec y, PetscBool *changed_y)
{
  PetscScalar    *yy;
  VecGetArray(y, &yy);  // new search direction and length

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = _system.region(n);
    if( region->type() != SemiconductorRegion ) continue;

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = *it;

      // we konw the fvm_node->local_offset() is psi in semiconductor region
      unsigned int local_offset = fvm_node->local_offset();

      if ( yy[local_offset+0] > 2.0*V )
        yy[local_offset+0] = 2.0*V;

      // the lattice temperature limit
      if ( yy[local_offset+3] > 50.0*K )
        yy[local_offset+3] = 50.0*K;
    }
  }

  VecRestoreArray(y, &yy);
  *changed_y = PETSC_TRUE;
  return;
}

#endif



/*------------------------------------------------------------------
 * Potential Newton Damping
 */
void DDM2Solver::potential_damping(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{

  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;

  VecGetArray(x, &xx);  // previous iterate value
  VecGetArray(y, &yy);  // new search direction and length
  VecGetArray(w, &ww);  // current candidate iterate

  PetscScalar dV_max = 0.0; // the max changes of psi
  const PetscScalar onePerCMC = 1.0*std::pow(cm,-3);

  // the lattice temperature limit
  // I think 50K under T_external is ok even for semiconductor cooling.
  PetscScalar TemperatureLimit = this->get_system().T_external() - 50*K;

  // we should find dV_max;
  // first, we find in local
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = _system.region(n);
    if( region->type() == SemiconductorRegion )
    {
      SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        const FVM_Node * fvm_node = *it;

        // we konw the fvm_node->local_offset() is psi in semiconductor region
        unsigned int local_offset = fvm_node->local_offset();
        dV_max = std::max(dV_max, std::abs(yy[local_offset]));

        // the lattice temperature limit
        if ( ww[local_offset+3] < TemperatureLimit )
          ww[local_offset+3] = TemperatureLimit;

        //prevent negative carrier density
        if ( ww[local_offset+1] < 0 )
          ww[local_offset+1] = onePerCMC;
        if ( ww[local_offset+2] < 0 )
          ww[local_offset+2] = onePerCMC;
      }
    }
  }

  // for parallel situation, we should find the dv_max in global.
  Parallel::max( dV_max );

  if( dV_max > 1e-6*V )
  {
    // compute logarithmic potential damping factor f;
    PetscScalar Vut = kb*this->get_system().T_external()/e * SolverSpecify::potential_update;
    PetscScalar f = log(1+dV_max/Vut)/(dV_max/Vut);

    // do newton damping here
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      const SimulationRegion * region = _system.region(n);
      switch ( region->type() )
      {
          case SemiconductorRegion :
          case InsulatorRegion :
          case ElectrodeRegion :
          case MetalRegion :
          {
            SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
            SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
            for(; it!=it_end; ++it)
            {
              const FVM_Node * fvm_node = *it;
              unsigned int local_offset = fvm_node->local_offset();
              ww[local_offset] = xx[local_offset] - f*yy[local_offset];
            }
            break;
          }
          default: break;
      }
    }

    // only the last processor do this
    if(Genius::is_last_processor())
    {
      for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
      {
        BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
        unsigned int array_offset = bc->array_offset();
        if(array_offset != invalid_uint)
          ww[array_offset] = xx[array_offset] - f*yy[array_offset];
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
void DDM2Solver::bank_rose_damping(Vec , Vec , Vec , PetscBool *changed_y, PetscBool *changed_w)
{
  *changed_y = PETSC_FALSE;
  *changed_w = PETSC_FALSE;

  return;
}



/*------------------------------------------------------------------
 * positive density Newton Damping
 */
void DDM2Solver::positive_density_damping(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{

  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;

  VecGetArray(x, &xx);  // previous iterate value
  VecGetArray(y, &yy);  // new search direction and length
  VecGetArray(w, &ww);  // current candidate iterate

  int changed_flag=0;
  const PetscScalar onePerCMC = 1.0*std::pow(cm,-3);

  // the lattice temperature limit
  // I think 50K under T_external is ok even for semiconductor cooling.
  PetscScalar TemperatureLimit = this->get_system().T_external() - 50*K;

  // do newton damping here
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = _system.region(n);
    if( region->type() == SemiconductorRegion )
    {
      SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        const FVM_Node * fvm_node = *it;
        unsigned int local_offset = fvm_node->local_offset();

        // the maximum potential update is limited to 1V
        if ( fabs(yy[local_offset]) > 1.0 )
        { ww[local_offset] = xx[local_offset] - std::sign(yy[local_offset])*1.0; changed_flag = 1; }

        // the lattice temperature limit
        if ( ww[local_offset+3] < TemperatureLimit )
        { ww[local_offset+3] = TemperatureLimit; changed_flag=1; }


        //prevent negative carrier density
        if ( ww[local_offset+1] < 0 )
        { ww[local_offset+1] = onePerCMC; changed_flag = 1;}
        if ( ww[local_offset+2] < 0 )
        { ww[local_offset+2] = onePerCMC; changed_flag = 1;}
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



void DDM2Solver::projection_positive_density_check(Vec x, Vec xo)
{
  PetscScalar    *xx;
  PetscScalar    *oo;

  VecGetArray(x, &xx);
  VecGetArray(xo, &oo);

  const PetscScalar onePerCMC = 1.0*std::pow(cm,-3);

  // the lattice temperature limit
  // I think 50K under T_external is ok even for semiconductor cooling.
  PetscScalar TemperatureLimit = this->get_system().T_external() - 50*K;

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = _system.region(n);
    if( region->type() == SemiconductorRegion )
    {
      SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        const FVM_Node * fvm_node = *it;
        unsigned int local_offset = fvm_node->local_offset();

        // lattice temperature limit
        if ( xx[local_offset+3] < TemperatureLimit ) xx[local_offset+3] = TemperatureLimit;

        //prevent negative carrier density
        if ( xx[local_offset+1] < 0 )
          xx[local_offset+1]= onePerCMC;
        if ( xx[local_offset+2] < 0 )
          xx[local_offset+2]= onePerCMC;
      }
    }
  }

  VecRestoreArray(x,&xx);
  VecRestoreArray(xo,&oo);
}



/*------------------------------------------------------------------
 * evaluate local truncation error
 */
PetscScalar DDM2Solver::LTE_norm()
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
          SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
          SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
          for(; it!=it_end; ++it)
          {
            const FVM_Node * fvm_node = *it;
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
        case ElectrodeRegion :
        case MetalRegion     :
        {
          SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
          SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
          for(; it!=it_end; ++it)
          {
            const FVM_Node * fvm_node = *it;
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


  // error estimate for each bc
  if(Genius::is_last_processor())
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      unsigned int array_offset = bc->array_offset();
      if( array_offset != invalid_uint )
        ll[array_offset] = 0;
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




void DDM2Solver::error_norm()
{
  PetscScalar    *xx;
  PetscScalar    *ff;

  // scatte global solution vector x to local vector lx
  // this is not necessary since it had already done in function evaluation
  //VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  //VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  // scatte global function vector f to local vector lf
  VecScatterBegin(scatter, f, lf, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, f, lf, INSERT_VALUES, SCATTER_FORWARD);



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

    SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      const FVM_Node * fvm_node = *it;
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
          case ElectrodeRegion :
          {
            potential_norm   += xx[offset]*xx[offset];
            temperature_norm += xx[offset+1]*xx[offset+1];

            poisson_norm        += ff[offset]*ff[offset];
            heat_equation_norm  += ff[offset+1]*ff[offset+1];
            break;
          }
          case MetalRegion :
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
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      const BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      unsigned int offset = bc->local_offset();
      if( offset != invalid_uint )
      {
        potential_norm += xx[offset]*xx[offset];
        electrode_norm += ff[offset]*ff[offset];
      }
    }

  // sum of variable value on all processors
  std::vector<PetscScalar> norm_buffer;

  norm_buffer.push_back(potential_norm);
  norm_buffer.push_back(electron_norm);
  norm_buffer.push_back(hole_norm);
  norm_buffer.push_back(temperature_norm);

  norm_buffer.push_back(poisson_norm);
  norm_buffer.push_back(elec_continuity_norm);
  norm_buffer.push_back(hole_continuity_norm);
  norm_buffer.push_back(heat_equation_norm);
  norm_buffer.push_back(electrode_norm);

  Parallel::sum(norm_buffer);


  // sqrt to get L2 norm
  potential_norm   = sqrt(norm_buffer[0]);
  electron_norm    = sqrt(norm_buffer[1]);
  hole_norm        = sqrt(norm_buffer[2]);
  temperature_norm = sqrt(norm_buffer[3]);

  poisson_norm        = sqrt(norm_buffer[4]);
  elec_continuity_norm = sqrt(norm_buffer[5]);
  hole_continuity_norm = sqrt(norm_buffer[6]);
  heat_equation_norm  = sqrt(norm_buffer[7]);
  electrode_norm      = sqrt(norm_buffer[8]);

  VecRestoreArray(lx, &xx);
  VecRestoreArray(lf, &ff);

}





///////////////////////////////////////////////////////////////
// provide function and jacobian evaluation for DDML1 solver //
///////////////////////////////////////////////////////////////





/*------------------------------------------------------------------
 * evaluate the residual of function f at x
 */
void DDM2Solver::build_petsc_sens_residual(Vec x, Vec r)
{
  START_LOG("DDM2Solver_Residual()", "DDM2Solver");

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

  // evaluate governing equations of DDML1 in all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Function(lxx, r, add_value_flag);
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate time derivative if necessary
  if(SolverSpecify::TimeDependent == true)
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DDM2_Time_Dependent_Function(lxx, r, add_value_flag);
    }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // process hanging node here
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Function_Hanging_Node(lxx, r, add_value_flag);
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // preprocess each bc
  VecAssemblyBegin(r);
  VecAssemblyEnd(r);
  std::vector<PetscInt> src_row,  dst_row,  clear_row;
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DDM2_Function_Preprocess(r, src_row, dst_row, clear_row);
  }
  //add source rows to destination rows, and clear rows
  PetscUtils::VecAddClearRow(r, src_row, dst_row, clear_row);
  add_value_flag = NOT_SET_VALUES;

  // evaluate governing equations of DDML1 for all the boundaries
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DDM2_Function(lxx, r, add_value_flag);
  }


#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);

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

  STOP_LOG("DDM2Solver_Residual()", "DDM2Solver");

}



/*------------------------------------------------------------------
 * evaluate the Jacobian J of function f at x
 */
void DDM2Solver::build_petsc_sens_jacobian(Vec x, Mat *, Mat *)
{
  START_LOG("DDM2Solver_Jacobian()", "DDM2Solver");

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

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate Jacobian matrix of time derivative if necessary
  if(SolverSpecify::TimeDependent == true)
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DDM2_Time_Dependent_Jacobian(lxx, &J, add_value_flag);
    }

  // process hanging node here
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Jacobian_Hanging_Node(lxx, &J, add_value_flag);
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // before first assemble, resereve none zero pattern for each boundary
  if( !jacobian_matrix_first_assemble )
  {
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      bc->DDM2_Jacobian_Reserve(&J, add_value_flag);
    }
  }



#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate Jacobian matrix of governing equations of DDML2 for all the boundaries
  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);

  // we do not allow zero insert/add to matrix
  if( !jacobian_matrix_first_assemble )
    genius_assert(!MatSetOption(J, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));

  std::vector<PetscInt> src_row,  dst_row,  clear_row;
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DDM2_Jacobian_Preprocess(&J, src_row, dst_row, clear_row);
  }

  //add source rows to destination rows
  PetscUtils::MatAddRowToRow(J, src_row, dst_row);
  // clear row
  PetscUtils::MatZeroRows(J, clear_row.size(), clear_row.empty() ? NULL : &clear_row[0], 0.0);
  add_value_flag = NOT_SET_VALUES;

  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DDM2_Jacobian(lxx, &J, add_value_flag);
  }


#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);

  // assembly the matrix
  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (J, MAT_FINAL_ASSEMBLY);

  //scaling the matrix
  MatDiagonalScale(J, L, PETSC_NULL);

  // we use the reciprocal of matrix diagonal as scaling value
  // this will take place at next call
  //MatGetDiagonal(J, L);
  //VecReciprocal(L);


  //MatView(J,PETSC_VIEWER_STDOUT_SELF );
  //getchar();

  if(!jacobian_matrix_first_assemble)
    jacobian_matrix_first_assemble = true;

  STOP_LOG("DDM2Solver_Jacobian()", "DDM2Solver");
}


void DDM2Solver::set_trace_electrode(BoundaryCondition *bc)
{
  // we needn't scatter again
  // scatte global solution vector x to local vector lx
  //VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  //VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  bc->DDM2_Electrode_Trace(lx, &J, pdI_pdx, pdF_pdV);
}

