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
#include "petsc_utils.h"

//#define DEBUG

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::cm;
using PhysicalUnit::K;
using PhysicalUnit::A;
using PhysicalUnit::V;

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

    // for all the bcs
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      if(bc->is_spice_electrode())
        bc->MixA_DDM2_Fill_Value(x, L);
      else
        bc->DDM2_Fill_Value(x, L);
    }

    spice_fill_value(x, L);

    VecAssemblyBegin(x);
    VecAssemblyBegin(L);

    VecAssemblyEnd(x);
    VecAssemblyEnd(L);
  }

  // do bc pre process
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DDM2_Pre_Process();
  }

  return MixASolverBase::pre_solve_process(load_solution);
}





/*------------------------------------------------------------------
 * the main solve routine. which is under the control of ngspice
 */
int MixA2Solver::solve()
{

  START_LOG("solve()", "MixA2Solver");

  switch( SolverSpecify::Type )
  {
      case SolverSpecify::OP :
      solve_dcop(); break;

      case SolverSpecify::DCSWEEP :
      solve_dcsweep(); break;

      case SolverSpecify::TRANSIENT :
      solve_transient(); break;

      default:
      {
        MESSAGE<<"ERROR: MixA2Solver does not support this solve type." << std::endl; RECORD();
        //genius_error();
      }
  }

  STOP_LOG("solve()", "MixA2Solver");

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


  // extra work, calculate the electric field
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);

    SimulationRegion::processor_node_iterator it = region->on_processor_nodes_begin();
    SimulationRegion::processor_node_iterator it_end = region->on_processor_nodes_end();
    for(; it!=it_end; ++it)
    {
      FVM_Node * fvm_node = *it;
      FVM_NodeData * node_data = fvm_node->node_data();
      node_data->E() = -fvm_node->gradient(POTENTIAL, true);
    }
  }


  // update bcs
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->MixA_DDM2_Update_Solution(lxx);
  }

  VecRestoreArray(lx, &lxx);

  _circuit->save_solution();

  return MixASolverBase::post_solve_process();
}


/*------------------------------------------------------------------
 * write the (intermediate) solution to each region
 */
void MixA2Solver::flush_system(Vec v)
{
  VecScatterBegin(scatter, v, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, v, lx, INSERT_VALUES, SCATTER_FORWARD);

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
void MixA2Solver::potential_damping(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{
  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;

  VecGetArray(x, &xx);  // previous iterate value
  VecGetArray(y, &yy);  // new search direction and length
  VecGetArray(w, &ww);  // current candidate iterate

  const PetscScalar onePerMC = 1.0e-6*std::pow(cm,-3);
  const PetscScalar T_external = this->get_system().T_external();
  // the lattice temperature limit
  // I think 50K under T_external is ok even for semiconductor cooling.
  PetscScalar TemperatureLimit = T_external - 50*K;

  PetscScalar dV_min =  std::numeric_limits<PetscScalar>::max(); // the min changes of psi
  PetscScalar dV_max = -std::numeric_limits<PetscScalar>::max(); // the max changes of psi
  PetscScalar dV;
  PetscScalar dV_comm;
  // we should find dV_max/dV_min;
  // first, we find in local
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
      dV_max = std::max(dV_max, yy[local_offset]);
      dV_min = std::min(dV_min, yy[local_offset]);
    }
  }

  // for parallel situation, we should find the dv_max/dV_min in global.
  Parallel::max( dV_max );
  Parallel::min( dV_min );

  if(dV_max < dV_min) goto potential_damping_end;

  // find the max/min in absolutely value
  if( std::abs(dV_max) < std::abs(dV_min) ) std::swap(dV_max, dV_min);

  // compute the dV and dV in common
  dV = std::abs(dV_max - dV_min);
  dV_comm = dV_max*dV_min > 0.0 ? dV_min : 0.0;

  if( dV > 1e-6*V )
  {
    // compute logarithmic potential damping factor f;
    PetscScalar Vut = kb*this->get_system().T_external()/e * SolverSpecify::potential_update;
    PetscScalar f = log(1+dV/Vut)/(dV/Vut);

    // do newton damping here
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
              ww[local_offset+0] = xx[local_offset+0] - (dV_comm + f*(yy[local_offset+0]-dV_comm));

              //prevent negative carrier density
              if ( ww[local_offset+1] < 0 )
              { ww[local_offset+1] = 1e-2*fabs(xx[local_offset+1]) + onePerMC;}
              if ( ww[local_offset+2] < 0 )
              { ww[local_offset+2] = 1e-2*fabs(xx[local_offset+2]) + onePerMC;}
              // the lattice temperature limit
              if ( ww[local_offset+3] < TemperatureLimit )
              { ww[local_offset+3] = TemperatureLimit; }
            }
            break;
          }
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
              ww[local_offset+0] = xx[local_offset+0] - (dV_comm + f*(yy[local_offset+0]-dV_comm));
            }
            break;
          }
          default: break;
      }
    }
  }

potential_damping_end:

  if(Genius::is_last_processor() && SolverSpecify::damping_spice)
  {
    for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
    {
      if(_circuit->is_voltage_node(n) && !_circuit->is_link_to_electrode(n))
      {
        unsigned int array_offset_x = _circuit->array_offset_x(n);
        PetscScalar dV_max = std::abs(yy[array_offset_x]);
        if (dV_max>1)
        {
          PetscScalar damp_factor = log(1+dV_max/20)/(dV_max/20);
          ww[array_offset_x] = xx[array_offset_x] - yy[array_offset_x]*damp_factor;
        }
      }

      if(_circuit->is_current_node(n) && !_circuit->is_link_to_electrode(n))
      {
        unsigned int array_offset_x = _circuit->array_offset_x(n);
        PetscScalar dI_max = std::abs(yy[array_offset_x]);
        if (dI_max>1)
        {
          PetscScalar damp_factor = log(1+dI_max/10)/(dI_max/10);
          ww[array_offset_x] = xx[array_offset_x] - yy[array_offset_x]*damp_factor;
        }
      }

    }
  }



  *changed_w = PETSC_TRUE;

  VecRestoreArray(x, &xx);
  VecRestoreArray(y, &yy);
  VecRestoreArray(w, &ww);

  return;
}



/*------------------------------------------------------------------
 * Bank-Rose Newton Damping
 */
void MixA2Solver::bank_rose_damping(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{
  return;
}



/*------------------------------------------------------------------
 * positive density Newton Damping
 */
void MixA2Solver::check_positive_density(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{

  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;

  VecGetArray(x, &xx);  // previous iterate value
  VecGetArray(y, &yy);  // new search direction and length
  VecGetArray(w, &ww);  // current candidate iterate

  int changed_flag=0;
  const PetscScalar onePerCMC = 1.0*std::pow(cm,-3);
  const PetscScalar T_external = this->get_system().T_external();

  // the lattice temperature limit
  // I think 50K under T_external is ok even for semiconductor cooling.
  PetscScalar TemperatureLimit = T_external - 50*K;

  // do newton damping here
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

      unsigned int local_offset = fvm_node->local_offset();

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

  const PetscScalar onePerCMC = 1.0*std::pow(cm,-3);
  const PetscScalar T_external = this->get_system().T_external();

  // the lattice temperature limit
  // I think 50K under T_external is ok even for semiconductor cooling.
  PetscScalar TemperatureLimit = T_external - 50*K;

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
 * test if BDF2 can be used for next time step
 */
bool MixA2Solver::BDF2_positive_defined() const
{
  const double r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
  const double a = 1.0/(r*(1-r));
  const double b = (1-r)/r;

  unsigned int failure_count=0;
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    const SimulationRegion * region = _system.region(n);
    if ( region->type() == SemiconductorRegion)
    {
      SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
      SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
      for(; it!=it_end; ++it)
      {
        const FVM_Node * fvm_node = *it;
        const FVM_NodeData * node_data = fvm_node->node_data();

        if(a*node_data->n() < b*node_data->n_last() ) failure_count++;
        if(a*node_data->p() < b*node_data->p_last() ) failure_count++;
        if(a*node_data->T() < b*node_data->T_last() ) failure_count++;
      }
    }
  }

  Parallel::sum(failure_count);
  return failure_count != 0;
}


/*------------------------------------------------------------------
 * evaluate local truncation error
 */
PetscReal MixA2Solver::LTE_norm()
{

  // time steps
  PetscReal hn  = SolverSpecify::dt;
  PetscReal hn1 = SolverSpecify::dt_last;
  PetscReal hn2 = SolverSpecify::dt_last_last;

  // relative error
  PetscReal eps_r = SolverSpecify::TS_rtol;
  // abs error
  PetscReal eps_a = SolverSpecify::TS_atol;
  PetscReal concentration = 5e22*std::pow(cm, -3);
  PetscReal temperature = 10000*K;


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
    if(SolverSpecify::BDF2_LowerOrder)
    {
      VecAXPY(xp, 1+hn/hn1, x_n);
      VecAXPY(xp, -hn/hn1, x_n1);
      VecAXPY(LTE, hn/(hn+hn1), x);
      VecAXPY(LTE, -hn/(hn+hn1), xp);
    }
    else
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
  }

  int N=0;
  PetscReal r;


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
            ll[local_offset+1] = ll[local_offset+1]/(eps_r*xx[local_offset+1]+eps_a*concentration); // current jn
            ll[local_offset+2] = ll[local_offset+2]/(eps_r*xx[local_offset+2]+eps_a*concentration); // current jp
            ll[local_offset+3] = ll[local_offset+3]/(eps_r*xx[local_offset+3]+eps_a*temperature); // temperature
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
            ll[local_offset+1] = ll[local_offset+1]/(eps_r*xx[local_offset+1]+eps_a*temperature); // temperature
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
      unsigned int array_offset_f = _circuit->array_offset_f(n);
      unsigned int array_offset_x = _circuit->array_offset_x(n);
      if(_circuit->is_voltage_node(n))
      {
        ll[array_offset_f] = ll[array_offset_f]/(eps_r*xx[array_offset_x]+eps_a*1000*V);
        N++;
      }
    }
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


void MixA2Solver::error_norm()
{
  PetscScalar    *xx;
  PetscScalar    *ff;

  // scatte global solution vector x to local vector lx
  // this is not necessary since it had already done in function evaluation
  //VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  //VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  // unscale the fcuntion
  VecPointwiseDivide(f, f, L);

  // scatte global function vector f to local vector lf
  VecScatterBegin(scatter, f, lf, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, f, lf, INSERT_VALUES, SCATTER_FORWARD);

  // scale the function vector back
  VecPointwiseMult(f, f, L);


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
          case ElectrodeRegion :
          case MetalRegion     :
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
    spice_norm = _circuit->ckt_residual_norm2()*A;
  Parallel::broadcast(spice_norm, Genius::last_processor_id());

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

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
  
  // evaluate governing equations of DDML2 in all the regions
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

  build_spice_function(lxx, r, add_value_flag);

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate governing equations of MixA2 for all the boundaries

  // preprocess each bc
  VecAssemblyBegin(r);
  VecAssemblyEnd(r);
  std::vector<PetscInt> src_row,  dst_row,  clear_row;
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    if(bc->is_spice_electrode())
      bc->MixA_DDM2_Function_Preprocess(lxx, r, src_row, dst_row, clear_row);
    else
      bc->DDM2_Function_Preprocess(lxx, r, src_row, dst_row, clear_row);
  }
  //add source rows to destination rows, and clear rows
  PetscUtils::VecAddClearRow(r, src_row, dst_row, clear_row);
  add_value_flag = NOT_SET_VALUES;

  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    if(bc->is_spice_electrode())
      bc->MixA_DDM2_Function(lxx, r, add_value_flag);
    else
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
  VecPointwiseMult(r, r, L);

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

  Jac->zero();

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

  // evaluate Jacobian matrix of governing equations of DDML2 in all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Jacobian(lxx, Jac, add_value_flag);
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate Jacobian matrix of time derivative if necessary
  if(SolverSpecify::TimeDependent == true)
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DDM2_Time_Dependent_Jacobian(lxx, Jac, add_value_flag);
    }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  build_spice_jacobian(lxx, Jac, add_value_flag);



#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  
  // assembly matrix
  Jac->close(false);


  std::vector<PetscInt> src_row,  dst_row,  clear_row;
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    if(bc->is_spice_electrode())
      bc->MixA_DDM2_Jacobian_Preprocess(lxx, Jac, src_row, dst_row, clear_row);
    else
      bc->DDM2_Jacobian_Preprocess(lxx, Jac, src_row, dst_row, clear_row);
  }
  
  //add source rows to destination rows
  Jac->add_row_to_row(src_row, dst_row);
  // clear row
  Jac->clear_row(clear_row);
  
  add_value_flag = NOT_SET_VALUES;

  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    if(bc->is_spice_electrode())
      bc->MixA_DDM2_Jacobian(lxx, Jac, add_value_flag);
    else
      bc->DDM2_Jacobian(lxx, Jac, add_value_flag);
  }


#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);

  Jac->close(true);


  //scaling the matrix
  MatDiagonalScale(J, L, PETSC_NULL);


  STOP_LOG("MixA2Solver_Jacobian()", "MixA2Solver");


}


