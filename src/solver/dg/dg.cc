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




#include "dg/dg.h"
#include "parallel.h"
#include "petsc_utils.h"
#include "mat_analysis.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::cm;



/*------------------------------------------------------------------
 * create nonlinear solver contex and adjust some parameters
 */
int DGSolver::create_solver()
{

  MESSAGE<< '\n' << "Density Gradient init..." << std::endl;
  RECORD();

  return DDMSolverBase::create_solver();

}




/*------------------------------------------------------------------
 * set initial value to solution vector and scaling vector
 */
int DGSolver::pre_solve_process(bool load_solution)
{
  if(load_solution)
  {
    // for all the regions
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DG_Fill_Value(x, L);
    }


    // for all the bcs
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      bc->DG_Fill_Value(x, L);
    }

    VecAssemblyBegin(x);
    VecAssemblyBegin(L);

    VecAssemblyEnd(x);
    VecAssemblyEnd(L);
  }


  // do bc pre process
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DG_Pre_Process();
  }

  return DDMSolverBase::pre_solve_process(load_solution);
}





/*------------------------------------------------------------------
 * wrap to each solve implimention
 */
int DGSolver::solve()
{

  START_LOG("solve()", "DGSolver");

  int ierr=0;

  switch( SolverSpecify::Type )
  {
  case SolverSpecify::EQUILIBRIUM :
    ierr=solve_equ();
    break;

  case SolverSpecify::STEADYSTATE:
    ierr=solve_steadystate();
    break;

  case SolverSpecify::DCSWEEP:
    ierr=solve_dcsweep();
    break;

  case SolverSpecify::OP:
    ierr=solve_op();
    break;

  case SolverSpecify::TRANSIENT:
    ierr=solve_transient();
    break;

  case SolverSpecify::TRACE:
    ierr=solve_iv_trace();
    break;

  default:
    MESSAGE<< '\n' << "DGSolver: Unsupported solve type.";
    RECORD();
    //genius_error();
    break;
  }

  STOP_LOG("solve()", "DGSolver");

  return ierr;
}




/*------------------------------------------------------------------
 * restore the solution to each region/boundary
 */
int DGSolver::post_solve_process()
{

  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  VecGetArray(lx, &lxx);

  //search for all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DG_Update_Solution(lxx);
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
    bc->DG_Update_Solution(lxx);
  }

  // do bc post process
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DG_Post_Process();
  }

  VecRestoreArray(lx, &lxx);

  return DDMSolverBase::post_solve_process();
}





/*------------------------------------------------------------------
 * load previous state into solution vector
 */
int DGSolver::diverged_recovery()
{
  // for all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DG_Fill_Value(x, L);
  }

  // for all the bcs
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DG_Fill_Value(x, L);
  }

  VecAssemblyBegin(x);
  VecAssemblyBegin(L);

  VecAssemblyEnd(x);
  VecAssemblyEnd(L);

  return 0;
}


/*------------------------------------------------------------------
 * Potential Newton Damping
 */
void DGSolver::potential_damping(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{

  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;


  VecGetArray(x, &xx);  // previous iterate value
  VecGetArray(y, &yy);  // new search direction and length
  VecGetArray(w, &ww);  // current candidate iterate


  PetscScalar dV_max = 0.0; // the max changes of psi
  const PetscScalar T = this->get_system().T_external();
  const PetscScalar onePerMC = 1.0e-6*std::pow(cm,-3);

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
      }
    }
  }

  // for parallel situation, we should find the dv_max in global.
  Parallel::max( dV_max );

  if( dV_max > 1e-6 )
  {
    // compute logarithmic potential damping factor f;
    PetscScalar Vut = kb*T/e * SolverSpecify::potential_update;
    PetscScalar f = log(1+dV_max/Vut)/(dV_max/Vut);

    // do newton damping here
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      const SimulationRegion * region = _system.region(n);
      switch ( region->type() )
      {
      case SemiconductorRegion :
        {
          unsigned int qn_offset = invalid_uint;
          unsigned int qp_offset = invalid_uint;
          if(region->advanced_model().QNEnabled && !region->advanced_model().QPEnabled) qn_offset=3;
          if(!region->advanced_model().QNEnabled && region->advanced_model().QPEnabled) qp_offset=3;
          if(region->advanced_model().QNEnabled && region->advanced_model().QPEnabled)  {qn_offset=3; qp_offset=4;}

          SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
          SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
          for(; it!=it_end; ++it)
          {
            const FVM_Node * fvm_node = *it;
            unsigned int local_offset = fvm_node->local_offset();
            ww[local_offset+0] = xx[local_offset+0] - f*yy[local_offset+0];

            //prevent negative carrier density
            if ( ww[local_offset+1] < 0 )
              ww[local_offset+1] = 1e-2*fabs(xx[local_offset+1]) + onePerMC;
            if ( ww[local_offset+2] < 0 )
              ww[local_offset+2] = 1e-2*fabs(xx[local_offset+2]) + onePerMC;

            if(qn_offset!=invalid_uint)
              ww[local_offset+qn_offset] = xx[local_offset+qn_offset] - f*yy[local_offset+qn_offset];
            if(qp_offset!=invalid_uint)
              ww[local_offset+qp_offset] = xx[local_offset+qp_offset] - f*yy[local_offset+qp_offset];
          }
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
            ww[local_offset+0] = xx[local_offset+0] - f*yy[local_offset+0];
          }
          break;
        }
      default: break;
      }
    }

    /*
    // only the last processor do this
    if(Genius::is_last_processor())
    {
      for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
      {
        BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
        unsigned int array_offset = bc->array_offset();
        if(array_offset != invalid_uint)
          yy[array_offset] *= f;
      }
    }
    */
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
void DGSolver::bank_rose_damping(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{
  return;
}



/*------------------------------------------------------------------
 * check for positive carrier density
 */
void DGSolver::check_positive_density(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{
  PetscScalar    *xx;
  PetscScalar    *ww;

  VecGetArray(x, &xx);
  VecGetArray(w, &ww);  // current candidate iterate

  int changed_flag=0;
  const PetscScalar T = this->get_system().T_external();
  const PetscScalar onePerCMC = 1.0*std::pow(cm,-3);

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

        //prevent negative carrier density
        if ( ww[local_offset+1] < 0 )
        { ww[local_offset+1] = 1e-2*fabs(xx[local_offset+1]) + onePerCMC; changed_flag++; }
        if ( ww[local_offset+2] < 0 )
        { ww[local_offset+2] = 1e-2*fabs(xx[local_offset+2]) + onePerCMC; changed_flag++; }

      }
    }
  }
  VecRestoreArray(x, &xx);
  VecRestoreArray(w, &ww);

  //synch changed_flag, if it is not zero, the vector is changed
  Parallel::sum( changed_flag );

  if(changed_flag)
  {
    *changed_w = PETSC_TRUE;
  }

  return;
}



void DGSolver::projection_positive_density_check(Vec x, Vec xo)
{
  PetscScalar    *xx;
  PetscScalar    *oo;

  VecGetArray(x, &xx);
  VecGetArray(xo, &oo);

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

        //prevent negative carrier density
        if ( xx[local_offset+1] < 0 )
          xx[local_offset+1]= fabs(0.01*oo[local_offset+1]);
        if ( xx[local_offset+2] < 0 )
          xx[local_offset+2]= fabs(0.01*oo[local_offset+2]);
      }
    }
  }

  VecRestoreArray(x,&xx);
  VecRestoreArray(xo,&oo);
}



/*------------------------------------------------------------------
 * test if BDF2 can be used for next time step
 */
bool DGSolver::BDF2_positive_defined() const
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
        if(failure_count) goto end;
      }
    }
  }
end:
  Parallel::sum(failure_count);
  return failure_count != 0;
}



/*------------------------------------------------------------------
 * evaluate local truncation error
 */
PetscReal DGSolver::LTE_norm()
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

  int N=0; //total variable number for LTE evaluation
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
        unsigned int qn_offset = invalid_uint;
        unsigned int qp_offset = invalid_uint;
        if(region->advanced_model().QNEnabled && !region->advanced_model().QPEnabled) qn_offset=3;
        if(!region->advanced_model().QNEnabled && region->advanced_model().QPEnabled) qp_offset=3;
        if(region->advanced_model().QNEnabled && region->advanced_model().QPEnabled)  {qn_offset=3; qp_offset=4;}

        SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
        SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
        for(; it!=it_end; ++it)
        {
          const FVM_Node * fvm_node = *it;
          unsigned int local_offset = fvm_node->local_offset();

          ll[local_offset+0] = 0;
          ll[local_offset+1] = std::abs(ll[local_offset+1])/(eps_r*std::abs(xx[local_offset+1])+eps_a*concentration);
          ll[local_offset+2] = std::abs(ll[local_offset+2])/(eps_r*std::abs(xx[local_offset+2])+eps_a*concentration);
          N += 2;

          if(qn_offset!=invalid_uint) ll[local_offset+qn_offset] = 0;
          if(qp_offset!=invalid_uint) ll[local_offset+qp_offset] = 0;
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

          ll[local_offset] = 0;
        }
        break;
      }
    case VacuumRegion:
      break;
    default:
      genius_error();
    }
  }


  // error estimate for each bc
  if(Genius::is_last_processor())
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      const BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
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



void DGSolver::error_norm()
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

  // scale the function vector
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
  elec_quantum_norm     = 0;
  hole_quantum_norm     = 0;

  poisson_norm              = 0;
  elec_continuity_norm      = 0;
  hole_continuity_norm      = 0;
  heat_equation_norm        = 0;
  elec_energy_equation_norm = 0;
  hole_energy_equation_norm = 0;
  elec_quantum_equation_norm= 0;
  hole_quantum_equation_norm= 0;
  electrode_norm            = 0;

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // only consider semiconductor region
    const SimulationRegion * region = _system.region(n);
    switch ( region->type() )
    {
    case SemiconductorRegion :
      {
        unsigned int qn_offset = invalid_uint;
        unsigned int qp_offset = invalid_uint;
        if(region->advanced_model().QNEnabled && !region->advanced_model().QPEnabled) qn_offset=3;
        if(!region->advanced_model().QNEnabled && region->advanced_model().QPEnabled) qp_offset=3;
        if(region->advanced_model().QNEnabled && region->advanced_model().QPEnabled)  {qn_offset=3; qp_offset=4;}

        SimulationRegion::const_processor_node_iterator it = region->on_processor_nodes_begin();
        SimulationRegion::const_processor_node_iterator it_end = region->on_processor_nodes_end();
        for(; it!=it_end; ++it)
        {
          const FVM_Node * fvm_node = *it;
          unsigned int offset = fvm_node->local_offset();

          potential_norm += xx[offset+0]*xx[offset+0];
          electron_norm  += xx[offset+1]*xx[offset+1];
          hole_norm      += xx[offset+2]*xx[offset+2];

          poisson_norm         += ff[offset+0]*ff[offset+0];
          elec_continuity_norm += ff[offset+1]*ff[offset+1];
          hole_continuity_norm += ff[offset+2]*ff[offset+2];

          if(qn_offset != invalid_uint)
          {
            elec_quantum_norm += xx[offset+qn_offset]*xx[offset+qn_offset];
            elec_quantum_equation_norm += ff[offset+qn_offset]*ff[offset+qn_offset];
          }
          if(qp_offset != invalid_uint)
          {
            hole_quantum_norm += xx[offset+qp_offset]*xx[offset+qp_offset];
            hole_quantum_equation_norm += ff[offset+qp_offset]*ff[offset+qp_offset];
          }
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
          unsigned int offset = fvm_node->local_offset();

          potential_norm += xx[offset]*xx[offset];
          poisson_norm   += ff[offset]*ff[offset];
        }
        break;
      }
    case VacuumRegion:
      break;
    default:
      genius_error(); //we should never reach here
    }
  }

  if(Genius::is_last_processor())
  {
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      const BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      unsigned int offset = bc->local_offset();
      if( offset != invalid_uint )
      {
        potential_norm += xx[offset]*xx[offset];

        PetscScalar scaling = 1.0;
        if(bc->is_electrode())
          scaling = bc->ext_circuit()->mna_scaling(SolverSpecify::dt);

        electrode_norm += ff[offset]*ff[offset]/(scaling*scaling+1e-6);
      }
    }
  }


  // sum of variable value on all processors
  std::vector<PetscScalar> norm_buffer;

  norm_buffer.push_back(potential_norm);
  norm_buffer.push_back(electron_norm);
  norm_buffer.push_back(hole_norm);
  norm_buffer.push_back(elec_quantum_norm);
  norm_buffer.push_back(hole_quantum_norm);

  norm_buffer.push_back(poisson_norm);
  norm_buffer.push_back(elec_continuity_norm);
  norm_buffer.push_back(hole_continuity_norm);
  norm_buffer.push_back(elec_quantum_equation_norm);
  norm_buffer.push_back(hole_quantum_equation_norm);
  norm_buffer.push_back(electrode_norm);

  Parallel::sum(norm_buffer);

  // sqrt to get L2 norm
  potential_norm = sqrt(norm_buffer[0]);
  electron_norm  = sqrt(norm_buffer[1]);
  hole_norm      = sqrt(norm_buffer[2]);
  elec_quantum_norm = sqrt(norm_buffer[3]);
  hole_quantum_norm = sqrt(norm_buffer[4]);

  poisson_norm         = sqrt(norm_buffer[5]);
  elec_continuity_norm = sqrt(norm_buffer[6]);
  hole_continuity_norm = sqrt(norm_buffer[7]);
  elec_quantum_equation_norm = sqrt(norm_buffer[8]);
  hole_quantum_equation_norm = sqrt(norm_buffer[9]);
  electrode_norm       = sqrt(norm_buffer[10]);


  VecRestoreArray(lx, &xx);
  VecRestoreArray(lf, &ff);

}




///////////////////////////////////////////////////////////////
// provide function and jacobian evaluation for QDDM solver  //
///////////////////////////////////////////////////////////////





/*------------------------------------------------------------------
 * evaluate the residual of function f at x
 */
void DGSolver::build_petsc_sens_residual(Vec x, Vec r)
{

  START_LOG("DGSolver_Residual()", "DGSolver");

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
    region->DG_Function(lxx, r, add_value_flag);
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate time derivative if necessary
  if(SolverSpecify::TimeDependent == true)
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DG_Time_Dependent_Function(lxx, r, add_value_flag);
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
    bc->DG_Function_Preprocess(lxx, r, src_row, dst_row, clear_row);
  }
  //add source rows to destination rows, and clear rows
  PetscUtils::VecAddClearRow(r, src_row, dst_row, clear_row);
  add_value_flag = NOT_SET_VALUES;

  // evaluate governing equations of DDML1 for all the boundaries
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DG_Function(lxx, r, add_value_flag);
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

  //VecView(r, PETSC_VIEWER_STDOUT_SELF);
  //getchar();

  STOP_LOG("DGSolver_Residual()", "DGSolver");

}



/*------------------------------------------------------------------
 * evaluate the Jacobian J of function f at x
 */
void DGSolver::build_petsc_sens_jacobian(Vec x, Mat *, Mat *)
{

  START_LOG("DGSolver_Jacobian()", "DGSolver");

  // scatte global solution vector x to local vector lx
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray(lx, &lxx);

  Jac->zero();

  START_LOG("DGSolver_Jacobian(R)", "DGSolver");

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

  // evaluate Jacobian matrix of governing equations of DDML1 in all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DG_Jacobian(lxx, Jac, add_value_flag);
  }


#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate Jacobian matrix of time derivative if necessary
  if(SolverSpecify::TimeDependent == true)
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DG_Time_Dependent_Jacobian(lxx, Jac, add_value_flag);
    }


  STOP_LOG("DGSolver_Jacobian(R)", "DGSolver");

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  START_LOG("DGSolver_Jacobian(B)", "DGSolver");
  // before first assemble, resereve none zero pattern for each boundary


#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  
  // assembly matrix
  Jac->close(false);


  std::vector<PetscInt> src_row,  dst_row,  clear_row;
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DG_Jacobian_Preprocess(lxx, Jac, src_row, dst_row, clear_row);
  }


  //add source rows to destination rows
  Jac->add_row_to_row(src_row, dst_row);
  // clear row
  Jac->clear_row(clear_row);


  add_value_flag = NOT_SET_VALUES;
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->DG_Jacobian(lxx, Jac, add_value_flag);
  }

  STOP_LOG("DGSolver_Jacobian(B)", "DGSolver");

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);

  Jac->close(true);

  //scaling the matrix
  MatDiagonalScale(J, L, PETSC_NULL);


  STOP_LOG("DGSolver_Jacobian()", "DGSolver");

}

void DGSolver::set_trace_electrode(BoundaryCondition *bc)
{
  genius_error();
}


