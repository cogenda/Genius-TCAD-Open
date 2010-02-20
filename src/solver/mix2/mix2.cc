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



#include "mix2/mix2.h"
#include "parallel.h"

#ifndef CYGWIN

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::cm;
using PhysicalUnit::K;


/*------------------------------------------------------------------
 * create nonlinear solver contex and adjust some parameters
 */
int Mix2Solver::create_solver()
{

  MESSAGE<< '\n' << "Mixed Simulation with DDM Level 2 init..." << std::endl;
  RECORD();

  return MixSolverBase::create_solver();
}




/*------------------------------------------------------------------
 * set initial value to solution vector and scaling vector
 */
int Mix2Solver::pre_solve_process(bool load_solution)
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
  }

  VecAssemblyBegin(x);
  VecAssemblyBegin(L);

  VecAssemblyEnd(x);
  VecAssemblyEnd(L);

  return MixSolverBase::pre_solve_process(load_solution);
}





/*------------------------------------------------------------------
 * the main solve routine. which is under the control of ngspice
 */
int Mix2Solver::solve()
{

  START_LOG("Mix2Solver_SNES()", "Mix2Solver");

  run_under_ngspice();

  STOP_LOG("Mix2Solver_SNES()", "Mix2Solver");

  return 0;
}




/*------------------------------------------------------------------
 * restore the solution to each region
 */
int Mix2Solver::post_solve_process()
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

  return MixSolverBase::post_solve_process();
}



/*------------------------------------------------------------------
 * get pdI/pdw, pdI/pdV and pdF/pdV for each electrode
 * here we use another G matrix instead of J matrix.
 * maybe we can do some fast computing with G
 */
int Mix2Solver::get_electrode_load()
{

  // scatte global solution vector x to local vector lx
  VecScatterBegin(scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (scatter, x, lx, INSERT_VALUES, SCATTER_FORWARD);

  PetscScalar *lxx;
  // get PetscScalar array contains solution from local solution vector lx
  VecGetArray(lx, &lxx);

  MatZeroEntries(G);

  // flag for indicate ADD_VALUES operator.
  InsertMode add_value_flag = NOT_SET_VALUES;

  // evaluate Jacobian matrix of governing equations of DDML2 in all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Jacobian(lxx, &G, add_value_flag);
  }

  // before first assemble, resereve none zero pattern for each boundary
  {
    if( !g_matrix_first_assemble )
      for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
      {
        BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
        bc->Mix_DDM2_Jacobian_Reserve(&G, add_value_flag);
      }
    g_matrix_first_assemble = true;
  }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
  {
    // assembly the matrix
    MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd  (G, MAT_FINAL_ASSEMBLY);

    // if scale matrix, then pdF_pdV should also be scaled!

    // scale the function vec
    PetscScalar *scale;

    // L is the scaling vector, the Jacobian evaluate function may dynamically update it.
    VecGetArray(L, &scale);

    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      if( !bc->is_electrode() ) continue;

      int pin = get_pin_index(bc->label());
      bc->Mix_DDM2_Electrode_Load(lx, &G, PINinfos[pin].I, PINconds[pin].pdI_pdV, PINconds[pin].pdI_pdw, PINconds[pin].pdF_pdV );

      PetscScalar *ff;

      // get function array
      VecGetArray(PINconds[pin].pdF_pdV, &ff);
      // scale it!
      for(unsigned int n=0; n<n_local_dofs; n++)
        ff[n] *= scale[n];

      VecRestoreArray(PINconds[pin].pdF_pdV, &ff);
    }

    // restore back

    VecRestoreArray(L, &scale);
  }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // evaluate Jacobian matrix of time derivative if necessary
  if(SolverSpecify::TimeDependent == true)
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      region->DDM2_Time_Dependent_Jacobian(lxx, &G, add_value_flag);
    }


  // evaluate Jacobian matrix of governing equations of Mixed type simulation of DDML2 for all the boundaries
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->Mix_DDM2_Jacobian(lxx, &G, add_value_flag);
  }


#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);

  // assembly the matrix
  MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (G, MAT_FINAL_ASSEMBLY);

  //scaling the matrix
  MatDiagonalScale(G, L, PETSC_NULL);

  return 0;
}




/*------------------------------------------------------------------
 * load previous state into solution vector
 */
int Mix2Solver::diverged_recovery()
{
  // for all the regions
  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    SimulationRegion * region = _system.region(n);
    region->DDM2_Fill_Value(x, L);
  }

  VecAssemblyBegin(x);
  VecAssemblyBegin(L);

  VecAssemblyEnd(x);
  VecAssemblyEnd(L);

  return 0;
}


/*------------------------------------------------------------------
 * Save steady-state solution as the previous solution data in spice
 */
int Mix2Solver::init_spice_data()
{

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // consider semiconductor region
    SimulationRegion * region = _system.region(n);
    switch( region->type() )
    {
    case  SemiconductorRegion:
      {
        SimulationRegion::node_iterator it = region->nodes_begin();
        SimulationRegion::node_iterator it_end = region->nodes_end();
        for(; it!=it_end; ++it)
        {
          FVM_Node * fvm_node = (*it).second;
          FVM_NodeData * node_data = fvm_node->node_data();

          //if this node NOT belongs to this processor, continue
          if( !fvm_node->on_processor() ) continue;

          node_data->CreateUserScalarValue("SPICE_psi");
          node_data->CreateUserScalarValue("SPICE_n");
          node_data->CreateUserScalarValue("SPICE_p");
          node_data->CreateUserScalarValue("SPICE_psi_last");
          node_data->CreateUserScalarValue("SPICE_n_last");
          node_data->CreateUserScalarValue("SPICE_p_last");
          node_data->UserScalarValue("SPICE_psi") = node_data->psi();
          node_data->UserScalarValue("SPICE_n")   = node_data->n();
          node_data->UserScalarValue("SPICE_p")   = node_data->p();
          node_data->UserScalarValue("SPICE_psi_last") = node_data->psi();
          node_data->UserScalarValue("SPICE_n_last")   = node_data->n();
          node_data->UserScalarValue("SPICE_p_last")   = node_data->p();

          node_data->CreateUserScalarValue("SPICE_T");
          node_data->UserScalarValue("SPICE_T") = node_data->T();
          node_data->CreateUserScalarValue("SPICE_T_last");
          node_data->UserScalarValue("SPICE_T_last") = node_data->T();
        }
        break;
      }
    case  InsulatorRegion:
    case  ConductorRegion:
      {
        SimulationRegion::node_iterator it = region->nodes_begin();
        SimulationRegion::node_iterator it_end = region->nodes_end();

        for(; it!=it_end; ++it)
        {
          FVM_Node * fvm_node = (*it).second;
          FVM_NodeData * node_data = fvm_node->node_data();

          //if this node NOT belongs to this processor, continue
          if( !fvm_node->on_processor() ) continue;

          node_data->CreateUserScalarValue("SPICE_psi");
          node_data->UserScalarValue("SPICE_psi") = node_data->psi();
          node_data->CreateUserScalarValue("SPICE_psi_last");
          node_data->UserScalarValue("SPICE_psi_last") = node_data->psi();

          node_data->CreateUserScalarValue("SPICE_T");
          node_data->UserScalarValue("SPICE_T") = node_data->T();
          node_data->CreateUserScalarValue("SPICE_T_last");
          node_data->UserScalarValue("SPICE_T_last") = node_data->T();
        }
        break;
      }
    default: break;
    }
  }

  return 0;
}

/*------------------------------------------------------------------
 * load solution data previously accepted by spice
 */
int Mix2Solver::load_spice_data()
{

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // consider semiconductor region
    SimulationRegion * region = _system.region(n);
    switch( region->type() )
    {
    case  SemiconductorRegion:
      {
        SimulationRegion::node_iterator it = region->nodes_begin();
        SimulationRegion::node_iterator it_end = region->nodes_end();
        for(; it!=it_end; ++it)
        {
          FVM_Node * fvm_node = (*it).second;
          FVM_NodeData * node_data = fvm_node->node_data();

          //if this node NOT belongs to this processor, continue
          if( !fvm_node->on_processor() ) continue;

          node_data->psi()       = node_data->UserScalarValue("SPICE_psi");
          node_data->n()         = node_data->UserScalarValue("SPICE_n");
          node_data->p()         = node_data->UserScalarValue("SPICE_p");
          node_data->psi_last()  = node_data->UserScalarValue("SPICE_psi_last");
          node_data->n_last()    = node_data->UserScalarValue("SPICE_n_last");
          node_data->p_last()    = node_data->UserScalarValue("SPICE_p_last");

          node_data->T()       = node_data->UserScalarValue("SPICE_T");
          node_data->T_last()  = node_data->UserScalarValue("SPICE_T_last");
        }
        break;
      }
    case  InsulatorRegion:
    case  ConductorRegion:
      {
        SimulationRegion::node_iterator it = region->nodes_begin();
        SimulationRegion::node_iterator it_end = region->nodes_end();

        for(; it!=it_end; ++it)
        {
          FVM_Node * fvm_node = (*it).second;
          FVM_NodeData * node_data = fvm_node->node_data();

          //if this node NOT belongs to this processor, continue
          if( !fvm_node->on_processor() ) continue;

          node_data->psi()      = node_data->UserScalarValue("SPICE_psi");
          node_data->psi_last() = node_data->UserScalarValue("SPICE_psi_last");

          node_data->T()      = node_data->UserScalarValue("SPICE_T");
          node_data->T()      = node_data->UserScalarValue("SPICE_T_last");
        }
        break;
      }
    default: break;
    }
  }

  return 0;

}


/*------------------------------------------------------------------
 * since spice accept the solution, save solution data
 */
int Mix2Solver::save_spice_data()
{

  for(unsigned int n=0; n<_system.n_regions(); n++)
  {
    // consider semiconductor region
    SimulationRegion * region = _system.region(n);
    switch( region->type() )
    {
    case  SemiconductorRegion:
      {
        SimulationRegion::node_iterator it = region->nodes_begin();
        SimulationRegion::node_iterator it_end = region->nodes_end();
        for(; it!=it_end; ++it)
        {
          FVM_Node * fvm_node = (*it).second;
          FVM_NodeData * node_data = fvm_node->node_data();

          //if this node NOT belongs to this processor, continue
          if( !fvm_node->on_processor() ) continue;

          node_data->UserScalarValue("SPICE_psi")      = node_data->psi();
          node_data->UserScalarValue("SPICE_n")        = node_data->n();
          node_data->UserScalarValue("SPICE_p")        = node_data->p();
          node_data->UserScalarValue("SPICE_psi_last") = node_data->psi_last();
          node_data->UserScalarValue("SPICE_n_last")   = node_data->n_last();
          node_data->UserScalarValue("SPICE_p_last")   = node_data->p_last();

          node_data->UserScalarValue("SPICE_T")        = node_data->T();
          node_data->UserScalarValue("SPICE_T_last")   = node_data->T_last();
        }
        break;
      }
    case  InsulatorRegion:
    case  ConductorRegion:
      {
        SimulationRegion::node_iterator it = region->nodes_begin();
        SimulationRegion::node_iterator it_end = region->nodes_end();

        for(; it!=it_end; ++it)
        {
          FVM_Node * fvm_node = (*it).second;
          FVM_NodeData * node_data = fvm_node->node_data();

          //if this node NOT belongs to this processor, continue
          if( !fvm_node->on_processor() ) continue;

          node_data->UserScalarValue("SPICE_psi")      = node_data->psi();
          node_data->UserScalarValue("SPICE_psi_last") = node_data->psi_last();

          node_data->UserScalarValue("SPICE_T")        = node_data->T();
          node_data->UserScalarValue("SPICE_T_last")   = node_data->T_last();
        }
        break;
      }
    default: break;
    }
  }

  return 0;

}


/*------------------------------------------------------------------
 * Potential Newton Damping
 */
void Mix2Solver::potential_damping(Vec x, Vec y, Vec w, PetscTruth *changed_y, PetscTruth *changed_w)
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
      if( region->type() != SemiconductorRegion ) continue;

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
void Mix2Solver::bank_rose_damping(Vec , Vec , Vec , PetscTruth *changed_y, PetscTruth *changed_w)
{
  *changed_y = PETSC_FALSE;
  *changed_w = PETSC_FALSE;

  return;
}



/*------------------------------------------------------------------
 * positive density Newton Damping
 */
void Mix2Solver::positive_density_damping(Vec x, Vec y, Vec w, PetscTruth *changed_y, PetscTruth *changed_w)
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
    if ( fabs(yy[local_offset]) > 1.0 )  { ww[local_offset] = xx[local_offset] - std::sign(yy[local_offset])*1.0; changed_flag = 1; }

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



void Mix2Solver::projection_positive_density_check(Vec x, Vec xo)
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



#if PETSC_VERSION_LE(2,3,3)
#include "private/snesimpl.h"
#endif
void Mix2Solver::error_norm()
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
  elec_continuity_norm       = 0;
  hole_continuity_norm       = 0;
  heat_equation_norm        = 0;
  elec_energy_equation_norm = 0;
  hole_energy_equation_norm = 0;


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

          poisson_norm        += ff[offset+0]*ff[offset+0];
          elec_continuity_norm += ff[offset+1]*ff[offset+1];
          hole_continuity_norm += ff[offset+2]*ff[offset+2];
          heat_equation_norm  += ff[offset+3]*ff[offset+3];
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

  VecRestoreArray(lx, &xx);
  VecRestoreArray(lf, &ff);

}






///////////////////////////////////////////////////////////////
// provide function and jacobian evaluation for DDML2 solver //
///////////////////////////////////////////////////////////////





/*------------------------------------------------------------------
 * evaluate the residual of function f at x
 */
void Mix2Solver::build_petsc_sens_residual(Vec x, Vec r)
{

  START_LOG("Mix2Solver_Residual()", "Mix2Solver");

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

  // evaluate governing equations of Mixed type simulation of DDML2 for all the boundaries
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->Mix_DDM2_Function(lxx, r, add_value_flag);
  }


#ifdef HAVE_FENV_H
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

  STOP_LOG("Mix2Solver_Residual()", "Mix2Solver");

}




/*------------------------------------------------------------------
 * evaluate the Jacobian J of function f at x
 */
void Mix2Solver::build_petsc_sens_jacobian(Vec x, Mat *, Mat *)
{

  START_LOG("Mix2Solver_Jacobian()", "Mix2Solver");

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

  // before first assemble, resereve none zero pattern for each boundary
  if( !jacobian_matrix_first_assemble )
  {
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      bc->Mix_DDM2_Jacobian_Reserve(&J, add_value_flag);
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

  // evaluate Jacobian matrix of governing equations of Mixed type simulation of DDML2 for all the boundaries
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
    bc->Mix_DDM2_Jacobian(lxx, &J, add_value_flag);
  }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // restore array back to Vec
  VecRestoreArray(lx, &lxx);

  // assembly the matrix
  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (J, MAT_FINAL_ASSEMBLY);

  //scaling the matrix
  MatDiagonalScale(J, L, PETSC_NULL);


  //MatView(J,PETSC_VIEWER_STDOUT_SELF );
  //getchar();


  STOP_LOG("Mix2Solver_Jacobian()", "Mix2Solver");

}

#endif //#ifndef CYGWIN
