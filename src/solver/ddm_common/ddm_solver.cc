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

//  $Id: ddm_solver.cc,v 1.11 2008/07/09 05:58:16 gdiso Exp $
#include <iomanip>
#include <stack>

#include "solver_specify.h"
#include "physical_unit.h"
#include "electrical_source.h"
#include "resistance_region.h"
#include "simulation_system.h"
#include "field_source.h"
#include "ddm_solver.h"
#include "parallel.h"
#include "MXMLUtil.h"


using PhysicalUnit::A;
using PhysicalUnit::V;
using PhysicalUnit::W;
using PhysicalUnit::C;
using PhysicalUnit::s;
using PhysicalUnit::um;

DDMSolverBase::DDMSolverBase(SimulationSystem & system): FVM_NonlinearSolver(system)
{
  // do clear
  potential_norm            = 0.0;
  electron_norm             = 0.0;
  hole_norm                 = 0.0;
  temperature_norm          = 0.0;
  elec_temperature_norm     = 0.0;
  hole_temperature_norm     = 0.0;
  elec_quantum_norm         = 0.0;
  hole_quantum_norm         = 0.0;

  poisson_norm              = 0.0;
  elec_continuity_norm      = 0.0;
  hole_continuity_norm      = 0.0;
  heat_equation_norm        = 0.0;
  elec_energy_equation_norm = 0.0;
  hole_energy_equation_norm = 0.0;
  elec_quantum_equation_norm= 0.0;
  hole_quantum_equation_norm= 0.0;
  electrode_norm            = 0.0;

  function_norm             = 0.0;
  functions_norm.resize(9, 0.0);
  nonlinear_iteration       = 0;
}

int DDMSolverBase::create_solver()
{
  set_nonlinear_solver_type ( SolverSpecify::NS );
  set_linear_solver_type    ( SolverSpecify::LS );
  set_preconditioner_type   ( SolverSpecify::PC );

  // must setup nonlinear contex here!
  setup_nonlinear_data();

  //NOTE Tolerances here only be set as a reference

  //abstol = 1e-15                  - absolute convergence tolerance
  //rtol   = 1e-14                  - relative convergence tolerance
  //stol   = 1e-9                   - convergence tolerance in terms of the norm of the change in the solution between steps
  SNESSetTolerances(snes, 1e-15, 1e-14, 1e-9, SolverSpecify::MaxIteration, 1000);

  // rtol   = SolverSpecify::ksp_rtol  - the relative convergence tolerance (relative decrease in the residual norm)
  // abstol = 1e-20                    - the absolute convergence tolerance (absolute size of the residual norm)
  KSPSetTolerances(ksp, SolverSpecify::ksp_rtol, 1e-20, PETSC_DEFAULT, std::max(50, std::min(1000, static_cast<int>(n_global_dofs/10))) );

  // user can do further adjusment from command line
  SNESSetFromOptions (snes);

  return FVM_NonlinearSolver::create_solver();
}



void DDMSolverBase::set_extra_matrix_nonzero_pattern()
{
  if(_system.get_bcs()!=NULL)
  {
    for(unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n )
    {
      const BoundaryCondition * bc = _system.get_bcs()->get_bc(n);
      std::map<FVM_Node *, std::pair<unsigned int, unsigned int> > node_extra_dofs;
      bc->DDM_extra_dofs(node_extra_dofs);

      std::map<FVM_Node *, std::pair<unsigned int, unsigned int> >::const_iterator it = node_extra_dofs.begin();
      for(; it != node_extra_dofs.begin(); it++)
      {
        const FVM_Node * fvm_node = it->first;
        const SimulationRegion * region = _system.region(fvm_node->subdomain_id());
        std::pair<unsigned int, unsigned int> dofs = it->second;

        unsigned int local_offset = fvm_node->local_offset();
        unsigned int local_node_dofs = this->node_dofs( region );

        for(unsigned int i=0; i<local_node_dofs; ++i)
        {
          n_nz[local_offset + i] += dofs.first;
          n_oz[local_offset + i] += dofs.second;
        }
      }
    }
  }
}


int DDMSolverBase::post_solve_process()
{
  mxml_node_t *eSolution = new_dom_solution_elem();
  const BoundaryConditionCollector * bcs = get_system().get_bcs();

  if ( Genius::processor_id() == 0)
  {
    mxml_node_t *eLabel = mxmlNewElement(eSolution, "label");
    switch (SolverSpecify::Type)
    {
        case SolverSpecify::TRANSIENT:
        {
          mxml_node_t *eTime = mxmlNewElement(eLabel, "time");
          mxmlAdd(eTime, MXML_ADD_AFTER, NULL, MXMLQVariant::makeQVFloat(SolverSpecify::clock / PhysicalUnit::s));
          break;
        }
        case SolverSpecify::DCSWEEP:
        case SolverSpecify::TRACE:
        {
          if( this->solver_type() == SolverSpecify::DDML1 ||
              this->solver_type() == SolverSpecify::DDML2 ||
              this->solver_type() == SolverSpecify::EBML3 )
          {
            if (!SolverSpecify::Electrode_VScan.empty())
            {
              mxml_node_t *eVolt = mxmlNewElement(eLabel, "voltage");
              std::string bcLabel = SolverSpecify::Electrode_VScan[0];
              double volt = bcs->get_bc(bcLabel)->ext_circuit()->Vapp()/PhysicalUnit::V;
              mxmlAdd(eVolt, MXML_ADD_AFTER, NULL, MXMLQVariant::makeQVFloat(volt));
            }
            else if (!SolverSpecify::Electrode_IScan.empty())
            {
              mxml_node_t *eCurr = mxmlNewElement(eLabel, "current");
              std::string bcLabel = SolverSpecify::Electrode_IScan[0];
              double curr = bcs->get_bc(bcLabel)->ext_circuit()->Iapp()/PhysicalUnit::A;
              mxmlAdd(eCurr, MXML_ADD_AFTER, NULL, MXMLQVariant::makeQVFloat(curr));
            }
          }
          break;
        }
        default:
        mxmlAdd(eLabel, MXML_ADD_AFTER, NULL, MXMLQVariant::makeQVString("solution"));
    }

    mxml_node_t *eTerm = bcs->get_dom_terminal_info();
    if (eTerm)
      mxmlAdd(eSolution, MXML_ADD_AFTER, NULL, eTerm);
  }

  return FVM_NonlinearSolver::post_solve_process();
}

int DDMSolverBase::destroy_solver()
{
  // clear nonlinear matrix/vector
  clear_nonlinear_data();

#if defined(HAVE_FENV_H)
  feclearexcept(FE_INVALID);
#endif

  return FVM_NonlinearSolver::destroy_solver();
}



/* ----------------------------------------------------------------------------
 * compute equilibrium state
 * all the stimulate source(s) are set to zero. time step set to inf
 */
int DDMSolverBase::solve_equ()
{
  MESSAGE<<"Compute equilibrium\n";
  RECORD();

  // set electrode with transient time 0 value of stimulate source(s)
  _system.get_electrical_source()->update ( 0 );
  // set all the electrode to ground
  _system.get_bcs()->all_electrode_ground();

  // not time dependent
  SolverSpecify::TimeDependent = false;
  SolverSpecify::dt = 1e100;
  SolverSpecify::clock = 0.0;

  // call pre_solve_process
  this->pre_solve_process();

  // here call Petsc to solve the nonlinear equations
  sens_solve();

  // call post_solve_process
  this->post_solve_process();

  // get the converged reason
  SNESConvergedReason reason;
  SNESGetConvergedReason ( snes, &reason );

  // linear solver iteration
  PetscInt lits;
  SNESGetLinearSolveIterations(snes, &lits);

  // print convergence/divergence reason
  {
    MESSAGE
    <<"--------------------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits << "\n\n\n";
    RECORD();
  }

  SolverSpecify::tran_histroy = false;

  return 0;
}



/* ----------------------------------------------------------------------------
 * compute steadystate
 * all the stimulate source(s) are set with transient time 0 value. time step set to inf
 */
int DDMSolverBase::solve_steadystate()
{
  MESSAGE<<"Compute steady-state\n";
  RECORD();

  // set electrode with transient time 0 value of stimulate source(s)
  _system.get_electrical_source()->update ( 0 );
  _system.get_field_source()->update ( 0 );

  // not time dependent
  SolverSpecify::TimeDependent = false;
  SolverSpecify::dt = 1e100;
  SolverSpecify::clock = 0.0;

  if(SolverSpecify::PseudoTimeMethod)
    snes_solve_pseudo_time_step();
  else
  {
    // call pre_solve_process
    this->pre_solve_process();

    // here call Petsc to solve the nonlinear equations
    sens_solve();

    // call post_solve_process
    this->post_solve_process();

    // get the converged reason
    SNESConvergedReason reason;
    SNESGetConvergedReason ( snes, &reason );

    // linear solver iteration
    PetscInt lits;
    SNESGetLinearSolveIterations(snes, &lits);

    // print convergence/divergence reason
    {
      MESSAGE <<"--------------------------------------------------------------------------------\n"
              <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits << "\n\n\n";
      RECORD();
    }
  }

  SolverSpecify::tran_histroy = false;

  return 0;
}



/* ----------------------------------------------------------------------------
 * compute dcsweep, sweep V or I for one electrode and get the device IV curve.
 * stimulate source(s) for other electrode are set with transient time 0 value.
 * time step set to inf
 */
int DDMSolverBase::solve_dcsweep()
{

  // set electrode with transient time 0 value of stimulate source(s)
  _system.get_electrical_source()->update ( 0 );

  // not time dependent
  SolverSpecify::TimeDependent = false;
  SolverSpecify::dt = 1e100;
  SolverSpecify::clock = 0.0;


  // output DC Scan information
  if ( SolverSpecify::Electrode_VScan.size() )
  {

    MESSAGE
    <<"DC voltage scan from "  <<SolverSpecify::VStart/PhysicalUnit::V
    <<" step "                 <<SolverSpecify::VStep/PhysicalUnit::V
    <<" to "                   <<SolverSpecify::VStop/PhysicalUnit::V
    <<'\n';
    RECORD();
  }
  else
  {
    MESSAGE
    <<"DC current scan from " << SolverSpecify::IStart/PhysicalUnit::A
    <<" step "                << SolverSpecify::IStep/PhysicalUnit::A
    <<" to "                  << SolverSpecify::IStop/PhysicalUnit::A
    <<'\n';
    RECORD();
  }


  PetscInt total_lits = 0;

  // voltage scan
  if ( SolverSpecify::Electrode_VScan.size() )
  {
    // the current vscan voltage
    PetscScalar Vscan = SolverSpecify::VStart;

    // the current vscan step
    PetscScalar VStep = SolverSpecify::VStep;

    // saved solutions and vscan values for solution projection.
    Vec xs1, xs2, xs3;
    PetscScalar Vs1=Vscan, Vs2=Vscan, Vs3=Vscan;
    std::stack<PetscScalar> V_retry;
    VecDuplicate ( x,&xs1 );
    VecDuplicate ( x,&xs2 );
    VecDuplicate ( x,&xs3 );

    // main loop
    for ( SolverSpecify::DC_Cycles=0;  (Vscan*SolverSpecify::VStep) <= SolverSpecify::VStop*SolverSpecify::VStep* ( 1.0+1e-7 ); )
    {
      // show current vscan value
      MESSAGE << "DC Scan: V("  << SolverSpecify::Electrode_VScan[0];
      for ( unsigned int i=1; i<SolverSpecify::Electrode_VScan.size(); i++ )
        MESSAGE << ", "  << SolverSpecify::Electrode_VScan[i];
      MESSAGE << ") = "  << Vscan/PhysicalUnit::V  <<" V" << '\n'
      <<"--------------------------------------------------------------------------------\n";
      RECORD();

      // set current vscan voltage to corresponding electrode
      _system.get_electrical_source()->assign_voltage_to ( SolverSpecify::Electrode_VScan, Vscan );
      _system.get_field_source()->update ( 0, SolverSpecify::SourceCoupled );

      // call pre_solve_process
      if ( SolverSpecify::DC_Cycles == 0 )
        this->pre_solve_process();
      else
        this->pre_solve_process ( false );

      // here call Petsc to solve the nonlinear equations
      sens_solve();

      // get the converged reason
      SNESConvergedReason reason;
      SNESGetConvergedReason ( snes,&reason );

      // linear solver iteration
      PetscInt lits;
      SNESGetLinearSolveIterations(snes, &lits);
      total_lits += lits;

      if ( reason>0 ) //ok, converged.
      {

        // call post_solve_process
        this->post_solve_process();

        SolverSpecify::DC_Cycles++;

        // save solution for linear/quadratic projection
        Vs3=Vs2;
        Vs2=Vs1;
        Vs1=Vscan;

        VecCopy ( xs2,xs3 );
        VecCopy ( xs1,xs2 );
        VecCopy ( x,xs1 );

        if ( V_retry.empty() )
        {
          // add vstep to current voltage
          Vscan += VStep;
        }
        else
        {
          // pop
          Vscan = V_retry.top();
          V_retry.pop();
        }

        if ( fabs ( Vscan-SolverSpecify::VStop ) <1e-10 )
          Vscan=SolverSpecify::VStop;

        // if v step small than VStepMax, mult by factor of 1.1
        if ( fabs ( VStep ) < fabs ( SolverSpecify::VStepMax ) )  VStep *= 1.1;


        // however, for last step, we force V equal to VStop
        if ( (Vscan*SolverSpecify::VStep) > SolverSpecify::VStop*SolverSpecify::VStep &&
             (Vscan*SolverSpecify::VStep) < ( SolverSpecify::VStop + VStep - 1e-10*VStep ) *SolverSpecify::VStep
           )
          Vscan = SolverSpecify::VStop;


        MESSAGE
        <<"--------------------------------------------------------------------------------\n"
        <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits << "\n\n\n";
        RECORD();
      }
      else // oh, diverged... reduce step and try again
      {

        if(reason == SNES_DIVERGED_LINEAR_SOLVE)
        {
          KSPConvergedReason ksp_reason;
          KSPGetConvergedReason ( ksp, &ksp_reason );
          MESSAGE <<"------> linear solver "<<KSPConvergedReasons[ksp_reason];
        }
        else
          MESSAGE <<"------> nonlinear solver "<<SNESConvergedReasons[reason];

        // failed in the first step, we didn't know how to set the scan bias
        if ( SolverSpecify::DC_Cycles == 0 )
        {
          MESSAGE <<". Failed in the first step.\n\n\n";
          RECORD();
          break;
        }

        if ( V_retry.size() >=8 )
        {
          MESSAGE <<". Too many failed steps, give up tring.\n\n\n";
          RECORD();
          break;
        }

        MESSAGE <<", do recovery...\n\n\n"; RECORD();

        // load previous result into solution vector
        this->diverged_recovery();

        // reduce step by a factor of 2
        V_retry.push ( Vscan );
        Vscan= ( Vscan+Vs1 ) /2.0;

      }

      if ( SolverSpecify::Predict )
      {
        PetscScalar hn = Vscan-Vs1;
        PetscScalar hn1 = Vs1-Vs2;
        PetscScalar hn2 = Vs2-Vs3;

        if ( SolverSpecify::DC_Cycles>=3 )
        {
          // quadradic projection
          PetscScalar cn=hn* ( hn+2*hn1+hn2 ) / ( hn1* ( hn1+hn2 ) );
          PetscScalar cn1=-hn* ( hn+hn1+hn2 ) / ( hn1*hn2 );
          PetscScalar cn2=hn* ( hn+hn1 ) / ( hn2* ( hn1+hn2 ) );

          VecAXPY ( x,cn,xs1 );
          VecAXPY ( x,cn1,xs2 );
          VecAXPY ( x,cn2,xs3 );
          this->projection_positive_density_check ( x,xs1 );
        }
        else if ( SolverSpecify::DC_Cycles>=2 )
        {
          // linear projection
          VecAXPY ( x, hn/hn1,xs1 );
          VecAXPY ( x,-hn/hn1,xs2 );
          this->projection_positive_density_check ( x,xs1 );
        }
      }
    }

    VecDestroy ( PetscDestroyObject(xs1) );
    VecDestroy ( PetscDestroyObject(xs2) );
    VecDestroy ( PetscDestroyObject(xs3) );

  }



  // current scan
  if ( SolverSpecify::Electrode_IScan.size() )
  {
    // iscan current
    PetscScalar Iscan = SolverSpecify::IStart;

    // iscan step
    PetscScalar IStep = SolverSpecify::IStep;

    // saved solutions and iscan values for solution projection.
    Vec xs1, xs2, xs3;
    PetscScalar Is1=Iscan, Is2=Iscan, Is3=Iscan;
    std::stack<PetscScalar> I_retry;
    VecDuplicate ( x,&xs1 );
    VecDuplicate ( x,&xs2 );
    VecDuplicate ( x,&xs3 );

    // main loop
    for ( SolverSpecify::DC_Cycles=0;  (Iscan*SolverSpecify::IStep) <= SolverSpecify::IStop*SolverSpecify::IStep* ( 1.0+1e-7 ); )
    {
      // show current iscan value
      MESSAGE << "DC Scan: I("  << SolverSpecify::Electrode_IScan[0];
      for ( unsigned int i=1; i<SolverSpecify::Electrode_IScan.size(); i++ )
        MESSAGE << ", "  << SolverSpecify::Electrode_IScan[i];
      MESSAGE << ") = "  << Iscan/PhysicalUnit::A  <<" A" << '\n'
      <<"--------------------------------------------------------------------------------\n";
      RECORD();

      // set iscan current to corresponding electrode
      _system.get_electrical_source()->assign_current_to ( SolverSpecify::Electrode_IScan, Iscan );
      _system.get_field_source()->update ( 0, SolverSpecify::SourceCoupled );

      // call pre_solve_process
      if ( SolverSpecify::DC_Cycles == 0 )
        this->pre_solve_process();
      else
        this->pre_solve_process ( false );

      sens_solve();
      // get the converged reason
      SNESConvergedReason reason;
      SNESGetConvergedReason ( snes,&reason );

      // linear solver iteration
      PetscInt lits;
      SNESGetLinearSolveIterations(snes, &lits);
      total_lits += lits;

      if ( reason>0 ) //ok, converged.
      {

        // call post_solve_process
        this->post_solve_process();

        SolverSpecify::DC_Cycles++;

        // save solution for linear/quadratic projection
        Is3=Is2;
        Is2=Is1;
        Is1=Iscan;

        VecCopy ( xs2,xs3 );
        VecCopy ( xs1,xs2 );
        VecCopy ( x,xs1 );

        if ( I_retry.empty() )
        {
          // add vstep to current voltage
          Iscan += IStep;
        }
        else
        {
          // pop
          Iscan = I_retry.top();
          I_retry.pop();
        }

        if ( fabs ( Iscan-SolverSpecify::IStop ) <1e-10 )
          Iscan=SolverSpecify::IStop;

        // if I step small than IStepMax, mult by factor of 1.1
        if ( fabs ( IStep ) < fabs ( SolverSpecify::IStepMax ) )  IStep *= 1.1;


        // however, for last step, we force I equal to IStop
        if ( (Iscan*SolverSpecify::IStep) > SolverSpecify::IStop*SolverSpecify::IStep &&
             (Iscan*SolverSpecify::IStep) < ( SolverSpecify::IStop + IStep - 1e-10*IStep ) *SolverSpecify::IStep
           )
          Iscan = SolverSpecify::IStop;

        MESSAGE
        <<"--------------------------------------------------------------------------------\n"
        <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits << "\n\n\n";
        RECORD();
      }
      else // oh, diverged... reduce step and try again
      {
        if(reason == SNES_DIVERGED_LINEAR_SOLVE)
        {
          KSPConvergedReason ksp_reason;
          KSPGetConvergedReason ( ksp, &ksp_reason );
          MESSAGE <<"------> linear solver "<<KSPConvergedReasons[ksp_reason];
        }
        else
          MESSAGE <<"------> nonlinear solver "<<SNESConvergedReasons[reason];

        if ( SolverSpecify::DC_Cycles == 0 )
        {
          MESSAGE <<". Failed in the first step.\n\n\n";
          RECORD();
          break;
        }
        if ( I_retry.size() >=8 )
        {
          MESSAGE <<". Too many failed steps, give up tring.\n\n\n";
          RECORD();
          break;
        }

        MESSAGE <<", do recovery...\n\n\n"; RECORD();

        // load previous result into solution vector
        this->diverged_recovery();

        // reduce step by a factor of 2
        I_retry.push ( Iscan );
        Iscan= ( Iscan+Is1 ) /2.0;


      }

      if ( SolverSpecify::Predict )
      {
        PetscScalar hn = Iscan-Is1;
        PetscScalar hn1 = Is1-Is2;
        PetscScalar hn2 = Is2-Is3;

        if ( SolverSpecify::DC_Cycles>=3 )
        {
          // quadradic projection
          PetscScalar cn=hn* ( hn+2*hn1+hn2 ) / ( hn1* ( hn1+hn2 ) );
          PetscScalar cn1=-hn* ( hn+hn1+hn2 ) / ( hn1*hn2 );
          PetscScalar cn2=hn* ( hn+hn1 ) / ( hn2* ( hn1+hn2 ) );

          VecAXPY ( x,cn,xs1 );
          VecAXPY ( x,cn1,xs2 );
          VecAXPY ( x,cn2,xs3 );
          this->projection_positive_density_check ( x,xs1 );
        }
        else if ( SolverSpecify::DC_Cycles>=2 )
        {
          // linear projection
          VecAXPY ( x, hn/hn1,xs1 );
          VecAXPY ( x,-hn/hn1,xs2 );
          this->projection_positive_density_check ( x,xs1 );
        }
      }
    }

    VecDestroy ( PetscDestroyObject(xs1) );
    VecDestroy ( PetscDestroyObject(xs2) );
    VecDestroy ( PetscDestroyObject(xs3) );
  }


  SolverSpecify::tran_histroy = false;

  return 0;
}


int DDMSolverBase::solve_op()
{
  // set electrode with transient time 0 value of stimulate source(s)
  _system.get_field_source()->update ( 0 );

  // not time dependent
  SolverSpecify::TimeDependent = false;
  SolverSpecify::dt = 1e100;
  SolverSpecify::clock = 0.0;

  _system.get_electrical_source()->save_bc_source_state();

  const double dtao_init_potential = SolverSpecify::PseudoTimeStepPotential;
  const double dtao_init_carrier =  SolverSpecify::PseudoTimeStepCarrier;
  const double dtao_init_metal =  SolverSpecify::PseudoTimeStepMetal;

  // set sigma to a pseudo value for CMOS device
  if(SolverSpecify::PseudoTimeMethod && SolverSpecify::PseudoTimeCMOS)
  {
    const double sigma = sqrt(SolverSpecify::PseudoTimeCMOSCap/SolverSpecify::PseudoTimeCMOSRes/SolverSpecify::PseudoTimeCMOSTime)/SolverSpecify::PseudoTimeCMOSLambda;

    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      if( region->type() == MetalRegion )
      {
        MetalSimulationRegion * resistance_region = dynamic_cast<MetalSimulationRegion *>(region);
        resistance_region->find_low_resistance_solderpad();
        if( !resistance_region->connect_to_low_resistance_solderpad() )
        {
          region->set_conductance(sigma);
        }
      }
    }
  }

  // time step counter
  int t_cycles=0;
  // diverged counter
  int diverged_retry=0;

  const double ramp_steps = _system.get_electrical_source()->steps_by_limiter(SolverSpecify::VStepMax, SolverSpecify::IStepMax, 0.0);
  double step = 1.0;
  double dstep = 1.0;
  // the main loop of transient solver.
  do
  {
    if(SolverSpecify::PseudoTimeMethod)
    {
      MESSAGE <<"Op Step PseudoTime "<< step <<" of "<< static_cast<int>(ramp_steps) << '\n'
              <<"--------------------------------------------------------------------------------\n";
      RECORD();
    }
    else
    {
      MESSAGE <<"Op Step "<< step <<" of "<< static_cast<int>(ramp_steps) << '\n'
          <<"--------------------------------------------------------------------------------\n";
      RECORD();
    }

    //update sources to current step
    _system.get_electrical_source()->rampup(step/ramp_steps, 0.0);

    //we do solve here!

    // call pre_solve_process
    if ( t_cycles == 0 )
      this->pre_solve_process();
    else
      this->pre_solve_process ( false );

    sens_solve();
    // get the converged reason
    SNESConvergedReason reason;
    SNESGetConvergedReason ( snes,&reason );

    // linear solver iteration
    PetscInt lits;
    SNESGetLinearSolveIterations(snes, &lits);

    //nonlinear solution diverged? try to do recovery
    if ( reason<0 )
    {
      // increase diverged_retry
      diverged_retry++;

      if ( diverged_retry >= 8 ) //failed 8 times, stop tring
      {
        MESSAGE<<"------> Too many failed steps, give up tring.\n\n\n"; RECORD();
        break;
      }

      if(reason == SNES_DIVERGED_LINEAR_SOLVE)
      {
        KSPConvergedReason ksp_reason;
        KSPGetConvergedReason ( ksp, &ksp_reason );
        MESSAGE <<"------> linear solver "<<KSPConvergedReasons[ksp_reason]<<", do recovery...\n\n\n"; RECORD();
      }
      else
      {
        MESSAGE <<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n"; RECORD();
      }

      // reduce time step by a factor of two, also set clock to next
      SolverSpecify::PseudoTimeStepPotential /= 2.0;
      SolverSpecify::PseudoTimeStepCarrier /= 2.0;
      SolverSpecify::PseudoTimeStepMetal /= 2.0;

      dstep /= 2.0;
      step -= dstep;

      // load previous result into solution vector
      this->diverged_recovery();
      goto Predict;
    }

    MESSAGE<<"--------------------------------------------------------------------------------\n"
           <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits << "\n\n\n";
    RECORD();


    // clear the counter
    diverged_retry = 0;

    // time step counter ++
    t_cycles++;

    // call post_solve_process
    this->post_solve_process();

    // set next pseudo time step
    if(dstep < 1.0)
      dstep = std::min(1.5*dstep, 1.0);

    step += dstep;

    //make sure we can terminat at ramp_final_time, the relative error should less than 1e-10.
    if ( step > ramp_steps && step < ( ramp_steps + (1 - 1e-10)*dstep) )
    {
      dstep -= step - ramp_steps;
      step = ramp_steps;
    }

    Predict: continue;

  }while ( step < ramp_steps +0.5*dstep );


  if(SolverSpecify::PseudoTimeMethod && SolverSpecify::OpToSteady)
  {
    MESSAGE <<"Drive device to steady state..."<<"\n\n\n";
    RECORD();

    _system.get_electrical_source()->update ( 0 );
    SolverSpecify::PseudoTimeTolRelax = 1e7;
    snes_solve_pseudo_time_step();
  }


  if(SolverSpecify::PseudoTimeMethod && SolverSpecify::PseudoTimeCMOS)
  {
    for(unsigned int n=0; n<_system.n_regions(); n++)
    {
      SimulationRegion * region = _system.region(n);
      if( region->type() == MetalRegion )
      {
        MetalSimulationRegion * resistance_region = dynamic_cast<MetalSimulationRegion *>(region);
        region->set_conductance(resistance_region->material()->basic->Conductance());
      }
    }

    if(SolverSpecify::OpToSteady)
    {
      MESSAGE <<"Recover metal region of CMOS device..."<<"\n\n\n";
      RECORD();

      SolverSpecify::PseudoTimeTolRelax = 1e8;
      snes_solve_pseudo_time_step();
    }
  }


  // restore pseudo parameters
  SolverSpecify::PseudoTimeStepPotential = dtao_init_potential;
  SolverSpecify::PseudoTimeStepCarrier = dtao_init_carrier;
  SolverSpecify::PseudoTimeStepMetal = dtao_init_metal;


  SolverSpecify::tran_histroy = false;

  return 0;
}




/**
 * create ksp solver for trace mode
 */
void DDMSolverBase::solve_iv_trace_begin()
{
  VecDuplicate(x, &pdI_pdx);
  VecDuplicate(x, &pdF_pdV);
  VecDuplicate(x, &pdx_pdV);

  // build special linear solver contex
  {
    PetscErrorCode ierr;

    ierr = KSPCreate(PETSC_COMM_WORLD, &kspc); genius_assert(!ierr);

    ierr = KSPGetPC(kspc, &pcc); genius_assert(!ierr);

    if(Genius::n_processors()>1)
    {
#if defined(PETSC_HAVE_SUPERLU_DIST) || defined(PETSC_HAVE_MUMPS)
      ierr = KSPSetType(kspc, KSPPREONLY); genius_assert(!ierr);
      ierr = PCSetType(pcc, PCLU); genius_assert(!ierr);
#ifdef PETSC_HAVE_MUMPS
      ierr = PCFactorSetMatSolverPackage (pcc, "mumps"); genius_assert(!ierr);
#else
      ierr = PCFactorSetMatSolverPackage (pcc, "superlu_dist"); genius_assert(!ierr);
#endif
#else
      // no parallel LU solver? we have to use krylov method for parallel!
      ierr = KSPSetType(kspc, KSPBCGS); genius_assert(!ierr);
      ierr = PCSetType(pcc, PCASM); genius_assert(!ierr);
#endif
    }
    else
    {
      ierr = KSPSetType(kspc, KSPPREONLY); genius_assert(!ierr);
      ierr = PCSetType(pcc, PCLU); genius_assert(!ierr);
#ifdef PETSC_HAVE_MUMPS
      ierr = PCFactorSetMatSolverPackage (pcc, "mumps"); genius_assert(!ierr);
#endif
    }


    ierr = KSPSetOperators(kspc, J, J, SAME_NONZERO_PATTERN); genius_assert(!ierr);
  }
}


/**
 * destroy ksp solver for trace mode
 */
void DDMSolverBase::solve_iv_trace_end()
{
  KSPDestroy(PetscDestroyObject(kspc));
  VecDestroy(PetscDestroyObject(pdI_pdx));
  VecDestroy(PetscDestroyObject(pdF_pdV));
  VecDestroy(PetscDestroyObject(pdx_pdV));
}



/* ----------------------------------------------------------------------------
 * DDMSolverBase::solve_iv_trace:  This function use continuation method to trace
 * IV curve automatically
 */
int DDMSolverBase::solve_iv_trace()
{
  int         error=0;
  int         first_step=1;
  int         slope_flag=0;

  const double PI = 3.14159265358979323846264338327950;
  const double degree = PI/180.0;

  PetscScalar I=0;
  PetscScalar dI_dV;

  std::string electrode_trace = SolverSpecify::Electrode_VScan[0];
  BoundaryCondition * bc_trace = _system.get_bcs()->get_bc(electrode_trace);
  PetscScalar R_bak = bc_trace->ext_circuit()->serial_resistance();

  // the current vscan voltage
  PetscScalar V = SolverSpecify::VStart;
  // set current vscan voltage to corresponding electrode
  _system.get_electrical_source()->assign_voltage_to ( electrode_trace, V );

  // the current vscan step
  PetscScalar VStep = SolverSpecify::VStep;
  PetscScalar VStepFactor = 1.0;

  // the current potential of vscan electrode
  PetscScalar Potential = bc_trace->ext_circuit()->potential();

  // r=dI/dV
  PetscScalar r;
  PetscScalar r_new;

  // load resistance
  PetscScalar Rload=0;
  PetscScalar Rload_new=0;

  // scaling factor
  PetscScalar Rref=PhysicalUnit::V/(1e-5*PhysicalUnit::A);

  // slope of load line
  PetscScalar slope=PI/2;
  PetscScalar slope_new;
  PetscScalar slope_chord;

  // set electrode with transient time 0 value of stimulate source(s)
  _system.get_electrical_source()->update ( 0 );
  _system.get_field_source()->update ( 0 );

  // not time dependent
  SolverSpecify::TimeDependent = false;
  SolverSpecify::dt = 1e100;
  SolverSpecify::clock = 0.0;

  solve_iv_trace_begin();

  // output TRACE information
  MESSAGE<<"IV automatically trace by continuation method\n"; RECORD();

  // initial condition check
  bc_trace->ext_circuit()->set_serial_resistance(Rload);

  // call pre_solve_process
  this->pre_solve_process();
  // here call Petsc to solve the nonlinear equations
  sens_solve();

  // get the converged reason
  SNESConvergedReason reason;
  SNESGetConvergedReason ( snes, &reason );

  if(reason<0)
  {
    MESSAGE<<"I can't get convergence even at initial point, need a better initial condition.\n\n"; RECORD();
    error = 1;
    goto trace_end;
  }

  // linear solver iteration
  PetscInt lits;
  SNESGetLinearSolveIterations(snes, &lits);

  MESSAGE
      <<"--------------------------------------------------------------------------------\n"
      <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits << "\n\n\n";
  RECORD();

  // call post_solve_process
  this->post_solve_process();

  I = bc_trace->ext_circuit()->current();
  Potential = bc_trace->ext_circuit()->potential();

  //loop here
  while(Potential*SolverSpecify::VStep < SolverSpecify::VStop*SolverSpecify::VStep && fabs(I)<SolverSpecify::IStop)
  {
    // for last Rload, increase VStep
    int recovery=0;
    do
    {
      V += VStep;
      bc_trace->ext_circuit()->Vapp() =  V;

      MESSAGE << "Trace "<< electrode_trace <<" for VTrace=" << V << "(V), V=" << Potential << "(V)\n"; RECORD();

      this->pre_solve_process(false);
      sens_solve();

      SNESGetConvergedReason(snes,&reason);
      if(reason<0)
      {
        MESSAGE<<"--------------------------------------------------------------------------------\n"
               <<"I can't get convergence at this step, do recovery...\n\n\n";
        RECORD();

        this->diverged_recovery();
        V -= VStep;
        VStep/=2;
        recovery++;
      }
      if(recovery>8)
      {
        MESSAGE<<"------>  Too many failed steps, give up tring.\n\n\n";RECORD();
        error = 1;
        goto trace_end;
      }
      SNESGetLinearSolveIterations(snes, &lits);
    }
    while(reason<0);

    MESSAGE
        <<"--------------------------------------------------------------------------------\n"
        <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits << "\n\n\n";
    RECORD();


    // for the new bias point, recompute slope
    // calculate the dynamic resistance of IV curve by different approximation
    this->set_trace_electrode(bc_trace);

    KSPSetOperators(kspc,J,J,SAME_NONZERO_PATTERN);
    KSPSolve(kspc, pdF_pdV, pdx_pdV); // KSPSolve(ksp, b, x)
    VecDot(pdI_pdx, pdx_pdV, &dI_dV);

    // compute tangent line of IV curve as well as load resistance
    {

      I = bc_trace->ext_circuit()->current();
      Parallel::sum(I);
      Potential = bc_trace->ext_circuit()->potential();

      //if( std::abs(I) > 1e-5*PhysicalUnit::A )
      //  Rref = PhysicalUnit::V/std::abs(I);

      if(first_step)
      {
        r = r_new = 1/dI_dV;
        Rload_new=Rref*Rref/r_new;
        slope=slope_new=atan(Rload_new);
      }
      else
      {
        r_new = 1/dI_dV;
        Rload_new=Rref*Rref/r_new;
        slope_new = atan(Rload_new);
        Rload = Rref*Rref/r;
        slope = atan(Rload);
      }

      /*
      MESSAGE<<"I="<<std::abs(I)/PhysicalUnit::A << " " <<"V="<< std::abs(Potential) <<std::endl;
      RECORD();
      MESSAGE<<"Rref="<<Rref << " " <<"V/I="<<std::abs(Potential/I)<<" "<<"r_new="<< r_new << " " <<"Rload_new="<< Rload_new<<std::endl;
      RECORD();
      MESSAGE<<"slope="<<slope/degree<< " " <<"slope_new="<< slope_new/degree << "  "<<"dslope="<<(slope_new-slope)/degree  <<std::endl;
      RECORD();
      */
    }

    // check the slope change
    if(!first_step)
    {
      PetscScalar angle=fabs(slope-slope_new);
      if(angle>90*degree) angle=PI-angle;

      if(angle<1*degree)       VStepFactor=2.0;     // slope change less than 1 degree
      else if(angle<5*degree)  VStepFactor=1.5;   // slope change less than 5 degree
      else if(angle<10*degree) VStepFactor=1.0;   // slope change less than 10 degree, but greater than 5  degree
      else if(angle<15*degree) VStepFactor=0.5;   // slope change less than 15 degree, but greater than 10 degree
      else                                        // slope greater than 15 degree, reject this solution
      {
        MESSAGE<<"Slope of IV curve changes too quickly, do recovery...\n\n";
        RECORD();
        this->diverged_recovery();
        V -= VStep;
        VStep/=2;
        Potential=0.0;
        continue;
      }
    }

    // ok, update solutions
    this->post_solve_process();

    //resolve turning point
    if(Rload*Rload_new<0) //the sign is changed
    {
      PetscScalar chord = (bc_trace->ext_circuit()->potential()-bc_trace->ext_circuit()->potential_old())/
                          (bc_trace->ext_circuit()->current()-bc_trace->ext_circuit()->current_old());
      slope_chord = atan(Rref*Rref/chord);

      //for vertical turning point
      if(fabs(slope)>60*degree && fabs(slope_new)>60*degree &&
          ((slope_chord*slope>0 && fabs(slope_chord)>fabs(slope)) || (slope_chord*slope_new>0 && fabs(slope_chord)>fabs(slope_new))))
      {
        //MESSAGE<<"vertical turning point"<<std::endl;
        VStep = - VStep;
      }
      //for horizontal turning point
      else if(fabs(slope)<30*degree && fabs(slope_new)<30*degree &&
          ((slope_chord*slope>0 && fabs(slope_chord)<fabs(slope)) || (slope_chord*slope_new>0 && fabs(slope_chord)<fabs(slope_new))))
      {
        //MESSAGE<<"horizontal turning point"<<std::endl;
      }
      else
      {
        MESSAGE<<"Trace Internal Error"<<std::endl; RECORD();
        genius_error();
      }
    }

    VStep *= VStepFactor;

    // since Rload will be changed, we should change bias voltage to meet the truncation point.
    V += I*(Rload_new-Rload);
    bc_trace->ext_circuit()->set_serial_resistance(Rload_new);
    bc_trace->ext_circuit()->Vapp() = V;

    Rload = Rload_new;
    slope = slope_new;
    r     = r_new;


    // limit the VStep
    {
      /*
      double theta = atan(dI_dV);
      if( std::abs(VStep)*cos(theta)*cos(theta) > SolverSpecify::VStepMax )
      {
        MESSAGE<<"VStep="<<VStep<<std::endl; RECORD();
        VStep = std::sign(VStep)*SolverSpecify::VStepMax/cos(theta)/cos(theta);
      }
      MESSAGE<<"VStep="<<VStep<<std::endl; RECORD();
     */
    }

    first_step=0;
  }



trace_end:
  bc_trace->ext_circuit()->set_serial_resistance(R_bak);
  solve_iv_trace_end();


  SolverSpecify::tran_histroy = false;

  return error;

}


/*----------------------------------------------------------------------------
 * transient simulation!
 */
int DDMSolverBase::solve_transient()
{

  // init aux vectors used in transient simulation
  VecDuplicate ( x, &x_n );
  VecDuplicate ( x, &x_n1 );
  VecDuplicate ( x, &x_n2 );
  VecDuplicate ( x, &xp );
  VecDuplicate ( x, &LTE );

  // time dependent
  SolverSpecify::TimeDependent = true;

  // if BDF2 scheme is used, we should set SolverSpecify::BDF2_LowerOrder flag to true
  if ( SolverSpecify::TS_type==SolverSpecify::BDF2 )
    SolverSpecify::BDF2_LowerOrder = true;

  // we have a previous dc solution
  if(!SolverSpecify::tran_histroy)
  {
    _system.get_electrical_source()->update ( SolverSpecify::TStart );
    for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); b++)
    {
      BoundaryCondition * bc = _system.get_bcs()->get_bc(b);
      if(bc && bc->is_electrode())
        bc->ext_circuit()->tran_op_init();
    }
  }

  // transient simulation clock
  SolverSpecify::clock = SolverSpecify::TStart + SolverSpecify::TStep;

  // for the first step, dt equals TStep
  SolverSpecify::dt = SolverSpecify::TStep;

  MESSAGE<<"Transient compute from "<<SolverSpecify::TStart/s*1e12
      <<" ps step "<<SolverSpecify::TStep/s*1e12
      <<" ps to "  <<SolverSpecify::TStop/s*1e12<<" ps"
  <<'\n';
  RECORD();

  // diverged counter
  int diverged_retry=0;

  // auto time step counter
  int autostep_retry=0;

  // time step counter
  SolverSpecify::T_Cycles=0;

  double dt_dynamic_factor = 1.0;

  // the main loop of transient solver.
  do
  {
    MESSAGE
    <<"t = "<<SolverSpecify::clock/s*1e12<<" ps, "<< "dt = " << SolverSpecify::dt/s*1e12 <<" ps"<< '\n'
    <<"--------------------------------------------------------------------------------\n";
    RECORD();

    //update sources to current clock
    _system.get_electrical_source()->update ( SolverSpecify::clock );
    _system.get_field_source()->update ( SolverSpecify::clock, SolverSpecify::SourceCoupled );
    //we do solve here!

    // call pre_solve_process
    if ( SolverSpecify::T_Cycles == 0 )
      this->pre_solve_process();
    else
      this->pre_solve_process ( false );

    sens_solve();
    // get the converged reason
    SNESConvergedReason reason;
    SNESGetConvergedReason ( snes,&reason );


    // linear solver iteration
    PetscInt lits;
    SNESGetLinearSolveIterations(snes, &lits);

    //nonlinear solution diverged? try to do recovery
    if ( reason<0 )
    {
      // increase diverged_retry
      diverged_retry++;

      if ( diverged_retry >= 8 ) //failed 8 times, stop tring
      {
        MESSAGE
        <<"------> Too many failed steps, give up tring.\n\n\n";
        RECORD();
        break;
      }


      if(reason == SNES_DIVERGED_LINEAR_SOLVE)
      {
        KSPConvergedReason ksp_reason;
        KSPGetConvergedReason ( ksp, &ksp_reason );
        MESSAGE <<"------> linear solver "<<KSPConvergedReasons[ksp_reason]<<", do recovery...\n\n\n"; RECORD();
      }
      else
      {
        MESSAGE <<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n"; RECORD();
      }

      // reduce time step by a factor of two, also set clock to next
      SolverSpecify::dt /= 2.0;
      SolverSpecify::clock -= SolverSpecify::dt;

      if ( SolverSpecify::clock < SolverSpecify::TStart )
        SolverSpecify::clock = SolverSpecify::TStart;

      // load previous result into solution vector
      this->diverged_recovery();
      goto Predict;
    }

    //ok, nonlinear solution converged.

    dt_dynamic_factor = 1.0;

    //do LTE estimation and auto time step control
    if ( SolverSpecify::AutoStep &&
         ( ( SolverSpecify::TS_type==SolverSpecify::BDF1 && SolverSpecify::T_Cycles>=2 ) ||
           ( SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::T_Cycles>=3 ) ) )
    {
      PetscReal r = this->LTE_norm() + 1e-10;

      if ( SolverSpecify::TS_type==SolverSpecify::BDF1 )
        r = std::pow ( r, PetscReal ( -1.0/2 ) );
      else if ( SolverSpecify::TS_type==SolverSpecify::BDF2 )
      {
        if(SolverSpecify::BDF2_LowerOrder)
          r = std::pow ( r, PetscReal ( -1.0/2 ) );
        else
          r = std::pow ( r, PetscReal ( -1.0/3 ) );
      }

      // when r<0.9, reject this solution
      if ( SolverSpecify::RejectStep && r<0.9 && SolverSpecify::dt > SolverSpecify::TStepMin )
      {
        // clear the counter
        diverged_retry = 0;
        autostep_retry++;

        MESSAGE<<"------> LTE too large, time step rejected...\n\n\n";
        RECORD();

        // reduce time step by a factor of 0.9*r
        SolverSpecify::clock -= SolverSpecify::dt;
        PetscScalar hn  = SolverSpecify::dt;           // here dt is the current time step
        SolverSpecify::dt *= 0.9*r;
        SolverSpecify::clock += SolverSpecify::dt;
        PetscScalar hn_new = SolverSpecify::dt;           // next time step

        // use linear interpolation to predict solution x at next time step
        VecScale ( x, hn_new/hn );
        VecAXPY ( x, 1-hn_new/hn,  x_n );
        this->projection_positive_density_check ( x, x_n );

        continue;
      }
      else      // accept this solution
      {
        // set next time step
        if( autostep_retry || diverged_retry)
        {
          if ( r > 1.0 )
            dt_dynamic_factor = 1.0;
          else
            dt_dynamic_factor = std::min(r, 0.9);
          autostep_retry = 0;
          diverged_retry = 0;
        }
        else
        {
          if ( r > 1.0 )
            dt_dynamic_factor = 1.0 + log10(r);
          else
            dt_dynamic_factor = std::min(r, 0.9);
        }
      }
    }
    else // auto time step control not used
    {
      // set next time step
      if ( fabs ( SolverSpecify::dt ) < fabs ( SolverSpecify::TStep ) )
        dt_dynamic_factor = 1.1;
    }


    MESSAGE
        <<"--------------------------------------------------------------------------------\n"
        <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits << "\n\n\n";
    RECORD();


    // call post_solve_process
    this->post_solve_process();

    // clear the counter
    diverged_retry = 0;

    // time step counter ++
    SolverSpecify::T_Cycles++;

    // save time step information
    SolverSpecify::dt_last_last = SolverSpecify::dt_last;
    SolverSpecify::dt_last = SolverSpecify::dt;

    // prepare for next time step
    SolverSpecify::dt *= dt_dynamic_factor;

    // limit the max time step by TStepMin/TStepMax
    if ( SolverSpecify::dt < SolverSpecify::TStepMin )
      SolverSpecify::dt = SolverSpecify::TStepMin;
    if ( SolverSpecify::dt > SolverSpecify::TStepMax )
      SolverSpecify::dt = SolverSpecify::TStepMax;

    // limit time step by changes of external source, i.e. max allowed changes of vsource
    SolverSpecify::dt = _system.get_electrical_source()->limit_dt(SolverSpecify::clock, SolverSpecify::dt, SolverSpecify::TStepMin, SolverSpecify::VStepMax, SolverSpecify::IStepMax);

    // limit time step by field source
    SolverSpecify::dt = _system.get_field_source()->limit_dt(SolverSpecify::clock, SolverSpecify::dt);

    // set clock to next time step
    SolverSpecify::clock += SolverSpecify::dt;

    // make sure we can terminat near TStop, the relative error should less than 1e-10.
    if ( SolverSpecify::clock > SolverSpecify::TStop && SolverSpecify::clock < ( SolverSpecify::TStop + SolverSpecify::dt - 1e-10*SolverSpecify::dt ) )
    {
      SolverSpecify::dt -= SolverSpecify::clock - SolverSpecify::TStop;
      SolverSpecify::clock = SolverSpecify::TStop;
    }


    //check if BDF2 can be used?
    if ( SolverSpecify::TS_type==SolverSpecify::BDF2 )
      SolverSpecify::BDF2_LowerOrder = this->BDF2_positive_defined();

    // use by auto step control and predict
    if( SolverSpecify::AutoStep  || SolverSpecify::Predict )
    {
      VecCopy ( x_n1, x_n2 );
      VecCopy ( x_n, x_n1 );
      VecCopy ( x, x_n );
    }

  Predict:

    // predict next solution
    if ( SolverSpecify::Predict )
    {
      PetscScalar hn  = SolverSpecify::dt;           // here dt is the next time step
      PetscScalar hn1 = SolverSpecify::dt_last;      // time step n-1
      PetscScalar hn2 = SolverSpecify::dt_last_last; // time step n-2

      if ( SolverSpecify::TS_type == SolverSpecify::BDF1 && SolverSpecify::T_Cycles>=2)
      {
        VecZeroEntries ( x );
        // use linear interpolation to predict solution x
        VecAXPY ( x, 1+hn/hn1, x_n );
        VecAXPY ( x, -hn/hn1,  x_n1 );
        this->projection_positive_density_check ( x, x_n );
      }
      if ( SolverSpecify::TS_type == SolverSpecify::BDF2 && SolverSpecify::T_Cycles>=3)
      {
        VecZeroEntries ( x );

        if(SolverSpecify::BDF2_LowerOrder)
        {
          // use linear interpolation to predict solution x
          VecAXPY ( x, 1+hn/hn1, x_n );
          VecAXPY ( x, -hn/hn1,  x_n1 );
          this->projection_positive_density_check ( x, x_n );
        }
        else
        {
          // use second order polynomial to predict solution x
          PetscScalar cn  = 1+hn* ( hn+2*hn1+hn2 ) / ( hn1* ( hn1+hn2 ) );
          PetscScalar cn1 = -hn* ( hn+hn1+hn2 ) / ( hn1*hn2 );
          PetscScalar cn2 = hn* ( hn+hn1 ) / ( hn2* ( hn1+hn2 ) );

          VecAXPY ( x, cn,  x_n );
          VecAXPY ( x, cn1, x_n1 );
          VecAXPY ( x, cn2, x_n2 );
          this->projection_positive_density_check ( x, x_n );
        }
      }
    }

  }
  while ( SolverSpecify::clock < SolverSpecify::TStop+0.5*SolverSpecify::dt );

  // free aux vectors
  VecDestroy ( PetscDestroyObject(x_n) );
  VecDestroy ( PetscDestroyObject(x_n1) );
  VecDestroy ( PetscDestroyObject(x_n2) );
  VecDestroy ( PetscDestroyObject(xp) );
  VecDestroy ( PetscDestroyObject(LTE) );


  SolverSpecify::tran_histroy = true;

  return 0;
}



int DDMSolverBase::snes_solve_pseudo_time_step()
{
  // diverged counter
  int diverged_retry=0;


  for(int k=1; k<=SolverSpecify::PseudoTimeSteps; k++)
  {
    MESSAGE <<"PseudoTime Step "<< k <<'\n'
            <<"--------------------------------------------------------------------------------\n";
    RECORD();

    // call pre_solve_process
    if ( SolverSpecify::T_Cycles == 0 )
      this->pre_solve_process();
    else
      this->pre_solve_process ( false );

    sens_solve();
    // get the converged reason
    SNESConvergedReason reason;
    SNESGetConvergedReason ( snes,&reason );

    // linear solver iteration
    PetscInt lits;
    SNESGetLinearSolveIterations(snes, &lits);

    if ( reason<0 )
    {
      // increase diverged_retry
      diverged_retry++;

      if ( diverged_retry >= 8 ) //failed 8 times, stop tring
      {
        MESSAGE<<"------> Too many failed steps, give up tring.\n\n\n"; RECORD();
        break;
      }


      if(reason == SNES_DIVERGED_LINEAR_SOLVE)
      {
        KSPConvergedReason ksp_reason;
        KSPGetConvergedReason ( ksp, &ksp_reason );
        MESSAGE <<"------> linear solver "<<KSPConvergedReasons[ksp_reason]<<", do recovery...\n\n\n"; RECORD();
      }
      else
      {
        MESSAGE <<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n"; RECORD();
      }
      SolverSpecify::PseudoTimeStepPotential /= 2.0;
      SolverSpecify::PseudoTimeStepCarrier /= 2.0;
      SolverSpecify::PseudoTimeStepMetal /= 2.0;
      this->diverged_recovery();
      continue;
    }

    MESSAGE <<"--------------------------------------------------------------------------------\n"
            <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits << "\n\n\n";
    RECORD();


    // check for convergence of pseudo time method
    if( this->pseudo_time_step_convergence_test() )
    {
      this->post_solve_process();
      break;
    }

    // call post_solve_process
    this->post_solve_process();

    if(diverged_retry == 0)
    {
      // set next pseudo time step
      if(SolverSpecify::PseudoTimeStepPotential < SolverSpecify::PseudoTimeStepMax)
        SolverSpecify::PseudoTimeStepPotential *= 2.0;
      if(SolverSpecify::PseudoTimeStepCarrier < SolverSpecify::PseudoTimeStepMax)
        SolverSpecify::PseudoTimeStepCarrier *= 2.0;
      if(SolverSpecify::PseudoTimeStepMetal < SolverSpecify::PseudoTimeStepMax)
        SolverSpecify::PseudoTimeStepMetal *= 2.0;
    }

    // do predict
  }


  return 0;
}




/*------------------------------------------------------------------
 * snes convergence criteria
 */
#if PETSC_VERSION_LE(3, 2, 0)
  #include "private/snesimpl.h"
  #define SNES_CONVERGED_SNORM_RELATIVE  SNES_CONVERGED_PNORM_RELATIVE
#else
  #include "petsc-private/snesimpl.h"
#endif
void DDMSolverBase::petsc_snes_convergence_test ( PetscInt its, PetscReal , PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason )
{
  // update error norm
  this->error_norm();

  *reason = SNES_CONVERGED_ITERATING;

  // the first iteration
  if ( !its )
  {
    snes->ttol = fnorm*snes->rtol;

    MESSAGE<<" "<<" n ";
    MESSAGE<<"| Eq(V) | "<<"| Eq(n) | "<<"| Eq(p) | ";
    MESSAGE<<"| Eq(T) | ";
    if(this->solver_type() == SolverSpecify::DENSITY_GRADIENT)
    {
      MESSAGE<<"|Eq(Qn)|  ";
      MESSAGE<<"|Eq(Qp)|  ";
    }
    else
    {
      MESSAGE<<"|Eq(Tn)|  ";
      MESSAGE<<"|Eq(Tp)|  ";
    }
    MESSAGE<<"|Eq(BC)|  ";
    MESSAGE<<"Lg(dx)"<<'\n';
    MESSAGE<<"--------------------------------------------------------------------------------\n";
    RECORD();

    functions_norm[0] = poisson_norm;
    functions_norm[1] = elec_continuity_norm;
    functions_norm[2] = hole_continuity_norm;
    functions_norm[3] = heat_equation_norm;
    functions_norm[4] = elec_energy_equation_norm;
    functions_norm[5] = hole_energy_equation_norm;
    functions_norm[6] = elec_quantum_equation_norm;
    functions_norm[7] = hole_quantum_equation_norm;
    functions_norm[8] = electrode_norm;
  }

  unsigned int dim = this->system().dim();
  double z_width = (dim == 2 ? 1.0*um : 1.0);
  if(dim == 2)
  {
    poisson_norm         *= z_width;
    elec_continuity_norm *= z_width;
    hole_continuity_norm *= z_width;
    heat_equation_norm   *= z_width;
    elec_energy_equation_norm *= z_width;
    hole_energy_equation_norm *= z_width;
    elec_quantum_equation_norm *= z_width;
    hole_quantum_equation_norm *= z_width;
  }

  double  toler_relax = SolverSpecify::toler_relax;
  if(its)
  {
    //toler_relax = std::max(SolverSpecify::toler_relax, 1.0/(pnorm+1e-6));
  }
  bool  poisson_conv               = poisson_norm               < toler_relax*SolverSpecify::poisson_abs_toler;
  bool  elec_continuity_conv       = elec_continuity_norm       < toler_relax*SolverSpecify::elec_continuity_abs_toler;
  bool  hole_continuity_conv       = hole_continuity_norm       < toler_relax*SolverSpecify::hole_continuity_abs_toler;
  bool  electrode_conv             = electrode_norm             < toler_relax*SolverSpecify::electrode_abs_toler;
  bool  heat_equation_conv         = heat_equation_norm         < toler_relax*SolverSpecify::heat_equation_abs_toler;
  bool  elec_energy_equation_conv  = elec_energy_equation_norm  < toler_relax*SolverSpecify::elec_energy_abs_toler;
  bool  hole_energy_equation_conv  = hole_energy_equation_norm  < toler_relax*SolverSpecify::hole_energy_abs_toler;
  bool  elec_quantum_equation_conv = elec_quantum_equation_norm < toler_relax*SolverSpecify::elec_quantum_abs_toler;
  bool  hole_quantum_equation_conv = hole_quantum_equation_norm < toler_relax*SolverSpecify::hole_quantum_abs_toler;

  bool conv = poisson_conv         &&
              elec_continuity_conv &&
              hole_continuity_conv &&
              electrode_conv       &&
              heat_equation_conv   &&
              elec_energy_equation_conv &&
              hole_energy_equation_conv &&
              elec_quantum_equation_conv &&
              hole_quantum_equation_conv;


  bool abs_conv =  poisson_norm               < SolverSpecify::poisson_abs_toler         &&
                   elec_continuity_norm       < SolverSpecify::elec_continuity_abs_toler &&
                   hole_continuity_norm       < SolverSpecify::hole_continuity_abs_toler &&
                   electrode_norm             < SolverSpecify::electrode_abs_toler       &&
                   heat_equation_norm         < SolverSpecify::heat_equation_abs_toler   &&
                   elec_energy_equation_norm  < SolverSpecify::elec_energy_abs_toler     &&
                   hole_energy_equation_norm  < SolverSpecify::hole_energy_abs_toler     &&
                   elec_quantum_equation_norm < SolverSpecify::elec_quantum_abs_toler    &&
                   hole_quantum_equation_norm < SolverSpecify::hole_quantum_abs_toler;

#ifdef WINDOWS
  MESSAGE.precision ( 1 );
#else
  MESSAGE.precision ( 2 );
#endif

  MESSAGE<< std::setw(3) << its << " " ;
  MESSAGE<< std::scientific;
  MESSAGE<< poisson_norm/C              << (poisson_conv ? "* " : "  ");
  MESSAGE<< elec_continuity_norm/A      << (elec_continuity_conv ? "* " : "  ");
  MESSAGE<< hole_continuity_norm/A      << (hole_continuity_conv ? "* " : "  ");
  MESSAGE<< heat_equation_norm/W        << (heat_equation_conv ? "* " : "  ");
  if(this->solver_type() == SolverSpecify::DENSITY_GRADIENT)
  {
    MESSAGE<< elec_quantum_equation_norm/W << (elec_quantum_equation_conv ? "* " : "  ");
    MESSAGE<< hole_quantum_equation_norm/W << (hole_quantum_equation_conv ? "* " : "  ");
  }
  else
  {
    MESSAGE<< elec_energy_equation_norm/W << (elec_energy_equation_conv ? "* " : "  ");
    MESSAGE<< hole_energy_equation_norm/W << (hole_energy_equation_conv ? "* " : "  ");
  }
  MESSAGE<< electrode_norm/A            << (electrode_conv ? "* " : "  ");
  MESSAGE<< std::fixed << std::setw(4) << (pnorm==0.0 ? -std::numeric_limits<PetscScalar>::infinity():log10(pnorm))
         << (pnorm < SolverSpecify::relative_toler ? "*" : " ") << "\n" ;
  RECORD();
  MESSAGE.precision ( 6 );
  MESSAGE<< std::scientific;


  // check for NaN (Not a Number)
  if ( fnorm != fnorm )
  {
    *reason = SNES_DIVERGED_FNORM_NAN;
  }
  else if ( snes->nfuncs >= snes->max_funcs )
  {
    *reason = SNES_DIVERGED_FUNCTION_COUNT;
  }

  if ( *reason == SNES_CONVERGED_ITERATING )
  {
    // special treatment in PseudoTimeMethod
    if(SolverSpecify::PseudoTimeMethod == true)
    {
      // in pseudo time method mode, should always do newton step at least once
      if(!its) goto end;

      if(pnorm <= SolverSpecify::PseudoTimeMethodRXTol &&
         (poisson_norm <= SolverSpecify::PseudoTimeMethodRFTol*functions_norm[0] )&&
         (elec_continuity_norm <= SolverSpecify::PseudoTimeMethodRFTol*functions_norm[1] )&&
         (hole_continuity_norm <= SolverSpecify::PseudoTimeMethodRFTol*functions_norm[2] )&&
         (heat_equation_norm <= SolverSpecify::PseudoTimeMethodRFTol*functions_norm[3] ) &&
         (elec_energy_equation_norm <= SolverSpecify::PseudoTimeMethodRFTol*functions_norm[4] ) &&
         (elec_energy_equation_norm <= SolverSpecify::PseudoTimeMethodRFTol*functions_norm[5] ) &&
         (electrode_norm <= SolverSpecify::PseudoTimeMethodRFTol*functions_norm[6] || electrode_conv ) )
      {
        *reason = SNES_CONVERGED_SNORM_RELATIVE;
        goto end;
      }
    }

    if( fnorm<snes->abstol )
    {
      *reason = SNES_CONVERGED_FNORM_ABS;
    }
    else
    {
      // check for absolute convergence
      if ( its && abs_conv )
      {
        *reason = SNES_CONVERGED_FNORM_ABS;
      }
      else if ( std::abs(fnorm-function_norm)/fnorm <= snes->rtol  && conv )
      {
        *reason = SNES_CONVERGED_FNORM_RELATIVE;
      }
      // check for relative convergence, should have at least one iteration here
      else if ( its && pnorm < SolverSpecify::relative_toler && conv)
      {
        *reason = SNES_CONVERGED_SNORM_RELATIVE;
      }
    }
  }

  end:

  // record function norm of this iteration
  function_norm = fnorm;
  // record iteration
  nonlinear_iteration = its;

  return;
}



/*------------------------------------------------------------------
 * ksp convergence criteria
 */
void DDMSolverBase::petsc_ksp_convergence_test ( PetscInt its, PetscReal rnorm, KSPConvergedReason* reason )
{

  PetscInt kspit = std::max(200, std::min(1000, static_cast<int>(n_global_dofs/10)));
  PetscScalar rtol = SolverSpecify::ksp_rtol;
  PetscScalar abstol = std::max ( SolverSpecify::ksp_atol_fnorm*function_norm, SolverSpecify::ksp_atol);

  if(its > static_cast<PetscInt>(0.3*kspit))
    abstol *= 1e1;
  if(its > static_cast<PetscInt>(0.5*kspit))
    abstol *= 1e2;
  if(its > static_cast<PetscInt>(0.7*kspit))
    abstol *= 1e3;

  KSPSetTolerances(ksp, rtol, abstol, 1e10, kspit );
  KSPDefaultConverged ( ksp, its, rnorm, reason, this );

  // stupid code, but KSPSetTolerances does NOT affect KSPDefaultConverged here!
  // seems KSPDefaultConverged cache the value...
  if(*reason == 0)
  {
    if(rnorm < abstol) *reason = KSP_CONVERGED_ATOL;
  }

  /*
   * it is impossible to change inner linear solver during snes iteration
  if( SolverSpecify::linear_solver_category(_linear_solver_type) == SolverSpecify::ITERATIVE )
  {
    ksp->rtol = PetscMax ( SolverSpecify::ksp_rtol, 1e-12*std::sqrt ( double ( n_global_dofs ) ) );
    ksp->abstol = PetscMax ( SolverSpecify::ksp_atol_fnorm*function_norm, SolverSpecify::ksp_atol*std::sqrt ( double ( n_global_dofs ) ) );
    KSPDefaultConverged ( ksp, its, rnorm, reason, this );
    if( *reason < 0 || (*reason == 0 && its >= ksp->max_it) )
    {
      _linear_solver_stack.push_back(_linear_solver_type);
      _linear_solver_type = SolverSpecify::LU;
      set_petsc_linear_solver_type(_linear_solver_type);
      ksp->reason = KSP_CONVERGED_ITERATING;
      ksp->its = 0;
      return;
    }
  }

  else if( SolverSpecify::linear_solver_category(_linear_solver_type) == SolverSpecify::DIRECT )
  {
    if(!_linear_solver_stack.empty())
    {
      _linear_solver_type = _linear_solver_stack.back();
      _linear_solver_stack.pop_back();
      set_petsc_linear_solver_type(_linear_solver_type);
      set_petsc_preconditioner_type(_preconditioner_type);
    }
  }
  */
}

