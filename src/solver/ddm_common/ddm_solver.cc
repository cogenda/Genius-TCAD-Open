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


#include <stack>

#include "solver_specify.h"
#include "physical_unit.h"
#include "electrical_source.h"
#include "field_source.h"
#include "ddm_solver.h"
#include "MXMLUtil.h"

int DDMSolverBase::create_solver()
{
  // must setup nonlinear contex here!
  setup_nonlinear_data();

  //NOTE Tolerances here only be set as a reference

  //abstol = 1e-12*n_global_dofs    - absolute convergence tolerance
  //rtol   = 1e-14                  - relative convergence tolerance
  //stol   = 1e-9                   - convergence tolerance in terms of the norm of the change in the solution between steps
  SNESSetTolerances(snes, 1e-12*n_global_dofs, 1e-14, 1e-9, SolverSpecify::MaxIteration, 1000);

  // rtol   = 1e-12*n_global_dofs  - the relative convergence tolerance (relative decrease in the residual norm)
  // abstol = 1e-20*n_global_dofs  - the absolute convergence tolerance (absolute size of the residual norm)
  KSPSetTolerances(ksp, 1e-12*n_global_dofs, 1e-20*n_global_dofs, PETSC_DEFAULT, std::max(50, static_cast<int>(n_global_dofs/10)));

  // user can do further adjusment from command line
  SNESSetFromOptions (snes);

  return FVM_NonlinearSolver::create_solver();
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
        mxmlAdd(eTime, MXML_ADD_AFTER, NULL, MXMLQVariant::makeQVFloat(SolverSpecify::clock));
        break;
      }
    case SolverSpecify::DCSWEEP:
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

  // print convergence/divergence reason
  {
    MESSAGE
    <<"-----------------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[reason]<<"\n\n\n";
    RECORD();
  }

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
  _system.get_sources()->update ( 0 );
  _system.get_field_source()->update ( 0 );

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

  // print convergence/divergence reason
  {
    MESSAGE
    <<"-----------------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[reason]<<"\n\n\n";
    RECORD();
  }

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
  _system.get_sources()->update ( 0 );
  _system.get_field_source()->update ( 0 );

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



  // voltage scan
  if ( SolverSpecify::Electrode_VScan.size() )
  {
    // the current vscan voltage
    PetscScalar Vscan = SolverSpecify::VStart;

    // the current vscan step
    PetscScalar VStep = SolverSpecify::VStep;

    // saved solutions and vscan values for solution projection.
    Vec xs1, xs2, xs3;
    PetscScalar Vs1, Vs2, Vs3;
    std::stack<PetscScalar> V_retry;
    VecDuplicate ( x,&xs1 );
    VecDuplicate ( x,&xs2 );
    VecDuplicate ( x,&xs3 );

    // main loop
    for ( SolverSpecify::DC_Cycles=0;  (Vscan*SolverSpecify::VStep) < SolverSpecify::VStop*SolverSpecify::VStep* ( 1.0+1e-7 ); )
    {
      // show current vscan value
      MESSAGE << "DC Scan: V("  << SolverSpecify::Electrode_VScan[0];
      for ( unsigned int i=1; i<SolverSpecify::Electrode_VScan.size(); i++ )
        MESSAGE << ", "  << SolverSpecify::Electrode_VScan[i];
      MESSAGE << ") = "  << Vscan/PhysicalUnit::V  <<" V" << '\n'
      <<"-----------------------------------------------------------------------------\n";
      RECORD();

      // set current vscan voltage to corresponding electrode
      _system.get_sources()->assign_voltage_to ( SolverSpecify::Electrode_VScan, Vscan );

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

      if ( reason>0 ) //ok, converged.
      {

        // call post_solve_process
        this->post_solve_process();

        SolverSpecify::DC_Cycles++;

        // save solution for linear/quadratic projection
        VecCopy ( xs2,xs3 ); Vs3=Vs2;
        VecCopy ( xs1,xs2 ); Vs2=Vs1;
        VecCopy ( x,xs1 );   Vs1=Vscan;

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
        <<"-----------------------------------------------------------------------------\n"
        <<"      "<<SNESConvergedReasons[reason]<<"\n\n\n";
        RECORD();
      }
      else // oh, diverged... reduce step and try again
      {

        if ( SolverSpecify::DC_Cycles == 0 )
        {
          MESSAGE <<"------> nonlinear solver " << SNESConvergedReasons[reason] <<". Failed in the first step.\n\n\n";
          RECORD();
          break;
        }
        if ( V_retry.size() >=8 )
        {
          MESSAGE <<"------> nonlinear solver " << SNESConvergedReasons[reason] <<". Too many failed steps, give up tring.\n\n\n";
          RECORD();
          break;
        }

        // load previous result into solution vector
        this->diverged_recovery();

        // reduce step by a factor of 2
        V_retry.push ( Vscan );
        Vscan= ( Vscan+Vs1 ) /2.0;

        MESSAGE <<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
        RECORD();
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

    VecDestroy ( xs1 );
    VecDestroy ( xs2 );
    VecDestroy ( xs3 );

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
    PetscScalar Is1, Is2, Is3;
    std::stack<PetscScalar> I_retry;
    VecDuplicate ( x,&xs1 );
    VecDuplicate ( x,&xs2 );
    VecDuplicate ( x,&xs3 );

    // main loop
    for ( SolverSpecify::DC_Cycles=0;  (Iscan*SolverSpecify::IStep) < SolverSpecify::IStop*SolverSpecify::IStep* ( 1.0+1e-7 ); )
    {
      // show current iscan value
      MESSAGE << "DC Scan: I("  << SolverSpecify::Electrode_IScan[0];
      for ( unsigned int i=1; i<SolverSpecify::Electrode_IScan.size(); i++ )
        MESSAGE << ", "  << SolverSpecify::Electrode_IScan[i];
      MESSAGE << ") = "  << Iscan/PhysicalUnit::A  <<" A" << '\n'
      <<"-----------------------------------------------------------------------------\n";
      RECORD();

      // set iscan current to corresponding electrode
      _system.get_sources()->assign_current_to ( SolverSpecify::Electrode_IScan, Iscan );

      // call pre_solve_process
      if ( SolverSpecify::DC_Cycles == 0 )
        this->pre_solve_process();
      else
        this->pre_solve_process ( false );

      sens_solve();
      // get the converged reason
      SNESConvergedReason reason;
      SNESGetConvergedReason ( snes,&reason );

      if ( reason>0 ) //ok, converged.
      {

        // call post_solve_process
        this->post_solve_process();

        SolverSpecify::DC_Cycles++;

        // save solution for linear/quadratic projection
        VecCopy ( xs2,xs3 ); Is3=Is2;
        VecCopy ( xs1,xs2 ); Is2=Is1;
        VecCopy ( x,xs1 );   Is1=Iscan;

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
        <<"-----------------------------------------------------------------------------\n"
        <<"      "<<SNESConvergedReasons[reason]<<"\n\n\n";
        RECORD();
      }
      else // oh, diverged... reduce step and try again
      {
        if ( SolverSpecify::DC_Cycles == 0 )
        {
          MESSAGE <<"------> nonlinear solver " << SNESConvergedReasons[reason] <<". Failed in the first step.\n\n\n";
          RECORD();
          break;
        }
        if ( I_retry.size() >=8 )
        {
          MESSAGE <<"------> nonlinear solver " << SNESConvergedReasons[reason] <<". Too many failed steps, give up tring.\n\n\n";
          RECORD();
          break;
        }

        // load previous result into solution vector
        this->diverged_recovery();

        // reduce step by a factor of 2
        I_retry.push ( Iscan );
        Iscan= ( Iscan+Is1 ) /2.0;

        MESSAGE <<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
        RECORD();
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

    VecDestroy ( xs1 );
    VecDestroy ( xs2 );
    VecDestroy ( xs3 );
  }

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

#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
    if(Genius::n_processors()>1)
    {
#if defined(PETSC_HAVE_LIBSUPERLU_DIST_2) || defined(PETSC_HAVE_MUMPS)
      ierr = KSPSetType(kspc, KSPPREONLY); genius_assert(!ierr);
      ierr = PCSetType(pcc, PCLU); genius_assert(!ierr);
#ifdef PETSC_HAVE_MUMPS
      ierr = PCFactorSetMatSolverPackage (pcc, MAT_SOLVER_MUMPS); genius_assert(!ierr);
#else
      ierr = PCFactorSetMatSolverPackage (pcc, MAT_SOLVER_SUPERLU_DIST); genius_assert(!ierr);
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
    }
#endif

#if PETSC_VERSION_LE(2,3,3)
    if(Genius::n_processors()>1)
    {
#ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
      ierr = KSPSetType(kspc, KSPPREONLY); genius_assert(!ierr);
      ierr = PCSetType(pcc, PCLU); genius_assert(!ierr);
#else
      // no parallel LU solver? we have to use krylov method for parallel!
      ierr = KSPSetType(kspc, KSPBCGS); genius_assert(!ierr);
      ierr = PCSetType(pcc, PCASM); genius_assert(!ierr);
#endif

    }
    else // we can use LU for sigle CPU
    {
      ierr = KSPSetType(kspc, KSPPREONLY); genius_assert(!ierr);
      ierr = PCSetType(pcc, PCLU); genius_assert(!ierr);
    }
#endif
    ierr = KSPSetOperators(kspc, J, J, SAME_NONZERO_PATTERN); genius_assert(!ierr);
  }
}


/**
 * destroy ksp solver for trace mode
 */
void DDMSolverBase::solve_iv_trace_end()
{
  KSPDestroy(kspc);
  VecDestroy(pdI_pdx);
  VecDestroy(pdF_pdV);
  VecDestroy(pdx_pdV);
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

  const PetscScalar PI = 3.14159265359;

  PetscScalar I=0;
  PetscScalar dI_dV;

  std::string electrode_trace = SolverSpecify::Electrode_VScan[0];
  BoundaryCondition * bc_trace = _system.get_bcs()->get_bc(electrode_trace);

  // the current vscan voltage
  PetscScalar V = SolverSpecify::VStart;
  // set current vscan voltage to corresponding electrode
  _system.get_sources()->assign_voltage_to ( electrode_trace, V );

  // the current vscan step
  PetscScalar VStep = SolverSpecify::VStep;

  // the current potential of vscan electrode
  PetscScalar Potential = bc_trace->ext_circuit()->potential();

  PetscScalar Rload=0;
  PetscScalar Rload_new;
  PetscScalar Rref;

  PetscScalar angle;
  PetscScalar slope=PI/2;
  PetscScalar slope_new;
  PetscScalar slope_chord;

  // set electrode with transient time 0 value of stimulate source(s)
  _system.get_sources()->update ( 0 );
  _system.get_field_source()->update ( 0 );

  // not time dependent
  SolverSpecify::TimeDependent = false;
  SolverSpecify::dt = 1e100;
  SolverSpecify::clock = 0.0;

  solve_iv_trace_begin();


  // output TRACE information
  MESSAGE<<"IV automatically trace by continuation method\n"; RECORD();

  // initial condition check
  bc_trace->ext_circuit()->R() = Rload;

  // call pre_solve_process
  this->pre_solve_process();
  // here call Petsc to solve the nonlinear equations
  sens_solve();

  // get the converged reason
  SNESConvergedReason reason;
  SNESGetConvergedReason ( snes, &reason );

  if(reason<0)
  {
    MESSAGE<<"I can't get convergence even at initial point, please give a good initial condition.\n\n"; RECORD();
    error = 1;
    goto trace_end;
  }

  // call post_solve_process
  this->post_solve_process();


  //loop here
  while(Potential*SolverSpecify::VStep < SolverSpecify::VStop*SolverSpecify::VStep && fabs(I)<SolverSpecify::IStop)
  {
    MESSAGE << "Trace for V(" << electrode_trace << ")=" << Potential << "(V)\n"; RECORD();

    // for last Rload, increase VStep
    int recovery=0;
    do
    {
      V += VStep;
      bc_trace->ext_circuit()->Vapp() =  V;

      this->pre_solve_process ( false );
      sens_solve();

      SNESGetConvergedReason(snes,&reason);
      if(reason<0)
      {
        MESSAGE<<"I can't get convergence at this step, do recovery...\n\n";RECORD();
        diverged_recovery();
        V -= VStep;
        VStep/=2;
        recovery++;
      }
      if(recovery>8)
      {
        MESSAGE<<"Too many failed steps, give up tring.\n\n\n";RECORD();
        error = 1;
        goto trace_end;
      }
    }
    while(reason<0);



    // for the new bias point, recompute slope
    // calculate the dynamic resistance of IV curve by different approximation
    this->set_trace_electrode(bc_trace);

    KSPSolve(kspc, pdF_pdV, pdx_pdV);
    VecDot(pdI_pdx, pdx_pdV, &dI_dV);

    PetscScalar r=1/dI_dV;


    if(first_step)
    {
      Rref = 1.0;
      slope=slope_new=atan(Rref/r);
    }
    else
    {
      // Rref, for scaling the dynamic resistance
      I = bc_trace->ext_circuit()->current_itering();
      Potential = bc_trace->ext_circuit()->potential_itering();

      if( fabs(1e6*I) < fabs(Potential) )
        Rref = 1e6;
      else
        Rref = Potential/I;

      slope_new = atan(Rref/r);
    }


    //resolve turning point
    if(tan(slope)*tan(slope_new)<0) //the sign is changed
    {
      PetscScalar chord = (Potential-bc_trace->ext_circuit()->potential_old())/(I-bc_trace->ext_circuit()->current_old());
      slope_chord = atan(Rref/chord);

      //for vertical turning point
      if(fabs(slope)/PI*180 >70 && fabs(slope_new)/PI*180 >70 &&
          ((slope_chord*slope>0 && fabs(slope_chord)>fabs(slope)) || (slope_chord*slope_new>0 && fabs(slope_chord)>fabs(slope_new))))
      {
        VStep*=-1.0;
      }
      //for horizontal turning point
      if( fabs(slope)/PI*180 <20 && fabs(slope_new)/PI*180 <20 &&
          ((slope_chord*slope>0 && fabs(slope_chord)<fabs(slope)) || (slope_chord*slope_new>0 && fabs(slope_chord)<fabs(slope_new))))
      {
        VStep*=1.0;
      }
    }

    // check the slope change
    angle=fabs(slope-slope_new);
    if(angle>PI/2) angle=PI-angle;

    if(angle<PI/36)      VStep*=1.5;   // slope change less than 5 degree
    else if(angle<PI/24) VStep=VStep;  // slope change less than 10 degree, but greater than 5  degree
    else if(angle<PI/12) VStep/=2;     // slope change less than 15 degree, but greater than 10 degree
    else                                               // slope greater than 15 degree, reject this solution
    {
      //printf("slope change old=%e new=%e angle=%e \n",slope/PI*180, slope_new/PI*180,  angle/PI*180);
      MESSAGE<<"Slope of IV curve changes too quickly, do recovery...\n\n";RECORD();
      VStep/=2;
    }

    // ok, update solutions
    this->post_solve_process();

    // since Rload will be changed, we should change bias voltage to meet the truncation point.

    Rload_new=Rref/r;
    V += I*(Rload_new-Rload);
    bc_trace->ext_circuit()->R() = Rload_new;
    bc_trace->ext_circuit()->Vapp() = V;

    Rload=Rload_new;
    slope = slope_new;
    first_step=0;

  }



trace_end:
  solve_iv_trace_end();

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

  // if BDF2 scheme is used, we should set SolverSpecify::BDF2_restart flag to true
  if ( SolverSpecify::TS_type==SolverSpecify::BDF2 )
    SolverSpecify::BDF2_restart = true;

  // transient simulation clock
  SolverSpecify::clock = SolverSpecify::TStart + SolverSpecify::TStep;

  // for the first step, dt equals TStep
  SolverSpecify::dt = SolverSpecify::TStep;

  MESSAGE<<"Transient compute from "<<SolverSpecify::TStart
  <<" ps step "<<SolverSpecify::TStep
  <<" ps to "  <<SolverSpecify::TStop<<" ps"
  <<'\n';
  RECORD();

  // diverged counter
  int diverged_retry=0;

  // time step counter
  SolverSpecify::T_Cycles=0;

  // the main loop of transient solver.
  do
  {
    MESSAGE
    <<"t = "<<SolverSpecify::clock<<" ps"<<'\n'
    <<"-----------------------------------------------------------------------------\n";
    RECORD();

    //update sources to current clock
    _system.get_sources()->update ( SolverSpecify::clock );
    _system.get_field_source()->update ( SolverSpecify::clock );
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


    //nonlinear solution diverged? try to do recovery
    if ( reason<0 )
    {
      if ( ++diverged_retry >= 8 ) //failed 8 times, stop tring
      {
        MESSAGE
        <<"------> Too many failed steps, give up tring.\n\n\n";
        RECORD();
        break;
      }
      // load previous result into solution vector
      this->diverged_recovery();

      MESSAGE
      <<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
      RECORD();

      // reduce time step by a factor of two, also set clock to next
      SolverSpecify::dt /= 2.0;
      SolverSpecify::clock -= SolverSpecify::dt;

      if ( SolverSpecify::clock < SolverSpecify::TStart )
        SolverSpecify::clock = SolverSpecify::TStart;

      continue;
    }

    //ok, nonlinear solution converged.

    // clear the counter
    diverged_retry = 0;

    //do LTE estimation and auto time step control
    if ( SolverSpecify::AutoStep &&
         ( ( SolverSpecify::TS_type==SolverSpecify::BDF1 && SolverSpecify::T_Cycles>=3 ) ||
           ( SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::T_Cycles>=4 ) ) )
    {
      PetscScalar r = this->LTE_norm();

      if ( SolverSpecify::TS_type==SolverSpecify::BDF1 )
        r = std::pow ( r, PetscScalar ( -1.0/2 ) );
      else if ( SolverSpecify::TS_type==SolverSpecify::BDF2 )
        r = std::pow ( r, PetscScalar ( -1.0/3 ) );

      // when r<0.9, reject this solution
      if ( r<0.9 )
      {
        this->diverged_recovery();

        MESSAGE
        <<"------> LTE too large, time step rejected...\n\n\n";
        RECORD();

        // reduce time step by a factor of r
        SolverSpecify::clock -= SolverSpecify::dt;
        SolverSpecify::dt *= r;
        SolverSpecify::clock += SolverSpecify::dt;

        continue;
      }
      else      // else, accept this solution
      {
        // save time step information
        SolverSpecify::dt_last_last = SolverSpecify::dt_last;
        SolverSpecify::dt_last = SolverSpecify::dt;
        // set next time step
        SolverSpecify::dt = SolverSpecify::dt * ( r > 1.1 ? 1.1 : r );

        // limit the max time step to TStepMax
        if ( SolverSpecify::dt > SolverSpecify::TStepMax )
          SolverSpecify::dt = SolverSpecify::TStepMax;
      }

    }
    else // auto time step control not used
    {
      // save time step information
      SolverSpecify::dt_last_last = SolverSpecify::dt_last;
      SolverSpecify::dt_last = SolverSpecify::dt;

      // set next time step
      if ( fabs ( SolverSpecify::dt ) < fabs ( SolverSpecify::TStep ) )
        SolverSpecify::dt *= 1.1;
    }

    MESSAGE
    <<"-----------------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[reason]<<"\n\n\n";
    RECORD();

    // time step counter ++
    SolverSpecify::T_Cycles++;

    // call post_solve_process
    this->post_solve_process();

    // prepare for next time step
    VecCopy ( x_n1, x_n2 );
    VecCopy ( x_n, x_n1 );
    VecCopy ( x, x_n );

    // predict next solution
    if ( SolverSpecify::Predict &&
         ( ( SolverSpecify::TS_type==SolverSpecify::BDF1 && SolverSpecify::T_Cycles>=3 ) ||
           ( SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::T_Cycles>=4 ) ) )
    {
      PetscScalar hn  = SolverSpecify::dt;           // here dt is the next time step
      PetscScalar hn1 = SolverSpecify::dt_last;      // time step n-1
      PetscScalar hn2 = SolverSpecify::dt_last_last; // time step n-2

      VecZeroEntries ( x );

      if ( SolverSpecify::TS_type == SolverSpecify::BDF1 )
      {
        // use linear interpolation to predict solution x
        VecAXPY ( x, 1+hn/hn1, x_n );
        VecAXPY ( x, -hn/hn1,  x_n1 );
        this->projection_positive_density_check ( x,x_n );
      }
      else if ( SolverSpecify::TS_type == SolverSpecify::BDF2 )
      {
        // use second order polynomial to predict solution x
        PetscScalar cn  = 1+hn* ( hn+2*hn1+hn2 ) / ( hn1* ( hn1+hn2 ) );
        PetscScalar cn1 = -hn* ( hn+hn1+hn2 ) / ( hn1*hn2 );
        PetscScalar cn2 = hn* ( hn+hn1 ) / ( hn2* ( hn1+hn2 ) );

        VecAXPY ( x, cn,  x_n );
        VecAXPY ( x, cn1, x_n1 );
        VecAXPY ( x, cn2, x_n2 );
        this->projection_positive_density_check ( x,x_n );
      }
    }

    // set clock to next time step
    SolverSpecify::clock += SolverSpecify::dt;

    //make sure we can terminat near TStop, the relative error should less than 1e-10.
    if ( SolverSpecify::clock > SolverSpecify::TStop && SolverSpecify::clock < ( SolverSpecify::TStop + SolverSpecify::dt - 1e-10*SolverSpecify::dt ) )
    {
      SolverSpecify::dt -= SolverSpecify::clock - SolverSpecify::TStop;
      SolverSpecify::clock = SolverSpecify::TStop;
    }

    //clear the first step flag of BDF2
    if ( SolverSpecify::TS_type==SolverSpecify::BDF2 )
      SolverSpecify::BDF2_restart = false;

  }
  while ( SolverSpecify::clock < SolverSpecify::TStop+0.5*SolverSpecify::dt );

  // free aux vectors
  VecDestroy ( x_n );
  VecDestroy ( x_n1 );
  VecDestroy ( x_n2 );
  VecDestroy ( xp );
  VecDestroy ( LTE );

  return 0;
}






/*------------------------------------------------------------------
 * snes convergence criteria
 */
#include "private/snesimpl.h"
void DDMSolverBase::petsc_snes_convergence_test ( PetscInt its, PetscReal , PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason )
{
  // update error norm
  this->error_norm();

  *reason = SNES_CONVERGED_ITERATING;

  // the first iteration
  if ( !its )
  {
    snes->ttol = fnorm*snes->rtol;

    MESSAGE<<" "<<"its\t"<<"| Eq(V) | "<<"| Eq(n) | "<<"| Eq(p) | "<<"| Eq(T) | "<<"|Eq(Tn)|  "<<"|Eq(Tp)|  "<<"|delta x|\n"
    <<"-----------------------------------------------------------------------------\n";
    RECORD();
  }

#ifdef CYGWIN
  MESSAGE.precision ( 1 );
#else
  MESSAGE.precision ( 2 );
#endif

  MESSAGE<< "  " << its << "\t" << std::scientific
  << poisson_norm << "  "
  << elec_continuity_norm << "  "
  << hole_continuity_norm << "  "
  << heat_equation_norm   << "  "
  << elec_energy_equation_norm << "  "
  << hole_energy_equation_norm << "  "
  << pnorm << "\n" ;
  RECORD();
  MESSAGE.precision ( 6 );

  // check for NaN (Not a Number)
  if ( fnorm != fnorm )
  {
    *reason = SNES_DIVERGED_FNORM_NAN;
  }
  // check for absolute convergence
  else if ( poisson_norm               < SolverSpecify::poisson_abs_toler         &&
            elec_continuity_norm       < SolverSpecify::elec_continuity_abs_toler &&
            hole_continuity_norm       < SolverSpecify::hole_continuity_abs_toler &&
            electrode_norm             < SolverSpecify::electrode_abs_toler       &&
            heat_equation_norm         < SolverSpecify::heat_equation_abs_toler   &&
            elec_energy_equation_norm  < SolverSpecify::elec_energy_abs_toler     &&
            hole_energy_equation_norm  < SolverSpecify::hole_energy_abs_toler )
  {
    *reason = SNES_CONVERGED_FNORM_ABS;
  }
  else if ( snes->nfuncs >= snes->max_funcs )
  {
    *reason = SNES_DIVERGED_FUNCTION_COUNT;
  }


  if ( *reason == SNES_CONVERGED_ITERATING )
  {
    if ( fnorm <= snes->ttol                                                                             &&
         poisson_norm              < SolverSpecify::toler_relax*SolverSpecify::poisson_abs_toler         &&
         elec_continuity_norm      < SolverSpecify::toler_relax*SolverSpecify::elec_continuity_abs_toler &&
         hole_continuity_norm      < SolverSpecify::toler_relax*SolverSpecify::hole_continuity_abs_toler &&
         electrode_norm            < SolverSpecify::toler_relax*SolverSpecify::electrode_abs_toler       &&
         heat_equation_norm        < SolverSpecify::toler_relax*SolverSpecify::heat_equation_abs_toler   &&
         elec_energy_equation_norm < SolverSpecify::toler_relax*SolverSpecify::elec_energy_abs_toler     &&
         hole_energy_equation_norm < SolverSpecify::toler_relax*SolverSpecify::hole_energy_abs_toler )
    {
      *reason = SNES_CONVERGED_FNORM_RELATIVE;
    }
    // check for relative convergence
    else if ( pnorm                     < SolverSpecify::relative_toler                                       &&
              poisson_norm              < SolverSpecify::toler_relax*SolverSpecify::poisson_abs_toler         &&
              elec_continuity_norm      < SolverSpecify::toler_relax*SolverSpecify::elec_continuity_abs_toler &&
              hole_continuity_norm      < SolverSpecify::toler_relax*SolverSpecify::hole_continuity_abs_toler &&
              electrode_norm            < SolverSpecify::toler_relax*SolverSpecify::electrode_abs_toler       &&
              heat_equation_norm        < SolverSpecify::toler_relax*SolverSpecify::heat_equation_abs_toler   &&
              elec_energy_equation_norm < SolverSpecify::toler_relax*SolverSpecify::elec_energy_abs_toler     &&
              hole_energy_equation_norm < SolverSpecify::toler_relax*SolverSpecify::hole_energy_abs_toler )
    {
      *reason = SNES_CONVERGED_PNORM_RELATIVE;
    }
  }

  // record function norm of this iteration
  function_norm = fnorm;

  return;
}



/*------------------------------------------------------------------
 * ksp convergence criteria
 */
#include "private/kspimpl.h"
void DDMSolverBase::petsc_ksp_convergence_test ( PetscInt its, PetscReal rnorm, KSPConvergedReason* reason )
{

  ksp->rtol = PetscMax ( SolverSpecify::ksp_rtol, 1e-12*std::sqrt ( double ( n_global_dofs ) ) );

  ksp->abstol = PetscMax ( SolverSpecify::ksp_atol_fnorm*function_norm, SolverSpecify::ksp_atol*std::sqrt ( double ( n_global_dofs ) ) );

  KSPDefaultConverged ( ksp, its, rnorm, reason, this );

}

