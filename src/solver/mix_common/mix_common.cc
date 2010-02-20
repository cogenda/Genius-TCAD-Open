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



#include "mix_common/mix_common.h"
#include "parallel.h"
#include "MXMLUtil.h"

#ifndef CYGWIN
#include <unistd.h>

using PhysicalUnit::A;
using PhysicalUnit::s;


static sigjmp_buf net_buf;

/**
 * This function gets called whenever the pipe broken.
 * See the signal(2) manpage for more information.
 */
inline void signal_handler(int signum)
{
  switch (signum)
  {
  case SIGPIPE:
    printf("Received signal SIGPIPE. It seems ngspice terminated.\n");
    siglongjmp(net_buf,1);
  default:
    break;
  }
}

int MixSolverBase::create_solver()
{

  // must setup nonlinear contex here!
  setup_nonlinear_data();

  //abstol = 1e-12*n_global_dofs    - absolute convergence tolerance
  //rtol   = 1e-14                  - relative convergence tolerance
  //stol   = 1e-9                   - convergence tolerance in terms of the norm of the change in the solution between steps
  SNESSetTolerances(snes, 1e-12*n_global_dofs, 1e-14, 1e-9, SolverSpecify::MaxIteration, 1000);

  // rtol   = 1e-12*n_global_dofs  - the relative convergence tolerance (relative decrease in the residual norm)
  // abstol = 1e-20*n_global_dofs  - the absolute convergence tolerance (absolute size of the residual norm)
  KSPSetTolerances(ksp, 1e-12*n_global_dofs, 1e-20*n_global_dofs, PETSC_DEFAULT, std::max(50, static_cast<int>(n_global_dofs/10)));

  // user can do further adjusment from command line
  SNESSetFromOptions (snes);

  // connect to NGSPICE
  connect_to_ngspice();

  // init (user defined) hook functions here
  // set this flag before hook initial, as mix simulation call dc before tran
  // so we set hook always using transient mode
  SolverSpecify::Type = SolverSpecify::TRANSIENT;

  return FVM_NonlinearSolver::create_solver();
}

int MixSolverBase::post_solve_process()
{
  mxml_node_t *eSolution = new_dom_solution_elem();
  const BoundaryConditionCollector * bcs = get_system().get_bcs();

  if ( Genius::processor_id() == 0)
  {
    mxml_node_t *eLabel = mxmlNewElement(eSolution, "label");

    //TODO add label according to the solve type (get from spice)
    mxmlAdd(eLabel, MXML_ADD_AFTER, NULL, MXMLQVariant::makeQVString("solution"));


    mxml_node_t *eTerm = bcs->get_dom_terminal_info();
    if (eTerm)
      mxmlAdd(eSolution, MXML_ADD_AFTER, NULL, eTerm);
  }

  return FVM_NonlinearSolver::post_solve_process();
}

int MixSolverBase::destroy_solver()
{

  // clear nonlinear matrix/vector
  clear_nonlinear_data();

  disconnet_to_ngspice();

  return FVM_NonlinearSolver::destroy_solver();
}

/* ----------------------------------------------------------------------------
 * NDEVServer:  This function set TCP/IP socket and connect to NGSPICE.
 */
int MixSolverBase::NDEVServer()
{
  // only processor 0  can call this function
  genius_assert(Genius::processor_id() == 0);

  char dotted_ip[15];    /* Buffer for converting the resolved address to a readable format. */
  struct sockaddr_in sa; /* Connection address. */

  listener = socket(PF_INET, SOCK_STREAM, IPPROTO_IP);
  if (listener < 0)
  {
    MESSAGE<<"Unable to create a listener socket: "<<strerror(errno)<<"\n";
    RECORD();
    return 1;
  }

  // Now bind the listener to a local address. This uses the same sockaddr_in structure as connect
  socklen_t sa_len;
  sa_len = sizeof(sa);
  memset(&sa,0,sa_len);
  sa.sin_family = AF_INET;
  sa.sin_port = htons(port);
  sa.sin_addr.s_addr = htonl(INADDR_ANY); /* Listen on all interfaces. */
  if (bind(listener, (sockaddr *)&sa, sa_len) < 0)
  {
    MESSAGE<<"Unable to bind to port "<<port<<": "<<strerror(errno)<<"\n";
    RECORD();
    return 1;
  }

  /**
   * Let the networking system know we're accepting
   * connections on this socket. Ask for a connection
   * queue of five clients. (If more than five clients
   * try to connect before we call accept, some will
   * be denied.)
   */
  if (listen(listener, 1) < 0)
  {
    MESSAGE<<"Unable to listen: "<<strerror(errno)<<std::endl;
    RECORD();
    return 1;
  }
  MESSAGE<<"Waiting for ngspice...";
  RECORD();
  client = accept(listener, (sockaddr *)&sa, &sa_len);

  /**
   *  We now have a live client. Print information
   *  about it and then send something over the wire.
   */
  inet_ntop(AF_INET, &sa.sin_addr, dotted_ip, 15);

  // socket buffer
  char sendbuf[128];

  /* Waiting for ngspice query data */
  recv(client, sendbuf, 128, MSG_FLAG);
  if(strcmp(sendbuf, NG_QUERY)) /* Query error */
  {
    MESSAGE<<"Query information error from ngspice.\n";
    RECORD();
    return 1;
  }
  MESSAGE<<"Received ngspice connection from "<<dotted_ip<<std::endl;
  RECORD();

  /* ACK to ngspice */
  sprintf(sendbuf, NDEV_REPLY);
  send(client, sendbuf, 128, 0);

  return 0;

}




int MixSolverBase::connect_to_ngspice()
{
  // set signal SIGPIPE call back function
  signal(SIGPIPE, signal_handler);

  //connect to ngspice, only process 0 can do it!
  if(Genius::processor_id() == 0)
  {
    if( NDEVServer() )
    {
      MESSAGE<<"I can't connect to ngspice, mix simulation stopped.\n";
      RECORD();
      genius_error();
    }

    //get ngspice terminal information
    recv(client, &(Deviceinfo), sizeof(sDeviceinfo), MSG_FLAG);
  }

  // synchronize ngspice terminal information with all the processors
  MPI_Bcast (&(Deviceinfo), sizeof(sDeviceinfo), MPI_BYTE, 0, PETSC_COMM_WORLD);

  MESSAGE<<"Ngspice demands device "<<Deviceinfo.NDEVname<<" to be calculated.\n";
  RECORD();

  for(int i=0; i<Deviceinfo.term; i++)
  {
    // get ngspice pin information
    if(Genius::processor_id() == 0)
      recv(client, &(PINinfos[i]), sizeof(sPINinfo), MSG_FLAG);

    // synchronize ngspice pin information with all the processors
    MPI_Bcast (&(PINinfos[i]), sizeof(sPINinfo), MPI_BYTE, 0, PETSC_COMM_WORLD);

    // the pin name is case insensitivity in ngspice
    BoundaryCondition * bc = _system.get_bcs()->get_bc_nocase(PINinfos[i].name);
    if( bc == NULL)
    {
      MESSAGE<<"I can't link ngspice terminal name "<<PINinfos[i].name<<" to Genius eletrode boundary.\n";
      RECORD();
      return 1;
    }

    // map pin name (use bc label as unique name) to pin index for later usage
    pin_name_to_pin_index_map[bc->label()] = i;

    MESSAGE<<"   Terminal "<<i<<": "<<PINinfos[i].name<<"\n";
    RECORD();
  }

  // build required PIN data structure
  // this should be done for all the processors
  for(int i=0; i<Deviceinfo.term; i++)
  {
    VecDuplicate(x, &PINconds[i].pdI_pdw);
    VecDuplicate(x, &PINconds[i].pdF_pdV);
    VecDuplicate(x, &PINconds[i].pdw_pdV);
  }

  // build special linear solver contex
  {
    PetscErrorCode ierr;

    // create the jacobian matrix
    ierr = MatCreate(PETSC_COMM_WORLD, &G); genius_assert(!ierr);

    ierr = MatSetSizes(G, n_local_dofs, n_local_dofs, n_global_dofs, n_global_dofs); genius_assert(!ierr);

#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
    if (Genius::n_processors()>1)
    {
      ierr = MatSetType(G, MATMPIAIJ); genius_assert(!ierr);
      ierr = MatMPIAIJSetPreallocation(G, 0, &n_nz[0], 0, &n_oz[0]); genius_assert(!ierr);
    }
    else
    {
      ierr = MatSetType(G, MATSEQAIJ); genius_assert(!ierr);
      ierr = MatSeqAIJSetPreallocation(G, 0, &n_nz[0]); genius_assert(!ierr);
    }
#endif

#if PETSC_VERSION_LE(2,3,3)
    if (Genius::n_processors()>1)
    {
 #ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
      ierr = MatSetType(G, MATSUPERLU_DIST); genius_assert(!ierr);
      ierr = PetscOptionsSetValue("-mat_superlu_dist_rowperm","NATURAL");
 #else
      ierr = MatSetType(G, MATMPIAIJ); genius_assert(!ierr);
 #endif
      ierr = MatMPIAIJSetPreallocation(G, 0, &n_nz[0], 0, &n_oz[0]); genius_assert(!ierr);
    }
    else
    {
      ierr = MatSetType(G, MATSEQAIJ); genius_assert(!ierr);
      ierr = MatSeqAIJSetPreallocation(G, 0, &n_nz[0]); genius_assert(!ierr);
    }
#endif

    // indicates when MatZeroRows() is called the zeroed entries are kept in the nonzero structure
#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
    ierr = MatSetOption(G, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE); genius_assert(!ierr);
#endif

#if PETSC_VERSION_LE(2,3,3)
    ierr = MatSetOption(G, MAT_KEEP_ZEROED_ROWS); genius_assert(!ierr);
#endif

    g_matrix_first_assemble = false;

    ierr = KSPCreate(PETSC_COMM_WORLD, &linear_lu); genius_assert(!ierr);

    ierr = KSPGetPC(linear_lu, &pc_direct); genius_assert(!ierr);

#if (PETSC_VERSION_GE(3,0,0) || defined(HAVE_PETSC_DEV))
    if(Genius::n_processors()>1)
    {
 #if defined(PETSC_HAVE_LIBSUPERLU_DIST_2) || defined(PETSC_HAVE_MUMPS)
      ierr = KSPSetType(linear_lu, KSPPREONLY); genius_assert(!ierr);
      ierr = PCSetType(pc_direct, PCLU); genius_assert(!ierr);
  #ifdef PETSC_HAVE_MUMPS
      ierr = PCFactorSetMatSolverPackage (pc_direct, MAT_SOLVER_MUMPS); genius_assert(!ierr);
  #else
      ierr = PCFactorSetMatSolverPackage (pc_direct, MAT_SOLVER_SUPERLU_DIST); genius_assert(!ierr);
  #endif
 #else
      // no parallel LU solver? we have to use krylov method for parallel!
      ierr = KSPSetType(linear_lu, KSPBCGS); genius_assert(!ierr);
      ierr = PCSetType(pc_direct, PCASM); genius_assert(!ierr);
 #endif
    }
    else
    {
      ierr = KSPSetType(linear_lu, KSPPREONLY); genius_assert(!ierr);
      ierr = PCSetType(pc_direct, PCLU); genius_assert(!ierr);
    }
#endif

#if PETSC_VERSION_LE(2,3,3)
    if(Genius::n_processors()>1)
    {
 #ifdef PETSC_HAVE_LIBSUPERLU_DIST_2
      ierr = KSPSetType(linear_lu, KSPPREONLY); genius_assert(!ierr);
      ierr = PCSetType(pc_direct, PCLU); genius_assert(!ierr);
 #else
      // no parallel LU solver? we have to use krylov method for parallel!
      ierr = KSPSetType(linear_lu, KSPBCGS); genius_assert(!ierr);
      ierr = PCSetType(pc_direct, PCASM); genius_assert(!ierr);
 #endif

    }
    else // we can use LU for sigle CPU
    {
      ierr = KSPSetType(linear_lu, KSPPREONLY); genius_assert(!ierr);
      ierr = PCSetType(pc_direct, PCLU); genius_assert(!ierr);
    }
#endif

    ierr = KSPSetOperators(linear_lu, G, G, SAME_NONZERO_PATTERN); genius_assert(!ierr);

  }

  return 0;
}




/* ----------------------------------------------------------------------------
 * This function receive ngspice device information.
 */
int MixSolverBase::run_under_ngspice()
{
  SolverSpecify::clock = -1e100;

  int bytes;


  this->init_spice_data();
  spice_dt_last=0.0;
  VecDuplicate(x, &xs_n);
  VecDuplicate(x, &xs_n1);

  for(;;)
  {
    if (sigsetjmp(net_buf, 1) == 1) break;

    //receive solver information
    if(Genius::processor_id() == 0)
      bytes=recv(client, &(CKTInfo), sizeof(sCKTinfo), MSG_FLAG);

    // synchronize solver information with all the processors
    MPI_Bcast (&(CKTInfo), sizeof(sCKTinfo), MPI_BYTE, 0, PETSC_COMM_WORLD);
    MPI_Bcast (&bytes, 1, MPI_INT, 0, PETSC_COMM_WORLD);

    //if ngspice terminated, stop loop.
    if( bytes < static_cast<int>(sizeof(sCKTinfo)) )
    {
      MESSAGE<<"Lost connection to remote ngspice.\n";
      RECORD();
      break;
    }

    switch(CKTInfo.DEV_CALL)
    {
    case   NDEV_LOAD                : this->DEV_LOAD();      break;
    case   NDEV_ACCEPT              : this->DEV_ACCEPT();    break;
    case   NDEV_CONVERGINCE_TEST    : this->DEV_CONV_TEST(); break;
    default : break;
    }
  }

  return 0;
}



/* ----------------------------------------------------------------------------
 * This function close connection to ngspice
 */
int MixSolverBase::disconnet_to_ngspice()
{
  // free work space

  for(int i=0; i<Deviceinfo.term; i++)
  {
    VecDestroy(PINconds[i].pdI_pdw);
    VecDestroy(PINconds[i].pdF_pdV);
    VecDestroy(PINconds[i].pdw_pdV);
  }

  if(Genius::processor_id() == 0)
  {
    close(client);
    close(listener);
  }

  // free extra ksp and matrix
  KSPDestroy(linear_lu);
  MatDestroy(G);

  return 0;
}



int  MixSolverBase::DEV_LOAD()
{
  int Converged = 0;

  //receive terminal voltage of device
  for(int i=0;i<Deviceinfo.term;i++)
  {
    if(Genius::processor_id() == 0)
      recv(client, &PINinfos[i], sizeof(sPINinfo), MSG_FLAG);

    // synchronize ngspice pin information with all the processors
    MPI_Bcast (&(PINinfos[i]), sizeof(sPINinfo), MPI_BYTE, 0, PETSC_COMM_WORLD);

    _system.get_bcs()->set_vapp_nocase(PINinfos[i].name, PINinfos[i].V);
  }

  // load previous solution data
  this->load_spice_data();
  SolverSpecify::dt_last = spice_dt_last;

  //transient calculation settings
  if(CKTInfo.CKTmode & MODETRAN)
  {
    SolverSpecify::TimeDependent = true;
    SolverSpecify::Type = SolverSpecify::TRANSIENT;
    SolverSpecify::TS_type=SolverSpecify::BDF2;

    //the first step of transient simulation
    if(CKTInfo.CKTmode & MODEINITTRAN)
    {
      MESSAGE<<"-----------------------------------------------------------------------------\n";
      MESSAGE<<"NGSPICE Transient Mode Start   Time = "<<CKTInfo.time<<"s  dt = "<<CKTInfo.dt<<"s\n";
      MESSAGE<<"-----------------------------------------------------------------------------\n\n";
      RECORD();
      SolverSpecify::BDF2_restart = true;
    }

    //the time matching flag
    else if(CKTInfo.CKTmode & MODEINITPRED)
    {
      MESSAGE<<"-----------------------------------------------------------------------------\n";
      MESSAGE<<"NGSPICE Transient Mode Update  Time = "<<CKTInfo.time<<"s  dt = "<<CKTInfo.dt<<"s\n";
      if(SolverSpecify::clock > CKTInfo.time*s)
      {
        SolverSpecify::BDF2_restart = true;
        MESSAGE<<"NGSPICE back trace...\n";
        RECORD();
      }
      else
      {
        SolverSpecify::BDF2_restart = false;
      }
      MESSAGE<<"-----------------------------------------------------------------------------\n\n";
      RECORD();
    }

    SolverSpecify::clock = CKTInfo.time*s;

    // call the real computation routine
    Converged = this->tran_solve();
  }

  //DC calculation settings
  else if(CKTInfo.CKTmode & MODEDC)
  {
    MESSAGE<<"-----------------------------------------------------------------------------\n";
    MESSAGE<<"NGSPICE DC Mode\n";
    MESSAGE<<"-----------------------------------------------------------------------------\n\n";
    RECORD();

    SolverSpecify::TimeDependent = false;
    SolverSpecify::Type = SolverSpecify::DCSWEEP;
    SolverSpecify::dt = 1e100;

    // call the real computation routine
    Converged = this->dc_solve();
  }

  if (!Converged)
    return 0;

  //get pdI/pdw, pdI/pdV and pdF/pdV for each electrode
  this->get_electrode_load();

  /*
   * compute conductance and rhs current source.
   */

  // we should call KSPSetOperators here, or the result is uncertain!
  // does some global data was changed by petsc?
  KSPSetOperators(linear_lu, G, G, SAME_NONZERO_PATTERN);

  for(int i=0; i<Deviceinfo.term; i++)
  {
    PetscScalar Ig=0;
    for(int j=0; j<Deviceinfo.term; j++)
    {
      KSPSolve(linear_lu, PINconds[j].pdF_pdV, PINconds[i].pdw_pdV);
      PetscScalar pdI_pdV;
      VecDot(PINconds[i].pdI_pdw, PINconds[i].pdw_pdV, &pdI_pdV);
      PINinfos[i].dI_dV[j] = pdI_pdV;
      PINinfos[i].dI_dV[j] += PINconds[i].pdI_pdV;
      PINinfos[i].dI_dV[j] /= A; //scale the element of dI/dV.
      Ig += PINinfos[i].dI_dV[j]*PINinfos[j].V;
    }
    PINinfos[i].I = Ig - PINinfos[i].I/A;
  }

  /*
   * sent conductance matrix and rhs current back to ngspice.
   */
  if(Genius::processor_id() == 0)
  {
    for(int i=0; i<Deviceinfo.term; i++)
    {
      send(client,&PINinfos[i],sizeof(sPINinfo),0);
    }
  }

  return 0;
}



int  MixSolverBase::DEV_ACCEPT()
{
  MESSAGE
    <<"-----------------------------------------------------------------------------\n"
    <<"NGSPICE accepted solution. Saving solution.\n"
    <<"-----------------------------------------------------------------------------\n\n";
  RECORD();
  this->save_spice_data();
  spice_dt_last = CKTInfo.dt_old*s;

  VecCopy(xs_n,xs_n1);
  VecCopy(x,xs_n);

  // call (user defined) hook function hook_post_solve_process
  hook_list()->post_solve();

  return 0;
}



int  MixSolverBase::DEV_CONV_TEST()
{
  // get the converged reason
  SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
  SNESGetConvergedReason(snes,&reason);

  if(Genius::processor_id() == 0)
  {
    CKTInfo.convergence_flag = reason;
    send(client,&(CKTInfo),sizeof(sCKTinfo),0);
  }

  return 0;
}



int MixSolverBase::tran_solve()
{
  PetscInt rework = 1;
  PetscInt max_rework = 64;
  PetscInt Converged = 0;

  int istep=0;

  SNESConvergedReason reason;
  do
  {
    reason = SNES_CONVERGED_ITERATING;
    for(int w=1; w<=rework; w++)
    {
      SolverSpecify::dt = CKTInfo.dt*s/rework;

      MESSAGE<<"Solving: ";
      for(int i=0; i<Deviceinfo.term; i++)
      {
        PetscScalar V_step = PINinfos[i].V_old+(PINinfos[i].V-PINinfos[i].V_old)*w/rework;
        _system.get_bcs()->set_vapp_nocase(PINinfos[i].name,V_step);
        MESSAGE<<PINinfos[i].name<<":"<<V_step<<"V ";
      }
      MESSAGE<<"with time step "<<SolverSpecify::dt/s<<"s"<<std::endl;
      RECORD();

      // call pre_solve_process
      if(istep++ == 0)
        this->pre_solve_process();
      else
        this->pre_solve_process(false);

      SNESSolve(snes,PETSC_NULL,x);

      SNESGetConvergedReason(snes,&reason);
      if(reason<0)
      {
        MESSAGE<<"------> Genius mixed solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
        RECORD();
        this->diverged_recovery();
        rework*=2;
        break;
      }
      // print convergence/divergence reason
      if(reason>0)
      {
        this->post_solve_process();
        MESSAGE
        <<"-----------------------------------------------------------------------------\n"
        <<"      "<<SNESConvergedReasons[reason]<<"\n\n\n";
        RECORD();
      }

    }

    if(reason>0)          {Converged = 1; break;}

    if(rework>max_rework) {
      // Converged = 0; break;
      MESSAGE<<"Genius does not converge with "<<max_rework<<" sub-steps. Quit...\n";
      RECORD();
      genius_error();
    }
    else
    {
      MESSAGE<<"Retry with "<<rework<<" sub-steps."<<std::endl;
    }

  }
  while(1);

  CKTInfo.convergence_flag = reason;
  send(client,&(CKTInfo),sizeof(sCKTinfo),0);

  return Converged;
}



int MixSolverBase::dc_solve()
{
  PetscInt rework = 1;
  PetscInt max_rework = 64;
  PetscInt Converged = 0;

  int istep=0;

  SNESConvergedReason reason;
  do
  {
    reason = SNES_CONVERGED_ITERATING;

    MESSAGE<<"Solving: ";
    for(int w=1;w<=rework;w++)
    {
      for(int i=0;i<Deviceinfo.term;i++)
      {
        PetscScalar V_step = PINinfos[i].V_old+(PINinfos[i].V-PINinfos[i].V_old)*w/rework;
        _system.get_bcs()->set_vapp_nocase(PINinfos[i].name,V_step);
        MESSAGE<<PINinfos[i].name<<":"<<V_step<<"V ";
      }
      MESSAGE<<std::endl;
      RECORD();

      // call pre_solve_process
      if(istep++ == 0)
      {
        if (SolverSpecify::dt_last>0)
        {
          this->pre_solve_process(false);
          PetscScalar hn  = SolverSpecify::dt;           // here dt is the next time step
          PetscScalar hn1 = SolverSpecify::dt_last;      // time step n-1
          VecAXPY(x, 1+hn/hn1, xs_n);
          VecAXPY(x, -hn/hn1,  xs_n1);
          this->projection_positive_density_check(x,xs_n);
        }
        else
          this->pre_solve_process();
      }
      else
        this->pre_solve_process(false);

      SNESSolve(snes,PETSC_NULL,x);

      SNESGetConvergedReason(snes,&reason);
      if(reason<0)
      {
        MESSAGE<<"------> Genius mixed solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
        RECORD();
        this->diverged_recovery();
        rework*=2;
        break;
      }
      // print convergence/divergence reason
      if(reason>0)
      {
        MESSAGE
        <<"-----------------------------------------------------------------------------\n"
        <<"      "<<SNESConvergedReasons[reason]<<"\n\n\n";
        RECORD();
      }
    }
    if(reason>0) {Converged = 1; break;}
    if(rework>max_rework) {}
    if(rework>max_rework) {
      // Converged = 0; break;
      MESSAGE<<"Genius does not converge with "<<max_rework<<" sub-steps. Reporting failure to spice...\n";
      RECORD();
      Converged = 0; break;
    }
    else
    {
      MESSAGE<<"Retry with "<<rework<<" sub-steps."<<std::endl;
    }
  }
  while(1);

  CKTInfo.convergence_flag = reason;
  send(client,&(CKTInfo),sizeof(sCKTinfo),0);
  return Converged;
}


/*------------------------------------------------------------------
 * snes convergence criteria
 */
#include "private/snesimpl.h"
void MixSolverBase::petsc_snes_convergence_test(PetscInt its, PetscReal , PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason)
{
  // update error norm
  this->error_norm();

  *reason = SNES_CONVERGED_ITERATING;

  // the first iteration
  if (!its)
  {
    snes->ttol = fnorm*snes->rtol;

    MESSAGE<<" "<<"its\t"<<"| Eq(V) | "<<"| Eq(n) | "<<"| Eq(p) | "<<"| Eq(T) | "<<"|Eq(Tn)|  "<<"|Eq(Tp)|  "<<"|delta x|\n"
                           <<"-----------------------------------------------------------------------------\n";
    RECORD();
  }

  MESSAGE.precision(2);
  MESSAGE<< "  " << its << "\t"
                         << poisson_norm << "  "
                         << elec_continuity_norm << "  "
                         << hole_continuity_norm << "  "
                         << heat_equation_norm  << "  "
                         << elec_energy_equation_norm << "  "
                         << hole_energy_equation_norm << "  "
                         << pnorm << "\n" ;
  RECORD();
  MESSAGE.precision(6);

  // check for NaN (Not a Number)
  if (fnorm != fnorm)
  {
    *reason = SNES_DIVERGED_FNORM_NAN;
  }
  // check for absolute convergence
  else if (poisson_norm              < SolverSpecify::poisson_abs_toler         &&
           elec_continuity_norm      < SolverSpecify::elec_continuity_abs_toler &&
           hole_continuity_norm      < SolverSpecify::hole_continuity_abs_toler &&
           heat_equation_norm        < SolverSpecify::heat_equation_abs_toler   &&
           elec_energy_equation_norm < SolverSpecify::elec_energy_abs_toler     &&
           hole_energy_equation_norm < SolverSpecify::hole_energy_abs_toler)
  {
    *reason = SNES_CONVERGED_FNORM_ABS;
  }
  else if (snes->nfuncs >= snes->max_funcs)
  {
    *reason = SNES_DIVERGED_FUNCTION_COUNT;
  }


  if (its && !*reason)
  {
    if (fnorm <= snes->ttol)
    {
      *reason = SNES_CONVERGED_FNORM_RELATIVE;
    }
    // check for relative convergence
    else if (pnorm                     < SolverSpecify::relative_toler                                       &&
             poisson_norm              < SolverSpecify::toler_relax*SolverSpecify::poisson_abs_toler         &&
             elec_continuity_norm      < SolverSpecify::toler_relax*SolverSpecify::elec_continuity_abs_toler &&
             hole_continuity_norm      < SolverSpecify::toler_relax*SolverSpecify::hole_continuity_abs_toler &&
             heat_equation_norm        < SolverSpecify::toler_relax*SolverSpecify::heat_equation_abs_toler   &&
             elec_energy_equation_norm < SolverSpecify::toler_relax*SolverSpecify::elec_energy_abs_toler     &&
             hole_energy_equation_norm < SolverSpecify::toler_relax*SolverSpecify::hole_energy_abs_toler)
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
void MixSolverBase::petsc_ksp_convergence_test(PetscInt its, PetscReal rnorm, KSPConvergedReason* reason)
{

  ksp->rtol = PetscMax(SolverSpecify::ksp_rtol, 1e-12*std::sqrt(double(n_global_dofs)));

  ksp->abstol = PetscMax(SolverSpecify::ksp_atol_fnorm*function_norm, SolverSpecify::ksp_atol*std::sqrt(double(n_global_dofs)));

  KSPDefaultConverged(ksp, its, rnorm, reason, this);

}

#endif //#ifndef CYGWIN
