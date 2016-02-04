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


#include <stack>
#include <iomanip>

#include "solver_specify.h"
#include "physical_unit.h"
#include "field_source.h"
#include "mixA_solver.h"
#include "spice_ckt_define.h"
#include "spice_ckt.h"
#include "parallel.h"
#include "petsc_utils.h"


using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::W;
using PhysicalUnit::C;
using PhysicalUnit::s;
using PhysicalUnit::um;

int MixASolverBase::create_solver()
{

  link_electrode_to_spice_node();

  set_nonlinear_solver_type ( SolverSpecify::NS );
  set_linear_solver_type    ( SolverSpecify::LS );
  set_preconditioner_type   ( SolverSpecify::PC );

  // must setup nonlinear contex here!
  setup_nonlinear_data();

  //abstol = 1e-12*n_global_dofs    - absolute convergence tolerance
  //rtol   = 1e-10                  - relative convergence tolerance
  //stol   = 1e-9                   - convergence tolerance in terms of the norm of the change in the solution between steps
  SNESSetTolerances(snes, 1e-12*n_global_dofs, 1e-10, 1e-9, SolverSpecify::MaxIteration, 1000);

  // rtol   = SolverSpecify::ksp_rtol  - the relative convergence tolerance (relative decrease in the residual norm)
  // abstol = 1e-20                    - the absolute convergence tolerance (absolute size of the residual norm)
  KSPSetTolerances(ksp, SolverSpecify::ksp_rtol, 1e-20, PETSC_DEFAULT, std::max(50, std::min(1000, static_cast<int>(n_global_dofs/10))) );

  // user can do further adjusment from command line
  SNESSetFromOptions (snes);

  return FVM_FlexNonlinearSolver::create_solver();
}


int MixASolverBase::destroy_solver()
{

  // clear nonlinear matrix/vector
  clear_nonlinear_data();

#if defined(HAVE_FENV_H)
  feclearexcept(FE_INVALID);
#endif

  return FVM_FlexNonlinearSolver::destroy_solver();
}






//FIXME the default LU solver of petsc has many problems with mixA solver
// one must use PCFactorSetShiftNonzero to avoid zero pivot -- however, nonlinear convergence maybe poor with this argument
// and PCFactorReorderForNonzeroDiagonal may cause a segment fault error.
// As a result, superlu, umfpack , bcgs_ilu are recommned for serial problems
// mumps or superlu_dist for parallel problems

void MixASolverBase::link_electrode_to_spice_node()
{
  _circuit->clear_electrode();
  for(unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(n);
    bc->set_spice_electrode(false);
  }

  const std::map<std::string, unsigned int> & electrode_to_spice_node = _circuit->get_electrode_info();
  std::map<std::string, unsigned int>::const_iterator it = electrode_to_spice_node.begin();
  for(; it!=electrode_to_spice_node.end(); ++it)
  {
    std::vector<BoundaryCondition *> bcs =  _system.get_bcs()->get_bcs_by_electrode_label_nocase(it->first);
    for(size_t b=0; b<bcs.size(); b++)
    {
      BoundaryCondition * bc = bcs[b];
      if(bc && bc->is_electrode())
      {
        bc->set_spice_electrode(true);
        _circuit->set_ckt_node_electrode_flag(it->second);
        _circuit->link_electrode(it->second, bc);
      }
    }
  }

  // if some electrode bc not linked to spice node?
  for(unsigned int n=0; n<_system.get_bcs()->n_bcs(); ++n)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(n);
    if( !bc->is_electrode() ) continue;  //skip non electrode
    if( bc->is_spice_electrode() ) continue; // skip electrode already link to spice
    if( bc->is_inter_connect_bc() && !bc->is_inter_connect_hub() ) continue; // skip electrode belongs to interconnect

    if(_circuit->get_spice_node_by_bc(const_cast<const BoundaryCondition *>(bc))==invalid_uint)
    {
      MESSAGE<<"Warning: Electrode "<<bc->label()<<" not linked to SPICE, set it to ground.\n";
      RECORD();
      bc->set_spice_electrode(true);
      _circuit->set_ckt_node_electrode_flag(0);
      _circuit->link_electrode(0, bc);
    }
  }

}



/*------------------------------------------------------------------
 * set the matrix nonzero pattern for spice circuit
 */
void MixASolverBase::set_extra_matrix_nonzero_pattern()
{

  DDMSolverBase::set_extra_matrix_nonzero_pattern();

  //spice ckt should know offset in petsc matrix and local scatter vec
  _circuit->set_offset(n_global_dofs-this->extra_dofs(),
                       local_index_array.size()-this->extra_dofs(),
                       n_local_dofs-this->extra_dofs());

#if 0 
  // extra dofs for spice matrix
  unsigned int n_extra_dofs = this->extra_dofs();
  unsigned int n_max_spice_nonzeros = _circuit->max_row_nonzeros();
  unsigned int n_electrode = _circuit->n_electrode();
  for(unsigned int n=0; n<n_extra_dofs; ++n)
  {
    // nnz and noz of spice circuit
    unsigned int on_processor_dofs = _system.get_bcs()->n_electrode_bcs() + n_max_spice_nonzeros + n_electrode; //slightly overkill
    unsigned int off_processor_dofs = 0;

    // only the last processor own this extra dofs
    if(Genius::is_last_processor())
    {
      n_nz[n_local_dofs - n_extra_dofs + n] += on_processor_dofs;
      n_oz[n_local_dofs - n_extra_dofs + n] += off_processor_dofs;
    }
  }

  // extra dofs for coupled electrode
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); ++b)
  {
    const BoundaryCondition * bc = _system.get_bcs()->get_bc(b);

    //this spice circuit node links to a electrode bc?
    if(bc->is_spice_electrode())
    {
      unsigned int spice_node_index = _circuit->get_spice_node_by_bc(bc);
      genius_assert(spice_node_index!=invalid_uint);

      // get the nodes belongs to this boundary condition
      const std::vector<const Node *> & bc_nodes = bc->nodes();

      // statistic neighbor information of boundary node
      std::vector<unsigned int> neighbors;
      for(unsigned int n=0; n<bc_nodes.size(); ++n  )
      {
        if(bc_nodes[n]->processor_id() == Genius::processor_id())
          neighbors.push_back( bc->n_node_neighbors(bc_nodes[n]) );
      }

      // statistic the on- and off- processor matrix bandwidth contributed by boundary node
      unsigned int on_processor_dofs  = 0;
      for(unsigned int n=0; n<bc_nodes.size(); ++n  )
      {
        if( bc_nodes[n]->processor_id() == Genius::processor_id() )
          on_processor_dofs += (neighbors[n]+1)*this->n_bc_node_dofs( bc ) ;
      }

      // total dofs equals to the sum of all the on_processor_dofs
      unsigned int total_dofs = on_processor_dofs;
      Parallel::sum(total_dofs);
      unsigned int off_processor_dofs = total_dofs - on_processor_dofs;

      // prevent overflow, this may be happened for very small problems.
      if ( on_processor_dofs  > n_local_dofs )
      {
        on_processor_dofs = n_local_dofs ;
      }

      if ( on_processor_dofs + off_processor_dofs > n_global_dofs )
      {
        off_processor_dofs = n_global_dofs - on_processor_dofs;
      }

      // set matrix entry for device to circuit node
      for(unsigned int nb=0; nb<bc_nodes.size(); ++nb  ) //for all the nodes on electrode boundary
      {
        const Node * bd_node = bc_nodes[nb];
        if( bd_node->processor_id()!=Genius::processor_id() )  continue;

        // get the region/fvm_node information
        BoundaryCondition::const_region_node_iterator  rnode_it     = bc->region_node_begin(bd_node);
        BoundaryCondition::const_region_node_iterator  end_rnode_it = bc->region_node_end(bd_node);
        for(; rnode_it!=end_rnode_it; ++rnode_it  )
        {
          const FVM_Node * bd_fvm_node = (*rnode_it).second.second;
          const SimulationRegion * region = (*rnode_it).second.first;
          unsigned int local_node_dofs = this->node_dofs( region );
          for(unsigned int i=0; i<local_node_dofs; ++i)
          {
            if(Genius::is_last_processor())
              n_nz[bd_fvm_node->local_offset() + i] += 1;
            else
              n_oz[bd_fvm_node->local_offset() + i] += 1;
          }
        }
      }

      // set matrix entry for circuit node to device variable
      // only the last processor own this extra dofs
      if(Genius::is_last_processor())
      {
        n_nz[_circuit->array_offset_f(spice_node_index)] += on_processor_dofs;
        n_oz[_circuit->array_offset_f(spice_node_index)] += off_processor_dofs;
      }

    }

  }
#endif

  // set the global/local offset of spice current conservation equation
  for(unsigned int b=0; b<_system.get_bcs()->n_bcs(); ++b)
  {
    BoundaryCondition * bc = _system.get_bcs()->get_bc(b);

    //is this electrode bc links to a spice circuit node?
    if(bc->is_spice_electrode())
    {
      unsigned int spice_node_index = _circuit->get_spice_node_by_bc(bc);
      bc->set_global_offset( _circuit->global_offset_f(spice_node_index) );
      bc->set_local_offset ( _circuit->local_offset_f(spice_node_index) );
    }
  }
}






/*------------------------------------------------------------------
 * return the extra dofs of spice circuit
 */
unsigned int MixASolverBase::extra_dofs() const
{
  return _circuit->n_ckt_nodes();
}



void MixASolverBase::spice_fill_value(Vec x, Vec L)
{
  if(Genius::is_last_processor())
  {
    _circuit->do_ic(SolverSpecify::UIC);
    _circuit->do_node_set(SolverSpecify::NodeSet);

    std::vector<int> ix;
    std::vector<PetscScalar> y;
    std::vector<PetscScalar> s;

    for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
    {
      ix.push_back(_circuit->global_offset_x(n));
      y.push_back(_circuit->rhs_old(n));
      s.push_back(1.0);
    }

    if( ix.size() )
    {
      VecSetValues(x, ix.size(), &ix[0], &y[0], INSERT_VALUES) ;
      VecSetValues(L, ix.size(), &ix[0], &s[0], INSERT_VALUES) ;
    }
  }
}






int MixASolverBase::pre_solve_process(bool load_solution)
{
  // spice voltage source/current source maybe changed, here we load solution vector again
  if(Genius::is_last_processor())
  {
    std::vector<int> ix;
    std::vector<PetscScalar> y;

    for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
    {
      ix.push_back(_circuit->global_offset_x(n));
      y.push_back(_circuit->rhs_old(n));
    }

    if( ix.size() )
    {
      VecSetValues(x, ix.size(), &ix[0], &y[0], INSERT_VALUES) ;
    }
  }

  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  return DDMSolverBase::pre_solve_process(load_solution);
}




void MixASolverBase::build_spice_function(PetscScalar *lxx, Vec f, InsertMode &add_value_flag)
{

  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // only the last processor do this
  if(Genius::is_last_processor())
  {
    // ask spice to build new rhs and matrix
    _circuit->circuit_load();
    //_circuit->print_ckt_matrix();

    std::vector<PetscInt> iy;
    std::vector<PetscScalar> y;
    for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
    {
      PetscInt global_index;
      Real value;
      _circuit->ckt_residual(n, global_index, value);
      iy.push_back(global_index);
      y.push_back(value);
    }

    if(iy.size())
      VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}


void MixASolverBase::build_spice_jacobian(PetscScalar *lxx, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // only the last processor do this
  if(Genius::is_last_processor())
  {
    // ckt load has been done in build_spice_function

    for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
    {
      std::vector<PetscInt> cols;
      std::vector<Real> _v;
      PetscInt global_row;

      _circuit->ckt_matrix_row(n, global_row, cols, _v);

      std::vector<PetscScalar> v;
      for(unsigned int e=0; e<_v.size(); ++e)
        v.push_back(_v[e]);
      jac->add_row(global_row, cols.size(), &cols[0], &v[0]);
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}



void MixASolverBase::dump_spice_matrix_petsc(const std::string &file) const
{
  // only the last processor do this
  if(Genius::is_last_processor())
  {
    Mat A;
    unsigned int N = _circuit->n_ckt_nodes();
    unsigned int n_max_spice_nonzeros = _circuit->max_row_nonzeros();

    MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, n_max_spice_nonzeros, PETSC_NULL, &A);
    for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
    {
      PetscInt row = n;
      std::vector<PetscInt> cols;
      std::vector<Real> _v;

      _circuit->ckt_matrix_row(n, cols, _v);

      std::vector<PetscScalar> v;
      for(unsigned int e=0; e<_v.size(); ++e)
        v.push_back(_v[e]);
      MatSetValues(A, 1, &row, cols.size(), &cols[0], &v[0], ADD_VALUES);
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);

    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_SELF, file.c_str(), FILE_MODE_WRITE, &viewer);
    //PetscViewerSetFormat(viewer,PetscViewerFormat format)
    MatView(A, viewer);
    PetscViewerDestroy(PetscDestroyObject(viewer));

    MatDestroy(PetscDestroyObject(A));
  }
}




void MixASolverBase::print_spice_node() const
{
  if(Genius::is_last_processor())
  {
    std::cout << std::left << std::setw(30) << "  Node" << std::setw(15) << "Voltage" << std::endl;
    std::cout << std::left << std::setw(30) << "  ----" << std::setw(15) << "-------" << std::endl;

    for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
    {
      std::cout << "  " << std::left << std::setw(30) << _circuit->ckt_node_name(n)  << std::setw(15) << _circuit->rhs_old(n) << std::endl;
    }

    std::cout<<std::endl;

    std::cout << std::left << std::setw(15) << "  state0" << std::setw(15) << "  state1" << std::setw(15) << "  state2" << std::endl;
    std::cout << std::left << std::setw(15) << "  ------" << std::setw(15) << "  ------" << std::setw(15) << "  ------" << std::endl;

    std::vector<double> state0;
    _circuit->get_state_vector(0, state0);

    std::vector<double> state1;
    _circuit->get_state_vector(1, state1);

    std::vector<double> state2;
    _circuit->get_state_vector(2, state2);

    for(unsigned int n=0; n<state0.size(); ++n)
    {
      std::cout << "  " << std::left
                << std::setw(15) << state0[n]
                << std::setw(15) << state1[n]
                << std::setw(15) << state2[n]
                << std::endl;
    }
    std::cout<<std::endl;
  }
}


void MixASolverBase::sens_line_search_post_check(Vec x, Vec y, Vec w, PetscBool *changed_y, PetscBool *changed_w)
{
#if 1
  //PetscScalar    *xx;
  //PetscScalar    *yy;
  PetscScalar    *ww;

  //VecGetArray(x, &xx);  // previous iterate value
  //VecGetArray(y, &yy);  // new search direction and length
  VecGetArray(w, &ww);  // current candidate iterate


  if(Genius::is_last_processor())
  {
    // insert x into spice rhs old
    std::vector<double> rhs;
    for(unsigned int n=0; n<_circuit->n_ckt_nodes(); ++n)
      rhs.push_back(ww[_circuit->array_offset_x(n)]);
    _circuit->update_rhs_old(rhs);
  }


  //VecRestoreArray(x, &xx);
  //VecRestoreArray(y, &yy);
  VecRestoreArray(w, &ww);
#endif

  FVM_FlexNonlinearSolver::sens_line_search_post_check(x, y, w, changed_y, changed_w);
}






/* ----------------------------------------------------------------------------
 * compute steadystate
 * all the stimulate source(s) are set with transient time 0 value. time step set to inf
 */
int MixASolverBase::solve_dcop(bool tran_op)
{
  MESSAGE<<"Compute dc operator\n";
  RECORD();

  // not time dependent
  SolverSpecify::TimeDependent = false;
  SolverSpecify::dt = 1e100;

  if(tran_op)
  {
    if(SolverSpecify::UIC)
    {
      // Set Initial Condition
      MESSAGE<<"Use Initial Condition\n";
      RECORD();

      if(Genius::is_last_processor())
      {
        _circuit->set_ckt_mode( MODEUIC | MODETRANOP | MODEINITJCT );
        _circuit->set_time(SolverSpecify::clock/PhysicalUnit::s);
      }
    }
    else
    {
      if(Genius::is_last_processor())
      {
        _circuit->set_ckt_mode( MODETRANOP | MODEINITJCT );
        _circuit->set_time(SolverSpecify::clock/PhysicalUnit::s);
      }
    }

    SolverSpecify::clock = SolverSpecify::TStart;
  }
  else
  {
    if(Genius::is_last_processor())
    {
      _circuit->init_dcop();
      _circuit->set_ckt_mode( MODEDCOP | MODEINITJCT );
    }

    SolverSpecify::clock = 0.0;
  }

  _system.get_field_source()->update(SolverSpecify::clock);



  // set gmin to a big value
  if(Genius::is_last_processor())
  {
    _circuit->ckt_set_gmin(SolverSpecify::GminInit);
  }


  int rampup_steps;
  std::map<std::string, double> vsrc_dc, isrc_dc;

  if(Genius::is_last_processor())
  {
    std::vector<std::string> vsrcs, isrcs;
    _circuit->get_voltage_sources(vsrcs);
    _circuit->get_current_sources(isrcs);

    double vabsmax=0, iabsmax=0;
    for(unsigned int n=0; n<vsrcs.size(); ++n)
    {
      double vdc = _circuit->get_voltage_from(vsrcs[n]);
      vsrc_dc[vsrcs[n]]= vdc;
      if( std::abs(vdc) > vabsmax )
        vabsmax = std::abs(vdc);
    }

    for(unsigned int n=0; n<isrcs.size(); ++n)
    {
      double idc =  _circuit->get_current_from(isrcs[n]);
      isrc_dc[isrcs[n]]= idc;
      if( std::abs(idc) > iabsmax )
        iabsmax = std::abs(idc);
    }

    rampup_steps = std::max( static_cast<int>(vabsmax/(SolverSpecify::RampUpVStep)),
                             static_cast<int>(iabsmax/(SolverSpecify::RampUpIStep)) );

    rampup_steps = std::max(1, rampup_steps);
  }

  if(SolverSpecify::RampUpSteps>0)
    rampup_steps = SolverSpecify::RampUpSteps;

  Parallel::broadcast(rampup_steps, Genius::last_processor_id() );

  MESSAGE<<"DC rampup process..."<<'\n'; RECORD();

  // saved solutions and vscan values for solution projection.
  Vec xs1, xs2, xs3;
  PetscScalar Vs1, Vs2, Vs3;
  VecDuplicate(x,&xs1);
  VecDuplicate(x,&xs2);
  VecDuplicate(x,&xs3);

  SolverSpecify::DC_Cycles=0;
  int retry=0;
  for(int step=1; step<=rampup_steps; ++step)
  {
    // show current rampup value
    MESSAGE << "DC rampup step "  << step << " of " << rampup_steps << '\n'
    <<"--------------------------------------------------------------------------------\n";
    RECORD();

    if(Genius::is_last_processor())
    {
      std::map<std::string, double>::iterator it;
      for(it = vsrc_dc.begin(); it!=vsrc_dc.end(); ++it)
        _circuit->set_voltage_to(it->first, it->second*step/rampup_steps);

      for(it = isrc_dc.begin(); it!=isrc_dc.end(); ++it)
        _circuit->set_current_to(it->first, it->second*step/rampup_steps);
    }


    // call pre_solve_process
    if( SolverSpecify::DC_Cycles == 0 )
      this->pre_solve_process();
    else
      this->pre_solve_process(false);

    // here call Petsc to solve the nonlinear equations
    snes_solve();

    // get the converged reason
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes,&reason);
    // linear solver iteration
    PetscInt lits;
    SNESGetLinearSolveIterations(snes, &lits);

    if( reason>0  ) //ok, converged.
    {
      // call post_solve_process
      this->post_solve_process();

      // save solution for linear/quadratic projection
      if (SolverSpecify::Predict)
      {
        VecCopy(xs2,xs3);
        Vs3=Vs2;
        VecCopy(xs1,xs2);
        Vs2=Vs1;
        VecCopy(x,xs1);
        Vs1=double(step)/rampup_steps;
      }

      // output op
      //print_spice_node();
      retry=0;
      SolverSpecify::DC_Cycles++;
      MESSAGE
      <<"--------------------------------------------------------------------------------\n"
      <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits << "\n\n\n";
      RECORD();
    }
    else // oh, diverged... reduce step and try again
    {

      if(retry++>=3)
      {
        MESSAGE <<"------> Too many failed steps, give up tring.\n\n\n";
        RECORD();
        break;
      }

      // load previous result into solution vector
      this->diverged_recovery();

      step *= 2;
      rampup_steps *= 2;
      step -= 2;

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
    }

    if (SolverSpecify::Predict)
    {
      PetscScalar hn = double(step)/rampup_steps-Vs1;
      PetscScalar hn1 = Vs1-Vs2;
      PetscScalar hn2 = Vs2-Vs3;
      if(SolverSpecify::DC_Cycles>=3)
      {
        // quadradic projection
        PetscScalar cn=hn*(hn+2*hn1+hn2)/(hn1*(hn1+hn2));
        PetscScalar cn1=-hn*(hn+hn1+hn2)/(hn1*hn2);
        PetscScalar cn2=hn*(hn+hn1)/(hn2*(hn1+hn2));

        VecAXPY(x,cn,xs1);
        VecAXPY(x,cn1,xs2);
        VecAXPY(x,cn2,xs3);
        this->projection_positive_density_check(x,xs1);
      }
      else if(SolverSpecify::DC_Cycles>=2)
      {
        // linear projection
        VecAXPY(x, hn/hn1,xs1);
        VecAXPY(x,-hn/hn1,xs2);
        this->projection_positive_density_check(x,xs1);
      }
    }
  }

  VecDestroy(PetscDestroyObject(xs1));
  VecDestroy(PetscDestroyObject(xs2));
  VecDestroy(PetscDestroyObject(xs3));


  // scale gmin back to user given value

  double gmin = SolverSpecify::GminInit;
  double gmin_user = SolverSpecify::Gmin;

  while( gmin > gmin_user )
  {
    gmin = std::max(gmin*1e-2, gmin_user);

    MESSAGE << "DC reduce gmin to "  << gmin << '\n'
            <<"--------------------------------------------------------------------------------\n";
    RECORD();

    if(Genius::is_last_processor())
    {
      _circuit->ckt_set_gmin(gmin);
    }

    // call pre_solve_process
    if( SolverSpecify::DC_Cycles == 0 )
      this->pre_solve_process();
    else
      this->pre_solve_process(false);

    // here call Petsc to solve the nonlinear equations
    snes_solve();

    // get the converged reason
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes,&reason);
    // linear solver iteration
    PetscInt lits;
    SNESGetLinearSolveIterations(snes, &lits);
    if( reason < 0  )
    {
      if(reason == SNES_DIVERGED_LINEAR_SOLVE)
      {
        KSPConvergedReason ksp_reason;
        KSPGetConvergedReason ( ksp, &ksp_reason );
        MESSAGE <<"------> linear solver "<<KSPConvergedReasons[ksp_reason]<<", stop reduce gmin\n\n\n"; RECORD();
      }
      else
      {
        MESSAGE <<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", stop reduce gmin\n\n\n"; RECORD();
      }

      break;
    }

    // output op
    //print_spice_node();

    MESSAGE
    <<"--------------------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[reason] << ", total linear iteration " << lits << "\n\n\n";
    RECORD();

    // call post_solve_process
    this->post_solve_process();
  }

  print_spice_node();


  SolverSpecify::tran_histroy = false;


  return 0;
}


int MixASolverBase::solve_dcsweep()
{
  _system.get_field_source()->update(0);

  // not time dependent
  SolverSpecify::TimeDependent = false;
  SolverSpecify::dt = 1e100;
  SolverSpecify::clock = 0.0;

  // set spice circuit
  if(Genius::is_last_processor())
  {
    _circuit->init_dctrcurv();
    _circuit->set_ckt_mode( MODEDCTRANCURVE | MODEINITJCT );
  }

  // output DC Scan information
  if(SolverSpecify::Electrode_VScan.size())
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
  if(SolverSpecify::Electrode_VScan.size())
  {
    std::string vsrc = SolverSpecify::Electrode_VScan[0];

    if(Genius::is_last_processor())
      _circuit->get_ckt_voltage_source(vsrc);

    // the current vscan voltage
    PetscScalar Vscan = SolverSpecify::VStart;
    // the current vscan step
    PetscScalar VStep = SolverSpecify::VStep;

    // saved solutions and vscan values for solution projection.
    Vec xs1, xs2, xs3;
    PetscScalar Vs1, Vs2, Vs3;
    std::stack<PetscScalar> V_retry;
    VecDuplicate(x,&xs1);
    VecDuplicate(x,&xs2);
    VecDuplicate(x,&xs3);

    // main loop
    for(SolverSpecify::DC_Cycles=0;  Vscan*SolverSpecify::VStep < SolverSpecify::VStop*SolverSpecify::VStep*(1.0+1e-7);)
    {
      // show current vscan value
      MESSAGE << "DC Scan: "  << vsrc << " = " << Vscan/PhysicalUnit::V  <<" V" << '\n'
      <<"--------------------------------------------------------------------------------\n";
      RECORD();

      // set current vscan voltage to corresponding electrode
      if(Genius::is_last_processor())
        _circuit->set_voltage_to(vsrc, Vscan);

      SolverSpecify::Electrode_VScan_Voltage = Vscan;

      // call pre_solve_process
      if( SolverSpecify::DC_Cycles == 0 )
        this->pre_solve_process();
      else
        this->pre_solve_process(false);

      // here call Petsc to solve the nonlinear equations
      snes_solve();

      // get the converged reason
      SNESConvergedReason reason;
      SNESGetConvergedReason(snes,&reason);
       // linear solver iteration
      PetscInt lits;
      SNESGetLinearSolveIterations(snes, &lits);
      if( reason>0  ) //ok, converged.
      {
        // call post_solve_process
        this->post_solve_process();

        SolverSpecify::DC_Cycles++;

        // save solution for linear/quadratic projection
        if (SolverSpecify::Predict)
        {
          VecCopy(xs2,xs3);
          Vs3=Vs2;
          VecCopy(xs1,xs2);
          Vs2=Vs1;
          VecCopy(x,xs1);
          Vs1=Vscan;
        }

        if(V_retry.empty())
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

        if(fabs(Vscan-SolverSpecify::VStop)<1e-10)
          Vscan=SolverSpecify::VStop;

        // if v step small than VStepMax, mult by factor of 1.1
        if( fabs(VStep) < fabs(SolverSpecify::VStepMax) )
          VStep *= 1.1;


        // however, for last step, we force V equal to VStop
        if( Vscan*SolverSpecify::VStep > SolverSpecify::VStop*SolverSpecify::VStep &&
            Vscan*SolverSpecify::VStep < (SolverSpecify::VStop + VStep - 1e-10*VStep)*SolverSpecify::VStep
          )
          Vscan = SolverSpecify::VStop;

        MESSAGE
        <<"--------------------------------------------------------------------------------\n"
        <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits <<"\n\n\n";
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
        V_retry.push(Vscan);
        Vscan=(Vscan+Vs1)/2.0;
      }

      if (SolverSpecify::Predict)
      {
        PetscScalar hn = Vscan-Vs1;
        PetscScalar hn1 = Vs1-Vs2;
        PetscScalar hn2 = Vs2-Vs3;
        if(SolverSpecify::DC_Cycles>=3)
        {
          // quadradic projection
          PetscScalar cn=hn*(hn+2*hn1+hn2)/(hn1*(hn1+hn2));
          PetscScalar cn1=-hn*(hn+hn1+hn2)/(hn1*hn2);
          PetscScalar cn2=hn*(hn+hn1)/(hn2*(hn1+hn2));

          VecAXPY(x,cn,xs1);
          VecAXPY(x,cn1,xs2);
          VecAXPY(x,cn2,xs3);
          this->projection_positive_density_check(x,xs1);
        }
        else if(SolverSpecify::DC_Cycles>=2)
        {
          // linear projection
          VecAXPY(x, hn/hn1,xs1);
          VecAXPY(x,-hn/hn1,xs2);
          this->projection_positive_density_check(x,xs1);
        }
      }
    }

    VecDestroy(PetscDestroyObject(xs1));
    VecDestroy(PetscDestroyObject(xs2));
    VecDestroy(PetscDestroyObject(xs3));

  }



  // current scan
  if(SolverSpecify::Electrode_IScan.size())
  {
    std::string isrc = SolverSpecify::Electrode_IScan[0];

    if(Genius::is_last_processor())
      _circuit->get_ckt_current_source(isrc);

    // iscan current
    PetscScalar Iscan = SolverSpecify::IStart;

    // iscan step
    PetscScalar IStep = SolverSpecify::IStep;

    // saved solutions and iscan values for solution projection.
    Vec xs1, xs2, xs3;
    PetscScalar Is1, Is2, Is3;
    std::stack<PetscScalar> I_retry;
    VecDuplicate(x,&xs1);
    VecDuplicate(x,&xs2);
    VecDuplicate(x,&xs3);

    // main loop
    for(SolverSpecify::DC_Cycles=0;  Iscan*SolverSpecify::IStep < SolverSpecify::IStop*SolverSpecify::IStep*(1.0+1e-7);)
    {
      // show current iscan value
      MESSAGE << "DC Scan: "  << isrc << " = " << Iscan/PhysicalUnit::A  <<" A" << '\n'
      <<"--------------------------------------------------------------------------------\n";

      // set current to corresponding electrode
      if(Genius::is_last_processor())
        _circuit->set_current_to(isrc, Iscan/PhysicalUnit::A);

      SolverSpecify::Electrode_IScan_Current = Iscan;

      // call pre_solve_process
      if( SolverSpecify::DC_Cycles == 0 )
        this->pre_solve_process();
      else
        this->pre_solve_process(false);

      // here call Petsc to solve the nonlinear equations
      snes_solve();

      // get the converged reason
      SNESConvergedReason reason;
      SNESGetConvergedReason(snes,&reason);
      // linear solver iteration
      PetscInt lits;
      SNESGetLinearSolveIterations(snes, &lits);

      if( reason>0  ) //ok, converged.
      {
        // call post_solve_process
        this->post_solve_process();

        SolverSpecify::DC_Cycles++;

        // save solution for linear/quadratic projection
        if (SolverSpecify::Predict)
        {
          VecCopy(xs2,xs3);
          Is3=Is2;
          VecCopy(xs1,xs2);
          Is2=Is1;
          VecCopy(x,xs1);
          Is1=Iscan;
        }

        if(I_retry.empty())
        {
          // add istep to current
          Iscan += IStep;
        }
        else
        {
          // pop
          Iscan = I_retry.top();
          I_retry.pop();
        }

        if(fabs(Iscan-SolverSpecify::IStop)<1e-10)
          Iscan=SolverSpecify::IStop;

        // if I step small than IStepMax, mult by factor of 1.1
        if( fabs(IStep) < fabs(SolverSpecify::IStepMax) )
          IStep *= 1.1;


        // however, for last step, we force I equal to IStop
        if( Iscan*SolverSpecify::IStep > SolverSpecify::IStop*SolverSpecify::IStep &&
            Iscan*SolverSpecify::IStep < (SolverSpecify::IStop + IStep - 1e-10*IStep)*SolverSpecify::IStep
          )
          Iscan = SolverSpecify::IStop;

        MESSAGE
        <<"--------------------------------------------------------------------------------\n"
        <<"      "<<SNESConvergedReasons[reason]<<", total linear iteration " << lits <<"\n\n\n";
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
        I_retry.push(Iscan);
        Iscan=(Iscan+Is1)/2.0;
      }

      if (SolverSpecify::Predict)
      {
        PetscScalar hn = Iscan-Is1;
        PetscScalar hn1 = Is1-Is2;
        PetscScalar hn2 = Is2-Is3;
        if(SolverSpecify::DC_Cycles>=3)
        {
          // quadradic projection
          PetscScalar cn=hn*(hn+2*hn1+hn2)/(hn1*(hn1+hn2));
          PetscScalar cn1=-hn*(hn+hn1+hn2)/(hn1*hn2);
          PetscScalar cn2=hn*(hn+hn1)/(hn2*(hn1+hn2));

          VecAXPY(x,cn,xs1);
          VecAXPY(x,cn1,xs2);
          VecAXPY(x,cn2,xs3);
          this->projection_positive_density_check(x,xs1);
        }
        else if(SolverSpecify::DC_Cycles>=2)
        {
          // linear projection
          VecAXPY(x, hn/hn1,xs1);
          VecAXPY(x,-hn/hn1,xs2);
          this->projection_positive_density_check(x,xs1);
        }
      }
    }

    VecDestroy(PetscDestroyObject(xs1));
    VecDestroy(PetscDestroyObject(xs2));
    VecDestroy(PetscDestroyObject(xs3));
  }


  SolverSpecify::tran_histroy = false;

  return 0;
}



/*----------------------------------------------------------------------------
 * transient simulation!
 */
int MixASolverBase::solve_transient()
{
    // diverged counter
  int diverged_retry=0;

  // auto time step counter
  int autostep_retry=0;

  // init aux vectors used in transient simulation
  VecDuplicate(x, &x_n);
  VecDuplicate(x, &x_n1);
  VecDuplicate(x, &x_n2);
  VecDuplicate(x, &xp);
  VecDuplicate(x, &LTE);

  // set spice circuit, does uic required?
  if(Genius::is_last_processor())
  {
    _circuit->init_dctran( SolverSpecify::UIC ? MODEUIC : 0);
  }


  // if we need do op
  if(SolverSpecify::tran_op)
  {
    solve_dcop(true);
  }


  // output op
  //print_spice_node();

  // time dependent
  SolverSpecify::TimeDependent = true;

  // if BDF2 scheme is used, we should set SolverSpecify::BDF2_LowerOrder flag to true
  if(SolverSpecify::TS_type==SolverSpecify::BDF2)
    SolverSpecify::BDF2_LowerOrder = true;

  // for the first step, dt equals to TStep
  SolverSpecify::dt = SolverSpecify::TStep;

  // transient simulation clock
  SolverSpecify::clock = SolverSpecify::TStart + SolverSpecify::TStep;

  MESSAGE<<"Transient compute from "<<SolverSpecify::TStart <<" ps"
         <<" to "  <<SolverSpecify::TStop<<" ps" <<'\n';
  RECORD();

  // time step counter
  SolverSpecify::T_Cycles=0;

  double dt_dynamic_factor = 1.0;

  // set spice circuit
  if(Genius::is_last_processor())
  {
    _circuit->set_modeinittran(SolverSpecify::dt/PhysicalUnit::s);
    // by default, spice use TRAPEZOIDAL for integration
    // however, it is not stable for sharp changes of circuit state (has oscillation).
    // (I don't know how spice keep this method works, maybe I use different time step strategy)
    // As a result, here I set GEAR (BDF) method for integration in spice, which is A and L stable.
    // And GEAR (BDF) is also the integration method used by Genius.
    _circuit->set_integrate_method(GEAR);
  }


  // the main loop of transient solver.
  do
  {
    MESSAGE
    <<"t = "<<SolverSpecify::clock/s*1e12<<" ps, "<< "dt = " << SolverSpecify::dt/s*1e12 <<" ps"<< '\n'
    <<"--------------------------------------------------------------------------------\n";
    RECORD();

    //print_spice_node();

    //update sources to current clock
    if(Genius::is_last_processor())
    {
      _circuit->set_time(SolverSpecify::clock/PhysicalUnit::s);
      _circuit->set_delta(SolverSpecify::dt/PhysicalUnit::s);
    }

    _system.get_field_source()->update(SolverSpecify::clock);

    //we do solve here!

    // call pre_solve_process
    if( SolverSpecify::T_Cycles == 0 )
      this->pre_solve_process();
    else
      this->pre_solve_process(false);

    // here call Petsc to solve the nonlinear equations
    snes_solve();

    //print_spice_node();

    // get the converged reason
    SNESConvergedReason reason;
    SNESGetConvergedReason(snes,&reason);
    // linear solver iteration
    PetscInt lits;
    SNESGetLinearSolveIterations(snes, &lits);
    //nonlinear solution diverged? try to do recovery
    if(reason<0)
    {
      // increase diverged_retry
      diverged_retry++;

      if(diverged_retry >= 8) //failed 8 times, stop tring
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

      if(SolverSpecify::clock < SolverSpecify::TStart)
        SolverSpecify::clock = SolverSpecify::TStart;

      // load previous result into solution vector
      this->diverged_recovery();
      continue;
    }

    //ok, nonlinear solution converged.

    dt_dynamic_factor = 1.0;

    //do LTE estimation and auto time step control
    if ( SolverSpecify::AutoStep &&
         ((SolverSpecify::TS_type==SolverSpecify::BDF1 && SolverSpecify::T_Cycles>=3) ||
          (SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::T_Cycles>=4) )  )
    {
      PetscScalar r = this->LTE_norm();

      if(SolverSpecify::TS_type==SolverSpecify::BDF1)
        r = std::pow(r, PetscScalar(-1.0/2));
      else if(SolverSpecify::TS_type==SolverSpecify::BDF2)
      {
        if(SolverSpecify::BDF2_LowerOrder)
          r = std::pow ( r, PetscScalar ( -1.0/2 ) );
        else
          r = std::pow ( r, PetscScalar ( -1.0/3 ) );
      }


      // when r<0.9, reject this solution
      if ( SolverSpecify::RejectStep && r<0.9 && SolverSpecify::dt > SolverSpecify::TStepMin )
      {
        // clear the counter
        diverged_retry = 0;
        autostep_retry++;

        MESSAGE <<"------> LTE too large, time step rejected...\n\n\n";
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

        if(SolverSpecify::T_Cycles==0 && Genius::is_last_processor())
          _circuit->set_ckt_mode( MODETRAN | MODEINITTRAN );

        continue;
      }
      else      // else, accept this solution
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

    if(Genius::is_last_processor() && SolverSpecify::T_Cycles==0)
      _circuit->prepare_ckt_state_first_time();

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
    if ( SolverSpecify::TStepMax >0 && SolverSpecify::dt > SolverSpecify::TStepMax )
      SolverSpecify::dt = SolverSpecify::TStepMax;

    // limit time step by field source
    SolverSpecify::dt = _system.get_field_source()->limit_dt(SolverSpecify::clock, SolverSpecify::dt, SolverSpecify::TStepMin);

    // set clock to next time step
    SolverSpecify::clock += SolverSpecify::dt;

    //make sure we can terminat near TStop, the relative error should less than 1e-10.
    if ( SolverSpecify::clock > SolverSpecify::TStop
      && SolverSpecify::clock < (SolverSpecify::TStop + SolverSpecify::dt - 1e-10*SolverSpecify::dt) )
    {
      SolverSpecify::dt -= SolverSpecify::clock - SolverSpecify::TStop;
      SolverSpecify::clock = SolverSpecify::TStop;
    }

    //check if BDF2 can be used?
    if(SolverSpecify::TS_type==SolverSpecify::BDF2)
      SolverSpecify::BDF2_LowerOrder = this->BDF2_positive_defined();

    // use by auto step control and predict
    if( SolverSpecify::AutoStep  || SolverSpecify::Predict )
    {
      VecCopy ( x_n1, x_n2 );
      VecCopy ( x_n, x_n1 );
      VecCopy ( x, x_n );
    }

    if(Genius::is_last_processor())
    {
      _circuit->rotate_state_vectors();
      _circuit->set_ckt_mode( MODETRAN | MODEINITPRED);
      //if(SolverSpecify::TS_type==SolverSpecify::BDF2 && !SolverSpecify::BDF2_LowerOrder)
        _circuit->set_time_order(2);
      //_circuit->set_integrate_method(1);
    }

  Predict:

    // predict next solution
    if ( SolverSpecify::Predict )
    {
      PetscScalar hn  = SolverSpecify::dt;           // here dt is the next time step
      PetscScalar hn1 = SolverSpecify::dt_last;      // time step n-1
      PetscScalar hn2 = SolverSpecify::dt_last_last; // time step n-2

      if ( SolverSpecify::TS_type == SolverSpecify::BDF1 && SolverSpecify::T_Cycles>=3 )
      {
        VecZeroEntries ( x );
        // use linear interpolation to predict solution x
        VecAXPY ( x, 1+hn/hn1, x_n );
        VecAXPY ( x, -hn/hn1,  x_n1 );
        this->projection_positive_density_check ( x, x_n );
      }
      else if ( SolverSpecify::TS_type == SolverSpecify::BDF2 && SolverSpecify::T_Cycles>=4)
      {
        VecZeroEntries ( x );
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
  while(SolverSpecify::clock < SolverSpecify::TStop+0.5*SolverSpecify::dt);


  // free aux vectors
  VecDestroy(PetscDestroyObject(x_n));
  VecDestroy(PetscDestroyObject(x_n1));
  VecDestroy(PetscDestroyObject(x_n2));
  VecDestroy(PetscDestroyObject(xp));
  VecDestroy(PetscDestroyObject(LTE));


  SolverSpecify::tran_histroy = true;

  return 0;
}


/*------------------------------------------------------------------
 * snes convergence criteria
 */
#if PETSC_VERSION_LT(3, 6, 0)
  #if PETSC_VERSION_LE(3, 2, 0)
    #include "private/snesimpl.h"
    #define SNES_CONVERGED_SNORM_RELATIVE  SNES_CONVERGED_PNORM_RELATIVE
  #else
    #include "petsc-private/snesimpl.h"
  #endif
#endif

#if PETSC_VERSION_GE(3, 6, 0)
  #include <petsc/private/snesimpl.h>
#endif
void MixASolverBase::petsc_snes_convergence_test(PetscInt its, PetscReal , PetscReal pnorm, PetscReal fnorm, SNESConvergedReason *reason)
{
  // update error norm
  this->error_norm();

  *reason = SNES_CONVERGED_ITERATING;

  // the first iteration
  if (!its)
  {
    snes->ttol = fnorm*snes->rtol;

    MESSAGE<<" "<<" n ";
    MESSAGE<<"| Eq(V) | "<<"| Eq(n) | "<<"| Eq(p) | ";
    MESSAGE<<"| Eq(T) | ";
    MESSAGE<<"|Eq(Tn)|  ";
    MESSAGE<<"|Eq(Tp)|  ";
    MESSAGE<<"| SPICE | ";
    MESSAGE<<"Lg(dx)"<<'\n';
    MESSAGE<<"--------------------------------------------------------------------------------\n";
    RECORD();
  }

  unsigned int dim = this->system().dim();
  double z_width = (dim == 2 ? 1.0*um : 1.0);

  double  toler_relax = SolverSpecify::toler_relax;
  if(its)
  {
    //toler_relax = std::max(SolverSpecify::toler_relax, 1/(pnorm+1e-6));
  }
  bool  poisson_conv              = poisson_norm*z_width              < toler_relax*SolverSpecify::poisson_abs_toler;
  bool  elec_continuity_conv      = elec_continuity_norm*z_width      < toler_relax*SolverSpecify::elec_continuity_abs_toler;
  bool  hole_continuity_conv      = hole_continuity_norm*z_width      < toler_relax*SolverSpecify::hole_continuity_abs_toler;
  bool  spice_conv                = spice_norm                        < toler_relax*SolverSpecify::spice_abs_toler;
  bool  heat_equation_conv        = heat_equation_norm*z_width        < toler_relax*SolverSpecify::heat_equation_abs_toler;
  bool  elec_energy_equation_conv = elec_energy_equation_norm*z_width < toler_relax*SolverSpecify::elec_energy_abs_toler;
  bool  hole_energy_equation_conv = hole_energy_equation_norm*z_width < toler_relax*SolverSpecify::hole_energy_abs_toler;


  bool div = pnorm > 1e5 ||
             poisson_norm > SolverSpecify::divergence_factor*SolverSpecify::poisson_abs_toler ||
             elec_continuity_norm > SolverSpecify::divergence_factor*SolverSpecify::elec_continuity_abs_toler ||
             hole_continuity_norm > SolverSpecify::divergence_factor*SolverSpecify::hole_continuity_abs_toler ||
             spice_norm       > SolverSpecify::divergence_factor*SolverSpecify::spice_abs_toler;


#ifdef WINDOWS
  MESSAGE.precision(1);
#else

  MESSAGE.precision(2);
#endif

  MESSAGE<< std::setw(3) << its << " " ;
  MESSAGE<< std::scientific;
  if(dim == 2)
  {
    MESSAGE<< poisson_norm*z_width/C << (poisson_conv ? "* " : "  ");
    MESSAGE<< elec_continuity_norm*z_width/A << (elec_continuity_conv ? "* " : "  ");
    MESSAGE<< hole_continuity_norm*z_width/A << (hole_continuity_conv ? "* " : "  ");
    MESSAGE<< heat_equation_norm*z_width/W   << (heat_equation_conv ? "* " : "  ");
    MESSAGE<< elec_energy_equation_norm*z_width/W << (elec_energy_equation_conv ? "* " : "  ");
    MESSAGE<< hole_energy_equation_norm*z_width/W << (hole_energy_equation_conv ? "* " : "  ");
  }

  if(dim == 3)
  {
    MESSAGE<< poisson_norm/C << (poisson_conv ? "* " : "  ");
    MESSAGE<< elec_continuity_norm/A << (elec_continuity_conv ? "* " : "  ");
    MESSAGE<< hole_continuity_norm/A << (hole_continuity_conv ? "* " : "  ");
    MESSAGE<< heat_equation_norm/W   << (heat_equation_conv ? "* " : "  ");
    MESSAGE<< elec_energy_equation_norm/W << (elec_energy_equation_conv ? "* " : "  ");
    MESSAGE<< hole_energy_equation_norm/W << (hole_energy_equation_conv ? "* " : "  ");
  }

  MESSAGE<< spice_norm/A << (spice_conv ? "* " : "  ");
  MESSAGE<< std::fixed << std::setw(4) << (pnorm==0.0 ? -std::numeric_limits<PetscScalar>::infinity():log10(pnorm))
    << (pnorm < SolverSpecify::relative_toler ? "*" : " ") << "\n" ;
  RECORD();
  MESSAGE.precision ( 6 );
  MESSAGE<< std::scientific;

  // check for NaN (Not a Number)
  if (fnorm != fnorm)
  {
    *reason = SNES_DIVERGED_FNORM_NAN;
  }
  else if (snes->nfuncs >= snes->max_funcs)
  {
    *reason = SNES_DIVERGED_FUNCTION_COUNT;
  }
  else if (its > 4 && div)
  {
    *reason = SNES_DIVERGED_FUNCTION_DOMAIN;
  }


  if ( *reason == SNES_CONVERGED_ITERATING )
  {
    if(its == 0 && fnorm<snes->abstol )
    {
      *reason = SNES_CONVERGED_FNORM_ABS;
    }
    else
    {
      // check for absolute convergence
      if ( its &&
           poisson_norm               < SolverSpecify::poisson_abs_toler         &&
           elec_continuity_norm       < SolverSpecify::elec_continuity_abs_toler &&
           hole_continuity_norm       < SolverSpecify::hole_continuity_abs_toler &&
           spice_norm                 < SolverSpecify::spice_abs_toler       &&
           heat_equation_norm         < SolverSpecify::heat_equation_abs_toler   &&
           elec_energy_equation_norm  < SolverSpecify::elec_energy_abs_toler     &&
           hole_energy_equation_norm  < SolverSpecify::hole_energy_abs_toler )
      {
        *reason = SNES_CONVERGED_FNORM_ABS;
      }
      else if ( std::abs(fnorm-function_norm)/(fnorm+function_norm) <=  SolverSpecify::snes_rtol  &&
                poisson_conv         &&
                elec_continuity_conv &&
                hole_continuity_conv &&
                spice_conv           &&
                heat_equation_conv   &&
                elec_energy_equation_conv &&
                hole_energy_equation_conv   )
      {
        *reason = SNES_CONVERGED_FNORM_RELATIVE;
      }
      // check for relative convergence, should have at least one iteration here
      else if ( its && pnorm < SolverSpecify::relative_toler &&
                poisson_conv         &&
                elec_continuity_conv &&
                hole_continuity_conv &&
                spice_conv           &&
                heat_equation_conv   &&
                elec_energy_equation_conv &&
                hole_energy_equation_conv   )
      {
        *reason = SNES_CONVERGED_SNORM_RELATIVE;
      }
      else if ( fnorm <=  SolverSpecify::spice_abs_toler  &&
                poisson_conv         &&
                elec_continuity_conv &&
                hole_continuity_conv &&
                spice_conv           &&
                heat_equation_conv   &&
                elec_energy_equation_conv &&
                hole_energy_equation_conv   )
      {
        *reason = SNES_CONVERGED_FNORM_ABS;
      }
    }

  }


  // update CTKMode
  // when the bit MODEINITJCT in CKTmode is set, nodeset/ic will be considered.
  // spice first stick the node solution to nodeset/ic value, until convergence is reached.
  // after that, we must cancle the MODEINITJCT flag to find the real solution without
  // the effect of nodeset/ic.
  //
  int ckt_mode;
  if(Genius::is_last_processor())
  {
    ckt_mode = _circuit->change_ckt_mode(int(*reason));
  }
  Parallel::broadcast(ckt_mode, Genius::last_processor_id());

  // if ckt_mode, continue newton iteration
  if( ckt_mode )
  {
    *reason = SNES_CONVERGED_ITERATING;
  }

  //in OP mode, should do newton iteration at least once
  if (!its && SolverSpecify::Type==SolverSpecify::OP)
  {
    *reason = SNES_CONVERGED_ITERATING;
  }


  // record function norm of this iteration
  function_norm = fnorm;
  // record iteration
  nonlinear_iteration = its;


  return;
}


/*------------------------------------------------------------------
 * ksp convergence criteria
 */
void MixASolverBase::petsc_ksp_convergence_test ( PetscInt its, PetscReal rnorm, KSPConvergedReason* reason )
{
  PetscInt kspit = std::max(200, std::min(1000, static_cast<int>(n_global_dofs/10)));

  PetscScalar rtol = SolverSpecify::ksp_rtol;
  PetscScalar abstol = std::max ( std::min(1e-3,SolverSpecify::ksp_atol_fnorm*function_norm), SolverSpecify::ksp_atol);
  if(its > static_cast<PetscInt>(0.3*kspit))
    abstol *= 1e1;
  if(its > static_cast<PetscInt>(0.5*kspit))
    abstol *= 1e2;
  if(its > static_cast<PetscInt>(0.7*kspit))
    abstol *= 1e3;
  if(its > static_cast<PetscInt>(0.9*kspit))
    abstol *= 1e4;

  KSPSetTolerances(ksp, rtol, abstol, 1e10, kspit );
#if PETSC_VERSION_GE(3,5,0)
  KSPConvergedDefault(ksp, its, rnorm, reason, this);
#else
  KSPDefaultConverged(ksp, its, rnorm, reason, this);
#endif
  // stupid code, but KSPSetTolerances does NOT affect KSPDefaultConverged here!
  // seems KSPDefaultConverged cache the value...
  if(*reason == 0)
  {
    if(rnorm < abstol) *reason = KSP_CONVERGED_ATOL;
  }

}


