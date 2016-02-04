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



#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_is.h"
#include "petsc_utils.h"
#include "parallel.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::eV;
using PhysicalUnit::A;
using PhysicalUnit::V;
using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::nm;
using PhysicalUnit::C;
using PhysicalUnit::hbar;
using PhysicalUnit::me;

/*---------------------------------------------------------------------
 * set scaling constant
 */
void InsulatorSemiconductorInterfaceBC::DG_Fill_Value(Vec , Vec L)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      if( region->type() != InsulatorRegion ) continue;

      const FVM_Node * fvm_node = (*rnode_it).second.second;
      VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
    }
  }
}

///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for Density Gradient solver
 */
void InsulatorSemiconductorInterfaceBC::DG_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
          // Insulator-Semiconductor interface at Semiconductor side, do nothing
          case SemiconductorRegion:  break;

          // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
          case InsulatorRegion:
          {
            // record the source row and dst row
            src_row.push_back(fvm_nodes[i]->global_offset());
            dst_row.push_back(fvm_nodes[0]->global_offset());
            clear_row.push_back(fvm_nodes[i]->global_offset());
            break;
          }
          case VacuumRegion:
          break;

          default: genius_error(); //we should never reach here
      }
    }

  }

}


/*---------------------------------------------------------------------
 * build function and its jacobian for Density Gradient solver
 */
void InsulatorSemiconductorInterfaceBC::DG_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Insulator-Semiconductor interface is processed here

  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  const PetscScalar T   = T_external();

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const InsulatorSimulationRegion * insulator_region = dynamic_cast<const InsulatorSimulationRegion *> ( _r2 );
  unsigned int qn_offset = semiconductor_region->dg_variable_offset(EQC);
  unsigned int qp_offset = semiconductor_region->dg_variable_offset(EQV);
  const PetscScalar QNFactor =  semiconductor_region->get_advanced_model()->QNFactor;
  const PetscScalar QPFactor =  semiconductor_region->get_advanced_model()->QPFactor;


  // buffer for Vec location
  std::vector<PetscInt> iy;
  iy.reserve(n_nodes());
  // buffer for Vec value
  std::vector<PetscScalar> y;
  y.reserve(n_nodes());

  const PetscScalar qf = this->scalar("qf");

  const PetscScalar  melec = insulator_region->material()->band->EffecElecMass(T);
  const PetscScalar  mhole = insulator_region->material()->band->EffecHoleMass(T);
  const PetscScalar  Affinity_ins  = insulator_region->material()->basic->Affinity(T);
  const PetscScalar  Affinity_semi = semiconductor_region->material()->basic->Affinity(T);
  const PetscScalar  Eg_ins  = insulator_region->material()->band->Eg(T);
  const PetscScalar  Eg_semi = semiconductor_region->material()->band->Eg(T);
  const PetscScalar  Phi_elec = std::abs(Affinity_semi - Affinity_ins);// elec barrier hight
  const PetscScalar  Phi_hole = std::abs(Affinity_ins+Eg_ins - Affinity_semi - Eg_semi);// hole barrier hight

  // characteristic penetration depth obtained from the Wentzel-Kramers-Brillouin (WKB) approximation
  const PetscScalar x_np = hbar/sqrt(2*melec*Phi_elec);
  const PetscScalar x_pp = hbar/sqrt(2*mhole*Phi_hole);

  const PetscScalar gn = insulator_region->material()->band->Gamman();
  const PetscScalar gp = insulator_region->material()->band->Gammap();
  const PetscScalar b_nox = gn*hbar*hbar/(6*e*melec);
  const PetscScalar b_pox = gp*hbar*hbar/(6*e*mhole);

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;


    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

    const FVM_Node * insulator_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * insulator_node_data = insulator_node->node_data();

    // Insulator-Semiconductor interface at Semiconductor side, do nothing

    PetscScalar V_semiconductor = x[semiconductor_node->local_offset()];
    PetscScalar n   =  x[semiconductor_node->local_offset()+1];                         // electron density
    PetscScalar p   =  x[semiconductor_node->local_offset()+2];                         // hole density

    const PetscScalar boundary_area = std::abs(semiconductor_node->outside_boundary_surface_area());

    //the WBK approx of Si/SiO2 interface

    if(qn_offset != invalid_uint)
    {
      PetscScalar grad_qpn = -QNFactor*b_nox/x_np*boundary_area;
      VecSetValue(f, semiconductor_node->global_offset()+qn_offset, grad_qpn, ADD_VALUES);
    }

    if(qp_offset != invalid_uint)
    {
      PetscScalar grad_qpp = QPFactor*b_pox/x_pp*boundary_area;
      VecSetValue(f, semiconductor_node->global_offset()+qp_offset, grad_qpp, ADD_VALUES);
    }


    // process interface fixed charge density
    VecSetValue(f, semiconductor_node->global_offset(), (qf+semiconductor_node_data->interface_charge())*boundary_area, ADD_VALUES);

    {
      // surface recombination
      semiconductor_region->material()->mapping(semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock);
      PetscScalar GSurf = - semiconductor_region->material()->band->R_Surf(p, n, T) * boundary_area;

      VecSetValue(f, semiconductor_node->global_offset()+1, GSurf, ADD_VALUES);
      VecSetValue(f, semiconductor_node->global_offset()+2, GSurf, ADD_VALUES);
    }

    if (semiconductor_region->get_advanced_model()->Trap)
    {
      // calculate interface trap occupancy
      PetscScalar ni = semiconductor_region->material()->band->nie(p, n, T);
      semiconductor_region->material()->trap->Calculate(false,p,n,ni,T);

      // contribution of trapped charge to Poisson's eqn
      PetscScalar TrappedC = semiconductor_region->material()->trap->Charge(false) * boundary_area;
      if (TrappedC !=0)
        VecSetValue(f, semiconductor_node->global_offset(), TrappedC, ADD_VALUES);

      // electron/hole capture rates and contribution to continuity equations
      PetscScalar TrapElec = semiconductor_region->material()->trap->ElectronTrapRate(false,n,ni,T) * boundary_area;
      PetscScalar TrapHole = semiconductor_region->material()->trap->HoleTrapRate    (false,p,ni,T) * boundary_area;
      if (TrapElec != 0)
        VecSetValue(f, semiconductor_node->global_offset()+1, -TrapElec, ADD_VALUES);
      if (TrapHole != 0)
        VecSetValue(f, semiconductor_node->global_offset()+2, -TrapHole, ADD_VALUES);
    }


    // Insulator-Semiconductor interface at Insulator side, we get the previous f for insulator region,
    // as well as surface charge density, plus to corresponding function of semiconductor. and then
    // force the potential equal to corresponding point in semiconductor


    iy.push_back(insulator_node->global_offset());

    // the governing equation of this fvm node

    // psi of this node
    PetscScalar V_insulator = x[insulator_node->local_offset()];


    // the phi of this node is equal to corresponding psi of semiconductor node
    // since psi should be continuous for the interface
    PetscScalar f_phi = V_insulator - V_semiconductor;

    y.push_back(f_phi);
  }


  // insert new value to src row
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y[0]), ADD_VALUES);

  add_value_flag = ADD_VALUES;

  // gate current
  _gate_current_function(x, f, add_value_flag);
}





/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for Density Gradient solver
 */
void InsulatorSemiconductorInterfaceBC::DG_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * insulator_node = get_region_fvm_node ( ( *node_it ), _r2 );

    // record the source row and dst row
    src_row.push_back(insulator_node->global_offset());

    // find the position matrix row will be add to
    // since we know only one ghost node exit, there is ghost_node_begin()
    dst_row.push_back(semiconductor_node->global_offset());

    clear_row.push_back(insulator_node->global_offset());
  }

}


/*---------------------------------------------------------------------
 * build function and its jacobian for Density Gradient solver
 */
void InsulatorSemiconductorInterfaceBC::DG_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{

  const PetscScalar T   = T_external();

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const InsulatorSimulationRegion * insulator_region = dynamic_cast<const InsulatorSimulationRegion *> ( _r2 );


  //the indepedent variable number, 3 for each node
  adtl::AutoDScalar::numdir = 3;

  // after that, set values to source rows
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

    const FVM_Node * insulator_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * insulator_node_data = insulator_node->node_data();

    PetscInt index[3] = {semiconductor_node->global_offset()+0, semiconductor_node->global_offset()+1, semiconductor_node->global_offset()+2};
    AutoDScalar V_semiconductor = x[semiconductor_node->local_offset()];        V_semiconductor.setADValue(0,1.0);
    AutoDScalar n   =  x[semiconductor_node->local_offset()+1];   n.setADValue(1, 1.0);              // electron density
    AutoDScalar p   =  x[semiconductor_node->local_offset()+2];   p.setADValue(2, 1.0);              // hole density

    PetscScalar boundary_area = std::abs(semiconductor_node->outside_boundary_surface_area());

    Material::MaterialSemiconductor *mt =  semiconductor_region->material();

    //synchronize with material database
    mt->set_ad_num(adtl::AutoDScalar::numdir);

    mt->mapping(semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock);

    {
      // surface recombination
      AutoDScalar GSurf = - mt->band->R_Surf(p, n, T) * boundary_area;

      jac->add_row(  index[1],  3,  &index[0],  GSurf.getADValue() );
      jac->add_row(  index[2],  3,  &index[0],  GSurf.getADValue() );
    }

    if (semiconductor_region->get_advanced_model()->Trap)
    {
      AutoDScalar ni = mt->band->nie(p, n, T);
      mt->trap->Calculate(false,p,n,ni,T);
      AutoDScalar TrappedC = mt->trap->ChargeAD(false) * boundary_area;
      jac->add_row(  index[0],  3,  &index[0],  TrappedC.getADValue() );

      AutoDScalar GElec = - mt->trap->ElectronTrapRate(false,n,ni,T) * boundary_area;
      AutoDScalar GHole = - mt->trap->HoleTrapRate    (false,p,ni,T) * boundary_area;
      jac->add_row(  index[1],  3,  &index[0],  GElec.getADValue() );
      jac->add_row(  index[2],  3,  &index[0],  GHole.getADValue() );
    }


    // phi of insulator node
    AutoDScalar  V_insulator = x[insulator_node->local_offset()]; V_insulator.setADValue(1, 1.0);

    // the phi of this node is equal to corresponding phi of semiconductor node
    AutoDScalar  f_phi = V_insulator - V_semiconductor;

    // set Jacobian of governing equation ff
    jac->add( insulator_node->global_offset(),  insulator_node->global_offset(),  f_phi.getADValue(1) );
    jac->add( insulator_node->global_offset(),  semiconductor_node->global_offset(),  f_phi.getADValue(0) );
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

  // gate current
  _gate_current_jacobian(x, jac, add_value_flag);
}



