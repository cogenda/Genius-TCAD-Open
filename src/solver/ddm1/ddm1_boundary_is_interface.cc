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

//  $Id: ddm1_boundary_is_interface.cc,v 1.5 2008/07/09 05:58:16 gdiso Exp $


#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_is.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::A;
using PhysicalUnit::V;
using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::C;


/*---------------------------------------------------------------------
 * set scaling constant
 */
void InsulatorSemiconductorInterfaceBC::DDM1_Fill_Value(Vec , Vec L)
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
 * do pre-process to function for DDM1 solver
 */
void InsulatorSemiconductorInterfaceBC::DDM1_Function_Preprocess(Vec f, std::vector<PetscInt> &src_row,
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
 * build function and its jacobian for DDML1 solver
 */
void InsulatorSemiconductorInterfaceBC::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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


  // buffer for Vec location
  std::vector<PetscInt> iy;
  iy.reserve(n_nodes());
  // buffer for Vec value
  std::vector<PetscScalar> y;
  y.reserve(n_nodes());

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


    // process interface fixed charge density
    PetscScalar boundary_area = semiconductor_node->outside_boundary_surface_area();
    VecSetValue(f, semiconductor_node->global_offset(), this->Qf()*boundary_area, ADD_VALUES);

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
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML1 solver
 */
void InsulatorSemiconductorInterfaceBC::DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

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
          // Insulator-Semiconductor interface at Semiconductor side, we should reserve entrance for later add operator
          case SemiconductorRegion:
          {
            // semiconductor region should be the first region
            genius_assert(i==0);

            // since we know only one ghost node exit, there is ghost_node_begin()
            FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
            const FVM_Node * ghost_fvm_node = (*gn_it).first;
            MatSetValue(*jac, fvm_nodes[i]->global_offset(), ghost_fvm_node->global_offset(), 0, ADD_VALUES);

            FVM_Node::fvm_neighbor_node_iterator  gnb_it = ghost_fvm_node->neighbor_node_begin();
            for(; gnb_it != ghost_fvm_node->neighbor_node_end(); ++gnb_it)
              MatSetValue(*jac, fvm_nodes[i]->global_offset(), (*gnb_it).second->global_offset(), 0, ADD_VALUES);

            break;
          }

          // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
          case InsulatorRegion:
          {
            genius_assert(i==1);

            // reserve for later operator
            MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), 0, ADD_VALUES);

            break;
          }
          case VacuumRegion:
          break;
          default: genius_error(); //we should never reach here
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void InsulatorSemiconductorInterfaceBC::DDM1_Jacobian_Preprocess(Mat *jac, std::vector<PetscInt> &src_row,
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
 * build function and its jacobian for DDML1 solver
 */
void InsulatorSemiconductorInterfaceBC::DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

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

    PetscScalar boundary_area = semiconductor_node->outside_boundary_surface_area();

    Material::MaterialSemiconductor *mt =  semiconductor_region->material();

    //synchronize with material database
    mt->set_ad_num(adtl::AutoDScalar::numdir);

    mt->mapping(semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock);

    {
      // surface recombination
      AutoDScalar GSurf = - mt->band->R_Surf(p, n, T) * boundary_area;

      MatSetValues(*jac, 1, &index[1], 3, &index[0], GSurf.getADValue(), ADD_VALUES);
      MatSetValues(*jac, 1, &index[2], 3, &index[0], GSurf.getADValue(), ADD_VALUES);
    }

    if (semiconductor_region->get_advanced_model()->Trap)
    {
      AutoDScalar ni = mt->band->nie(p, n, T);
      mt->trap->Calculate(false,p,n,ni,T);
      AutoDScalar TrappedC = mt->trap->ChargeAD(false) * boundary_area;
      MatSetValues(*jac, 1, &index[0], 3, &index[0], TrappedC.getADValue(), ADD_VALUES);

      AutoDScalar GElec = - mt->trap->ElectronTrapRate(false,n,ni,T) * boundary_area;
      AutoDScalar GHole = - mt->trap->HoleTrapRate    (false,p,ni,T) * boundary_area;
      MatSetValues(*jac, 1, &index[1], 3, &index[0], GElec.getADValue(), ADD_VALUES);
      MatSetValues(*jac, 1, &index[2], 3, &index[0], GHole.getADValue(), ADD_VALUES);
    }


    // phi of insulator node
    AutoDScalar  V_insulator = x[insulator_node->local_offset()]; V_insulator.setADValue(1, 1.0);

    // the phi of this node is equal to corresponding phi of semiconductor node
    AutoDScalar  f_phi = V_insulator - V_semiconductor;

    // set Jacobian of governing equation ff
    MatSetValue(*jac, insulator_node->global_offset(), insulator_node->global_offset(), f_phi.getADValue(1), ADD_VALUES);
    MatSetValue(*jac, insulator_node->global_offset(), semiconductor_node->global_offset(), f_phi.getADValue(0), ADD_VALUES);
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}



void InsulatorSemiconductorInterfaceBC::DDM1_Post_Process()
{
  // calculate gate curret
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;
  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const InsulatorSimulationRegion * insulator_region = dynamic_cast<const InsulatorSimulationRegion *> ( _r2 );

  // do nothing when HotCarrierInjection flag is false
  if( !semiconductor_region->advanced_model().HotCarrierInjection ) return;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();

  PetscScalar Qf=0;

  std::multimap<const Node *, NearestPoint>::const_iterator node_it = _node_nearest_point_map.begin();
  for(; node_it != _node_nearest_point_map.end(); ++node_it)
  {
    // fvm_node at semiconductor side
    const FVM_Node * semiconductor_node = get_region_fvm_node(node_it->first, semiconductor_region);    assert(semiconductor_node);
    assert(semiconductor_node->on_processor());

    PetscScalar Affinity_semi = semiconductor_node->node_data()->affinity();
    PetscScalar Eg_semi = semiconductor_node->node_data()->Eg();
    PetscScalar V_semi = semiconductor_node->node_data()->psi();
    PetscScalar n_semi = semiconductor_node->node_data()->n();
    PetscScalar p_semi = semiconductor_node->node_data()->p();

    // injection end, at gate side
    SimulationRegion * region = node_it->second.region;
    BoundaryCondition * bc = node_it->second.bc;
    const Elem * bc_elem = node_it->second.elem;
    const unsigned int bc_elem_face_index = node_it->second.side;
    const Point & point_injection = node_it->second.p;

    AutoPtr<Elem> bc_elem_face = bc_elem->build_side(bc_elem_face_index, false);
    std::vector<PetscScalar> psi_elem;
    for(unsigned int n=0; n<bc_elem_face->n_nodes(); ++n)
    {
      const Node * node = bc_elem_face->get_node(n); assert(node->on_local());
      const FVM_Node * fvm_node = region->region_fvm_node(node);
      assert(fvm_node->on_local());
      psi_elem.push_back(fvm_node->node_data()->psi());
    }
    PetscScalar V_gate = bc_elem_face->interpolation(psi_elem, point_injection);

    // the electrical field in insulator
    double injection_distance = (*(semiconductor_node->root_node())-point_injection).size();
    PetscScalar E_insulator = (V_gate - V_semi)/injection_distance;

    //std::cout<<V_gate<<" "<<V_semi<<" "<<injection_distance/nm<<std::endl;

    const VectorValue<Real> & norm = semiconductor_node->norm(); // norm to insulator interface
    VectorValue<Real> E = -semiconductor_node->gradient(POTENTIAL, false);

    PetscScalar   E_eff_n = (E - (E*norm)*norm).size();//
    PetscScalar   E_eff_p = (E - (E*norm)*norm).size();//

    PetscScalar phi_barrier_n = insulator_region->material()->band->HCI_Barrier_n(Affinity_semi, Eg_semi, injection_distance, E_insulator);
    PetscScalar phi_barrier_p = insulator_region->material()->band->HCI_Barrier_p(Affinity_semi, Eg_semi, injection_distance, E_insulator);

    // possibility in insulator
    PetscScalar P_insulator_n = insulator_region->material()->band->HCI_Probability_Insulator_n( injection_distance, E_insulator);
    PetscScalar P_insulator_p = insulator_region->material()->band->HCI_Probability_Insulator_p( injection_distance, E_insulator );

    PetscScalar Jn_HCI = n_semi*P_insulator_n*semiconductor_region->material()->band->HCI_Integral_Fiegna_n(phi_barrier_n, E_eff_n);
    PetscScalar Jp_HCI = p_semi*P_insulator_p*semiconductor_region->material()->band->HCI_Integral_Fiegna_p(phi_barrier_p, E_eff_p);

    PetscScalar J_FN = insulator_region->material()->band->J_FN_Tunneling(E_insulator);

    PetscScalar In = Jn_HCI*semiconductor_node->outside_boundary_surface_area()*current_scale;
    PetscScalar Ip = Jp_HCI*semiconductor_node->outside_boundary_surface_area()*current_scale;
    PetscScalar I_FN = (E_insulator > 0 ? -J_FN : J_FN)*semiconductor_node->outside_boundary_surface_area()*current_scale;

    //std::cout<<E_eff_n/(V/cm)  <<" " <<phi_barrier_n<<" "<< In/A<<std::endl;

    switch ( bc->bc_type() )
    {
        case  ChargedContact :
        {
          BoundaryCondition * charge_integral_bc = bc->inter_connect_hub();
          // only do this when transient simulation is on
          if(SolverSpecify::TimeDependent == true)
          {
            charge_integral_bc->Qf() += (Ip - In + I_FN)*SolverSpecify::dt;
            Qf=charge_integral_bc->Qf();
          }
          break;
        }

        case GateContact     :
        {
          bc->ext_circuit()->Iapp() += (Ip - In + I_FN);
          break;
        }
    }
  }
  //std::cout<<"Q="<<Qf/C<<std::endl;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}

