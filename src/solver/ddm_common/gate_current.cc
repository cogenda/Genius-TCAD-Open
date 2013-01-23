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

#include <numeric>

#include "simulation_system.h"
#include "semiconductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_is.h"
#include "boundary_condition_collector.h"
#include "petsc_utils.h"
#include "parallel.h"

#define DEBUG



using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::A;
using PhysicalUnit::V;
using PhysicalUnit::cm;
using PhysicalUnit::um;
using PhysicalUnit::C;


void InsulatorSemiconductorInterfaceBC::_gate_current()
{

  // calculate gate current
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;
  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const InsulatorSimulationRegion * insulator_region = dynamic_cast<const InsulatorSimulationRegion *> ( _r2 );

    // do nothing when HotCarrierInjection/FNTunneling flag is false
  if( !semiconductor_region->advanced_model().HotCarrierInjection &&
      !semiconductor_region->advanced_model().DIRTunneling &&
      !semiconductor_region->advanced_model().FNTunneling)
    return;


  I_to_BC_Prev.clear();
  I_to_BC_Prev = I_to_BC;
  I_to_BC.clear();

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();
  PetscScalar I=0.0;

  // statistic current density
  std::vector<PetscScalar> J_CBET_Buffer;
  std::vector<PetscScalar> J_VBHT_Buffer;
  std::vector<PetscScalar> J_VBET_Buffer;

  std::multimap<const Node *, NearestPoint>::const_iterator node_it = _node_nearest_point_map.begin();
  for(; node_it != _node_nearest_point_map.end(); ++node_it)
  {
    // injection end
    SimulationRegion * region = node_it->second.region;
    BoundaryCondition * bc = node_it->second.bc;
    const Elem * bc_elem = node_it->second.elem;
    const unsigned int bc_elem_face_index = node_it->second.side;
    const Point & point_injection = node_it->second.p;

    AutoPtr<Elem> bc_elem_face = bc_elem->build_side(bc_elem_face_index, false);
    std::vector<PetscScalar> psi_elem;
    std::vector<PetscScalar> Ec_elem;
    std::vector<PetscScalar> Ev_elem;
    std::vector<PetscScalar> qfn_elem;
    std::vector<PetscScalar> qfp_elem;
    for(unsigned int n=0; n<bc_elem_face->n_nodes(); ++n)
    {
      const Node * node = bc_elem_face->get_node(n); assert(node->on_local());
      const FVM_Node * fvm_node = region->region_fvm_node(node);
      assert(fvm_node->on_local());
      psi_elem.push_back(fvm_node->node_data()->psi());
      Ec_elem.push_back(fvm_node->node_data()->Ec());
      Ev_elem.push_back(fvm_node->node_data()->Ev());
      qfn_elem.push_back(fvm_node->node_data()->qFn());
      qfp_elem.push_back(fvm_node->node_data()->qFp());
    }
    PetscScalar V_gate = bc_elem_face->interpolation(psi_elem, point_injection);
    PetscScalar Ec_gate = bc_elem_face->interpolation(Ec_elem, point_injection);
    PetscScalar Ev_gate = bc_elem_face->interpolation(Ev_elem, point_injection);
    PetscScalar Efn_gate = bc_elem_face->interpolation(qfn_elem, point_injection);
    PetscScalar Efp_gate = bc_elem_face->interpolation(qfp_elem, point_injection);
    // Ec/Ev is meaning less for metal
    /*
    if(region->type() == MetalRegion ||  region->type() == ElectrodeRegion )
    {
      Ec_gate = Efn_gate - 5.0;
      Ev_gate = Efp_gate + 5.0;
    }
    */

    // fvm_node at semiconductor side
    const FVM_Node * semiconductor_node = get_region_fvm_node(node_it->first, semiconductor_region);    assert(semiconductor_node);
    assert(semiconductor_node->on_processor());

    PetscScalar T_semi = semiconductor_node->node_data()->T();
    PetscScalar Eg_semi = semiconductor_node->node_data()->Eg();
    PetscScalar V_semi = semiconductor_node->node_data()->psi();
    PetscScalar n_semi = semiconductor_node->node_data()->n();
    PetscScalar p_semi = semiconductor_node->node_data()->p();
    PetscScalar Ec_semi = semiconductor_node->node_data()->Ec();
    PetscScalar Ev_semi = semiconductor_node->node_data()->Ev();
    PetscScalar Affinity_semi = semiconductor_node->node_data()->affinity();
    PetscScalar Efn_semi = semiconductor_node->node_data()->qFn();
    PetscScalar Efp_semi = semiconductor_node->node_data()->qFp();
    PetscScalar me_semi = semiconductor_region->material()->band->EffecElecMass(T_semi);
    PetscScalar mh_semi = semiconductor_region->material()->band->EffecHoleMass(T_semi);

    // barraier (of conduction band)
    const FVM_Node * insulator_node = get_region_fvm_node(node_it->first, insulator_region);    assert(insulator_node);
    assert(insulator_node->on_processor());
    PetscScalar Affinity_ins = insulator_node->node_data()->affinity();
    PetscScalar Eg_ins = insulator_node->node_data()->Eg();
    PetscScalar Bc_semi = -e*V_semi - Affinity_ins;
    PetscScalar Bc_gate = -e*V_gate - Affinity_ins;
    PetscScalar Bv_semi = -e*V_semi - Affinity_ins - Eg_ins;
    PetscScalar Bv_gate = -e*V_gate - Affinity_ins - Eg_ins;

    // the electrical field in insulator
    double t = (*(semiconductor_node->root_node())-point_injection).size();
    PetscScalar E_insulator = (V_gate - V_semi)/t;

    PetscScalar I_DIR = 0.0;
    if( semiconductor_region->advanced_model().DIRTunneling )
    {
      PetscScalar J_CBET = insulator_region->material()->band->J_CBET_Tunneling(me_semi, T_semi, Efn_semi, Efn_gate, Ec_semi, Ec_gate, Bc_semi, Bc_gate, t);
      PetscScalar J_VBHT = insulator_region->material()->band->J_VBHT_Tunneling(mh_semi, T_semi, Efp_semi, Efp_gate, Ev_semi, Ev_gate, Bv_semi, Bv_gate, t);
      PetscScalar J_VBET = insulator_region->material()->band->J_VBET_Tunneling(0.5*(me_semi+mh_semi), T_semi, Efn_semi, Efn_gate, Ec_semi, Ec_gate, Ev_semi, Ev_gate, Bc_semi, Bc_gate, t);
      I_DIR = (J_VBHT-J_CBET-J_VBET)*semiconductor_node->outside_boundary_surface_area()*current_scale;

      J_CBET_Buffer.push_back(J_CBET);
      J_VBHT_Buffer.push_back(J_VBHT);
      J_VBET_Buffer.push_back(J_VBET);
    }


    PetscScalar I_FN = 0.0;
    if( semiconductor_region->advanced_model().FNTunneling )
    {
      if(E_insulator > 0 ) //from semiconductor to gate
      {
        PetscScalar J_FN = insulator_region->material()->band->J_FN_Tunneling(E_insulator, 1.0);
        I_FN = -J_FN*semiconductor_node->outside_boundary_surface_area()*current_scale;
        //std::cout<<"FN " << J_FN/(A/cm/cm) << std::endl;
      }
      else
      {
        PetscScalar J_FN = insulator_region->material()->band->J_FN_Tunneling(-E_insulator, 1.0);
        I_FN = J_FN*semiconductor_node->outside_boundary_surface_area()*current_scale;
        //std::cout<<"FN " << J_FN/(A/cm/cm) << std::endl;
      }
    }


    PetscScalar In = 0.0, Ip=0.0;
    if( semiconductor_region->advanced_model().HotCarrierInjection )
    {
      const VectorValue<Real> & norm = semiconductor_node->norm(); // norm to insulator interface
      VectorValue<Real> E = -semiconductor_node->gradient(POTENTIAL, false);

      PetscScalar   E_eff_n = (E - (E*norm)*norm).size();//
      PetscScalar   E_eff_p = (E - (E*norm)*norm).size();//

      PetscScalar phi_barrier_n = insulator_region->material()->band->HCI_Barrier_n(Affinity_semi, Eg_semi, t, E_insulator);
      PetscScalar phi_barrier_p = insulator_region->material()->band->HCI_Barrier_p(Affinity_semi, Eg_semi, t, E_insulator);

      // possibility in insulator
      PetscScalar P_insulator_n = insulator_region->material()->band->HCI_Probability_Insulator_n( t, E_insulator);
      PetscScalar P_insulator_p = insulator_region->material()->band->HCI_Probability_Insulator_p( t, E_insulator );

      PetscScalar Jn_HCI = n_semi*P_insulator_n*semiconductor_region->material()->band->HCI_Integral_Fiegna_n(phi_barrier_n, E_eff_n);
      PetscScalar Jp_HCI = p_semi*P_insulator_p*semiconductor_region->material()->band->HCI_Integral_Fiegna_p(phi_barrier_p, E_eff_p);

      In = Jn_HCI*semiconductor_node->outside_boundary_surface_area()*current_scale;
      Ip = Jp_HCI*semiconductor_node->outside_boundary_surface_area()*current_scale;
    }

    switch ( bc->bc_type() )
    {
      case  ChargedContact :
      {
        I_to_BC[bc->label()] += Ip - In + I_FN + I_DIR;
        break;
      }

      case GateContact     :
      {
        // gate current, positive direction: flow into gate electrode
        bc->ext_circuit()->current() += -(Ip - In + I_FN + I_DIR);
        break;
      }
    }

    I += -(Ip - In + I_FN + I_DIR);
  }


  // sync current/ charge

  this->current() = I;
  Parallel::sum(this->current());

  BoundaryConditionCollector * bcs = system().get_bcs();
  for(unsigned int n=0; n<bcs->n_bcs(); ++n)
  {
    BoundaryCondition * bc = bcs->get_bc(n);
    if( bc->bc_type() == ChargedContact)
    {
      PetscScalar Ic = 0;
      if(I_to_BC.find(bc->label())!=I_to_BC.end())
        Ic = I_to_BC.find(bc->label())->second;
      Parallel::sum(Ic);

      if(SolverSpecify::TimeDependent == true)
      {
        PetscScalar Ic_Prev = 0;
        if(I_to_BC_Prev.find(bc->label())!=I_to_BC_Prev.end())
          Ic_Prev = I_to_BC_Prev.find(bc->label())->second;
        Parallel::sum(Ic_Prev);

        BoundaryCondition * charge_integral_bc = bc->inter_connect_hub();
        charge_integral_bc->scalar("qf") += 0.5*(Ic_Prev+Ic)*SolverSpecify::dt;
      }

      bc->current() += Ic;
    }
  }

  Parallel::allgather(J_CBET_Buffer);
  Parallel::allgather(J_VBHT_Buffer);
  Parallel::allgather(J_VBET_Buffer);

  if(J_CBET_Buffer.empty())
    scalar("J_CBET") = 0.0;
  else
    scalar("J_CBET") = std::accumulate(J_CBET_Buffer.begin(), J_CBET_Buffer.end(), 0.0)/(J_CBET_Buffer.size());

  if(J_VBHT_Buffer.empty())
    scalar("J_VBHT") = 0.0;
  else
    scalar("J_VBHT") = std::accumulate(J_VBHT_Buffer.begin(), J_VBHT_Buffer.end(), 0.0)/(J_VBHT_Buffer.size());

  if(J_VBET_Buffer.empty())
    scalar("J_VBET") = 0.0;
  else
    scalar("J_VBET") = std::accumulate(J_VBET_Buffer.begin(), J_VBET_Buffer.end(), 0.0)/(J_VBET_Buffer.size());

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}




void InsulatorSemiconductorInterfaceBC::DDM_extra_dofs(std::map<FVM_Node *, std::pair<unsigned int, unsigned int> >& dofs) const
{
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const InsulatorSimulationRegion * insulator_region = dynamic_cast<const InsulatorSimulationRegion *> ( _r2 );

  std::multimap<const Node *, NearestPoint>::const_iterator node_it = _node_nearest_point_map.begin();
  for(; node_it != _node_nearest_point_map.end(); ++node_it)
  {
    // injection begin
    FVM_Node * semiconductor_node = get_region_fvm_node(node_it->first, semiconductor_region);

    // injection end
    SimulationRegion * region = node_it->second.region;
    BoundaryCondition * bc = node_it->second.bc;
    const Elem * bc_elem = node_it->second.elem;
    const unsigned int bc_elem_face_index = node_it->second.side;
    const unsigned int nv = region->type() == SemiconductorRegion ? 3 : 1;

    AutoPtr<Elem> bc_elem_face = bc_elem->build_side(bc_elem_face_index, false);
    for(unsigned int n=0; n<bc_elem_face->n_nodes(); ++n)
    {
      Node * node = bc_elem_face->get_node(n); assert(node->on_local());
      FVM_Node * fvm_node = region->region_fvm_node(node);
      if(fvm_node->on_processor())
        dofs[semiconductor_node].first += nv;
      else
        dofs[semiconductor_node].second += nv;

      dofs[fvm_node].first += 3;
    }
  }
}



#define __F_SELF_CONSISTANCE__

void InsulatorSemiconductorInterfaceBC::_gate_current_function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const InsulatorSimulationRegion * insulator_region = dynamic_cast<const InsulatorSimulationRegion *> ( _r2 );

  if(!semiconductor_region->advanced_model().TunnelingSelfConsistently) return;


  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }


  // statistic current density
  std::vector<PetscScalar> J_CBET_Buffer;
  std::vector<PetscScalar> J_VBHT_Buffer;
  std::vector<PetscScalar> J_VBET_Buffer;

  std::multimap<const Node *, NearestPoint>::const_iterator node_it = _node_nearest_point_map.begin();
  for(; node_it != _node_nearest_point_map.end(); ++node_it)
  {
    // injection end
    SimulationRegion * region = node_it->second.region;
    BoundaryCondition * bc = node_it->second.bc;
    const Elem * bc_elem = node_it->second.elem;
    const unsigned int bc_elem_face_index = node_it->second.side;
    const Point & point_injection = node_it->second.p;
    AutoPtr<Elem> bc_elem_face = bc_elem->build_side(bc_elem_face_index, false);

    std::vector<PetscScalar> psi_elem;
    std::vector<PetscScalar> Ec_elem;
    std::vector<PetscScalar> Ev_elem;
    std::vector<PetscScalar> qfn_elem;
    std::vector<PetscScalar> qfp_elem;
    for(unsigned int v=0; v<bc_elem_face->n_nodes(); ++v)
    {
      const Node * node = bc_elem_face->get_node(v); assert(node->on_local());
      const FVM_Node * fvm_node = region->region_fvm_node(node);
      const FVM_NodeData * node_data = fvm_node->node_data();

      const unsigned int local_offset = fvm_node->local_offset();
      if(region->type() == SemiconductorRegion)
      {
        PetscScalar V = x[local_offset+0];
#ifdef __F_SELF_CONSISTANCE__
        PetscScalar n = x[local_offset+1];
        PetscScalar p = x[local_offset+2];
#else
        PetscScalar n = node_data->n();//x[local_offset+1];
        PetscScalar p = node_data->p();//x[local_offset+2];
#endif
        psi_elem.push_back(V);
        Ec_elem.push_back(-(e*V + node_data->affinity()) );
        Ev_elem.push_back(-(e*V + node_data->affinity() + node_data->Eg() ));
        qfn_elem.push_back(-(e*V + node_data->affinity()) + log(fabs(n/node_data->Nc()))*kb*node_data->T());
        qfp_elem.push_back(-(e*V + node_data->affinity() + node_data->Eg()) - log(fabs(p/node_data->Nv()))*kb*node_data->T());
      }
      else
      {
        PetscScalar V = x[local_offset];

        psi_elem.push_back(V);
        Ec_elem.push_back(-(e*V + node_data->affinity()) );
        Ev_elem.push_back(-(e*V + node_data->affinity()) );
        qfn_elem.push_back(-(e*V + node_data->affinity()) );
        qfp_elem.push_back(-(e*V + node_data->affinity()) );
      }
    }
    PetscScalar V_gate = bc_elem_face->interpolation(psi_elem, point_injection);
    PetscScalar Ec_gate = bc_elem_face->interpolation(Ec_elem, point_injection);
    PetscScalar Ev_gate = bc_elem_face->interpolation(Ev_elem, point_injection);
    PetscScalar Efn_gate = bc_elem_face->interpolation(qfn_elem, point_injection);
    PetscScalar Efp_gate = bc_elem_face->interpolation(qfp_elem, point_injection);

    // fvm_node at semiconductor side
    const FVM_Node * semiconductor_node = get_region_fvm_node(node_it->first, semiconductor_region);    assert(semiconductor_node);
    assert(semiconductor_node->on_processor());
    const unsigned int local_offset = semiconductor_node->local_offset();
    const unsigned int global_offset = semiconductor_node->global_offset();

    PetscScalar T_semi = semiconductor_node->node_data()->T();
    PetscScalar Eg_semi = semiconductor_node->node_data()->Eg();
    PetscScalar V_semi = x[local_offset+0];
#ifdef __F_SELF_CONSISTANCE__
    PetscScalar n_semi = x[local_offset+1];
    PetscScalar p_semi = x[local_offset+2];
#else
    PetscScalar n_semi = semiconductor_node->node_data()->n();//x[local_offset+1];
    PetscScalar p_semi = semiconductor_node->node_data()->p();//x[local_offset+2];
#endif
    PetscScalar Affinity_semi = semiconductor_node->node_data()->affinity();
    PetscScalar Ec_semi = -(e*V_semi + Affinity_semi);
    PetscScalar Ev_semi = Ec_semi - Eg_semi;
    PetscScalar Efn_semi = -(e*V_semi + Affinity_semi) + log(fabs(n_semi/semiconductor_node->node_data()->Nc()))*kb*T_semi;
    PetscScalar Efp_semi = -(e*V_semi + Affinity_semi + Eg_semi) - log(fabs(p_semi/semiconductor_node->node_data()->Nv()))*kb*T_semi;
    PetscScalar me_semi = semiconductor_region->material()->band->EffecElecMass(T_semi);
    PetscScalar mh_semi = semiconductor_region->material()->band->EffecHoleMass(T_semi);

    // barraier (of conduction band)
    const FVM_Node * insulator_node = get_region_fvm_node(node_it->first, insulator_region);    assert(insulator_node);
    assert(insulator_node->on_processor());
    PetscScalar Affinity_ins = insulator_node->node_data()->affinity();
    PetscScalar Eg_ins = insulator_node->node_data()->Eg();
    PetscScalar Bc_semi = -e*V_semi - Affinity_ins;
    PetscScalar Bc_gate = -e*V_gate - Affinity_ins;
    PetscScalar Bv_semi = -e*V_semi - Affinity_ins - Eg_ins;
    PetscScalar Bv_gate = -e*V_gate - Affinity_ins - Eg_ins;

    // the electrical field in insulator
    double t = (*(semiconductor_node->root_node())-point_injection).size();
    PetscScalar E_insulator = (V_gate - V_semi)/t;

    PetscScalar In_DIR=0.0, Ip_DIR=0.0;
    {
      PetscScalar J_CBET = insulator_region->material()->band->J_CBET_Tunneling(me_semi, T_semi, Efn_semi, Efn_gate, Ec_semi, Ec_gate, Bc_semi, Bc_gate, t);
      PetscScalar J_VBHT = insulator_region->material()->band->J_VBHT_Tunneling(mh_semi, T_semi, Efp_semi, Efp_gate, Ev_semi, Ev_gate, Bv_semi, Bv_gate, t);
      PetscScalar J_VBET = insulator_region->material()->band->J_VBET_Tunneling(0.5*(me_semi+mh_semi), T_semi, Efn_semi, Efn_gate, Ec_semi, Ec_gate, Ev_semi, Ev_gate, Bc_semi, Bc_gate, t);

      In_DIR = (J_CBET+J_VBET)*semiconductor_node->outside_boundary_surface_area();
      Ip_DIR = (J_VBHT)*semiconductor_node->outside_boundary_surface_area();

      J_CBET_Buffer.push_back(J_CBET);
      J_VBHT_Buffer.push_back(J_VBHT);
      J_VBET_Buffer.push_back(J_VBET);
    }

    //std::cout<<"In_DIR " << In_DIR << std::endl;
    //std::cout<<"Ip_DIR " << Ip_DIR << std::endl;

    VecSetValue(f, global_offset+1, -In_DIR, ADD_VALUES);
    VecSetValue(f, global_offset+2, -Ip_DIR, ADD_VALUES);

    unsigned int n_piece = bc_elem_face->n_nodes();
    for(unsigned int v=0; v<bc_elem_face->n_nodes(); ++v)
    {
      const Node * node = bc_elem_face->get_node(v); assert(node->on_local());
      const FVM_Node * fvm_node = region->region_fvm_node(node);

      const unsigned int global_offset = fvm_node->global_offset();
      if(region->type() == SemiconductorRegion)
      {
        VecSetValue(f, global_offset+1, In_DIR/n_piece, ADD_VALUES);
        VecSetValue(f, global_offset+2, Ip_DIR/n_piece, ADD_VALUES);
      }
      else
      {
        // current in meral region has a negative sign
        VecSetValue(f, global_offset, -(Ip_DIR+In_DIR)/n_piece, ADD_VALUES);
      }
    }

  }


  Parallel::allgather(J_CBET_Buffer);
  Parallel::allgather(J_VBHT_Buffer);
  Parallel::allgather(J_VBET_Buffer);

  if(J_CBET_Buffer.empty())
    scalar("J_CBET") = 0.0;
  else
    scalar("J_CBET") = std::accumulate(J_CBET_Buffer.begin(), J_CBET_Buffer.end(), 0.0)/(J_CBET_Buffer.size());

  if(J_VBHT_Buffer.empty())
    scalar("J_VBHT") = 0.0;
  else
    scalar("J_VBHT") = std::accumulate(J_VBHT_Buffer.begin(), J_VBHT_Buffer.end(), 0.0)/(J_VBHT_Buffer.size());

  if(J_VBET_Buffer.empty())
    scalar("J_VBET") = 0.0;
  else
    scalar("J_VBET") = std::accumulate(J_VBET_Buffer.begin(), J_VBET_Buffer.end(), 0.0)/(J_VBET_Buffer.size());


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



void InsulatorSemiconductorInterfaceBC::_gate_current_jacobian_reserve(Mat *jac, InsertMode &add_value_flag)
{
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const InsulatorSimulationRegion * insulator_region = dynamic_cast<const InsulatorSimulationRegion *> ( _r2 );

  if(!semiconductor_region->advanced_model().TunnelingSelfConsistently) return;

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  std::multimap<const Node *, NearestPoint>::const_iterator node_it = _node_nearest_point_map.begin();
  for(; node_it != _node_nearest_point_map.end(); ++node_it)
  {
    // injection begin
    FVM_Node * semiconductor_node = get_region_fvm_node(node_it->first, semiconductor_region);

    // injection end
    SimulationRegion * region = node_it->second.region;
    BoundaryCondition * bc = node_it->second.bc;
    const Elem * bc_elem = node_it->second.elem;
    const unsigned int bc_elem_face_index = node_it->second.side;
    const unsigned int nv = region->type() == SemiconductorRegion ? 3 : 1;

    AutoPtr<Elem> bc_elem_face = bc_elem->build_side(bc_elem_face_index, false);
    for(unsigned int n=0; n<bc_elem_face->n_nodes(); ++n)
    {
      Node * node = bc_elem_face->get_node(n); assert(node->on_local());
      FVM_Node * fvm_node = region->region_fvm_node(node);

      for(unsigned int v=0; v<nv; v++)
      {
        MatSetValue(*jac, semiconductor_node->global_offset()+0, fvm_node->global_offset()+v, 0, ADD_VALUES);
        MatSetValue(*jac, semiconductor_node->global_offset()+1, fvm_node->global_offset()+v, 0, ADD_VALUES);
        MatSetValue(*jac, semiconductor_node->global_offset()+2, fvm_node->global_offset()+v, 0, ADD_VALUES);

        MatSetValue(*jac, fvm_node->global_offset()+v, semiconductor_node->global_offset()+0, 0, ADD_VALUES);
        MatSetValue(*jac, fvm_node->global_offset()+v, semiconductor_node->global_offset()+1, 0, ADD_VALUES);
        MatSetValue(*jac, fvm_node->global_offset()+v, semiconductor_node->global_offset()+2, 0, ADD_VALUES);
      }
    }
  }


  
  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


#define __J_SELF_CONSISTANCE__

void InsulatorSemiconductorInterfaceBC::_gate_current_jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  const PetscScalar T   = T_external();

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const InsulatorSimulationRegion * insulator_region = dynamic_cast<const InsulatorSimulationRegion *> ( _r2 );

  if(!semiconductor_region->advanced_model().TunnelingSelfConsistently) return;

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  std::multimap<const Node *, NearestPoint>::const_iterator node_it = _node_nearest_point_map.begin();
  for(; node_it != _node_nearest_point_map.end(); ++node_it)
  {
    // injection end
    SimulationRegion * region = node_it->second.region;
    BoundaryCondition * bc = node_it->second.bc;
    const Elem * bc_elem = node_it->second.elem;
    const unsigned int bc_elem_face_index = node_it->second.side;
    const Point & point_injection = node_it->second.p;
    AutoPtr<Elem> bc_elem_face = bc_elem->build_side(bc_elem_face_index, false);

    if(region->type() == SemiconductorRegion)
      adtl::AutoDScalar::numdir = 3 + 3*bc_elem_face->n_nodes();
    else
      adtl::AutoDScalar::numdir = 3 + bc_elem_face->n_nodes();

    std::vector<PetscInt> col_index(adtl::AutoDScalar::numdir);

    std::vector<AutoDScalar> psi_elem;
    std::vector<AutoDScalar> Ec_elem;
    std::vector<AutoDScalar> Ev_elem;
    std::vector<AutoDScalar> qfn_elem;
    std::vector<AutoDScalar> qfp_elem;
    for(unsigned int v=0; v<bc_elem_face->n_nodes(); ++v)
    {
      const Node * node = bc_elem_face->get_node(v); assert(node->on_local());
      const FVM_Node * fvm_node = region->region_fvm_node(node);
      const FVM_NodeData * node_data = fvm_node->node_data();

      const unsigned int local_offset = fvm_node->local_offset();
      const unsigned int global_offset = fvm_node->global_offset();
      if(region->type() == SemiconductorRegion)
      {
        AutoDScalar V = x[local_offset+0]; V.setADValue(3 + 3*v +0, 1.0);
#ifdef  __J_SELF_CONSISTANCE__
        AutoDScalar n = x[local_offset+1]; n.setADValue(3 + 3*v +1, 1.0);
        AutoDScalar p = x[local_offset+2]; p.setADValue(3 + 3*v +2, 1.0);
#else
        AutoDScalar n = x[local_offset+1]; //n.setADValue(3 + 3*v +1, 1.0);
        AutoDScalar p = x[local_offset+2]; //p.setADValue(3 + 3*v +2, 1.0);
#endif
        col_index[3 + 3*v +0] = global_offset + 0;
        col_index[3 + 3*v +1] = global_offset + 1;
        col_index[3 + 3*v +2] = global_offset + 2;

        psi_elem.push_back(V);
        Ec_elem.push_back(-(e*V + node_data->affinity()) );
        Ev_elem.push_back(-(e*V + node_data->affinity() + node_data->Eg() ));
        qfn_elem.push_back(-(e*V + node_data->affinity()) + log(fabs(n/node_data->Nc()))*kb*node_data->T());
        qfp_elem.push_back(-(e*V + node_data->affinity() + node_data->Eg()) - log(fabs(p/node_data->Nv()))*kb*node_data->T());
      }
      else
      {
        AutoDScalar V = x[local_offset];  V.setADValue(3 + v, 1.0);
        col_index[3 + v] = global_offset;

        psi_elem.push_back(V);
        Ec_elem.push_back(-(e*V + node_data->affinity()) );
        Ev_elem.push_back(-(e*V + node_data->affinity()) );
        qfn_elem.push_back(-(e*V + node_data->affinity()) );
        qfp_elem.push_back(-(e*V + node_data->affinity()) );
      }
    }

    AutoDScalar V_gate = bc_elem_face->interpolation(psi_elem, point_injection);
    AutoDScalar Ec_gate = bc_elem_face->interpolation(Ec_elem, point_injection);
    AutoDScalar Ev_gate = bc_elem_face->interpolation(Ev_elem, point_injection);
    AutoDScalar Efn_gate = bc_elem_face->interpolation(qfn_elem, point_injection);
    AutoDScalar Efp_gate = bc_elem_face->interpolation(qfp_elem, point_injection);


    // fvm_node at semiconductor side
    const FVM_Node * semiconductor_node = get_region_fvm_node(node_it->first, semiconductor_region);    assert(semiconductor_node);
    assert(semiconductor_node->on_processor());
    const unsigned int local_offset = semiconductor_node->local_offset();
    const unsigned int global_offset = semiconductor_node->global_offset();
    col_index[0] = global_offset + 0;
    col_index[1] = global_offset + 1;
    col_index[2] = global_offset + 2;

    PetscScalar T_semi = semiconductor_node->node_data()->T();
    PetscScalar Eg_semi = semiconductor_node->node_data()->Eg();
    AutoDScalar V_semi = x[local_offset+0];  V_semi.setADValue(0, 1.0);
#ifdef  __J_SELF_CONSISTANCE__
    AutoDScalar n_semi = x[local_offset+1];  n_semi.setADValue(1, 1.0);
    AutoDScalar p_semi = x[local_offset+2];  p_semi.setADValue(2, 1.0);
#else
    PetscScalar n_semi = semiconductor_node->node_data()->n();//x[local_offset+1];
    PetscScalar p_semi = semiconductor_node->node_data()->p();//x[local_offset+2];
#endif
    PetscScalar Affinity_semi = semiconductor_node->node_data()->affinity();
    AutoDScalar Ec_semi = -(e*V_semi + Affinity_semi);
    AutoDScalar Ev_semi = Ec_semi - Eg_semi;
    AutoDScalar Efn_semi = -(e*V_semi + Affinity_semi) + log(fabs(n_semi/semiconductor_node->node_data()->Nc()))*kb*T_semi;
    AutoDScalar Efp_semi = -(e*V_semi + Affinity_semi + Eg_semi) - log(fabs(p_semi/semiconductor_node->node_data()->Nv()))*kb*T_semi;
    PetscScalar me_semi = semiconductor_region->material()->band->EffecElecMass(T_semi);
    PetscScalar mh_semi = semiconductor_region->material()->band->EffecHoleMass(T_semi);

    // barraier (of conduction band)
    const FVM_Node * insulator_node = get_region_fvm_node(node_it->first, insulator_region);    assert(insulator_node);
    assert(insulator_node->on_processor());
    PetscScalar Affinity_ins = insulator_node->node_data()->affinity();
    PetscScalar Eg_ins = insulator_node->node_data()->Eg();
    AutoDScalar Bc_semi = -e*V_semi - Affinity_ins;
    AutoDScalar Bc_gate = -e*V_gate - Affinity_ins;
    AutoDScalar Bv_semi = -e*V_semi - Affinity_ins - Eg_ins;
    AutoDScalar Bv_gate = -e*V_gate - Affinity_ins - Eg_ins;



    // the electrical field in insulator
    double t = (*(semiconductor_node->root_node())-point_injection).size();
    AutoDScalar E_insulator = (V_gate - V_semi)/t;


    AutoDScalar In_DIR=0.0, Ip_DIR=0.0;
    {
      AutoDScalar J_CBET = insulator_region->material()->band->J_CBET_Tunneling(me_semi, T_semi, Efn_semi, Efn_gate, Ec_semi, Ec_gate, Bc_semi, Bc_gate, t);
      AutoDScalar J_VBHT = insulator_region->material()->band->J_VBHT_Tunneling(mh_semi, T_semi, Efp_semi, Efp_gate, Ev_semi, Ev_gate, Bv_semi, Bv_gate, t);
      AutoDScalar J_VBET = insulator_region->material()->band->J_VBET_Tunneling(0.5*(me_semi+mh_semi), T_semi, Efn_semi, Efn_gate, Ec_semi, Ec_gate, Ev_semi, Ev_gate, Bc_semi, Bc_gate, t);

      In_DIR = (J_CBET+J_VBET)*semiconductor_node->outside_boundary_surface_area();
      Ip_DIR = (J_VBHT)*semiconductor_node->outside_boundary_surface_area();
    }

    //std::cout<<"In_DIR " << In_DIR << std::endl;
    //std::cout<<"Ip_DIR " << Ip_DIR << std::endl;

    PetscInt fn_index = global_offset + 1;
    MatSetValues(*jac, 1, &fn_index, adtl::AutoDScalar::numdir, &col_index[0], (-In_DIR).getADValue(), ADD_VALUES);

    PetscInt fp_index = global_offset + 2;
    MatSetValues(*jac, 1, &fp_index, adtl::AutoDScalar::numdir, &col_index[0], (-Ip_DIR).getADValue(), ADD_VALUES);

    unsigned int n_piece = bc_elem_face->n_nodes();
    for(unsigned int v=0; v<bc_elem_face->n_nodes(); ++v)
    {
      const Node * node = bc_elem_face->get_node(v); assert(node->on_local());
      const FVM_Node * fvm_node = region->region_fvm_node(node);

      const unsigned int global_offset = fvm_node->global_offset();
      if(region->type() == SemiconductorRegion)
      {
        PetscInt fn_index = global_offset + 1;
        MatSetValues(*jac, 1, &fn_index, adtl::AutoDScalar::numdir, &col_index[0], (In_DIR/n_piece).getADValue(), ADD_VALUES);

        PetscInt fp_index = global_offset + 2;
        MatSetValues(*jac, 1, &fp_index, adtl::AutoDScalar::numdir, &col_index[0], (Ip_DIR/n_piece).getADValue(), ADD_VALUES);
      }
      else
      {
        // current in meral region has a negative sign
        PetscInt f_index = global_offset;
        AutoDScalar I = -(In_DIR+Ip_DIR)/n_piece;
        MatSetValues(*jac, 1, &f_index, adtl::AutoDScalar::numdir, &col_index[0], I.getADValue(), ADD_VALUES);
      }
    }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif


  }
}

