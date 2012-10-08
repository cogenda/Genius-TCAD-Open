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
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_is.h"
#include "boundary_condition_collector.h"
#include "petsc_utils.h"
#include "parallel.h"


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


      if(node_it == _node_nearest_point_map.begin())
      {
        /*
        std::cout<<"Ec  " << Ec_semi <<" " << Ec_gate << std::endl;
        std::cout<<"Efn " << Efn_semi<<" " << Efn_gate << std::endl;
        std::cout<<"Bc  " << Bc_semi <<" " << Bc_gate << std::endl<< std::endl;


        std::cout<<"Ev  " << Ev_semi <<" " << Ev_gate << std::endl;
        std::cout<<"Efp " << Efp_semi<<" " << Efp_gate << std::endl;
        std::cout<<"Bv  " << Bv_semi <<" " << Bv_gate << std::endl<< std::endl;

        std::cout<< "J_CBET = " <<  J_CBET/(A/cm/cm) <<  " "
                 << "J_VBHT = " <<  J_VBHT/(A/cm/cm) <<  " "
                 << "J_VBET = " <<  J_VBET/(A/cm/cm) << std::endl;
        */
      }

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
