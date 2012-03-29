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
  if( !semiconductor_region->advanced_model().HotCarrierInjection && !semiconductor_region->advanced_model().FNTunneling) return;


  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();


  std::map<BoundaryCondition *, PetscScalar> I_to_FGate;
  std::map<BoundaryCondition *, PetscScalar> Qf_to_FGate;
  PetscScalar I=0.0;

  std::multimap<const Node *, NearestPoint>::const_iterator node_it = _node_nearest_point_map.begin();
  for(; node_it != _node_nearest_point_map.end(); ++node_it)
  {
    // fvm_node at semiconductor side
    const FVM_Node * semiconductor_node = get_region_fvm_node(node_it->first, semiconductor_region);    assert(semiconductor_node);
    assert(semiconductor_node->on_processor());

    PetscScalar T_semi = semiconductor_node->node_data()->T();
    PetscScalar Affinity_semi = semiconductor_node->node_data()->affinity();
    PetscScalar Eg_semi = semiconductor_node->node_data()->Eg();
    PetscScalar V_semi = semiconductor_node->node_data()->psi();
    PetscScalar n_semi = semiconductor_node->node_data()->n();
    PetscScalar p_semi = semiconductor_node->node_data()->p();
    PetscScalar Ec_semi = semiconductor_node->node_data()->Ec();
    PetscScalar Efn_semi = semiconductor_node->node_data()->qFn();


    // barraier (of conduction band)
    const FVM_Node * insulator_node = get_region_fvm_node(node_it->first, insulator_region);    assert(insulator_node);
    assert(insulator_node->on_processor());
    PetscScalar Affinity_ins = insulator_node->node_data()->affinity();

    PetscScalar Eb = Affinity_semi-Affinity_ins;


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

    PetscScalar I_DIR = 0.0;
    if( semiconductor_region->advanced_model().DIRTunneling )
    {
      PetscScalar alpha = std::max(1e-4, std::min(1.0, (V_gate - V_semi)/((Ec_semi+Eb-Efn_semi)/e)));
      PetscScalar J_DIR = insulator_region->material()->band->J_FN_Tunneling(E_insulator, alpha);
      I_DIR = (E_insulator > 0 ? -J_DIR : J_DIR)*semiconductor_node->outside_boundary_surface_area()*current_scale;
    }

    PetscScalar I_FN = 0.0;
    if( semiconductor_region->advanced_model().FNTunneling )
    {
      PetscScalar J_FN = insulator_region->material()->band->J_FN_Tunneling(E_insulator, 1.0);
      I_FN = (E_insulator > 0 ? -J_FN : J_FN)*semiconductor_node->outside_boundary_surface_area()*current_scale;
    }

    PetscScalar In = 0.0, Ip=0.0;
    if( semiconductor_region->advanced_model().HotCarrierInjection )
    {
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

      In = Jn_HCI*semiconductor_node->outside_boundary_surface_area()*current_scale;
      Ip = Jp_HCI*semiconductor_node->outside_boundary_surface_area()*current_scale;
    }

    switch ( bc->bc_type() )
    {
      case  ChargedContact :
      {
        I_to_FGate[bc] += Ip - In + I_FN + I_DIR;

        BoundaryCondition * charge_integral_bc = bc->inter_connect_hub();
        // only do this when transient simulation is on
        if(SolverSpecify::TimeDependent == true)
        {
          Qf_to_FGate[charge_integral_bc] += (Ip - In + I_FN + I_DIR)*SolverSpecify::dt;
        }
        break;
      }

      case GateContact     :
      {
        bc->ext_circuit()->Iapp() += (Ip - In + I_FN + I_DIR);
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
      if(I_to_FGate.find(bc)!=I_to_FGate.end())
        Ic = I_to_FGate.find(bc)->second;
      Parallel::sum(Ic);
      bc->current() += Ic;

      BoundaryCondition * charge_integral_bc = bc->inter_connect_hub();
      PetscScalar Qf = 0;
      if(Qf_to_FGate.find(charge_integral_bc)!=Qf_to_FGate.end())
        Qf = Qf_to_FGate.find(charge_integral_bc)->second;
      Parallel::sum(Qf);
      charge_integral_bc->scalar("qf") += Qf;
    }
  }




#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}

