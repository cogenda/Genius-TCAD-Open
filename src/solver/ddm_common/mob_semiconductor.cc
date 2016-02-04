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

#include "elem.h"
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "solver_specify.h"
#include "parallel.h"
#include "jflux1.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::V;
using PhysicalUnit::cm;

#define DEBUG


void SemiconductorSimulationRegion::Mob_Evaluation()
{

  const double Imin  = 1*V/cm * 1*std::pow(cm, -3); // 1V/cm * 1 carrier per cm-3

  // clear mun and mup
  for(local_node_iterator node_it = on_local_nodes_begin(); node_it!=on_local_nodes_end(); ++node_it)
  {
    FVM_Node * fvm_node = *node_it;
    FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data!=NULL);
    node_data->mun() = 0.0;
    node_data->mup() = 0.0;
  }

  std::map< FVM_Node *, double> MunWeightMap;
  std::map< FVM_Node *, double> MupWeightMap;


  bool  highfield_mob   = this->highfield_mobility() && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM;

  const_element_iterator it = elements_begin();
  const_element_iterator it_end = elements_end();
  for(; it!=it_end; ++it)
  {
    const Elem * elem = *it;
    bool insulator_interface_elem = is_elem_on_insulator_interface(elem);
    bool mos_channel_elem = is_elem_in_mos_channel(elem);
    bool truncation =  SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationAlways ||
                       (SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationBoundary && (elem->on_boundary() || elem->on_interface())) ;


    VectorValue<PetscScalar> E;
    VectorValue<PetscScalar> Jnv;
    VectorValue<PetscScalar> Jpv;


    // E field parallel to current flow
    PetscScalar Epn=0;
    PetscScalar Epp=0;

    // E field vertical to current flow
    PetscScalar Etn=0;
    PetscScalar Etp=0;

    // evaluate E field parallel and vertical to current flow
    if(highfield_mob)
    {
      // build the gradient of psi and fermi potential in this cell.
      // which are the vector of electric field and current density.
      std::vector<PetscScalar> psi_vertex(elem->n_nodes());
      std::vector<PetscScalar> phin_vertex(elem->n_nodes());
      std::vector<PetscScalar> phip_vertex(elem->n_nodes());

      for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
      {
        const FVM_Node * fvm_node = elem->get_fvm_node(nd);
        const FVM_NodeData * fvm_node_data = fvm_node->node_data();

        PetscScalar V  =  fvm_node_data->psi();
        PetscScalar n  =  std::abs(fvm_node_data->n()) + 1.0*std::pow(cm, -3);
        PetscScalar p  =  std::abs(fvm_node_data->p()) + 1.0*std::pow(cm, -3);
        PetscScalar T  =  fvm_node_data->T();

        psi_vertex[nd] = V;
        //fermi potential
        phin_vertex[nd] = V - kb*T/e*log(n/fvm_node_data->ni());
        phip_vertex[nd] = V + kb*T/e*log(p/fvm_node_data->ni());
      }

      // compute the gradient
      E   = - elem->gradient(psi_vertex);  // E = - grad(psi)
      Jnv = - elem->gradient(phin_vertex); // we only need the direction of Jnv, here Jnv = - gradient of Fn
      Jpv = - elem->gradient(phip_vertex); // Jpv = - gradient of Fp
    }

    if(highfield_mob)
    {
      // for elem on insulator interface, we will do special treatment to electrical field
      if(get_advanced_model()->ESurface && insulator_interface_elem)
      {
        // get all the sides on insulator interface
        std::vector<unsigned int> sides;
        std::vector<SimulationRegion *> regions;
        elem_on_insulator_interface(elem, sides, regions);

        VectorValue<PetscScalar> E_insul(0,0,0);
        unsigned int side_insul=invalid_uint;
        SimulationRegion * region_insul;
        // find the neighbor element which has max E field
        for(unsigned int ne=0; ne<sides.size(); ++ne)
        {
          const Elem * elem_neighbor = elem->neighbor(sides[ne]);
          std::vector<PetscScalar> psi_vertex_neighbor;
          for(unsigned int nd=0; nd<elem_neighbor->n_nodes(); ++nd)
          {
            const FVM_Node * fvm_node_neighbor = elem_neighbor->get_fvm_node(nd);
            psi_vertex_neighbor.push_back(fvm_node_neighbor->node_data()->psi());
          }
          VectorValue<PetscScalar> E_neighbor = - elem_neighbor->gradient(psi_vertex_neighbor);
          if(E_neighbor.size()>=E_insul.size())
          {
            E_insul = E_neighbor;
            side_insul = sides[ne];
            region_insul = regions[ne];
          }
        }
        if(side_insul!=invalid_uint)
        {
          // interface normal, point to semiconductor side
          Point _norm = - elem->outside_unit_normal(side_insul);

          // stupid code... we can not dot point with VectorValue<PetscScalar> yet.
          VectorValue<PetscScalar> norm(_norm(0), _norm(1), _norm(2));
          // effective electric fields in vertical
          PetscScalar ZETAN = mt->mob->ZETAN();
          PetscScalar ETAN  = mt->mob->ETAN();
          PetscScalar ZETAP = mt->mob->ZETAP();
          PetscScalar ETAP  = mt->mob->ETAP();
          PetscScalar E_eff_v_n = ZETAN*(E*norm) + ETAN*((region_insul->get_eps()/this->get_eps())*E_insul*norm-E*norm);
          PetscScalar E_eff_v_p = ZETAP*(E*norm) + ETAP*((region_insul->get_eps()/this->get_eps())*E_insul*norm-E*norm);
          // effective electric fields in parallel
          VectorValue<PetscScalar> E_eff_p = E - norm*(E*norm);

          // E field parallel to current flow
          Epn = E_eff_p.size();
          Epp = E_eff_p.size();

          // E field vertical to current flow
          Etn = std::max(0.0,  E_eff_v_n);
          Etp = std::max(0.0, -E_eff_v_p);
        }
      }
      else
      {
        if(get_advanced_model()->Mob_Force == ModelSpecify::EQF)
        {
          // E field parallel to current flow
          Epn = Jnv.size();
          Epp = Jpv.size();

          if(mos_channel_elem)
          {
            // E field vertical to current flow
            Etn = (E.cross(Jnv.unit(true))).size();
            Etp = (E.cross(Jpv.unit(true))).size();
          }
        }

        if(get_advanced_model()->Mob_Force == ModelSpecify::EJ)
        {
          // E field parallel to current flow
          Epn = std::max(E.dot(Jnv.unit(true)), 0.0);
          Epp = std::max(E.dot(Jpv.unit(true)), 0.0);

          if(mos_channel_elem)
          {
            // E field vertical to current flow
            Etn = (E.cross(Jnv.unit(true))).size();
            Etp = (E.cross(Jpv.unit(true))).size();
          }
        }
      }
    }

    // process conservation terms: laplace operator of poisson's equation and div operator of continuation equation
    // search for all the Edge this cell own
    for(unsigned int ne=0; ne<elem->n_edges(); ++ne )
    {
      std::pair<unsigned int, unsigned int> edge_nodes;
      elem->nodes_on_edge(ne, edge_nodes);

      // the length of this edge
      const double length = elem->edge_length(ne);

      // fvm_node of node1
      FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes.first);
      // fvm_node of node2
      FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes.second);

      // fvm_node_data of node1
      FVM_NodeData * n1_data =  fvm_n1->node_data();
      // fvm_node_data of node2
      FVM_NodeData * n2_data =  fvm_n2->node_data();

      // partial area associated with this edge
      double partial_area = elem->partial_area_with_edge(ne);
      double truncated_partial_area =  partial_area;
      if(truncation)
      {
        // use truncated partial area to avoid negative area due to bad mesh elem
        truncated_partial_area =  this->truncated_partial_area(elem, ne);
      }

      {

        //for node 1 of the edge
        mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

        PetscScalar V1    =  n1_data->psi();       // electrostatic potential
        PetscScalar n1    =  n1_data->n();           // electron density
        PetscScalar p1    =  n1_data->p();           // hole density
        PetscScalar T1    =  n1_data->T();           // lattice temperature

        //for node 2 of the edge
        mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

        PetscScalar V2    =  n2_data->psi();       // electrostatic potential
        PetscScalar n2    =  n2_data->n();           // electron density
        PetscScalar p2    =  n2_data->p();           // hole density
        PetscScalar T2    =  n2_data->T();           // lattice temperature

        PetscScalar mun1;  // electron mobility
        PetscScalar mup1;  // hole mobility
        PetscScalar mun2;   // electron mobility
        PetscScalar mup2;   // hole mobility

        if(highfield_mob)
        {
          if ( get_advanced_model()->Mob_Force == ModelSpecify::ESimple && !insulator_interface_elem )
          {
            PetscScalar Epn = std::abs((V1-V2)/length);
            PetscScalar Epp = std::abs((V2-V1)/length);
            PetscScalar Et = 0;
            if(mos_channel_elem)
            {
              Point _dir = (*fvm_n1->root_node() - *fvm_n2->root_node()).unit();
              VectorValue<PetscScalar> dir(_dir(0), _dir(1), _dir(2));
              Et = (E - dir*(E*dir)).size();
            }

            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T1, Epn, Et, T1);
            mup1 = mt->mob->HoleMob(p1, n1, T2, Epp, Et, T2);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T1, Epn, Et, T1);
            mup2 = mt->mob->HoleMob(p2, n2, T2, Epp, Et, T2);
          }
          else // ModelSpecify::EJ || ModelSpecify::EQF
          {
            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T1, Epn, Etn, T1);
            mup1 = mt->mob->HoleMob(p1, n1, T1, Epp, Etp, T1);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T2, Epn, Etn, T2);
            mup2 = mt->mob->HoleMob(p2, n2, T2, Epp, Etp, T2);
          }
        }
        else
        {
          mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
          mun1 = mt->mob->ElecMob(p1, n1, T1, 0, 0, T1);
          mup1 = mt->mob->HoleMob(p1, n1, T1, 0, 0, T1);

          mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
          mun2 = mt->mob->ElecMob(p2, n2, T2, 0, 0, T2);
          mup2 = mt->mob->HoleMob(p2, n2, T2, 0, 0, T2);
        }

        PetscScalar In = std::max(Imin, highfield_mob ? std::abs(In_dd(kb*T1/e,(V1-V2),n1,n2,length)) : 1.0);
        PetscScalar Ip = std::max(Imin, highfield_mob ? std::abs(Ip_dd(kb*T2/e,(V1-V2),p1,p2,length)) : 1.0);


        n1_data->mun() += mun1*truncated_partial_area*In;
        n1_data->mup() += mup1*truncated_partial_area*Ip;

        n2_data->mun() += mun2*truncated_partial_area*In;
        n2_data->mup() += mup2*truncated_partial_area*Ip;

        MunWeightMap[fvm_n1] += truncated_partial_area*In;
        MunWeightMap[fvm_n2] += truncated_partial_area*In;

        MupWeightMap[fvm_n1] += truncated_partial_area*Ip;
        MupWeightMap[fvm_n2] += truncated_partial_area*Ip;
      }
    }// end of scan all edges of the cell

  }// end of scan all the cell

  for(processor_node_iterator node_it = on_processor_nodes_begin(); node_it!=on_processor_nodes_end(); ++node_it)
  {
    FVM_Node * fvm_node = *node_it;
    FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data!=NULL);
    double surface_area_mun = MunWeightMap[fvm_node];
    double surface_area_mup = MupWeightMap[fvm_node];
    node_data->mun() /= surface_area_mun;
    node_data->mup() /= surface_area_mup;
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}


void SemiconductorSimulationRegion::Mob_Evaluation( std::vector< std::pair<unsigned int, unsigned int> > &edge,
    std::vector< std::pair<double, double> > & mob,
    std::vector< double > & weight) const
{

  bool  highfield_mob   = this->highfield_mobility();

  typedef std::pair< std::pair< double, double > , double > EdgeMob;
  typedef std::map< std::pair<unsigned int, unsigned int>,  std::vector< EdgeMob > > WeightedMobMap;
  WeightedMobMap weighted_mob_map;


  const_element_iterator it = elements_begin();
  const_element_iterator it_end = elements_end();
  for(; it!=it_end; ++it)
  {
    const Elem * elem = *it;
    bool insulator_interface_elem = is_elem_on_insulator_interface(elem);
    bool mos_channel_elem = is_elem_in_mos_channel(elem);
    bool truncation =  SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationAlways ||
                       (SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationBoundary && (elem->on_boundary() || elem->on_interface())) ;


    VectorValue<PetscScalar> E;
    VectorValue<PetscScalar> Jnv;
    VectorValue<PetscScalar> Jpv;


    // E field parallel to current flow
    PetscScalar Epn=0;
    PetscScalar Epp=0;

    // E field vertical to current flow
    PetscScalar Etn=0;
    PetscScalar Etp=0;

    // evaluate E field parallel and vertical to current flow
    if(highfield_mob)
    {
      // build the gradient of psi and fermi potential in this cell.
      // which are the vector of electric field and current density.
      std::vector<PetscScalar> psi_vertex(elem->n_nodes());
      std::vector<PetscScalar> phin_vertex(elem->n_nodes());
      std::vector<PetscScalar> phip_vertex(elem->n_nodes());

      for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
      {
        const FVM_Node * fvm_node = elem->get_fvm_node(nd);
        const FVM_NodeData * fvm_node_data = fvm_node->node_data();

        PetscScalar V  =  fvm_node_data->psi();
        PetscScalar n  =  std::abs(fvm_node_data->n()) + 1.0*std::pow(cm, -3);
        PetscScalar p  =  std::abs(fvm_node_data->p()) + 1.0*std::pow(cm, -3);
        PetscScalar T  =  fvm_node_data->T();

        psi_vertex[nd] = V;
        //fermi potential
        phin_vertex[nd] = V - kb*T/e*log(n/fvm_node_data->ni());
        phip_vertex[nd] = V + kb*T/e*log(p/fvm_node_data->ni());
      }

      // compute the gradient
      E   = - elem->gradient(psi_vertex);  // E = - grad(psi)
      Jnv = - elem->gradient(phin_vertex); // we only need the direction of Jnv, here Jnv = - gradient of Fn
      Jpv = - elem->gradient(phip_vertex); // Jpv = - gradient of Fp
    }

    if(highfield_mob)
    {
      // for elem on insulator interface, we will do special treatment to electrical field
      if(get_advanced_model()->ESurface && insulator_interface_elem)
      {
        // get all the sides on insulator interface
        std::vector<unsigned int> sides;
        std::vector<SimulationRegion *> regions;
        elem_on_insulator_interface(elem, sides, regions);

        VectorValue<PetscScalar> E_insul(0,0,0);
        unsigned int side_insul=invalid_uint;
        SimulationRegion * region_insul;
        // find the neighbor element which has max E field
        for(unsigned int ne=0; ne<sides.size(); ++ne)
        {
          const Elem * elem_neighbor = elem->neighbor(sides[ne]);
          std::vector<PetscScalar> psi_vertex_neighbor;
          for(unsigned int nd=0; nd<elem_neighbor->n_nodes(); ++nd)
          {
            const FVM_Node * fvm_node_neighbor = elem_neighbor->get_fvm_node(nd);
            psi_vertex_neighbor.push_back(fvm_node_neighbor->node_data()->psi());
          }
          VectorValue<PetscScalar> E_neighbor = - elem_neighbor->gradient(psi_vertex_neighbor);
          if(E_neighbor.size()>=E_insul.size())
          {
            E_insul = E_neighbor;
            side_insul = sides[ne];
            region_insul = regions[ne];
          }
        }
        if(side_insul!=invalid_uint)
        {
          // interface normal, point to semiconductor side
          Point _norm = - elem->outside_unit_normal(side_insul);

          // stupid code... we can not dot point with VectorValue<PetscScalar> yet.
          VectorValue<PetscScalar> norm(_norm(0), _norm(1), _norm(2));
          // effective electric fields in vertical
          PetscScalar ZETAN = mt->mob->ZETAN();
          PetscScalar ETAN  = mt->mob->ETAN();
          PetscScalar ZETAP = mt->mob->ZETAP();
          PetscScalar ETAP  = mt->mob->ETAP();
          PetscScalar E_eff_v_n = ZETAN*(E*norm) + ETAN*((region_insul->get_eps()/this->get_eps())*E_insul*norm-E*norm);
          PetscScalar E_eff_v_p = ZETAP*(E*norm) + ETAP*((region_insul->get_eps()/this->get_eps())*E_insul*norm-E*norm);
          // effective electric fields in parallel
          VectorValue<PetscScalar> E_eff_p = E - norm*(E*norm);

          // E field parallel to current flow
          Epn = E_eff_p.size();
          Epp = E_eff_p.size();

          // E field vertical to current flow
          Etn = std::max(0.0,  E_eff_v_n);
          Etp = std::max(0.0, -E_eff_v_p);
        }
      }
      else
      {
        if(get_advanced_model()->Mob_Force == ModelSpecify::EQF)
        {
          // E field parallel to current flow
          Epn = Jnv.size();
          Epp = Jpv.size();

          if(mos_channel_elem)
          {
            // E field vertical to current flow
            Etn = (E.cross(Jnv.unit(true))).size();
            Etp = (E.cross(Jpv.unit(true))).size();
          }
        }

        if(get_advanced_model()->Mob_Force == ModelSpecify::EJ)
        {
          // E field parallel to current flow
          Epn = std::max(E.dot(Jnv.unit(true)), 0.0);
          Epp = std::max(E.dot(Jpv.unit(true)), 0.0);

          if(mos_channel_elem)
          {
            // E field vertical to current flow
            Etn = (E.cross(Jnv.unit(true))).size();
            Etp = (E.cross(Jpv.unit(true))).size();
          }
        }
      }
    }

    // process conservation terms: laplace operator of poisson's equation and div operator of continuation equation
    // search for all the Edge this cell own
    for(unsigned int ne=0; ne<elem->n_edges(); ++ne )
    {
      std::pair<unsigned int, unsigned int> edge_nodes;
      elem->nodes_on_edge(ne, edge_nodes);

      // the length of this edge
      const double length = elem->edge_length(ne);

      // fvm_node of node1
      const FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes.first);
      // fvm_node of node2
      const FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes.second);

      if(fvm_n1->root_node()->id() >  fvm_n2->root_node()->id() ) std::swap(fvm_n1, fvm_n2);
      if(!fvm_n1->on_processor()) continue;

      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data();
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data();

      // partial area associated with this edge
      double partial_area = elem->partial_area_with_edge(ne);
      double truncated_partial_area =  partial_area;
      if(truncation)
      {
        // use truncated partial area to avoid negative area due to bad mesh elem
        truncated_partial_area =  this->truncated_partial_area(elem, ne);
      }

      {

        //for node 1 of the edge
        mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

        PetscScalar V1    =  n1_data->psi();       // electrostatic potential
        PetscScalar n1    =  n1_data->n();           // electron density
        PetscScalar p1    =  n1_data->p();           // hole density
        PetscScalar T1    =  n1_data->T();           // lattice temperature

        //for node 2 of the edge
        mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

        PetscScalar V2    =  n2_data->psi();       // electrostatic potential
        PetscScalar n2    =  n2_data->n();           // electron density
        PetscScalar p2    =  n2_data->p();           // hole density
        PetscScalar T2    =  n2_data->T();           // lattice temperature

        PetscScalar mun1;  // electron mobility
        PetscScalar mup1;  // hole mobility
        PetscScalar mun2;   // electron mobility
        PetscScalar mup2;   // hole mobility

        if(highfield_mob)
        {
          if ( get_advanced_model()->Mob_Force == ModelSpecify::ESimple && !insulator_interface_elem )
          {
            PetscScalar Epn = std::abs((V1-V2)/length);
            PetscScalar Epp = std::abs((V2-V1)/length);
            PetscScalar Et = 0;
            if(mos_channel_elem)
            {
              Point _dir = (*fvm_n1->root_node() - *fvm_n2->root_node()).unit();
              VectorValue<PetscScalar> dir(_dir(0), _dir(1), _dir(2));
              Et = (E - dir*(E*dir)).size();
            }

            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T1, Epn, Et, T1);
            mup1 = mt->mob->HoleMob(p1, n1, T2, Epp, Et, T2);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T1, Epn, Et, T1);
            mup2 = mt->mob->HoleMob(p2, n2, T2, Epp, Et, T2);
          }
          else // ModelSpecify::EJ || ModelSpecify::EQF
          {
            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T1, Epn, Etn, T1);
            mup1 = mt->mob->HoleMob(p1, n1, T1, Epp, Etp, T1);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T2, Epn, Etn, T2);
            mup2 = mt->mob->HoleMob(p2, n2, T2, Epp, Etp, T2);
          }
        }
        else
        {
          mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
          mun1 = mt->mob->ElecMob(p1, n1, T1, 0, 0, T1);
          mup1 = mt->mob->HoleMob(p1, n1, T1, 0, 0, T1);

          mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
          mun2 = mt->mob->ElecMob(p2, n2, T2, 0, 0, T2);
          mup2 = mt->mob->HoleMob(p2, n2, T2, 0, 0, T2);
        }

        PetscScalar mun = 0.5*(mun1+mun2); // the electron mobility at the mid point of the edge, use linear interpolation
        PetscScalar mup = 0.5*(mup1+mup2); // the hole mobility at the mid point of the edge, use linear interpolation

        unsigned int node1_id = fvm_n1->root_node()->id();
        unsigned int node2_id = fvm_n2->root_node()->id();
        std::pair<unsigned int, unsigned int> edge_key = std::make_pair(node1_id, node2_id);

        EdgeMob edge_mob;
        edge_mob.first.first = mun;
        edge_mob.first.second = mup;
        edge_mob.second = truncated_partial_area;
        weighted_mob_map[edge_key].push_back(edge_mob);

      }
    }// end of scan all edges of the cell

  }// end of scan all the cell

  // reduce to
  std::map< std::pair<unsigned int, unsigned int>,  EdgeMob > mob_map;
  WeightedMobMap::const_iterator weighted_mob_map_it = weighted_mob_map.begin();
  for(; weighted_mob_map_it != weighted_mob_map.end(); ++weighted_mob_map_it)
  {
    const std::vector< EdgeMob > & weighted_mob = weighted_mob_map_it->second;
    EdgeMob edge_mob;
    double & elec_mob = edge_mob.first.first  = 0.0;
    double & hole_mob = edge_mob.first.second =0.0;
    double & total_weight = edge_mob.second = 1e-10;
    for(unsigned int n=0; n<weighted_mob.size(); ++n)
    {
      double weight = weighted_mob[n].second;
      elec_mob += weighted_mob[n].first.first*weight;
      hole_mob += weighted_mob[n].first.second*weight;
      total_weight += weight;
    }
    elec_mob /= total_weight;
    hole_mob /= total_weight;

    mob_map.insert( std::make_pair(weighted_mob_map_it->first, edge_mob) );
  }
  weighted_mob_map.clear();

  //gather
  std::vector<unsigned int> key_array;
  std::vector<double> value_array;

  std::map< std::pair<unsigned int, unsigned int>,  EdgeMob >::const_iterator  mob_map_it = mob_map.begin();
  for(; mob_map_it != mob_map.end(); ++mob_map_it)
  {
    key_array.push_back( mob_map_it->first.first);
    key_array.push_back( mob_map_it->first.second );
    value_array.push_back( mob_map_it->second.first.first );
    value_array.push_back( mob_map_it->second.first.second );
    value_array.push_back( mob_map_it->second.second );
  }

  Parallel::allgather(key_array);
  Parallel::allgather(value_array);

  unsigned int size = key_array.size()/2;
  for(unsigned int n=0; n<size; ++n)
  {
    std::pair<unsigned int, unsigned int> key = std::make_pair(key_array[2*n], key_array[2*n+1]);
    EdgeMob edge_mob;
    edge_mob.first.first  = value_array[3*n];
    edge_mob.first.second = value_array[3*n+1];
    edge_mob.second = value_array[3*n+2];
    mob_map.insert( std::make_pair(key, edge_mob) );
  }
  key_array.clear();
  value_array.clear();

  // output
  for(mob_map_it = mob_map.begin(); mob_map_it != mob_map.end(); ++mob_map_it)
  {
    edge.push_back( mob_map_it->first );
    mob.push_back ( mob_map_it->second.first );
    weight.push_back ( mob_map_it->second.second );
  }
}
