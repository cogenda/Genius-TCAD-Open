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

#include "log.h"
#include "jflux2.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::um;
using PhysicalUnit::cm;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


void SemiconductorSimulationRegion::DDM2_Fill_Value(Vec x, Vec L)
{
  std::vector<int> ix;
  std::vector<PetscScalar> y;
  std::vector<PetscScalar> s;

  ix.reserve(4*this->n_node());
  y.reserve(4*this->n_node());
  s.reserve(4*this->n_node());

  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();

    mt->mapping(fvm_node->root_node(), node_data, 0.0);

    // the first variable, psi
    ix.push_back(fvm_node->global_offset()+0);
    y.push_back(node_data->psi());
    s.push_back(1.0/(node_data->eps()*fvm_node->volume()));

    // the second variable, n
    ix.push_back(fvm_node->global_offset()+1);
    y.push_back(node_data->n());
    s.push_back(1.0/fvm_node->volume());

    // the third variable, p
    ix.push_back(fvm_node->global_offset()+2);
    y.push_back(node_data->p());
    s.push_back(1.0/fvm_node->volume());

    // the forth variable, lattice temperature
    ix.push_back(fvm_node->global_offset()+3);
    y.push_back(node_data->T());
    s.push_back(1.0/fvm_node->volume());
  }

  if( ix.size() )
  {
    VecSetValues(x, ix.size(), &ix[0], &y[0], INSERT_VALUES) ;
    VecSetValues(L, ix.size(), &ix[0], &s[0], INSERT_VALUES) ;
  }
}



/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void SemiconductorSimulationRegion::DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec first!
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // set local buf here
  std::vector<int>          iy;
  std::vector<PetscScalar>  y;

  // slightly overkill -- the HEX8 element has 12 edges, each edge has 2 node
  iy.reserve(4*(24*this->n_cell()+this->n_node()));
  y.reserve(4*(24*this->n_cell()+this->n_node()));


  bool  highfield_mob   = highfield_mobility() && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM;

  // first, search all the element in this region and process "cell" related terms
  // note, they are all local element, thus must be processed

  const_element_iterator it = elements_begin();
  const_element_iterator it_end = elements_end();
  for(unsigned int nelem=0 ; it!=it_end; ++it, ++nelem)
  {
    const Elem * elem = *it;

    FVM_CellData * elem_data = this->get_region_elem_data(nelem);

    bool insulator_interface_elem = is_elem_on_insulator_interface(elem);
    bool mos_channel_elem = is_elem_in_mos_channel(elem);
    bool truncation =  SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationAlways ||
        (SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationBoundary && (elem->on_boundary() || elem->on_interface())) ;
    // build the gradient of psi and fermi potential in this cell.
    // which are the vector of electric field and current density.

    VectorValue<PetscScalar> E;
    VectorValue<PetscScalar> Jnv;
    VectorValue<PetscScalar> Jpv;

    std::vector<PetscScalar> Jn_edge; //store all the edge Jn
    std::vector<PetscScalar> Jp_edge; //store all the edge Jp

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

        PetscScalar V;  // electrostatic potential
        PetscScalar n;  // electron density
        PetscScalar p;  // hole density
        PetscScalar Vt  = kb*fvm_node_data->T()/e;

        if(get_advanced_model()->HighFieldMobilitySelfConsistently)
        {
          // use current iteration value for potential, n and p still use previous solution value
          // this is the best method balanced stable and accurate.
          V  =  x[fvm_node->local_offset()+0];
          n  =  fabs(x[fvm_node->local_offset()+1]) + fvm_node_data->ni()*1e-2;
          p  =  fabs(x[fvm_node->local_offset()+2]) + fvm_node_data->ni()*1e-2;
        }
        else
        {
          // n and p use previous solution value
          V  =  x[fvm_node->local_offset()+0];
          n  =  fvm_node_data->n() + 1.0*std::pow(cm, -3);
          p  =  fvm_node_data->p() + 1.0*std::pow(cm, -3);
        }

        psi_vertex[nd] = V;
        //fermi potential
        phin_vertex[nd] = V - Vt*log(n/fvm_node_data->ni());
        phip_vertex[nd] = V + Vt*log(p/fvm_node_data->ni());
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
        unsigned int side_insul;
        SimulationRegion * region_insul;
        // find the neighbor element which has max E field
        for(unsigned int ne=0; ne<sides.size(); ++ne)
        {
          const Elem * elem_neighbor = elem->neighbor(sides[ne]);
          std::vector<PetscScalar> psi_vertex_neighbor;
          for(unsigned int nd=0; nd<elem_neighbor->n_nodes(); ++nd)
          {
            const FVM_Node * fvm_node_neighbor = elem_neighbor->get_fvm_node(nd);
            psi_vertex_neighbor.push_back(x[fvm_node_neighbor->local_offset()+0]);
          }
          VectorValue<PetscScalar> E_neighbor = - elem_neighbor->gradient(psi_vertex_neighbor);
          if(E_neighbor.size()>=E_insul.size())
          {
            E_insul = E_neighbor;
            side_insul = sides[ne];
            region_insul = regions[ne];
          }
        }
        // interface normal, point to semiconductor side
        Point norm = - elem->outside_unit_normal(side_insul);
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
        //Epn = std::max(E_eff_p.dot(Jnv.unit()), 0.0);
        //Epp = std::max(E_eff_p.dot(Jpv.unit()), 0.0);
        Epn = E_eff_p.size();
        Epp = E_eff_p.size();

        // E field vertical to current flow
        Etn = std::max(0.0,  E_eff_v_n);
        Etp = std::max(0.0, -E_eff_v_p);
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

    // process \nabla psi and S-G current along the cell's edge
    // search for all the edges this cell own
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

      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data();
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data();

      // partial area associated with this edge
      double partial_area = elem->partial_area_with_edge(ne);
      double partial_volume = elem->partial_volume_with_edge(ne);

      double truncated_partial_area =  partial_area;
      double truncated_partial_volume =  partial_volume;
      if(truncation)
      {
        // use truncated partial area to avoid negative area due to bad mesh elem
        truncated_partial_area =  this->truncated_partial_area(elem, ne);
        truncated_partial_volume =  elem->partial_volume_with_edge_truncated(ne);
      }

      const unsigned int n1_local_offset = fvm_n1->local_offset();
      const unsigned int n2_local_offset = fvm_n2->local_offset();

      // build governing equation of DDML2
      {

        //for node 1 of the edge
        mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

        PetscScalar V1   =  x[n1_local_offset+0];                  // electrostatic potential
        PetscScalar n1   =  x[n1_local_offset+1];                  // electron density
        PetscScalar p1   =  x[n1_local_offset+2];                  // hole density
        PetscScalar T1   =  x[n1_local_offset+3];                  // lattice temperature

        // NOTE: Here Ec1, Ev1 are not the conduction/valence band energy.
        // They are here for the calculation of effective driving field for electrons and holes
        // They differ from the conduction/valence band energy by the term with log(Nc), which
        // takes care of the change effective DOS.
        // Ec/Ev should not be used except when its difference between two nodes.
        // The same comment applies to Ec2/Ev2.
        PetscScalar Ec1 =  -(e*V1 + n1_data->affinity() + kb*T1*log(mt->band->nie(p1, n1, T1)));
        PetscScalar Ev1 =  -(e*V1 + n1_data->affinity() - kb*T1*log(mt->band->nie(p1, n1, T1)));
        if(get_advanced_model()->Fermi)
        {
          Ec1 = Ec1 - kb*T1*log(gamma_f(fabs(n1)/n1_data->Nc()));
          Ev1 = Ev1 + kb*T1*log(gamma_f(fabs(p1)/n1_data->Nv()));
        }

        PetscScalar eps1 =  n1_data->eps();                        // eps
        PetscScalar Eg1  =  mt->band->Eg(T1);
        PetscScalar kap1 =  mt->thermal->HeatConduction(T1);


        //for node 2 of the edge
        mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

        PetscScalar V2   =  x[n2_local_offset+0];                   // electrostatic potential
        PetscScalar n2   =  x[n2_local_offset+1];                   // electron density
        PetscScalar p2   =  x[n2_local_offset+2];                   // hole density
        PetscScalar T2   =  x[n2_local_offset+3];                   // lattice temperature


        PetscScalar Ec2 =  -(e*V2 + n2_data->affinity() + kb*T2*log(mt->band->nie(p2, n2, T2)));
        PetscScalar Ev2 =  -(e*V2 + n2_data->affinity() - kb*T2*log(mt->band->nie(p2, n2, T2)));
        if(get_advanced_model()->Fermi)
        {
          Ec2 = Ec2 - kb*T2*log(gamma_f(fabs(n2)/n2_data->Nc()));
          Ev2 = Ev2 + kb*T2*log(gamma_f(fabs(p2)/n2_data->Nv()));
        }

        PetscScalar eps2 =  n2_data->eps();                         // eps
        PetscScalar Eg2  =  mt->band->Eg(T2);
        PetscScalar kap2 =  mt->thermal->HeatConduction(T2);

        PetscScalar mun1;  // electron mobility
        PetscScalar mup1;  // hole mobility
        PetscScalar mun2;   // electron mobility
        PetscScalar mup2;   // hole mobility

        if(highfield_mob)
        {
          if (get_advanced_model()->Mob_Force == ModelSpecify::ESimple && !insulator_interface_elem )
          {
            Point dir = (*fvm_n1->root_node() - *fvm_n2->root_node()).unit();
            PetscScalar Ep = fabs((V2-V1)/length);
            PetscScalar Et = 0;
            if(mos_channel_elem)
              Et = (E - dir*(E*dir)).size();

            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T1, Ep, Et, T1);
            mup1 = mt->mob->HoleMob(p1, n1, T1, Ep, Et, T1);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T2, Ep, Et, T2);
            mup2 = mt->mob->HoleMob(p2, n2, T2, Ep, Et, T2);
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
        PetscScalar eps = 0.5*(eps1+eps2); // eps at mid point of the edge
        PetscScalar kap = 0.5*(kap1+kap2); // kapa at mid point of the edge

        // S-G current along the edge
        PetscScalar Jn =  mun*In_lt(kb,e,(Ec1-Ec2)/e,n1,n2,0.5*(T1+T2),T2-T1,length);
        PetscScalar Jp =  mup*Ip_lt(kb,e,(Ev1-Ev2)/e,p1,p2,0.5*(T1+T2),T2-T1,length);

        Jn_edge.push_back(Jn);
        Jp_edge.push_back(Jp);

        // joule heating
        PetscScalar H = 0.5*(V1-V2)*(Jn + Jp);

        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {
          // poisson's equation
          iy.push_back( fvm_n1->global_offset()+0 );
          y.push_back ( eps*(V2 - V1)/length*partial_area );

          // continuity equation of electron
          iy.push_back( fvm_n1->global_offset()+1 );
          y.push_back ( Jn*truncated_partial_area );

          // continuity equation of hole
          iy.push_back( fvm_n1->global_offset()+2 );
          y.push_back ( - Jp*truncated_partial_area );

          // heat transport equation
          iy.push_back( fvm_n1->global_offset()+3 );
          y.push_back ( kap*(T2 - T1)/length*partial_area + H*truncated_partial_area);

        }

        // for node 2.
        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {
          // poisson's equation
          iy.push_back( fvm_n2->global_offset()+0 );
          y.push_back ( eps*(V1 - V2)/length*partial_area );

          // continuity equation of electron
          iy.push_back( fvm_n2->global_offset()+1 );
          y.push_back ( - Jn*truncated_partial_area );

          // continuity equation of hole
          iy.push_back( fvm_n2->global_offset()+2 );
          y.push_back ( Jp*truncated_partial_area );

          // heat transport equation
          iy.push_back( fvm_n2->global_offset()+3 );
          y.push_back ( kap*(T1 - T2)/length*partial_area + H*truncated_partial_area);
        }

        if (get_advanced_model()->BandBandTunneling && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          PetscScalar GBTBT1 = mt->band->BB_Tunneling(T1, E.size());
          PetscScalar GBTBT2 = mt->band->BB_Tunneling(T2, E.size());

          if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            iy.push_back( fvm_n1->global_offset() + 1);
            y.push_back ( 0.5*GBTBT1*truncated_partial_volume );

            iy.push_back( fvm_n1->global_offset() + 2);
            y.push_back ( 0.5*GBTBT1*truncated_partial_volume );
          }

          if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            iy.push_back( fvm_n2->global_offset() + 1);
            y.push_back ( 0.5*GBTBT2*truncated_partial_volume );

            iy.push_back( fvm_n2->global_offset() + 2);
            y.push_back ( 0.5*GBTBT2*truncated_partial_volume );
          }
        }

        if (get_advanced_model()->ImpactIonization && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          // consider impact-ionization
          PetscScalar IIn,IIp,GIIn,GIIp;
          PetscScalar T,Eg;

          T  = std::min(T1,T2);
          Eg = 0.5* ( Eg1 + Eg2 );

          VectorValue<PetscScalar> ev = (elem->point(edge_nodes.second) - elem->point(edge_nodes.first));
          PetscScalar riin1 = 0.5 + 0.5* (ev.unit()).dot(Jnv.unit(true));
          PetscScalar riin2 = 1.0 - riin1;
          PetscScalar riip2 = 0.5 + 0.5* (ev.unit()).dot(Jpv.unit(true));
          PetscScalar riip1 = 1.0 - riip2;

          switch (get_advanced_model()->II_Force)
          {
          case ModelSpecify::IIForce_EdotJ:
            Epn = std::max(E.dot(Jnv.unit(true)), 0.0);
            Epp = std::max(E.dot(Jpv.unit(true)), 0.0);
            IIn = mt->gen->ElecGenRate(T,Epn,Eg);
            IIp = mt->gen->HoleGenRate(T,Epp,Eg);
            break;
          case ModelSpecify::EVector:
            IIn = mt->gen->ElecGenRate(T,E.size(),Eg);
            IIp = mt->gen->HoleGenRate(T,E.size(),Eg);
            break;
          case ModelSpecify::ESide:
            IIn = mt->gen->ElecGenRate(T,fabs(Ec2-Ec1)/e/length,Eg);
            IIp = mt->gen->HoleGenRate(T,fabs(Ev2-Ev1)/e/length,Eg);
            break;
          case ModelSpecify::GradQf:
            IIn = mt->gen->ElecGenRate(T,Jnv.size(),Eg);
            IIp = mt->gen->HoleGenRate(T,Jpv.size(),Eg);
            break;
          default:
            {
              MESSAGE<<"ERROR: Unsupported Impact Ionization Type."<<std::endl; RECORD();
              genius_error();
            }
          }
          GIIn = IIn * fabs(Jn)/e;
          GIIp = IIp * fabs(Jp)/e;

          if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            iy.push_back( fvm_n1->global_offset() + 1);
            y.push_back ( (riin1*GIIn+riip1*GIIp)*truncated_partial_volume );

            iy.push_back( fvm_n1->global_offset() + 2);
            y.push_back ( (riin1*GIIn+riip1*GIIp)*truncated_partial_volume );
          }

          if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            iy.push_back( fvm_n2->global_offset() + 1);
            y.push_back ( (riin2*GIIn+riip2*GIIp)*truncated_partial_volume );

            iy.push_back( fvm_n2->global_offset() + 2);
            y.push_back ( (riin2*GIIn+riip2*GIIp)*truncated_partial_volume );
          }
        }

      }

    }

    // the average cell electron/hole current density vector
    elem_data->Jn() = -elem->reconstruct_vector(Jn_edge);
    elem_data->Jp() =  elem->reconstruct_vector(Jp_edge);

  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // process node related terms
  // including \rho of poisson's equation and recombination term of continuation equation
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();

    const unsigned int local_offset  = fvm_node->local_offset();
    const unsigned int global_offset = fvm_node->global_offset();


    //PetscScalar V   =  x[local_offset+0];                       // electrostatic potential
    PetscScalar n   =  x[local_offset+1];                         // electron density
    PetscScalar p   =  x[local_offset+2];                         // hole density
    PetscScalar T   =  x[local_offset+3];                         // lattice temperature

    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);      // map this node and its data to material database
    PetscScalar R   = mt->band->Recomb(p, n, T)*fvm_node->volume();         // the recombination term
    PetscScalar rho = e*(node_data->Net_doping() + p - n)*fvm_node->volume(); // the charge density
    PetscScalar HR  = R*(node_data->Eg()+3*kb*T);                             // heat due to carrier recombination

    // consider carrier generation
    PetscScalar Field_G = node_data->Field_G()*fvm_node->volume();
    PetscScalar OptQ = node_data->OptQ()*fvm_node->volume();

    iy.push_back(global_offset+0);                                // save index in the buffer
    iy.push_back(global_offset+1);
    iy.push_back(global_offset+2);
    iy.push_back(global_offset+3);

    y.push_back( rho );                                                       // save value in the buffer
    y.push_back( Field_G - R  + node_data->EIn());
    y.push_back( Field_G - R  + node_data->HIn());
    y.push_back( HR + OptQ );

    if (get_advanced_model()->Trap)
    {
      // consider charge trapping in semiconductor bulk (bulk_flag=true)

      // call the Trap MPI to calculate trap occupancy using the local carrier densities and lattice temperature
      PetscScalar ni = mt->band->nie(p, n, T);
      mt->trap->Calculate(true,p,n,ni,T);

      // calculate the contribution of trapped charge to Poisson's equation
      PetscScalar TrappedC = mt->trap->Charge(true) * fvm_node->volume();
      if (TrappedC !=0)
      {
        iy.push_back(fvm_node->global_offset());
        y.push_back(TrappedC);
      }

      // calculate the rates of electron and hole capture
      PetscScalar TrapElec = mt->trap->ElectronTrapRate(true,n,ni,T) * fvm_node->volume();
      PetscScalar TrapHole = mt->trap->HoleTrapRate    (true,p,ni,T) * fvm_node->volume();
      // contribution to the contribution to continuity equations
      if (TrapElec != 0)
      {
        iy.push_back(fvm_node->global_offset()+1);
        // we lose carrier when electron get trapped, therefore negative contribution
        y.push_back(-TrapElec);
      }
      if (TrapHole != 0)
      {
        iy.push_back(fvm_node->global_offset()+2);
        y.push_back(-TrapHole);
      }

      PetscScalar EcEi = 0.5*node_data->Eg() - kb*T*log(node_data->Nc()/node_data->Nv());
      PetscScalar EiEv = 0.5*node_data->Eg() + kb*T*log(node_data->Nc()/node_data->Nv());
      PetscScalar H = mt->trap->TrapHeat(true,p,n,ni,T,T,T,EcEi,EiEv);
      iy.push_back(fvm_node->global_offset()+3);
      y.push_back( H*fvm_node->volume() );

    }
  }


  // add into petsc vector, we should prevent zero length vector add here.
  if(iy.size())  VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);

  // after the first scan, every nodes are updated.
  // however, boundary condition should be processed later.

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}






/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 * AD is fully used here
 */
void SemiconductorSimulationRegion::DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // note, we will use ADD_VALUES to set values of matrix J
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  bool  highfield_mob   = highfield_mobility() && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM;

  // search all the element in this region.
  // note, they are all local element, thus must be processed

  const_element_iterator it = elements_begin();
  const_element_iterator it_end = elements_end();
  for(; it!=it_end; ++it)
  {
    const Elem * elem = *it;
    bool insulator_interface_elem = is_elem_on_insulator_interface(elem);
    bool mos_channel_elem = is_elem_in_mos_channel(elem);
    bool truncation =  SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationAlways ||
        (SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationBoundary && (elem->on_boundary() || elem->on_interface())) ;

    //the indepedent variable number, 4*n_nodes
    adtl::AutoDScalar::numdir = 4*elem->n_nodes();

    //synchronize with material database
    mt->set_ad_num(adtl::AutoDScalar::numdir);

    // indicate the column position of the variables in the matrix
    std::vector<PetscInt> cell_col;
    cell_col.reserve(4*elem->n_nodes());
    for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
    {
      const FVM_Node * fvm_node = elem->get_fvm_node(nd);
      cell_col.push_back( fvm_node->global_offset()+0 );
      cell_col.push_back( fvm_node->global_offset()+1 );
      cell_col.push_back( fvm_node->global_offset()+2 );
      cell_col.push_back( fvm_node->global_offset()+3 );
    }

    // first, we build the gradient of psi and fermi potential in this cell.
    VectorValue<AutoDScalar> E;
    VectorValue<AutoDScalar> Jnv;
    VectorValue<AutoDScalar> Jpv;

    // E field parallel to current flow
    AutoDScalar Epn=0;
    AutoDScalar Epp=0;

    // E field vertical to current flow
    AutoDScalar Etn=0;
    AutoDScalar Etp=0;

    // evaluate E field parallel and vertical to current flow
    if(highfield_mob)
    {
      // which are the vector of electric field and current density.
      // here use type AutoDScalar, we should make sure the order of independent variable keeps the
      // same all the time
      std::vector<AutoDScalar> psi_vertex;
      std::vector<AutoDScalar> phin_vertex;
      std::vector<AutoDScalar> phip_vertex;

      for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
      {
        const FVM_Node * fvm_node = elem->get_fvm_node(nd);
        const FVM_NodeData * fvm_node_data = fvm_node->node_data();

        AutoDScalar V;               // electrostatic potential
        AutoDScalar n;               // electron density
        AutoDScalar p;               // hole density
        PetscScalar Vt  = kb*fvm_node_data->T()/e;

        if(get_advanced_model()->HighFieldMobilitySelfConsistently)
        {
          // use current iteration value for potential, n and p still use previous solution value
          // this is the best method balanced stable and accurate.
          V  =  x[fvm_node->local_offset()+0];   V.setADValue(4*nd+0, 1.0);
          n  =  x[fvm_node->local_offset()+1];   n.setADValue(4*nd+1, 1.0);
          p  =  x[fvm_node->local_offset()+2];   p.setADValue(4*nd+2, 1.0);
          n  +=  fvm_node_data->ni()*1e-2;
          p  +=  fvm_node_data->ni()*1e-2;
        }
        else
        {
          // n and p use previous solution value
          V  =  x[fvm_node->local_offset()+0];   V.setADValue(4*nd+0, 1.0);
          n  =  fvm_node_data->n() + 1.0*std::pow(cm, -3);
          p  =  fvm_node_data->p() + 1.0*std::pow(cm, -3);
        }

        psi_vertex.push_back  ( V );
        //fermi potential
        phin_vertex.push_back ( V - Vt*log(n/fvm_node_data->ni()) );
        phip_vertex.push_back ( V + Vt*log(p/fvm_node_data->ni()) );
      }

      // compute the gradient
      E   = - elem->gradient(psi_vertex);  // E = - grad(psi)
      Jnv = - elem->gradient(phin_vertex); // we only need the direction of Jnv, here Jnv = - gradient of Fn
      Jpv = - elem->gradient(phip_vertex); // the same as Jnv
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

        VectorValue<AutoDScalar> E_insul(0.0,0.0,0.0);
        const Elem * elem_insul;
        unsigned int side_insul;
        SimulationRegion * region_insul;

        // find the neighbor element which has max E field
        for(unsigned int ne=0; ne<sides.size(); ++ne)
        {
          const Elem * elem_neighbor = elem->neighbor(sides[ne]);

          std::vector<AutoDScalar> psi_vertex_neighbor;
          for(unsigned int nd=0; nd<elem_neighbor->n_nodes(); ++nd)
          {
            const FVM_Node * fvm_node_neighbor = elem_neighbor->get_fvm_node(nd);
            AutoDScalar V_neighbor = x[fvm_node_neighbor->local_offset()+0];
            V_neighbor.setADValue(4*elem->n_nodes()+nd, 1.0);
            psi_vertex_neighbor.push_back(V_neighbor);
          }
          VectorValue<AutoDScalar> E_neighbor = - elem_neighbor->gradient(psi_vertex_neighbor);
          if(E_neighbor.size()>=E_insul.size())
          {
            E_insul = E_neighbor;
            elem_insul = elem_neighbor;
            side_insul = sides[ne];
            region_insul = regions[ne];
          }
        }

        // we need more AD variable
        adtl::AutoDScalar::numdir += 1*elem_insul->n_nodes();
        mt->set_ad_num(adtl::AutoDScalar::numdir);
        for(unsigned int nd=0; nd<elem_insul->n_nodes(); ++nd)
        {
          const FVM_Node * fvm_node = elem_insul->get_fvm_node(nd);
          cell_col.push_back( fvm_node->global_offset()+0 );
        }

        // interface normal, point to semiconductor side
        Point _norm = - elem->outside_unit_normal(side_insul);
        // stupid code... we can not dot point with VectorValue<AutoDScalar> yet.
        VectorValue<AutoDScalar> norm(_norm(0), _norm(1), _norm(2));

        // effective electric fields in vertical
        PetscScalar ZETAN = mt->mob->ZETAN();
        PetscScalar ETAN  = mt->mob->ETAN();
        PetscScalar ZETAP = mt->mob->ZETAP();
        PetscScalar ETAP  = mt->mob->ETAP();
        AutoDScalar E_eff_v_n = ZETAN*(E*norm) + ETAN*((region_insul->get_eps()/this->get_eps())*E_insul*norm-E*norm);
        AutoDScalar E_eff_v_p = ZETAP*(E*norm) + ETAP*((region_insul->get_eps()/this->get_eps())*E_insul*norm-E*norm);
        // effective electric fields in parallel
        VectorValue<AutoDScalar> E_eff_p = E - norm*(E*norm);

        // E field parallel to current flow
        //Epn = adtl::fmax(E_eff_p.dot(Jnv.unit()), 0.0);
        //Epp = adtl::fmax(E_eff_p.dot(Jpv.unit()), 0.0);
        Epn = E_eff_p.size();
        Epp = E_eff_p.size();

        // E field vertical to current flow
        Etn = adtl::fmax(0.0,  E_eff_v_n);
        Etp = adtl::fmax(0.0, -E_eff_v_p);
      }
      else // elem NOT on insulator interface
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
          Epn = adtl::fmax(E.dot(Jnv.unit(true)), 0.0);
          Epp = adtl::fmax(E.dot(Jpv.unit(true)), 0.0);

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

      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data();
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data();

      // partial area associated with this edge
      double partial_area = elem->partial_area_with_edge(ne);
      double partial_volume = elem->partial_volume_with_edge(ne);

      double truncated_partial_area =  partial_area;
      double truncated_partial_volume =  partial_volume;
      if(truncation)
      {
        // use truncated partial area to avoid negative area due to bad mesh elem
        truncated_partial_area =  this->truncated_partial_area(elem, ne);
        truncated_partial_volume =  elem->partial_volume_with_edge_truncated(ne);
      }

      const unsigned int n1_local_offset = fvm_n1->local_offset();
      const unsigned int n2_local_offset = fvm_n2->local_offset();

      // the row position of variables in the matrix
      PetscInt row[8];
      for(unsigned int i=0; i<4; ++i) row[i]   = fvm_n1->global_offset()+i;
      for(unsigned int i=0; i<4; ++i) row[i+4] = fvm_n2->global_offset()+i;


      // here we use AD again. Can we hand write it for more efficient?
      {

        //for node 1 of the edge
        mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

        AutoDScalar V1   =  x[n1_local_offset+0];       V1.setADValue(4*edge_nodes.first+0, 1.0);           // electrostatic potential
        AutoDScalar n1   =  x[n1_local_offset+1];       n1.setADValue(4*edge_nodes.first+1, 1.0);           // electron density
        AutoDScalar p1   =  x[n1_local_offset+2];       p1.setADValue(4*edge_nodes.first+2, 1.0);           // hole density
        AutoDScalar T1   =  x[n1_local_offset+3];       T1.setADValue(4*edge_nodes.first+3, 1.0);           // lattice temperature

        AutoDScalar Ec1 =  -(e*V1 + n1_data->affinity() + kb*T1*log(mt->band->nie(p1, n1, T1)));//conduct band energy level
        AutoDScalar Ev1 =  -(e*V1 + n1_data->affinity() - kb*T1*log(mt->band->nie(p1, n1, T1)));//valence band energy level
        if(get_advanced_model()->Fermi)
        {
          Ec1 = Ec1 - kb*T1*log(gamma_f(fabs(n1)/n1_data->Nc()));
          Ev1 = Ev1 + kb*T1*log(gamma_f(fabs(p1)/n1_data->Nv()));
        }
        PetscScalar eps1 =  n1_data->eps();                        // eps
        AutoDScalar Eg1  =  mt->band->Eg(T1);
        AutoDScalar kap1 =  mt->thermal->HeatConduction(T1);


        //for node 2 of the edge
        mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

        AutoDScalar V2   =  x[n2_local_offset+0];       V2.setADValue(4*edge_nodes.second+0, 1.0);             // electrostatic potential
        AutoDScalar n2   =  x[n2_local_offset+1];       n2.setADValue(4*edge_nodes.second+1, 1.0);             // electron density
        AutoDScalar p2   =  x[n2_local_offset+2];       p2.setADValue(4*edge_nodes.second+2, 1.0);             // hole density
        AutoDScalar T2   =  x[n2_local_offset+3];       T2.setADValue(4*edge_nodes.second+3, 1.0);             // hole density

        AutoDScalar Ec2 =  -(e*V2 + n2_data->affinity() + kb*T2*log(mt->band->nie(p2, n2, T2)));//conduct band energy level
        AutoDScalar Ev2 =  -(e*V2 + n2_data->affinity() - kb*T2*log(mt->band->nie(p2, n2, T2)));//valence band energy level
        if(get_advanced_model()->Fermi)
        {
          Ec2 = Ec2 - kb*T2*log(gamma_f(fabs(n2)/n2_data->Nc()));
          Ev2 = Ev2 + kb*T2*log(gamma_f(fabs(p2)/n2_data->Nv()));
        }
        PetscScalar eps2 =  n2_data->eps();                         // eps
        AutoDScalar Eg2  = mt->band->Eg(T2);
        AutoDScalar kap2 =  mt->thermal->HeatConduction(T2);


        AutoDScalar mun1;   // electron mobility for node 1 of the edge
        AutoDScalar mup1;   // hole mobility for node 1 of the edge
        AutoDScalar mun2;   // electron mobility  for node 2 of the edge
        AutoDScalar mup2;   // hole mobility for node 2 of the edge

        if(highfield_mob)
        {
          if ( get_advanced_model()->Mob_Force == ModelSpecify::ESimple && !insulator_interface_elem )
          {
            Point _dir = (*fvm_n1->root_node() - *fvm_n2->root_node()).unit();
            VectorValue<AutoDScalar> dir(_dir(0), _dir(1), _dir(2));
            AutoDScalar Ep = fabs((V2-V1)/length);
            AutoDScalar Et = 0;//
            if(mos_channel_elem)
              Et = (E - (E*dir)*dir).size();

            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T1, Ep, Et, T1);
            mup1 = mt->mob->HoleMob(p1, n1, T1, Ep, Et, T1);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T1, Ep, Et, T2);
            mup2 = mt->mob->HoleMob(p2, n2, T1, Ep, Et, T2);
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

        AutoDScalar mun = 0.5*(mun1+mun2);  // the electron mobility at the mid point of the edge, use linear interpolation
        AutoDScalar mup = 0.5*(mup1+mup2);  // the hole mobility at the mid point of the edge, use linear interpolation


        PetscScalar eps = 0.5*(eps1+eps2); // eps at mid point of the edge
        AutoDScalar kap = 0.5*(kap1+kap2); // kapa at mid point of the edge

        // S-G current along the edge
        AutoDScalar Jn =  mun*In_lt(kb,e,(Ec1-Ec2)/e,n1,n2,0.5*(T1+T2),T2-T1,length);
        AutoDScalar Jp =  mup*Ip_lt(kb,e,(Ev1-Ev2)/e,p1,p2,0.5*(T1+T2),T2-T1,length);

        // joule heating
        AutoDScalar H = 0.5*(V1-V2)*(Jn + Jp);

#if defined(HAVE_FENV_H) && defined(DEBUG)
        genius_assert( !fetestexcept(FE_INVALID) );
#endif

        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {
          AutoDScalar ff1 = ( eps*(V2 - V1)/length*partial_area );

          AutoDScalar ff2 = ( Jn*truncated_partial_area );

          AutoDScalar ff3 = ( - Jp*truncated_partial_area );

          AutoDScalar ff4 = ( kap*(T2 - T1)/length*partial_area + H*truncated_partial_area);

          // general coding always has some overkill... bypass it.
          MatSetValues(*jac, 1, &row[0], cell_col.size(), &cell_col[0], ff1.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[1], cell_col.size(), &cell_col[0], ff2.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[2], cell_col.size(), &cell_col[0], ff3.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[3], cell_col.size(), &cell_col[0], ff4.getADValue(), ADD_VALUES);
        }

        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {
          AutoDScalar ff1 = ( eps*(V1 - V2)/length*partial_area );

          AutoDScalar ff2 = ( - Jn*truncated_partial_area );

          AutoDScalar ff3 = ( Jp*truncated_partial_area );

          AutoDScalar ff4 = ( kap*(T1 - T2)/length*partial_area + H*truncated_partial_area);

          MatSetValues(*jac, 1, &row[4], cell_col.size(), &cell_col[0], ff1.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[5], cell_col.size(), &cell_col[0], ff2.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[6], cell_col.size(), &cell_col[0], ff3.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[7], cell_col.size(), &cell_col[0], ff4.getADValue(), ADD_VALUES);
        }

        if (get_advanced_model()->BandBandTunneling && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          AutoDScalar GBTBT1 = mt->band->BB_Tunneling(T1, E.size());
          AutoDScalar GBTBT2 = mt->band->BB_Tunneling(T2, E.size());

          if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            AutoDScalar continuity = 0.5*GBTBT1*truncated_partial_volume;
            MatSetValues(*jac, 1, &row[1], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row[2], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
          }

          if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            AutoDScalar continuity = 0.5*GBTBT2*truncated_partial_volume;
            MatSetValues(*jac, 1, &row[5], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row[6], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
          }
        }

        if (get_advanced_model()->ImpactIonization && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          // consider impact-ionization
          AutoDScalar IIn,IIp,GIIn,GIIp;
          AutoDScalar T,Eg;

          // FIXME should use weighted carrier temperature.
          T  = std::min(T1,T2);
          Eg = 0.5* ( Eg1 + Eg2 );

          VectorValue<PetscScalar> ev0 = (elem->point(edge_nodes.second) - elem->point(edge_nodes.first));
          VectorValue<AutoDScalar> ev;
          ev(0)=ev0(0); ev(1)=ev0(1); ev(2)=ev0(2);
          AutoDScalar riin1 = 0.5 + 0.5* (ev.unit()).dot(Jnv.unit(true));
          AutoDScalar riin2 = 1.0 - riin1;
          AutoDScalar riip2 = 0.5 + 0.5* (ev.unit()).dot(Jpv.unit(true));
          AutoDScalar riip1 = 1.0 - riip2;

          switch (get_advanced_model()->II_Force)
          {
          case ModelSpecify::IIForce_EdotJ:
            Epn = adtl::fmax(E.dot(Jnv.unit(true)), 0.0);
            Epp = adtl::fmax(E.dot(Jpv.unit(true)), 0.0);
            IIn = mt->gen->ElecGenRate(T,Epn,Eg);
            IIp = mt->gen->HoleGenRate(T,Epp,Eg);
            break;
          case ModelSpecify::EVector:
            IIn = mt->gen->ElecGenRate(T,E.size(),Eg);
            IIp = mt->gen->HoleGenRate(T,E.size(),Eg);
            break;
          case ModelSpecify::ESide:
            IIn = mt->gen->ElecGenRate(T,fabs(Ec2-Ec1)/e/length,Eg);
            IIp = mt->gen->HoleGenRate(T,fabs(Ev2-Ev1)/e/length,Eg);
            break;
          case ModelSpecify::GradQf:
            IIn = mt->gen->ElecGenRate(T,Jnv.size(),Eg);
            IIp = mt->gen->HoleGenRate(T,Jpv.size(),Eg);
            break;
          default:
            {
              MESSAGE<<"ERROR: Unsupported Impact Ionization Type."<<std::endl; RECORD();
              genius_error();
            }
          }
          GIIn = IIn * fabs(Jn)/e;
          GIIp = IIp * fabs(Jp)/e;

          if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            AutoDScalar electron_continuity = (riin1*GIIn+riip1*GIIp)*truncated_partial_volume ;
            AutoDScalar hole_continuity     = (riin1*GIIn+riip1*GIIp)*truncated_partial_volume ;
            MatSetValues(*jac, 1, &row[1], cell_col.size(), &cell_col[0], electron_continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row[2], cell_col.size(), &cell_col[0], hole_continuity.getADValue(), ADD_VALUES);
          }

          if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            AutoDScalar electron_continuity = (riin2*GIIn+riip2*GIIp)*truncated_partial_volume ;
            AutoDScalar hole_continuity     = (riin2*GIIn+riip2*GIIp)*truncated_partial_volume ;
            MatSetValues(*jac, 1, &row[5], cell_col.size(), &cell_col[0], electron_continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row[6], cell_col.size(), &cell_col[0], hole_continuity.getADValue(), ADD_VALUES);
          }
        }

      }
    }// end of scan all edges of the cell

  }// end of scan all the cell

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // process node related terms
  // including \rho of poisson's equation and recombination term of continuation equation

  //the indepedent variable number, 4 for each node
  adtl::AutoDScalar::numdir = 4;

  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();

    PetscInt index[4] = {fvm_node->global_offset()+0, fvm_node->global_offset()+1,
                         fvm_node->global_offset()+2, fvm_node->global_offset()+3};

    AutoDScalar V   =  x[fvm_node->local_offset()+0];   V.setADValue(0, 1.0);              // psi
    AutoDScalar n   =  x[fvm_node->local_offset()+1];   n.setADValue(1, 1.0);              // electron density
    AutoDScalar p   =  x[fvm_node->local_offset()+2];   p.setADValue(2, 1.0);              // hole density
    AutoDScalar T   =  x[fvm_node->local_offset()+3];   T.setADValue(3, 1.0);              // hole density

    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);                   // map this node and its data to material database
    AutoDScalar R   = mt->band->Recomb(p, n, T)*fvm_node->volume();                      // the recombination term
    AutoDScalar rho = e*(node_data->Net_doping() + p - n)*fvm_node->volume();              // the charge density
    AutoDScalar HR  = R*(node_data->Eg()+3*kb*T);                                          // heat due to carrier recombination

    // ADD to Jacobian matrix,
    MatSetValues(*jac, 1, &index[0], 4, &index[0], rho.getADValue(), ADD_VALUES);
    MatSetValues(*jac, 1, &index[1], 4, &index[0], (-R).getADValue(),   ADD_VALUES);
    MatSetValues(*jac, 1, &index[2], 4, &index[0], (-R).getADValue(),   ADD_VALUES);
    MatSetValues(*jac, 1, &index[3], 4, &index[0], HR.getADValue(),  ADD_VALUES);

    if (get_advanced_model()->Trap)
    {
      AutoDScalar ni = mt->band->nie(p, n, T);
      mt->trap->Calculate(true,p,n,ni,T);

      AutoDScalar TrappedC = mt->trap->ChargeAD(true) * fvm_node->volume();
      MatSetValues(*jac, 1, &index[0], 4, &index[0], TrappedC.getADValue(), ADD_VALUES);

      AutoDScalar GElec = - mt->trap->ElectronTrapRate(true,n,ni,T) * fvm_node->volume();
      AutoDScalar GHole = - mt->trap->HoleTrapRate    (true,p,ni,T) * fvm_node->volume();

      MatSetValues(*jac, 1, &index[1], 4, &index[0], GElec.getADValue(), ADD_VALUES);
      MatSetValues(*jac, 1, &index[2], 4, &index[0], GHole.getADValue(), ADD_VALUES);

      AutoDScalar EcEi = 0.5*node_data->Eg() - kb*T*log(node_data->Nc()/node_data->Nv());
      AutoDScalar EiEv = 0.5*node_data->Eg() + kb*T*log(node_data->Nc()/node_data->Nv());
      AutoDScalar H = mt->trap->TrapHeat(true,p,n,ni,T,T,T,EcEi,EiEv);
      MatSetValues(*jac, 1, &index[3], 4, &index[0], (H*fvm_node->volume()).getADValue(),   ADD_VALUES);

    }

  }


  // boundary condition should be processed later!

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}




void SemiconductorSimulationRegion::DDM2_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec first!
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // set local buf here
  std::vector<int>          iy;
  std::vector<PetscScalar>  y;
  iy.reserve(4*this->n_node());
  y.reserve(4*this->n_node());

  const double r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);

  // process node related terms
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();

    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

    //PetscScalar V   =  x[fvm_node->local_offset()+0];                         // electrostatic potential
    PetscScalar n   =  x[fvm_node->local_offset()+1];                         // electron density
    PetscScalar p   =  x[fvm_node->local_offset()+2];                         // hole density
    PetscScalar T   =  x[fvm_node->local_offset()+3];                         // lattice temperature

    PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(T);

    // process \partial t

    iy.push_back(fvm_node->global_offset()+1);                                // save index in the buffer
    iy.push_back(fvm_node->global_offset()+2);
    iy.push_back(fvm_node->global_offset()+3);

    //second order
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false)
    {
      PetscScalar Tn_BDF2 = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                            / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
      PetscScalar Tp_BDF2 = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*node_data->p() + (1-r)/r*node_data->p_last())
                            / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
      PetscScalar TT_BDF2 = -((2-r)/(1-r)*T - 1.0/(r*(1-r))*node_data->T() + (1-r)/r*node_data->T_last())*node_data->density()*HeatCapacity
                            / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();

      y.push_back( Tn_BDF2 );
      y.push_back( Tp_BDF2 );
      y.push_back( TT_BDF2 );
    }
    else //first order
    {
      PetscScalar Tn_BDF1 = -(n - node_data->n())/SolverSpecify::dt*fvm_node->volume();
      PetscScalar Tp_BDF1 = -(p - node_data->p())/SolverSpecify::dt*fvm_node->volume();
      PetscScalar TT_BDF1 = -(T - node_data->T())*node_data->density()*HeatCapacity/SolverSpecify::dt*fvm_node->volume();
      y.push_back( Tn_BDF1 );
      y.push_back( Tp_BDF1 );
      y.push_back( TT_BDF1 );
    }
  }


  // add into petsc vector, we should prevent zero length vector add here.
  if(iy.size())  VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}




void SemiconductorSimulationRegion::DDM2_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of matrix J
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  //the indepedent variable number, 4 for each node
  adtl::AutoDScalar::numdir = 4;
  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  const double r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);

  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();

    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

    PetscInt index[4] = {fvm_node->global_offset()+0, fvm_node->global_offset()+1, fvm_node->global_offset()+2, fvm_node->global_offset()+3};

    AutoDScalar V   =  x[fvm_node->local_offset()+0];   V.setADValue(0, 1.0);              // psi
    AutoDScalar n   =  x[fvm_node->local_offset()+1];   n.setADValue(1, 1.0);              // electron density
    AutoDScalar p   =  x[fvm_node->local_offset()+2];   p.setADValue(2, 1.0);              // hole density
    AutoDScalar T   =  x[fvm_node->local_offset()+3];   T.setADValue(3, 1.0);              // lattice temperature

    AutoDScalar HeatCapacity =  mt->thermal->HeatCapacity(T);

    // process \partial t

    //second order
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false)
    {
      AutoDScalar Tn_BDF2 = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                       / (SolverSpecify::dt_last+SolverSpecify::dt)*fvm_node->volume();
      AutoDScalar Tp_BDF2 = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*node_data->p() + (1-r)/r*node_data->p_last())
                       / (SolverSpecify::dt_last+SolverSpecify::dt)*fvm_node->volume();
      AutoDScalar TT_BDF2 = -((2-r)/(1-r)*T - 1.0/(r*(1-r))*node_data->T() + (1-r)/r*node_data->T_last())*node_data->density()*HeatCapacity
                       / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();

      // ADD to Jacobian matrix,
      MatSetValues(*jac, 1, &index[1], 4, &index[0], Tn_BDF2.getADValue(), ADD_VALUES);
      MatSetValues(*jac, 1, &index[2], 4, &index[0], Tp_BDF2.getADValue(), ADD_VALUES);
      MatSetValues(*jac, 1, &index[3], 4, &index[0], TT_BDF2.getADValue(), ADD_VALUES);
    }
    else //first order
    {
      AutoDScalar Tn_BDF1 = -(n - node_data->n())/SolverSpecify::dt*fvm_node->volume();
      AutoDScalar Tp_BDF1 = -(p - node_data->p())/SolverSpecify::dt*fvm_node->volume();
      AutoDScalar TT_BDF1 = -(T - node_data->T())*node_data->density()*HeatCapacity/SolverSpecify::dt*fvm_node->volume();

      // ADD to Jacobian matrix,
      MatSetValues(*jac, 1, &index[1], 4, &index[0], Tn_BDF1.getADValue(), ADD_VALUES);
      MatSetValues(*jac, 1, &index[2], 4, &index[0], Tp_BDF1.getADValue(), ADD_VALUES);
      MatSetValues(*jac, 1, &index[3], 4, &index[0], TT_BDF1.getADValue(), ADD_VALUES);
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}


void SemiconductorSimulationRegion::DDM2_Pseudo_Time_Step_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec first!
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }


  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    //if( fvm_node->boundary_id() == BoundaryInfo::invalid_id ) continue;

    const FVM_NodeData * node_data = fvm_node->node_data();

    const unsigned int local_offset  = fvm_node->local_offset();
    const unsigned int global_offset = fvm_node->global_offset();

    PetscScalar V   =  x[local_offset];                         // phi
    PetscScalar f_V = -node_data->eps()*(V-node_data->psi())/SolverSpecify::PseudoTimeStepPotential*fvm_node->volume();
    VecSetValue(f, global_offset, f_V, ADD_VALUES);
  }

  for(node_it = on_processor_nodes_begin(); node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    //if( fvm_node->boundary_id() != BoundaryInfo::invalid_id ) continue;

    const FVM_NodeData * node_data = fvm_node->node_data();

    const unsigned int local_offset  = fvm_node->local_offset();
    const unsigned int global_offset = fvm_node->global_offset();

    PetscScalar n   =  x[local_offset+1];                         // electron density
    PetscScalar p   =  x[local_offset+2];                         // hole density

    PetscScalar f_n   =  -(n-node_data->n())/SolverSpecify::PseudoTimeStepCarrier*fvm_node->volume();
    VecSetValue(f, global_offset+1, f_n, ADD_VALUES);

    PetscScalar f_p   =  -(p-node_data->p())/SolverSpecify::PseudoTimeStepCarrier*fvm_node->volume();
    VecSetValue(f, global_offset+2, f_p, ADD_VALUES);
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}


void SemiconductorSimulationRegion::DDM2_Pseudo_Time_Step_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // note, we will use ADD_VALUES to set values of matrix J
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  //the indepedent variable number, 1 for each node
  adtl::AutoDScalar::numdir = 1;
  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);


  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();

  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    //if( fvm_node->boundary_id() == BoundaryInfo::invalid_id ) continue;

    const FVM_NodeData * node_data = fvm_node->node_data();

    const unsigned int local_offset  = fvm_node->local_offset();
    const unsigned int global_offset = fvm_node->global_offset();

    AutoDScalar V(x[local_offset]);   V.setADValue(0, 1.0);              // psi
    AutoDScalar f_V = -node_data->eps()*(V-node_data->psi())/SolverSpecify::PseudoTimeStepPotential*fvm_node->volume();
    MatSetValue(*jac, global_offset, global_offset, f_V.getADValue(0), ADD_VALUES);
  }


  for(node_it = on_processor_nodes_begin(); node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    //if( fvm_node->boundary_id() != BoundaryInfo::invalid_id ) continue;

    const FVM_NodeData * node_data = fvm_node->node_data();

    const unsigned int local_offset  = fvm_node->local_offset();
    const unsigned int global_offset = fvm_node->global_offset();

    AutoDScalar n(x[local_offset+1]);   n.setADValue(0, 1.0);              // electron density
    AutoDScalar p(x[local_offset+2]);   p.setADValue(0, 1.0);              // hole density

    AutoDScalar f_n = -(n-node_data->n())/SolverSpecify::PseudoTimeStepCarrier*fvm_node->volume();
    MatSetValue(*jac, global_offset+1, global_offset+1, f_n.getADValue(0), ADD_VALUES);

    AutoDScalar f_p = -(p-node_data->p())/SolverSpecify::PseudoTimeStepCarrier*fvm_node->volume();
    MatSetValue(*jac, global_offset+2, global_offset+2, f_p.getADValue(0), ADD_VALUES);
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}



int SemiconductorSimulationRegion::DDM2_Pseudo_Time_Step_Convergence_Test(PetscScalar * x)
{
  int unconverged_node_potential = 0;
  int unconverged_node_n = 0;
  int unconverged_node_p = 0;

  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();

  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    //if( fvm_node->boundary_id() == BoundaryInfo::invalid_id ) continue;

    const FVM_NodeData * node_data = fvm_node->node_data();

    const unsigned int local_offset  = fvm_node->local_offset();
    const unsigned int global_offset = fvm_node->global_offset();

    PetscScalar V   =  x[local_offset];                         // phi

    PetscScalar fV_abs = std::abs(-node_data->eps()*(V-node_data->psi())/SolverSpecify::PseudoTimeStepPotential);
    PetscScalar V_rel  = std::abs(node_data->eps()*(V-node_data->psi()))/(std::abs(V) + 1e-30);

    if( fV_abs > SolverSpecify::PseudoTimeTolRelax*SolverSpecify::poisson_abs_toler && V_rel > SolverSpecify::relative_toler )
      unconverged_node_potential++;
  }


  for(node_it = on_processor_nodes_begin(); node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    //if( fvm_node->boundary_id() != BoundaryInfo::invalid_id ) continue;

    const FVM_NodeData * node_data = fvm_node->node_data();

    const unsigned int local_offset  = fvm_node->local_offset();

    PetscScalar n   =  x[local_offset+1];                         // electron density
    PetscScalar p   =  x[local_offset+2];                         // hole density

    PetscScalar fn_abs   = std::abs(-(n-node_data->n())/SolverSpecify::PseudoTimeStepCarrier);
    PetscScalar n_rel    = std::abs((n-node_data->n())/node_data->n());

    PetscScalar fp_abs   = std::abs(-(p-node_data->p())/SolverSpecify::PseudoTimeStepCarrier);
    PetscScalar p_rel    = std::abs((p-node_data->p())/node_data->p());

    if( fn_abs > SolverSpecify::PseudoTimeTolRelax*SolverSpecify::elec_continuity_abs_toler && n_rel  > SolverSpecify::relative_toler )
      unconverged_node_n++;

    if(fp_abs > SolverSpecify::PseudoTimeTolRelax*SolverSpecify::hole_continuity_abs_toler && p_rel  > SolverSpecify::relative_toler)
      unconverged_node_p++;
  }

  return unconverged_node_potential + unconverged_node_n + unconverged_node_p;
}


void SemiconductorSimulationRegion::DDM2_Update_Solution(PetscScalar *lxx)
{

  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;

    PetscScalar V          = lxx[fvm_node->local_offset()+0];
    PetscScalar n          = lxx[fvm_node->local_offset()+1];
    PetscScalar p          = lxx[fvm_node->local_offset()+2];
    PetscScalar T          = lxx[fvm_node->local_offset()+3];

    FVM_NodeData * node_data = fvm_node->node_data();


    //update psi
    node_data->psi_last() =  node_data->psi();
    node_data->psi()      =  V;

    // electron density
    node_data->n_last()   = node_data->n();
    node_data->n()        = n;

    // hole density
    node_data->p_last()   = node_data->p();
    node_data->p()        =  p;

    // lattice temperature
    node_data->T_last()   = node_data->T();
    node_data->T()        =  T;

    // update buffered parameters dependent on temperature
    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);
    node_data->Eg() = mt->band->Eg(T) - mt->band->EgNarrow(p, n, T);
    node_data->Ec() = -(e*V + node_data->affinity() + mt->band->EgNarrowToEc(p, n, T));
    node_data->Ev() = -(e*V + node_data->affinity() - mt->band->EgNarrowToEv(p, n, T) + mt->band->Eg(T));

    if(get_advanced_model()->Fermi)
    {
      node_data->qFn() = -(e*V + node_data->affinity()) + inv_fermi_half(fabs(n/node_data->Nc()))*kb*T;
      node_data->qFp() = -(e*V + node_data->affinity() + mt->band->Eg(T)) - inv_fermi_half(fabs(p/node_data->Nv()))*kb*T;
    }
    else
    {
      node_data->qFn() = -(e*V + node_data->affinity()) + log(fabs(n/node_data->Nc()))*kb*T;
      node_data->qFp() = -(e*V + node_data->affinity() + mt->band->Eg(T)) - log(fabs(p/node_data->Nv()))*kb*T;
    }

    node_data->Nc() = mt->band->Nc(T);
    node_data->Nv() = mt->band->Nv(T);

    node_data->Recomb() = mt->band->Recomb(p, n, T);
    node_data->Recomb_Dir() = mt->band->R_Direct(p, n, T);
    node_data->Recomb_SRH() = mt->band->R_SHR(p, n, T);
    node_data->Recomb_Auger() = mt->band->R_Auger(p, n, T);

    if (get_advanced_model()->Trap)
    {
      // update traps
      mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);
      mt->trap->Update(true, p, n, node_data->ni(), T_external());
      mt->trap->Update(false, p, n, node_data->ni(), T_external());
    }

  }

  // addtional work: compute electrical field for all the cell.
  // Since this value is only used for reference.
  // It can be done simply by weighted average of cell's electrical field
  for(unsigned int n=0; n<n_cell(); ++n)
  {
    const Elem * elem = this->get_region_elem(n);
    FVM_CellData * elem_data = this->get_region_elem_data(n);

    std::vector<PetscScalar> psi_vertex;

    for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
    {
      const FVM_Node * fvm_node = elem->get_fvm_node(nd);
      const FVM_NodeData * fvm_node_data = fvm_node->node_data();
      psi_vertex.push_back  ( fvm_node_data->psi() );
    }
    // compute the gradient in the cell
    elem_data->E()  = - elem->gradient(psi_vertex);  // E = - grad(psi)
  }


}




