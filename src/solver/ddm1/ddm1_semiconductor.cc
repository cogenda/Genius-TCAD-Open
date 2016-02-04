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

//  $Id: ddm1_semiconductor.cc,v 1.26 2008/07/09 05:58:16 gdiso Exp $

#include "elem.h"
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "solver_specify.h"
#include "log.h"

#include "jflux1.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::um;
using PhysicalUnit::cm;
using PhysicalUnit::us;
#define DEBUG


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


void SemiconductorSimulationRegion::DDM1_Fill_Value(Vec x, Vec L)
{
  std::vector<int> ix;
  std::vector<PetscScalar> y;
  std::vector<PetscScalar> s;

  ix.reserve(3*this->n_node());
  y.reserve(3*this->n_node());
  s.reserve(3*this->n_node());

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
  }

  if( ix.size() )
  {
    VecSetValues(x, ix.size(), &ix[0], &y[0], INSERT_VALUES) ;
    VecSetValues(L, ix.size(), &ix[0], &s[0], INSERT_VALUES) ;
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void SemiconductorSimulationRegion::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec first!
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // set local buf here

  // buffer for flux
  std::vector<PetscInt>          iflux;
  std::vector<PetscScalar>       flux;
  // slightly overkill -- the HEX8 element has 12 edges, each edge has 2 node
  iflux.reserve(3*(24*this->n_cell()));
  flux.reserve(3*(24*this->n_cell()));

  // buffer for source item
  std::vector<PetscInt>          isource;
  std::vector<PetscScalar>       source;
  isource.reserve(3*this->n_node());
  source.reserve(3*this->n_node());

  // buffer for band band tunneling
  std::vector<PetscInt>          ibbt;
  std::vector<PetscScalar>       bbt;
  if (get_advanced_model()->BandBandTunneling && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
  {
    // slightly overkill -- the HEX8 element has 12 edges, each edge has 2 node
    ibbt.reserve(3*(24*this->n_cell()));
    bbt.reserve(3*(24*this->n_cell()));
  }

  // buffer for impact ionization
  std::vector<PetscInt>          iii;
  std::vector<PetscScalar>       ii;
  if (get_advanced_model()->ImpactIonization && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
  {
    // slightly overkill -- the HEX8 element has 12 edges, each edge has 2 node
    iii.reserve(3*(24*this->n_cell()));
    ii.reserve(3*(24*this->n_cell()));

    processor_node_iterator node_it = on_processor_nodes_begin();
    processor_node_iterator node_it_end = on_processor_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      FVM_Node * fvm_node = *node_it;
      FVM_NodeData * node_data = fvm_node->node_data();

      node_data->ImpactIonization() = 0.0;
    }
  }

  //common used variable
  const PetscScalar T   = T_external();
  const PetscScalar Vt  = kb*T/e;
  bool  highfield_mob   = highfield_mobility() && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM;

  // precompute S-G current on each edge
  std::vector<PetscScalar> Jn_edge_buffer;
  std::vector<PetscScalar> Jp_edge_buffer;
  {
    Jn_edge_buffer.reserve(n_edge());
    Jp_edge_buffer.reserve(n_edge());

    // search all the edges of this region
    const_edge_iterator it = edges_begin();
    const_edge_iterator it_end = edges_end();
    for(; it!=it_end; ++it)
    {
      // fvm_node of node1
      const FVM_Node * fvm_n1 = (*it).first;
      // fvm_node of node2
      const FVM_Node * fvm_n2 = (*it).second;

      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data();
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data();

      const unsigned int n1_local_offset = fvm_n1->local_offset();
      const unsigned int n2_local_offset = fvm_n2->local_offset();

      const double length = fvm_n1->distance(fvm_n2);

      // build S-G current along edge

      //for node 1 of the edge
      mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

      const PetscScalar V1   =  x[n1_local_offset+0];                  // electrostatic potential
      const PetscScalar n1   =  x[n1_local_offset+1];                  // electron density
      const PetscScalar p1   =  x[n1_local_offset+2];                  // hole density

      // NOTE: Here Ec1, Ev1 are not the conduction/valence band energy.
      // They are here for the calculation of effective driving field for electrons and holes
      // They differ from the conduction/valence band energy by the term with kb*T*log(Nc or Nv), which
      // takes care of the change effective DOS.
      // Ec/Ev should not be used except when its difference between two nodes.
      // The same comment applies to Ec2/Ev2.
      PetscScalar Ec1 =  -(e*V1 + n1_data->affinity() - n1_data->dEcStrain() + mt->band->EgNarrowToEc(p1, n1, T) + kb*T*log(n1_data->Nc()));
      PetscScalar Ev1 =  -(e*V1 + n1_data->affinity() - n1_data->dEvStrain() - mt->band->EgNarrowToEv(p1, n1, T) - kb*T*log(n1_data->Nv()) + mt->band->Eg(T));
      if(get_advanced_model()->Fermi)
      {
        Ec1 = Ec1 - kb*T*log(gamma_f(fabs(n1)/n1_data->Nc()));
        Ev1 = Ev1 + kb*T*log(gamma_f(fabs(p1)/n1_data->Nv()));
      }
      const PetscScalar eps1 =  n1_data->eps();

      //for node 2 of the edge
      mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

      const PetscScalar V2   =  x[n2_local_offset+0];                   // electrostatic potential
      const PetscScalar n2   =  x[n2_local_offset+1];                   // electron density
      const PetscScalar p2   =  x[n2_local_offset+2];                   // hole density

      PetscScalar Ec2 =  -(e*V2 + n2_data->affinity() - n2_data->dEcStrain() + mt->band->EgNarrowToEc(p2, n2, T) + kb*T*log(n2_data->Nc()));
      PetscScalar Ev2 =  -(e*V2 + n2_data->affinity() - n2_data->dEvStrain() - mt->band->EgNarrowToEv(p2, n2, T) - kb*T*log(n2_data->Nv()) + mt->band->Eg(T));
      if(get_advanced_model()->Fermi)
      {
        Ec2 = Ec2 - kb*T*log(gamma_f(fabs(n2)/n2_data->Nc()));
        Ev2 = Ev2 + kb*T*log(gamma_f(fabs(p2)/n2_data->Nv()));
      }
      const PetscScalar eps2 =  n2_data->eps();

      // S-G current along the edge
      Jn_edge_buffer.push_back( In_dd(Vt,(Ec2-Ec1)/e,n1,n2,length) );
      Jp_edge_buffer.push_back( Ip_dd(Vt,(Ev2-Ev1)/e,p1,p2,length) );


      // poisson's equation

      PetscScalar eps = 0.5*(eps1+eps2);

      // "flux" from node 2 to node 1
      PetscScalar f =  eps*fvm_n1->cv_surface_area(fvm_n2)*(V2 - V1)/fvm_n1->distance(fvm_n2) ;

      // ignore thoese ghost nodes
      if( fvm_n1->on_processor() )
      {
        iflux.push_back(fvm_n1->global_offset());
        flux.push_back(f);
      }

      if( fvm_n2->on_processor() )
      {
        iflux.push_back(fvm_n2->global_offset());
        flux.push_back(-f);
      }
    }
  }

  // then, search all the element in this region and process "cell" related terms
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
        (SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationBoundary && is_elem_touch_boundary(elem)) ;

    std::vector<PetscScalar> Jn_edge_cell; //store all the edge Jn
    std::vector<PetscScalar> Jp_edge_cell; //store all the edge Jp

    // E field parallel to current flow
    PetscScalar Epn=0;
    PetscScalar Epp=0;

    // E field vertical to current flow
    PetscScalar Etn=0;
    PetscScalar Etp=0;


    VectorValue<PetscScalar> E;
    VectorValue<PetscScalar> Jnv;
    VectorValue<PetscScalar> Jpv;

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

        if(get_advanced_model()->HighFieldMobilitySelfConsistently)
        {
          double truc = get_advanced_model()->QuasiFermiCarrierTruc;
          // use values in the current iteration
          V  =  x[fvm_node->local_offset()+0];
          n  =  std::max(x[fvm_node->local_offset()+1], truc*fvm_node_data->ni());
          p  =  std::max(x[fvm_node->local_offset()+2], truc*fvm_node_data->ni());
        }
        else
        {
          // n and p will use previous solution value
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
      E = - elem->gradient(psi_vertex);  // E = - grad(psi)
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
        const Elem * elem_insul=elem->neighbor(sides[0]);
        SimulationRegion * region_insul=regions[0];

        std::vector<PetscScalar> psi_vertex_neighbor;
        for(unsigned int nd=0; nd<elem_insul->n_nodes(); ++nd)
        {
          const FVM_Node * fvm_node_neighbor = elem_insul->get_fvm_node(nd);
          PetscScalar V_neighbor = x[fvm_node_neighbor->local_offset()+0];
          psi_vertex_neighbor.push_back(V_neighbor);
        }

        VectorValue<PetscScalar> E_insul    = - elem_insul->gradient(psi_vertex_neighbor);

        // interface normal, point to semiconductor side
        Point _norm = - elem->outside_unit_normal(sides[0]);
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
        //Epn = std::max(E_eff_p.dot(Jnv.unit()), 0.0);
        //Epp = std::max(E_eff_p.dot(Jpv.unit()), 0.0);
        Epn = E_eff_p.size();
        Epp = E_eff_p.size();

        // E field vertical to current flow
        Etn = std::max(0.0,  E_eff_v_n);
        Etp = std::max(0.0, -E_eff_v_p);
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

      const unsigned int edge_index = this->elem_edge_index(elem, ne);

      const double length = elem->edge_length(ne);                         // the length of this edge

      FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes.first);   // fvm_node of node1
      FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes.second);  // fvm_node of node2

      double partial_area = elem->partial_area_with_edge(ne);        // partial area associated with this edge
      double partial_volume = elem->partial_volume_with_edge(ne);    // partial volume associated with this edge
      double truncated_partial_area =  partial_area;
      double truncated_partial_volume =  partial_volume;
      if(truncation)
      {
        // use truncated partial area to avoid negative area due to bad mesh elem
        truncated_partial_area =  this->truncated_partial_area(elem, ne);
        truncated_partial_volume =  elem->partial_volume_with_edge_truncated(ne);
      }


      bool inverse = fvm_n1->root_node()->id() > fvm_n2->root_node()->id();

      FVM_NodeData * n1_data = fvm_n1->node_data();            // fvm_node_data of node1
      FVM_NodeData * n2_data = fvm_n2->node_data();            // fvm_node_data of node2

      const unsigned int n1_local_offset  = fvm_n1->local_offset();
      const unsigned int n2_local_offset  = fvm_n2->local_offset();
      const unsigned int n1_global_offset = fvm_n1->global_offset();
      const unsigned int n2_global_offset = fvm_n2->global_offset();

      // build governing equation of DDML1
      {

        //for node 1 of the edge
        const PetscScalar V1   =  x[n1_local_offset+0];                  // electrostatic potential
        const PetscScalar n1   =  x[n1_local_offset+1];                  // electron density
        const PetscScalar p1   =  x[n1_local_offset+2];                  // hole density

        //for node 2 of the edge
        const PetscScalar V2   =  x[n2_local_offset+0];                   // electrostatic potential
        const PetscScalar n2   =  x[n2_local_offset+1];                   // electron density
        const PetscScalar p2   =  x[n2_local_offset+2];                   // hole density

        PetscScalar mun1;  // electron mobility
        PetscScalar mup1;  // hole mobility
        PetscScalar mun2;   // electron mobility
        PetscScalar mup2;   // hole mobility

        if(highfield_mob)
        {
          if(get_advanced_model()->ESurface && insulator_interface_elem)
          {
            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T, Epn, Etn, T);
            mup1 = mt->mob->HoleMob(p1, n1, T, Epp, Etp, T);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T, Epn, Etn, T);
            mup2 = mt->mob->HoleMob(p2, n2, T, Epp, Etp, T);
          }
          else
          {
            // high field mobility
            if (get_advanced_model()->Mob_Force == ModelSpecify::ESimple)
            {
              PetscScalar Epn = std::abs((V1-V2)/length);
              PetscScalar Epp = std::abs((V1-V2)/length);
              //PetscScalar Epn = std::max(0.0, (inverse ? -1 : 1) *(V1-V2)/length);
              //PetscScalar Epp = std::max(0.0, (inverse ? -1 : 1) *(V2-V1)/length);
              PetscScalar Et = 0;
              if(mos_channel_elem)
              {
                Point _dir = (*fvm_n1->root_node() - *fvm_n2->root_node()).unit();
                VectorValue<PetscScalar> dir(_dir(0), _dir(1), _dir(2));
                Et = (E - dir*(E*dir)).size();
              }

              mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
              mun1 = mt->mob->ElecMob(p1, n1, T, Epn, Et, T);
              mup1 = mt->mob->HoleMob(p1, n1, T, Epp, Et, T);

              mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
              mun2 = mt->mob->ElecMob(p2, n2, T, Epn, Et, T);
              mup2 = mt->mob->HoleMob(p2, n2, T, Epp, Et, T);
            }
            else// ModelSpecify::EJ || ModelSpecify::EQF
            {
              mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
              mun1 = mt->mob->ElecMob(p1, n1, T, Epn, Etn, T);
              mup1 = mt->mob->HoleMob(p1, n1, T, Epp, Etp, T);

              mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
              mun2 = mt->mob->ElecMob(p2, n2, T, Epn, Etn, T);
              mup2 = mt->mob->HoleMob(p2, n2, T, Epp, Etp, T);
            }
          }
        }
        else // low field mobility
        {
          mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
          mun1 = mt->mob->ElecMob(p1, n1, T, 0, 0, T);
          mup1 = mt->mob->HoleMob(p1, n1, T, 0, 0, T);

          mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
          mun2 = mt->mob->ElecMob(p2, n2, T, 0, 0, T);
          mup2 = mt->mob->HoleMob(p2, n2, T, 0, 0, T);
        }


        const PetscScalar mun = 0.5*(mun1+mun2); // the electron mobility at the mid point of the edge, use linear interpolation
        const PetscScalar mup = 0.5*(mup1+mup2); // the hole mobility at the mid point of the edge, use linear interpolation

        const PetscScalar eps = 0.5*(n1_data->eps()+n2_data->eps()); // eps at mid point of the edge

        // S-G current along the edge, use precomputed value
        PetscScalar Jn =  mun*(inverse ? -Jn_edge_buffer[edge_index] : Jn_edge_buffer[edge_index]);
        PetscScalar Jp =  mup*(inverse ? -Jp_edge_buffer[edge_index] : Jp_edge_buffer[edge_index]);

        Jn_edge_cell.push_back(Jn);
        Jp_edge_cell.push_back(Jp);


        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->on_processor() )
        {
          // poisson's equation
          //iflux.push_back( n1_global_offset );
          //flux.push_back ( eps*(V2 - V1)/length*partial_area );

          // continuity equation of electron
          iflux.push_back( n1_global_offset+1 );
          flux.push_back ( Jn*truncated_partial_area );

          // continuity equation of hole
          iflux.push_back( n1_global_offset+2 );
          flux.push_back ( - Jp*truncated_partial_area );
        }

        // for node 2.
        if( fvm_n2->on_processor() )
        {
          // poisson's equation
          //iflux.push_back( n2_global_offset );
          //flux.push_back ( -eps*(V2 - V1)/length*partial_area );

          // continuity equation of electron
          iflux.push_back( n2_global_offset+1);
          flux.push_back ( -Jn*truncated_partial_area );

          // continuity equation of hole
          iflux.push_back( n2_global_offset+2);
          flux.push_back ( Jp*truncated_partial_area );
        }

        if (get_advanced_model()->BandBandTunneling && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {

          PetscScalar GBTBT1 = mt->band->BB_Tunneling(T, E.size());
          PetscScalar GBTBT2 = mt->band->BB_Tunneling(T, E.size());

          if( fvm_n1->on_processor() )
          {
            // continuity equation
            ibbt.push_back( n1_global_offset + 1);
            bbt.push_back ( 0.5*GBTBT1*truncated_partial_volume );

            ibbt.push_back( n1_global_offset + 2);
            bbt.push_back ( 0.5*GBTBT1*truncated_partial_volume );
          }

          if( fvm_n2->on_processor() )
          {
            // continuity equation
            ibbt.push_back( n2_global_offset + 1);
            bbt.push_back ( 0.5*GBTBT2*truncated_partial_volume );

            ibbt.push_back( n2_global_offset + 2);
            bbt.push_back ( 0.5*GBTBT2*truncated_partial_volume );
          }
        }

        if (get_advanced_model()->ImpactIonization && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          // consider impact-ionization
          PetscScalar v = std::max(0.0, partial_volume);
          PetscScalar IIn,IIp,GIIn,GIIp;
          PetscScalar Eg = 0.5* ( n1_data->Eg() + n2_data->Eg() );

          VectorValue<Real> ev = (elem->point(edge_nodes.second) - elem->point(edge_nodes.first));
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
              IIn = mt->gen->ElecGenRate(T,fabs((V2-V1)/length),Eg);
              IIp = mt->gen->HoleGenRate(T,fabs((V2-V1)/length),Eg);
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

          if( fvm_n1->on_processor() )
          {
            // continuity equation
            iii.push_back( n1_global_offset + 1);
            ii.push_back ( (riin1*GIIn+riip1*GIIp)*truncated_partial_volume );

            iii.push_back( n1_global_offset + 2);
            ii.push_back ( (riin1*GIIn+riip1*GIIp)*truncated_partial_volume );

            n1_data->ImpactIonization() += (riin1*GIIn+riip1*GIIp)*truncated_partial_volume/fvm_n1->volume();
          }

          if( fvm_n2->on_processor() )
          {
            // continuity equation
            iii.push_back( n2_global_offset + 1);
            ii.push_back ( (riin2*GIIn+riip2*GIIp)*truncated_partial_volume );

            iii.push_back( n2_global_offset + 2);
            ii.push_back ( (riin2*GIIn+riip2*GIIp)*truncated_partial_volume );

            n2_data->ImpactIonization() += (riin2*GIIn+riip2*GIIp)*truncated_partial_volume/fvm_n2->volume();
          }
        }
      }
    }
    // the average cell electron/hole current density vector
    elem_data->Jn() = -elem->reconstruct_vector(Jn_edge_cell);
    elem_data->Jp() =  elem->reconstruct_vector(Jp_edge_cell);

  }

  // add into petsc vector, we should prevent zero length vector add here.
  if(iflux.size())    VecSetValues(f, iflux.size(), &iflux[0], &flux[0], ADD_VALUES);
  if(ibbt.size())     VecSetValues(f, ibbt.size(), &ibbt[0], &bbt[0], ADD_VALUES);
  if(iii.size())      VecSetValues(f, iii.size(), &iii[0], &ii[0], ADD_VALUES);

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

    PetscScalar V   =  x[local_offset+0];                         // electrostatic potential
    PetscScalar n   =  x[local_offset+1];                         // electron density
    PetscScalar p   =  x[local_offset+2];                         // hole density

    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);      // map this node and its data to material database

    PetscScalar R   = - mt->band->Recomb(p, n, T)*fvm_node->volume();         // the recombination term

    PetscScalar doping = node_data->Net_doping();
    if(get_advanced_model()->IncompleteIonization)
      doping = mt->band->Nd_II(n, T, get_advanced_model()->Fermi) - mt->band->Na_II(p, T, get_advanced_model()->Fermi);
    PetscScalar rho = e*( doping + p - n)*fvm_node->volume(); // the charge density


    // consider carrier generation
    PetscScalar Field_G = node_data->Field_G()*fvm_node->volume();

    isource.push_back(global_offset+0);                                // save index in the buffer
    isource.push_back(global_offset+1);
    isource.push_back(global_offset+2);
    source.push_back( rho );                                                       // save value in the buffer
    source.push_back( R + Field_G + node_data->EIn());
    source.push_back( R + Field_G + node_data->HIn());


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
        isource.push_back(fvm_node->global_offset());
        source.push_back(TrappedC);
      }

      // calculate the rates of electron and hole capture
      PetscScalar TrapElec = mt->trap->ElectronTrapRate(true,n,ni,T) * fvm_node->volume();
      PetscScalar TrapHole = mt->trap->HoleTrapRate    (true,p,ni,T) * fvm_node->volume();
      // contribution to the contribution to continuity equations
      if (TrapElec != 0)
      {
        isource.push_back(global_offset+1);
        // we lose carrier when electron get trapped, therefore negative contribution
        source.push_back(-TrapElec);
      }
      if (TrapHole != 0)
      {
        isource.push_back(global_offset+2);
        source.push_back(-TrapHole);
      }
    }
  }


  // add into petsc vector, we should prevent zero length vector add here.
  if(isource.size())  VecSetValues(f, isource.size(), &isource[0], &source[0], ADD_VALUES);


  // after the first scan, every nodes are updated.
  // however, boundary condition should be processed later.

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}






/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 * AD is fully used here
 */
void SemiconductorSimulationRegion::DDM1_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{

  //common used variable
  const PetscScalar T   = T_external();
  const PetscScalar Vt  = kb*T/e;
  bool  highfield_mob   = highfield_mobility() && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM;

  // precompute S-G current on each edge
  std::vector<AutoDScalar> Jn_edge_buffer;
  std::vector<AutoDScalar> Jp_edge_buffer;
  {
    Jn_edge_buffer.reserve(n_edge());
    Jp_edge_buffer.reserve(n_edge());

    //the indepedent variable number, 2 nodes * 3 variables per edge
    adtl::AutoDScalar::numdir = 6;

    //synchronize with material database
    mt->set_ad_num(adtl::AutoDScalar::numdir);

    // search all the edges of this region
    const_edge_iterator it = edges_begin();
    const_edge_iterator it_end = edges_end();
    for(; it!=it_end; ++it)
    {
      // fvm_node of node1
      const FVM_Node * fvm_n1 = (*it).first;
      // fvm_node of node2
      const FVM_Node * fvm_n2 = (*it).second;

      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data();
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data();

      const unsigned int n1_local_offset = fvm_n1->local_offset();
      const unsigned int n2_local_offset = fvm_n2->local_offset();

      const double length = fvm_n1->distance(fvm_n2);

      // build S-G current along edge


      //for node 1 of the edge
      mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

      AutoDScalar V1   =  x[n1_local_offset+0];   V1.setADValue(0, 1.0);               // electrostatic potential
      AutoDScalar n1   =  x[n1_local_offset+1];   n1.setADValue(1, 1.0);               // electron density
      AutoDScalar p1   =  x[n1_local_offset+2];   p1.setADValue(2, 1.0);               // hole density

      // NOTE: Here Ec1, Ev1 are not the conduction/valence band energy.
      // They are here for the calculation of effective driving field for electrons and holes
      // They differ from the conduction/valence band energy by the term with kb*T*log(Nc or Nv), which
      // takes care of the change effective DOS.
      // Ec/Ev should not be used except when its difference between two nodes.
      // The same comment applies to Ec2/Ev2.
      AutoDScalar Ec1 =  -(e*V1 + n1_data->affinity() - n1_data->dEcStrain() + mt->band->EgNarrowToEc(p1, n1, T) + kb*T*log(n1_data->Nc()));
      AutoDScalar Ev1 =  -(e*V1 + n1_data->affinity() - n1_data->dEvStrain() - mt->band->EgNarrowToEv(p1, n1, T) - kb*T*log(n1_data->Nv()) + mt->band->Eg(T));
      if(get_advanced_model()->Fermi)
      {
        Ec1 = Ec1 - kb*T*log(gamma_f(fabs(n1)/n1_data->Nc()));
        Ev1 = Ev1 + kb*T*log(gamma_f(fabs(p1)/n1_data->Nv()));
      }
      const PetscScalar eps1 =  n1_data->eps();

      //for node 2 of the edge
      mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

      AutoDScalar V2   =  x[n2_local_offset+0];   V2.setADValue(3, 1.0);                // electrostatic potential
      AutoDScalar n2   =  x[n2_local_offset+1];   n2.setADValue(4, 1.0);                // electron density
      AutoDScalar p2   =  x[n2_local_offset+2];   p2.setADValue(5, 1.0);                // hole density

      AutoDScalar Ec2 =  -(e*V2 + n2_data->affinity() - n2_data->dEcStrain() + mt->band->EgNarrowToEc(p2, n2, T) + kb*T*log(n2_data->Nc()));
      AutoDScalar Ev2 =  -(e*V2 + n2_data->affinity() - n2_data->dEvStrain() - mt->band->EgNarrowToEv(p2, n2, T) - kb*T*log(n2_data->Nv()) + mt->band->Eg(T));
      if(get_advanced_model()->Fermi)
      {
        Ec2 = Ec2 - kb*T*log(gamma_f(fabs(n2)/n2_data->Nc()));
        Ev2 = Ev2 + kb*T*log(gamma_f(fabs(p2)/n2_data->Nv()));
      }
      const PetscScalar eps2 =  n2_data->eps();

      // S-G current along the edge
      Jn_edge_buffer.push_back( In_dd(Vt,(Ec2-Ec1)/e,n1,n2,length) );
      Jp_edge_buffer.push_back( Ip_dd(Vt,(Ev2-Ev1)/e,p1,p2,length) );

      // poisson's equation

      const PetscScalar eps = 0.5*(eps1+eps2);
      AutoDScalar f_phi =  eps*fvm_n1->cv_surface_area(fvm_n2)*(V2 - V1)/length ;

      PetscInt row[2],col[2];
      row[0] = col[0] = fvm_n1->global_offset();
      row[1] = col[1] = fvm_n2->global_offset();

      // ignore thoese ghost nodes
      if( fvm_n1->on_processor() )
      {
        jac->add( row[0],  col[0],  f_phi.getADValue(0) );
        jac->add( row[0],  col[1],  f_phi.getADValue(3) );
      }

      if( fvm_n2->on_processor() )
      {
        jac->add( row[1],  col[0],  -f_phi.getADValue(0) );
        jac->add( row[1],  col[1],  -f_phi.getADValue(3) );
      }

    }
  }

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
        (SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationBoundary && is_elem_touch_boundary(elem)) ;


    //the indepedent variable number, 3*n_nodes
    adtl::AutoDScalar::numdir = 3*elem->n_nodes();

    //synchronize with material database
    mt->set_ad_num(adtl::AutoDScalar::numdir);

    // indicate the column position of the variables in the matrix
    std::vector<PetscInt> cell_col;
    cell_col.reserve(4*elem->n_nodes());
    for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
    {
      const FVM_Node * fvm_node = elem->get_fvm_node(nd);
      const unsigned int global_offset = fvm_node->global_offset();
      cell_col.push_back(global_offset+0);
      cell_col.push_back(global_offset+1);
      cell_col.push_back(global_offset+2);
    }


    // first, we build the gradient of psi and fermi potential in this cell.
    VectorValue<AutoDScalar> E;
    VectorValue<AutoDScalar> Jnv;
    VectorValue<AutoDScalar> Jpv;

    // E field parallel to current flow
    AutoDScalar Epn(0);
    AutoDScalar Epp(0);

    // E field vertical to current flow
    AutoDScalar Etn(0);
    AutoDScalar Etp(0);

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

        if(get_advanced_model()->HighFieldMobilitySelfConsistently)
        {
          double truc = get_advanced_model()->QuasiFermiCarrierTruc;
          // use values in the current iteration
          V  =  x[fvm_node->local_offset()+0];   V.setADValue(3*nd+0, 1.0);
          n  =  std::max(x[fvm_node->local_offset()+1], truc*fvm_node_data->ni());
          p  =  std::max(x[fvm_node->local_offset()+2], truc*fvm_node_data->ni());

          if(x[fvm_node->local_offset()+1] > truc*fvm_node_data->ni())
            n.setADValue(3*nd+1, 1.0);

          if(x[fvm_node->local_offset()+2] > truc*fvm_node_data->ni())
            p.setADValue(3*nd+2, 1.0);
        }
        else
        {
          // n and p use previous solution value
          V  =  x[fvm_node->local_offset()+0];   V.setADValue(3*nd+0, 1.0);
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
        genius_assert(!sides.empty());
        const Elem * elem_insul=elem->neighbor(sides[0]);
        SimulationRegion * region_insul=regions[0];

        std::vector<AutoDScalar> psi_vertex_neighbor;
        for(unsigned int nd=0; nd<elem_insul->n_nodes(); ++nd)
        {
          const FVM_Node * fvm_node_neighbor = elem_insul->get_fvm_node(nd);
          AutoDScalar V_neighbor = x[fvm_node_neighbor->local_offset()+0];
          V_neighbor.setADValue(3*elem->n_nodes()+nd, 1.0);
          psi_vertex_neighbor.push_back(V_neighbor);
        }

        VectorValue<AutoDScalar> E_insul    = - elem_insul->gradient(psi_vertex_neighbor);

        // we need more AD variable
        adtl::AutoDScalar::numdir += 1*elem_insul->n_nodes();
        mt->set_ad_num(adtl::AutoDScalar::numdir);
        for(unsigned int nd=0; nd<elem_insul->n_nodes(); ++nd)
        {
          const FVM_Node * fvm_node = elem_insul->get_fvm_node(nd);
          cell_col.push_back( fvm_node->global_offset()+0 );
        }

        // interface normal, point to semiconductor side
        Point _norm = - elem->outside_unit_normal(sides[0]);
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
        //Epn = adtl::fmax(E_eff_p.dot(Jnv.unit(true)), 0.0);
        //Epp = adtl::fmax(E_eff_p.dot(Jpv.unit(true)), 0.0);
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

      const unsigned int edge_index = this->elem_edge_index(elem, ne);

      // the length of this edge
      const double length = elem->edge_length(ne);

      const FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes.first);   // fvm_node of node1
      const FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes.second);  // fvm_node of node2

      double partial_area = elem->partial_area_with_edge(ne);        // partial area associated with this edge
      double partial_volume = elem->partial_volume_with_edge(ne);    // partial volume associated with this edge
      double truncated_partial_area =  partial_area;
      double truncated_partial_volume =  partial_volume;
      if(truncation)
      {
        // use truncated partial area to avoid negative area due to bad mesh elem
        truncated_partial_area =  this->truncated_partial_area(elem, ne);
        truncated_partial_volume =  elem->partial_volume_with_edge_truncated(ne);
      }

      bool inverse = fvm_n1->root_node()->id() > fvm_n2->root_node()->id();       // find the correct order


      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data();
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data();

      const unsigned int n1_local_offset = fvm_n1->local_offset();
      const unsigned int n2_local_offset = fvm_n2->local_offset();
      const unsigned int n1_global_offset = fvm_n1->global_offset();
      const unsigned int n2_global_offset = fvm_n2->global_offset();

      // the row position of variables in the matrix
      PetscInt row[6];
      for(int i=0; i<3; ++i) row[i]   = n1_global_offset+i;
      for(int i=0; i<3; ++i) row[i+3] = n2_global_offset+i;

      // here we use AD again. Can we hand write it for more efficient?
      {
        AutoDScalar V1(x[n1_local_offset+0]);       V1.setADValue(3*edge_nodes.first+0, 1.0);           // electrostatic potential
        AutoDScalar n1(x[n1_local_offset+1]);       n1.setADValue(3*edge_nodes.first+1, 1.0);           // electron density
        AutoDScalar p1(x[n1_local_offset+2]);       p1.setADValue(3*edge_nodes.first+2, 1.0);           // hole density

        AutoDScalar V2(x[n2_local_offset+0]);       V2.setADValue(3*edge_nodes.second+0, 1.0);          // electrostatic potential
        AutoDScalar n2(x[n2_local_offset+1]);       n2.setADValue(3*edge_nodes.second+1, 1.0);          // electron density
        AutoDScalar p2(x[n2_local_offset+2]);       p2.setADValue(3*edge_nodes.second+2, 1.0);          // hole density

        AutoDScalar mun1;   // electron mobility for node 1 of the edge
        AutoDScalar mup1;   // hole mobility for node 1 of the edge
        AutoDScalar mun2;   // electron mobility  for node 2 of the edge
        AutoDScalar mup2;   // hole mobility for node 2 of the edge

        if(highfield_mob)
        {


          if(get_advanced_model()->ESurface && insulator_interface_elem)
          {
            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T, Epn, Etn, T);
            mup1 = mt->mob->HoleMob(p1, n1, T, Epp, Etp, T);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T, Epn, Etn, T);
            mup2 = mt->mob->HoleMob(p2, n2, T, Epp, Etp, T);
          }
          else
          {
            // high field mobility
            if (get_advanced_model()->Mob_Force == ModelSpecify::ESimple)
            {
              Point _dir = (*fvm_n1->root_node() - *fvm_n2->root_node()).unit();
              VectorValue<AutoDScalar> dir(_dir(0), _dir(1), _dir(2));
              AutoDScalar Epn = adtl::fabs((V1-V2)/length);
              AutoDScalar Epp = adtl::fabs((V2-V1)/length);
              //AutoDScalar Epn = adtl::fmax(0.0, (inverse ? -1 : 1) *(V1-V2)/length);
              //AutoDScalar Epp = adtl::fmax(0.0, (inverse ? -1 : 1) *(V2-V1)/length);
              AutoDScalar Et = 0;//
              if(mos_channel_elem)
                Et = (E - (E*dir)*dir).size();

              mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
              mun1 = mt->mob->ElecMob(p1, n1, T, Epn, Et, T);
              mup1 = mt->mob->HoleMob(p1, n1, T, Epp, Et, T);

              mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
              mun2 = mt->mob->ElecMob(p2, n2, T, Epn, Et, T);
              mup2 = mt->mob->HoleMob(p2, n2, T, Epp, Et, T);
            }
            else // ModelSpecify::EJ || ModelSpecify::EQF
            {
              mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
              mun1 = mt->mob->ElecMob(p1, n1, T, Epn, Etn, T);
              mup1 = mt->mob->HoleMob(p1, n1, T, Epp, Etp, T);

              mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
              mun2 = mt->mob->ElecMob(p2, n2, T, Epn, Etn, T);
              mup2 = mt->mob->HoleMob(p2, n2, T, Epp, Etp, T);
            }
          }
        }
        else // low field mobility
        {
          mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
          mun1 = mt->mob->ElecMob(p1, n1, T, 0, 0, T);
          mup1 = mt->mob->HoleMob(p1, n1, T, 0, 0, T);

          mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
          mun2 = mt->mob->ElecMob(p2, n2, T, 0, 0, T);
          mup2 = mt->mob->HoleMob(p2, n2, T, 0, 0, T);
        }


        AutoDScalar mun = 0.5*(mun1+mun2);  // the electron mobility at the mid point of the edge, use linear interpolation
        AutoDScalar mup = 0.5*(mup1+mup2);  // the hole mobility at the mid point of the edge, use linear interpolation

        // S-G current along the edge
        const AutoDScalar & Jn_edge = Jn_edge_buffer[edge_index];
        const AutoDScalar & Jp_edge = Jp_edge_buffer[edge_index];

        // shift AD value since they have different location
        unsigned int order[6];
        if(inverse)
        {
          order[0]= 3*edge_nodes.second+0;
          order[1]= 3*edge_nodes.second+1;
          order[2]= 3*edge_nodes.second+2;
          order[3]= 3*edge_nodes.first+0;
          order[4]= 3*edge_nodes.first+1;
          order[5]= 3*edge_nodes.first+2;
        }
        else
        {
          order[0]= 3*edge_nodes.first+0;
          order[1]= 3*edge_nodes.first+1;
          order[2]= 3*edge_nodes.first+2;
          order[3]= 3*edge_nodes.second+0;
          order[4]= 3*edge_nodes.second+1;
          order[5]= 3*edge_nodes.second+2;
        }

        AutoDScalar Jn = (inverse ? -1.0 : 1.0)*mun*AutoDScalar(Jn_edge, order, 6);
        AutoDScalar Jp = (inverse ? -1.0 : 1.0)*mup*AutoDScalar(Jp_edge, order, 6);

        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->on_processor() )
        {
          // flux on edge
          AutoDScalar f_Jn  =  Jn*truncated_partial_area ;
          AutoDScalar f_Jp  = -Jp*truncated_partial_area;
          // general coding always has some overkill... bypass it.
          jac->add_row(  row[1],  cell_col.size(),  &cell_col[0],  f_Jn.getADValue() );
          jac->add_row(  row[2],  cell_col.size(),  &cell_col[0],  f_Jp.getADValue() );
        }

        if( fvm_n2->on_processor() )
        {
          // flux on edge
          AutoDScalar f_Jn  = -Jn*truncated_partial_area ;
          AutoDScalar f_Jp  =  Jp*truncated_partial_area;
          jac->add_row(  row[4],  cell_col.size(),  &cell_col[0],  f_Jn.getADValue() );
          jac->add_row(  row[5],  cell_col.size(),  &cell_col[0],  f_Jp.getADValue() );
        }

        // BandBandTunneling && ImpactIonization

        if (get_advanced_model()->BandBandTunneling && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          AutoDScalar GBTBT1 = mt->band->BB_Tunneling(T, E.size());
          AutoDScalar GBTBT2 = mt->band->BB_Tunneling(T, E.size());

          if( fvm_n1->on_processor() )
          {
            // continuity equation
            AutoDScalar continuity = 0.5*GBTBT1*truncated_partial_volume;
            jac->add_row(  row[1],  cell_col.size(),  &cell_col[0],  continuity.getADValue() );
            jac->add_row(  row[2],  cell_col.size(),  &cell_col[0],  continuity.getADValue() );
          }

          if( fvm_n2->on_processor() )
          {
            // continuity equation
            AutoDScalar continuity = 0.5*GBTBT2*truncated_partial_volume;
            jac->add_row(  row[4],  cell_col.size(),  &cell_col[0],  continuity.getADValue() );
            jac->add_row(  row[5],  cell_col.size(),  &cell_col[0],  continuity.getADValue() );
          }
        }

        if (get_advanced_model()->ImpactIonization && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          // consider impact-ionization
           AutoDScalar IIn,IIp,GIIn,GIIp;
          PetscScalar Eg = 0.5* ( n1_data->Eg() + n2_data->Eg() );

          // FIXME should use weighted carrier temperature.

          VectorValue<Real> ev0 = (elem->point(edge_nodes.second) - elem->point(edge_nodes.first));
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
              IIn = mt->gen->ElecGenRate(T,fabs((V2-V1)/length),Eg);
              IIp = mt->gen->HoleGenRate(T,fabs((V2-V1)/length),Eg);
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

          if( fvm_n1->on_processor() )
          {
            // continuity equation
            AutoDScalar electron_continuity = (riin1*GIIn+riip1*GIIp)*truncated_partial_volume ;
            AutoDScalar hole_continuity     = (riin1*GIIn+riip1*GIIp)*truncated_partial_volume ;
            jac->add_row(  row[1],  cell_col.size(),  &cell_col[0],  electron_continuity.getADValue() );
            jac->add_row(  row[2],  cell_col.size(),  &cell_col[0],  hole_continuity.getADValue() );
          }

          if( fvm_n2->on_processor() )
          {
            // continuity equation
            AutoDScalar electron_continuity = (riin2*GIIn+riip2*GIIp)*truncated_partial_volume ;
            AutoDScalar hole_continuity     = (riin2*GIIn+riip2*GIIp)*truncated_partial_volume ;
            jac->add_row(  row[4],  cell_col.size(),  &cell_col[0],  electron_continuity.getADValue() );
            jac->add_row(  row[5],  cell_col.size(),  &cell_col[0],  hole_continuity.getADValue() );
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

  //the indepedent variable number, 3 for each node
  adtl::AutoDScalar::numdir = 3;

  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;

    const unsigned int local_offset = fvm_node->local_offset();
    const unsigned int global_offset = fvm_node->global_offset();
    const FVM_NodeData * node_data = fvm_node->node_data();

    PetscInt index[3] = {global_offset+0, global_offset+1, global_offset+2};

    AutoDScalar V(x[local_offset+0]);   V.setADValue(0, 1.0);              // psi
    AutoDScalar n(x[local_offset+1]);   n.setADValue(1, 1.0);              // electron density
    AutoDScalar p(x[local_offset+2]);   p.setADValue(2, 1.0);              // hole density

    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);                   // map this node and its data to material database

    AutoDScalar R   = - mt->band->Recomb(p, n, T)*fvm_node->volume();                      // the recombination term

    AutoDScalar doping = node_data->Net_doping();
    if(get_advanced_model()->IncompleteIonization)
      doping = mt->band->Nd_II(n, T, get_advanced_model()->Fermi) - mt->band->Na_II(p, T, get_advanced_model()->Fermi);
    AutoDScalar rho = e*( doping + p - n)*fvm_node->volume(); // the charge density



    // ADD to Jacobian matrix,
    jac->add_row(  index[0],  3,  &index[0],  rho.getADValue() );
    jac->add_row(  index[1],  3,  &index[0],  R.getADValue() );
    jac->add_row(  index[2],  3,  &index[0],  R.getADValue() );

    if (get_advanced_model()->Trap)
    {
      AutoDScalar ni = mt->band->nie(p, n, T);
      mt->trap->Calculate(true,p,n,ni,T);

      AutoDScalar TrappedC = mt->trap->ChargeAD(true) * fvm_node->volume();
      jac->add_row(  index[0],  3,  &index[0],  TrappedC.getADValue() );

      AutoDScalar GElec = - mt->trap->ElectronTrapRate(true,n,ni,T) * fvm_node->volume();
      AutoDScalar GHole = - mt->trap->HoleTrapRate    (true,p,ni,T) * fvm_node->volume();

      jac->add_row(  index[1],  3,  &index[0],  GElec.getADValue() );
      jac->add_row(  index[2],  3,  &index[0],  GHole.getADValue() );
    }
  }


  // boundary condition should be processed later!

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}




void SemiconductorSimulationRegion::DDM1_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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
  iy.reserve(2*this->n_node());
  y.reserve(2*this->n_node());

  const double r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);

  // process node related terms
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();

    PetscScalar n   =  x[fvm_node->local_offset()+1];                         // electron density
    PetscScalar p   =  x[fvm_node->local_offset()+2];                         // hole density

    // process \partial t

    iy.push_back(fvm_node->global_offset()+1);                                // save index in the buffer
    iy.push_back(fvm_node->global_offset()+2);

    //second order
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false)
    {
      PetscScalar Tn_BDF2 = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                       / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
      PetscScalar Tp_BDF2 = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*node_data->p() + (1-r)/r*node_data->p_last())
                       / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();

      y.push_back( Tn_BDF2 );
      y.push_back( Tp_BDF2 );
    }
    else //first order
    {
      PetscScalar Tn_BDF1 = -(n - node_data->n())/SolverSpecify::dt*fvm_node->volume();
      PetscScalar Tp_BDF1 = -(p - node_data->p())/SolverSpecify::dt*fvm_node->volume();
      y.push_back( Tn_BDF1 );
      y.push_back( Tp_BDF1 );
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




void SemiconductorSimulationRegion::DDM1_Time_Dependent_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{

  //the indepedent variable number, 1 for each node
  adtl::AutoDScalar::numdir = 1;
  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  const double r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);

  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();

    PetscInt index[2] = {fvm_node->global_offset()+1, fvm_node->global_offset()+2};

    AutoDScalar n(x[fvm_node->local_offset()+1]);   n.setADValue(0, 1.0);              // electron density
    AutoDScalar p(x[fvm_node->local_offset()+2]);   p.setADValue(0, 1.0);              // hole density

    // process \partial t

    //second order
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false)
    {
      AutoDScalar Tn_BDF2 = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                       / (SolverSpecify::dt_last+SolverSpecify::dt)*fvm_node->volume();
      AutoDScalar Tp_BDF2 = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*node_data->p() + (1-r)/r*node_data->p_last())
                       / (SolverSpecify::dt_last+SolverSpecify::dt)*fvm_node->volume();
      // ADD to Jacobian matrix
      jac->add( index[0],  index[0],  Tn_BDF2.getADValue(0) );
      jac->add( index[1],  index[1],  Tp_BDF2.getADValue(0) );
    }
    else //first order
    {
      AutoDScalar Tn_BDF1 = -(n - node_data->n())/SolverSpecify::dt*fvm_node->volume();
      AutoDScalar Tp_BDF1 = -(p - node_data->p())/SolverSpecify::dt*fvm_node->volume();
      // ADD to Jacobian matrix
      jac->add( index[0],  index[0],  Tn_BDF1.getADValue(0) );
      jac->add( index[1],  index[1],  Tp_BDF1.getADValue(0) );
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}


void SemiconductorSimulationRegion::DDM1_Pseudo_Time_Step_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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


void SemiconductorSimulationRegion::DDM1_Pseudo_Time_Step_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
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
    jac->add( global_offset,  global_offset,  f_V.getADValue(0) );
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
    jac->add( global_offset+1,  global_offset+1,  f_n.getADValue(0) );

    AutoDScalar f_p = -(p-node_data->p())/SolverSpecify::PseudoTimeStepCarrier*fvm_node->volume();
    jac->add( global_offset+2,  global_offset+2,  f_p.getADValue(0) );
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}


int SemiconductorSimulationRegion::DDM1_Pseudo_Time_Step_Convergence_Test(PetscScalar * x)
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


void SemiconductorSimulationRegion::DDM1_Update_Solution(PetscScalar *lxx)
{
  //common used variable
  const PetscScalar T   = T_external();
  const PetscScalar Vt  = kb*T/e;

  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;

    PetscScalar V          = lxx[fvm_node->local_offset()+0];
    PetscScalar n          = lxx[fvm_node->local_offset()+1];
    PetscScalar p          = lxx[fvm_node->local_offset()+2];

    FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data!=NULL);
    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

    //update psi
    node_data->psi_old()  =  node_data->psi_last();
    node_data->psi_last() =  node_data->psi();
    node_data->psi()      =  V;

    // electron density
    node_data->n_last()   = node_data->n();
    node_data->n()        =  n;

    // hole density
    node_data->p_last()   = node_data->p();
    node_data->p()        =  p;


    node_data->Eg() = mt->band->Eg(T) - mt->band->EgNarrow(p, n, T);
    node_data->Ec() = -(e*V + node_data->affinity() + mt->band->EgNarrowToEc(p, n, T));
    node_data->Ev() = -(e*V + node_data->affinity() + mt->band->Eg(T) - mt->band->EgNarrowToEv(p, n, T) );

    if(get_advanced_model()->Fermi)
    {
      node_data->qFn() = node_data->Ec() + inv_fermi_half(fabs(n/node_data->Nc()))*kb*T;
      node_data->qFp() = node_data->Ev() - inv_fermi_half(fabs(p/node_data->Nv()))*kb*T;
    }
    else
    {
      node_data->qFn() = node_data->Ec() + log(fabs(n/node_data->Nc()))*kb*T;
      node_data->qFp() = node_data->Ev() - log(fabs(p/node_data->Nv()))*kb*T;
    }

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
    psi_vertex.reserve(elem->n_nodes());

    for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
    {
      const FVM_Node * fvm_node = elem->get_fvm_node(nd);
      const FVM_NodeData * fvm_node_data = fvm_node->node_data();
      psi_vertex.push_back  ( fvm_node_data->psi() );
    }
    // compute the gradient in the cell
    elem_data->E()  = - elem->gradient(psi_vertex);  // E = - grad(psi)
  }

  // calculate mobility on node
  Mob_Evaluation();

}






