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

#include "jflux1.h"
#include "jflux2.h"
#include "jflux3.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::um;
using PhysicalUnit::cm;


/*-------------------------------------------------------------------------------------
 * get nodal variable number, which depedent on EBM Level
 */
unsigned int SemiconductorSimulationRegion::ebm_n_variables() const
{
  switch(get_advanced_model()->EB_Level)
  {
      case ModelSpecify::NONE : return 3; // 3 basic dofs, psi, n and p
      case ModelSpecify::Tl   :
      case ModelSpecify::Tn   :
      case ModelSpecify::Tp   : return 4; // 4 dofs for one of the extra Tl, Tn and Tp equation
      case ModelSpecify::TnTp :
      case ModelSpecify::TnTl :
      case ModelSpecify::TpTl : return 5; // 5 dofs for Tn/Tp, Tn/Tl or Tp/Tl pair
      case ModelSpecify::ALL  : return 6; // 6 dofs for all of the extra Tl, Tn and Tp equations
  }
  return 0; // prevent compiler warning.
}


/*-------------------------------------------------------------------------------------
 * get offset of nodal variable, which depedent on EBM Level
 */
unsigned int SemiconductorSimulationRegion::ebm_variable_offset(SolutionVariable var) const
{
  switch(var)
  {
      case POTENTIAL     : return 0; // psi is always at offset 0
      case ELECTRON      : return 1; // followed by electron density
      case HOLE          : return 2; // and hole density
      case TEMPERATURE   :           // if lattice temperature is required?
      switch(get_advanced_model()->EB_Level)
      {
          case ModelSpecify::NONE :
          case ModelSpecify::Tn   :
          case ModelSpecify::Tp   :
          case ModelSpecify::TnTp : return invalid_uint; // no lattice temperature

          case ModelSpecify::Tl   :
          case ModelSpecify::TnTl :
          case ModelSpecify::TpTl :
          case ModelSpecify::ALL  : return 3;            // lattice temperature at offset 3
      }
      case E_TEMP :                  // if electron temperature is required?
      switch(get_advanced_model()->EB_Level)
      {
          case ModelSpecify::NONE :
          case ModelSpecify::Tl   :
          case ModelSpecify::Tp   :
          case ModelSpecify::TpTl : return invalid_uint; // no electron temperature
          case ModelSpecify::TnTp :
          case ModelSpecify::Tn   : return 3;            // electron temperature at offset 3
          case ModelSpecify::TnTl :
          case ModelSpecify::ALL  : return 4;            // electron temperature at offset 4
      }
      case H_TEMP :                  // if hole temperature is required?
      switch(get_advanced_model()->EB_Level)
      {
          case ModelSpecify::NONE :
          case ModelSpecify::Tl   :
          case ModelSpecify::Tn   :
          case ModelSpecify::TnTl : return invalid_uint; // no hole temperature

          case ModelSpecify::Tp   : return 3;            // hole temperature at offset 3
          case ModelSpecify::TnTp :
          case ModelSpecify::TpTl : return 4;            // hole temperature at offset 4
          case ModelSpecify::ALL  : return 5;            // hole temperature at offset 5
      }
      default : return invalid_uint;
  }
}


/*-------------------------------------------------------------------------------------
 * filling solution data from FVM_NodeData into petsc vector of L3 EBM.
 */
void SemiconductorSimulationRegion::EBM3_Fill_Value(Vec x, Vec L)
{
  // find the node variable offset
  unsigned int node_psi_offset = ebm_variable_offset(POTENTIAL);
  unsigned int node_n_offset   = ebm_variable_offset(ELECTRON);
  unsigned int node_p_offset   = ebm_variable_offset(HOLE);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);
  unsigned int node_Tn_offset  = ebm_variable_offset(E_TEMP);
  unsigned int node_Tp_offset  = ebm_variable_offset(H_TEMP);

  // data buffer
  std::vector<int> ix;
  std::vector<PetscScalar> y;
  std::vector<PetscScalar> s;

  // reserve menory for data buffer
  ix.reserve(6*this->n_node());
  y.reserve(6*this->n_node());
  s.reserve(6*this->n_node());

  // for all the on processor node, insert value to petsc vector
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();

    mt->mapping(fvm_node->root_node(), node_data, 0.0);

    // the first variable, psi
    ix.push_back(fvm_node->global_offset()+node_psi_offset);
    y.push_back(node_data->psi());
    s.push_back(1.0/(node_data->eps()*fvm_node->volume()));

    // the second variable, n
    ix.push_back(fvm_node->global_offset()+node_n_offset);
    y.push_back(node_data->n());
    s.push_back(1.0/fvm_node->volume());

    // the third variable, p
    ix.push_back(fvm_node->global_offset()+node_p_offset);
    y.push_back(node_data->p());
    s.push_back(1.0/fvm_node->volume());

    // for extra Tl temperature equations
    if(get_advanced_model()->enable_Tl())
    {
      ix.push_back(fvm_node->global_offset()+node_Tl_offset);
      y.push_back(node_data->T());
      s.push_back(1.0/fvm_node->volume());
    }

    // for extra Tn temperature equations, use n*Tn as indepedent variable
    if(get_advanced_model()->enable_Tn())
    {
      ix.push_back(fvm_node->global_offset()+node_Tn_offset);
      y.push_back(node_data->n()*node_data->Tn());
      s.push_back(1.0/fvm_node->volume());
    }

    // for extra Tp temperature equations, use p*Tp as indepedent variable
    if(get_advanced_model()->enable_Tp())
    {
      ix.push_back(fvm_node->global_offset()+node_Tp_offset);
      y.push_back(node_data->p()*node_data->Tp());
      s.push_back(1.0/fvm_node->volume());
    }

  }

  // call petsc VecSetValues routine to insert bufferred value
  if( ix.size() )
  {
    VecSetValues(x, ix.size(), &ix[0], &y[0], INSERT_VALUES) ;
    VecSetValues(L, ix.size(), &ix[0], &s[0], INSERT_VALUES) ;
  }
}




///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////



/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void SemiconductorSimulationRegion::EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // find the node variable offset
  unsigned int node_psi_offset = ebm_variable_offset(POTENTIAL);
  unsigned int node_n_offset   = ebm_variable_offset(ELECTRON);
  unsigned int node_p_offset   = ebm_variable_offset(HOLE);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);
  unsigned int node_Tn_offset  = ebm_variable_offset(E_TEMP);
  unsigned int node_Tp_offset  = ebm_variable_offset(H_TEMP);

  // set current and joule heating level
  unsigned int Jn_level = get_advanced_model()->Jn_level();
  unsigned int Jp_level = get_advanced_model()->Jp_level();
  unsigned int Hn_level = get_advanced_model()->Hn_level();
  unsigned int Hp_level = get_advanced_model()->Hp_level();

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

  iy.reserve(6*n_node());
  y.reserve(6*n_node());

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
          // use values in the current iteration
          V  =  x[fvm_node->local_offset() + node_psi_offset];
          n  =  x[fvm_node->local_offset()+node_n_offset] + fvm_node_data->ni()*1e-2;
          p  =  x[fvm_node->local_offset()+node_p_offset] + fvm_node_data->ni()*1e-2;
        }
        else
        {
          // use previous solution value
          V  =  x[fvm_node->local_offset() + node_psi_offset];
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

      // the length of this edge
      const double length = elem->edge_length(ne);

      // fvm_node of node1
      const FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes.first);
      // fvm_node of node2
      const FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes.second);

      const FVM_NodeData * n1_data = fvm_n1->node_data();  genius_assert(n1_data);            // fvm_node_data of node1
      const FVM_NodeData * n2_data = fvm_n2->node_data();  genius_assert(n2_data);            // fvm_node_data of node2

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

      unsigned int n1_local_offset = fvm_n1->local_offset();
      unsigned int n2_local_offset = fvm_n2->local_offset();

      // build governing equation of EBM
      {

        //for node 1 of the edge
        mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

        PetscScalar V1   =  x[n1_local_offset + node_psi_offset];                // electrostatic potential
        PetscScalar n1   =  x[n1_local_offset + node_n_offset];                  // electron density
        PetscScalar p1   =  x[n1_local_offset + node_p_offset];                  // hole density

        PetscScalar T1 = T_external();
        PetscScalar Tn1= T_external();
        PetscScalar Tp1= T_external();

        // lattice temperature if required
        if(get_advanced_model()->enable_Tl())
          T1 =  x[n1_local_offset + node_Tl_offset];

        // electron temperature if required
        if(get_advanced_model()->enable_Tn())
          Tn1 = x[n1_local_offset + node_Tn_offset]/n1;

        // hole temperature if required
        if(get_advanced_model()->enable_Tp())
          Tp1 = x[n1_local_offset + node_Tp_offset]/p1;

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
        PetscScalar kap1 =  mt->thermal->HeatConduction(T1);
        PetscScalar Eg1= mt->band->Eg(T1);


        //for node 2 of the edge
        mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

        PetscScalar V2   =  x[n2_local_offset + node_psi_offset];                // electrostatic potential
        PetscScalar n2   =  x[n2_local_offset + node_n_offset];                  // electron density
        PetscScalar p2   =  x[n2_local_offset + node_p_offset];                  // hole density

        PetscScalar T2 = T_external();
        PetscScalar Tn2= T_external();
        PetscScalar Tp2= T_external();

        // lattice temperature if required
        if(get_advanced_model()->enable_Tl())
          T2 =  x[n2_local_offset + node_Tl_offset];

        // electron temperature if required
        if(get_advanced_model()->enable_Tn())
          Tn2 = x[n2_local_offset + node_Tn_offset]/n2;

        // hole temperature if required
        if(get_advanced_model()->enable_Tp())
          Tp2 = x[n2_local_offset + node_Tp_offset]/p2;

        PetscScalar Ec2 =  -(e*V2 + n2_data->affinity() + kb*T2*log(mt->band->nie(p2, n2, T2)));
        PetscScalar Ev2 =  -(e*V2 + n2_data->affinity() - kb*T2*log(mt->band->nie(p2, n2, T2)));
        if(get_advanced_model()->Fermi)
        {
          Ec2 = Ec2 - kb*T2*log(gamma_f(fabs(n2)/n2_data->Nc()));
          Ev2 = Ev2 + kb*T2*log(gamma_f(fabs(p2)/n2_data->Nv()));
        }

        PetscScalar eps2 =  n2_data->eps();                         // eps
        PetscScalar kap2 =  mt->thermal->HeatConduction(T2);
        PetscScalar Eg2= mt->band->Eg(T2);


        PetscScalar mun1;  // electron mobility
        PetscScalar mup1;  // hole mobility
        PetscScalar mun2;  // electron mobility
        PetscScalar mup2;  // hole mobility

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
            mun1 = mt->mob->ElecMob(p1, n1, T1, Ep, Et, Tn1);
            mup1 = mt->mob->HoleMob(p1, n1, T1, Ep, Et, Tp1);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T2, Ep, Et, Tn2);
            mup2 = mt->mob->HoleMob(p2, n2, T2, Ep, Et, Tp2);
          }
          else
          {
            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T1, Epn, Etn, Tn1);
            mup1 = mt->mob->HoleMob(p1, n1, T1, Epp, Etp, Tp1);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T2, Epn, Etn, Tn2);
            mup2 = mt->mob->HoleMob(p2, n2, T2, Epp, Etp, Tp2);
          }
        }
        else
        {
          mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
          mun1 = mt->mob->ElecMob(p1, n1, T1, 0, 0, Tn1);
          mup1 = mt->mob->HoleMob(p1, n1, T1, 0, 0, Tp1);

          mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
          mun2 = mt->mob->ElecMob(p2, n2, T2, 0, 0, Tn2);
          mup2 = mt->mob->HoleMob(p2, n2, T2, 0, 0, Tp2);
        }


        PetscScalar mun = 0.5*(mun1+mun2); // the electron mobility at the mid point of the edge, use linear interpolation
        PetscScalar mup = 0.5*(mup1+mup2); // the hole mobility at the mid point of the edge, use linear interpolation
        PetscScalar eps = 0.5*(eps1+eps2); // eps at mid point of the edge
        PetscScalar kap = 0.5*(kap1+kap2); // kapa at mid point of the edge

        // S-G current along the edge, call different SG scheme selected by EBM level
        PetscScalar Jn, Jp, Sn=0, Sp=0;
        switch(Jn_level)
        {
            case 1:
            Jn =  mun*In_dd(kb*T1/e, (Ec2-Ec1)/e, n1, n2, length);
            break;
            case 2:
            Jn =  mun*In_lt(kb, e, (Ec1-Ec2)/e, n1, n2, 0.5*(T1+T2), T2-T1, length);
            break;
            case 3:
            Jn =  mun*In_eb(kb, e, -Ec1/e, -Ec2/e, n1, n2, Tn1, Tn2, length);
            Sn =  mun*Sn_eb(kb, e, -Ec1/e, -Ec2/e, n1, n2, Tn1, Tn2, length);
            break;
        }


        switch(Jp_level)
        {
            case 1:
            Jp =  mup*Ip_dd(kb*T2/e, (Ev2-Ev1)/e, p1, p2, length);
            break;
            case 2:
            Jp =  mup*Ip_lt(kb, e, (Ev1-Ev2)/e, p1, p2, 0.5*(T1+T2), T2-T1, length);
            break;
            case 3:
            Jp =  mup*Ip_eb(kb, e, -Ev1/e, -Ev2/e, p1, p2, Tp1, Tp2, length);
            Sp =  mup*Sp_eb(kb, e, -Ev1/e, -Ev2/e, p1, p2, Tp1, Tp2, length);
            break;
        }


        // joule heating
        PetscScalar H=0, Hn=0, Hp=0;

        switch(Hn_level)
        {
            case  0  : break;                         // no heat equation
            case  1  : H += 0.5*(V1-V2)*(Jn); break;  // use JdotE as heating source to lattice
            case  2  : Hn = 0.5*(V1-V2)*(Jn); break;  // use JdotE as heating source to electron system
        }
        switch(Hp_level)
        {
            case  0  : break;                         // no heat equation
            case  1  : H += 0.5*(V1-V2)*(Jp); break;  // use JdotE as heating source to lattice
            case  2  : Hp = 0.5*(V1-V2)*(Jp); break;  // use JdotE as heating source to hole system
        }

        Jn_edge.push_back(Jn);
        Jp_edge.push_back(Jp);


        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {

          // poisson's equation
          iy.push_back( fvm_n1->global_offset() + node_psi_offset );
          y.push_back ( eps*(V2 - V1)/length*partial_area );

          // continuity equation of electron
          iy.push_back( fvm_n1->global_offset() + node_n_offset );
          y.push_back ( Jn*truncated_partial_area );

          // continuity equation of hole
          iy.push_back( fvm_n1->global_offset() + node_p_offset );
          y.push_back ( - Jp*truncated_partial_area );


          // heat transport equation if required
          if(get_advanced_model()->enable_Tl())
          {
            iy.push_back( fvm_n1->global_offset() + node_Tl_offset );
            y.push_back ( kap*(T2 - T1)/length*partial_area + H*truncated_partial_area);
          }


          // energy balance equation for electron if required
          if(get_advanced_model()->enable_Tn())
          {
            iy.push_back( fvm_n1->global_offset() + node_Tn_offset );
            y.push_back ( -Sn*truncated_partial_area + Hn*truncated_partial_area);
          }


          // energy balance equation for hole if required
          if(get_advanced_model()->enable_Tp())
          {
            iy.push_back( fvm_n1->global_offset() + node_Tp_offset );
            y.push_back ( -Sp*truncated_partial_area + Hp*truncated_partial_area);
          }

        }

        // for node 2.
        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {

          // poisson's equation
          iy.push_back( fvm_n2->global_offset() + node_psi_offset );
          y.push_back ( eps*(V1 - V2)/length*partial_area );

          // continuity equation of electron
          iy.push_back( fvm_n2->global_offset() + node_n_offset );
          y.push_back ( - Jn*truncated_partial_area );

          // continuity equation of hole
          iy.push_back( fvm_n2->global_offset() + node_p_offset );
          y.push_back ( Jp*truncated_partial_area );


          // heat transport equation if required
          if(get_advanced_model()->enable_Tl())
          {
            iy.push_back( fvm_n2->global_offset() + node_Tl_offset );
            y.push_back ( kap*(T1 - T2)/length*partial_area + H*truncated_partial_area);
          }


          // energy balance equation for electron if required
          if(get_advanced_model()->enable_Tn())
          {
            iy.push_back( fvm_n2->global_offset() + node_Tn_offset );
            y.push_back ( Sn*truncated_partial_area + Hn*truncated_partial_area);
          }


          // energy balance equation for hole if required
          if(get_advanced_model()->enable_Tp())
          {
            iy.push_back( fvm_n2->global_offset() + node_Tp_offset );
            y.push_back ( Sp*truncated_partial_area + Hp*truncated_partial_area);
          }

        }

        if (get_advanced_model()->BandBandTunneling && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          PetscScalar GBTBT1 = mt->band->BB_Tunneling(T1, E.size());
          PetscScalar GBTBT2 = mt->band->BB_Tunneling(T2, E.size());

          if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            iy.push_back( fvm_n1->global_offset() + node_n_offset );
            y.push_back ( 0.5*GBTBT1*truncated_partial_volume );

            iy.push_back( fvm_n1->global_offset() + node_p_offset );
            y.push_back ( 0.5*GBTBT1*truncated_partial_volume );
          }

          if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            iy.push_back( fvm_n2->global_offset() + node_n_offset );
            y.push_back ( 0.5*GBTBT2*truncated_partial_volume );

            iy.push_back( fvm_n2->global_offset() + node_p_offset );
            y.push_back ( 0.5*GBTBT2*truncated_partial_volume );
          }
        }


        if (get_advanced_model()->ImpactIonization && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          // consider impact-ionization

          PetscScalar IIn,IIp,GIIn,GIIp;
          PetscScalar T,Tn,Tp,Eg;

          // FIXME should use weighted carrier temperature.
          T  = std::min(T1,T2);
          Eg = 0.5* ( Eg1 + Eg2 );
          Tn = std::min(Tn1,Tn2);
          Tp = std::min(Tp1,Tp2);

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
              case ModelSpecify::TempII:
              IIn = mt->gen->ElecGenRateEBM (Tn,T,Eg);
              IIp = mt->gen->HoleGenRateEBM (Tp,T,Eg);
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
            iy.push_back( fvm_n1->global_offset() + node_n_offset );
            y.push_back ( (riin1*GIIn+riip1*GIIp)*truncated_partial_volume );

            iy.push_back( fvm_n1->global_offset() + node_p_offset );
            y.push_back ( (riin1*GIIn+riip1*GIIp)*truncated_partial_volume );

            if (get_advanced_model()->enable_Tn())
            {
              Hn = - (Eg+1.5*kb*Tp) * riin1*GIIn + 1.5*kb*Tn * riip1*GIIp;
              iy.push_back(fvm_n1->global_offset()+node_Tn_offset);
              y.push_back( Hn*truncated_partial_volume );
            }
            if (get_advanced_model()->enable_Tp())
            {
              Hp = - (Eg+1.5*kb*Tn) * riip1*GIIp + 1.5*kb*Tp * riin1*GIIn;
              iy.push_back(fvm_n1->global_offset()+node_Tp_offset);
              y.push_back( Hp*truncated_partial_volume );
            }
          }

          if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            iy.push_back( fvm_n2->global_offset() + node_n_offset );
            y.push_back ( (riin2*GIIn+riip2*GIIp)*truncated_partial_volume );

            iy.push_back( fvm_n2->global_offset() + node_p_offset );
            y.push_back ( (riin2*GIIn+riip2*GIIp)*truncated_partial_volume );

            if (get_advanced_model()->enable_Tn())
            {
              Hn = - (Eg+1.5*kb*Tp) * riin2*GIIn + 1.5*kb*Tn * riip2*GIIp;
              iy.push_back(fvm_n2->global_offset()+node_Tn_offset);
              y.push_back( Hn*truncated_partial_volume );
            }
            if (get_advanced_model()->enable_Tp())
            {
              Hp = - (Eg+1.5*kb*Tn) * riip2*GIIp + 1.5*kb*Tp * riin2*GIIn;
              iy.push_back(fvm_n2->global_offset()+node_Tp_offset);
              y.push_back( Hp*truncated_partial_volume );
            }
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


    PetscScalar n   =  x[fvm_node->local_offset() + node_n_offset];           // electron density
    PetscScalar p   =  x[fvm_node->local_offset() + node_p_offset];           // hole density

    PetscScalar T  = T_external();
    PetscScalar Tn = T_external();
    PetscScalar Tp = T_external();

    // lattice temperature if required
    if(get_advanced_model()->enable_Tl())
      T =  x[fvm_node->local_offset() + node_Tl_offset];

    // electron temperature if required
    if(get_advanced_model()->enable_Tn())
      Tn = x[fvm_node->local_offset() + node_Tn_offset]/n;

    // hole temperature if required
    if(get_advanced_model()->enable_Tp())
      Tp = x[fvm_node->local_offset() + node_Tp_offset]/p;

    // map this node and its data to material database
    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

    // the charge density
    PetscScalar rho = e*(node_data->Net_doping() + p - n)*fvm_node->volume();
    iy.push_back(fvm_node->global_offset()+node_psi_offset);
    y.push_back( rho );

    // the recombination term
    PetscScalar R_SHR  = mt->band->R_SHR(p,n,T);
    PetscScalar R_AUG_N = mt->band->R_Auger_N(p,n,T);
    PetscScalar R_AUG_P = mt->band->R_Auger_P(p,n,T);
    PetscScalar R_DIR  = mt->band->R_Direct(p,n,T);

    PetscScalar R   = (R_SHR + R_AUG_N + R_AUG_P + R_DIR)*fvm_node->volume();
    // consider carrier generation
    PetscScalar Field_G = node_data->Field_G()*fvm_node->volume();
    PetscScalar OptQ = node_data->OptQ()*fvm_node->volume();

    iy.push_back(fvm_node->global_offset()+node_n_offset);
    iy.push_back(fvm_node->global_offset()+node_p_offset);
    y.push_back( Field_G - R  + node_data->EIn());
    y.push_back( Field_G - R  + node_data->HIn());

    // process heat consume due to R/G and collision
    PetscScalar H=0, Hn=0, Hp=0;
    PetscScalar Eg = mt->band->Eg(T);
    PetscScalar tao_en = mt->band->ElecEnergyRelaxTime(Tn, T);
    PetscScalar tao_ep = mt->band->HoleEnergyRelaxTime(Tp, T);

    // lattice temperature if required
    if(get_advanced_model()->enable_Tl())
      H  = R_SHR*(Eg+1.5*kb*Tn+1.5*kb*Tp) + OptQ;

    // electron temperature if required
    if(get_advanced_model()->enable_Tn())
    {
      Hn = R_AUG_N*(Eg+1.5*kb*Tp) - R_AUG_P * 1.5*kb*Tn - 1.5*kb*Tn*(R_SHR+R_DIR) - 1.5*kb*n*(Tn-T)/tao_en;
      H += 1.5*kb*n*(Tn-T)/tao_en;
    }

    // hole temperature if required
    if(get_advanced_model()->enable_Tp())
    {
      Hp = R_AUG_P*(Eg+1.5*kb*Tn) - R_AUG_N * 1.5*kb*Tp - 1.5*kb*Tp*(R_SHR+R_DIR) - 1.5*kb*p*(Tp-T)/tao_ep;
      H += 1.5*kb*p*(Tp-T)/tao_ep;
    }

    // save extra equation to data buffer if required
    if(get_advanced_model()->enable_Tl())
    {
      iy.push_back(fvm_node->global_offset()+node_Tl_offset);
      y.push_back( H*fvm_node->volume() );
    }

    if(get_advanced_model()->enable_Tn())
    {
      iy.push_back(fvm_node->global_offset()+node_Tn_offset);
      y.push_back( Hn*fvm_node->volume() );
    }

    if(get_advanced_model()->enable_Tp())
    {
      iy.push_back(fvm_node->global_offset()+node_Tp_offset);
      y.push_back( Hp*fvm_node->volume() );
    }


    if (get_advanced_model()->Trap)
    {
      // consider charge trapping in semiconductor bulk (bulk_flag=true)

      // call the Trap MPI to calculate trap occupancy using the local carrier densities and lattice temperature
      PetscScalar ni = mt->band->nie(p, n, T);
      mt->trap->Calculate(true,p,n,ni,T);

      // calculate the contribution of trapped charge to Poisson's equation
      PetscScalar TrappedC = mt->trap->Charge(true);
      if (TrappedC !=0)
      {
        iy.push_back(fvm_node->global_offset()+node_psi_offset);
        y.push_back(TrappedC * fvm_node->volume());
      }

      // calculate the rates of electron and hole capture
      PetscScalar TrapElec = mt->trap->ElectronTrapRate(true,n,ni,T);
      PetscScalar TrapHole = mt->trap->HoleTrapRate    (true,p,ni,T);

      iy.push_back(fvm_node->global_offset()+node_n_offset);
      iy.push_back(fvm_node->global_offset()+node_p_offset);
      y.push_back( - TrapElec * fvm_node->volume());
      y.push_back( - TrapHole * fvm_node->volume());

      if(get_advanced_model()->enable_Tn())
      {
        Hn = - 1.5 * kb*Tn * TrapElec;
        iy.push_back(fvm_node->global_offset()+node_Tn_offset);
        y.push_back( Hn*fvm_node->volume() );
      }

      if(get_advanced_model()->enable_Tp())
      {
        Hp = - 1.5 * kb*Tp * TrapHole;
        iy.push_back(fvm_node->global_offset()+node_Tp_offset);
        y.push_back( Hp*fvm_node->volume() );
      }

      if(get_advanced_model()->enable_Tl())
      {
        PetscScalar EcEi = 0.5*Eg - kb*T*log(node_data->Nc()/node_data->Nv());
        PetscScalar EiEv = 0.5*Eg + kb*T*log(node_data->Nc()/node_data->Nv());
        H = mt->trap->TrapHeat(true,p,n,ni,Tp,Tn,T,EcEi,EiEv);
        iy.push_back(fvm_node->global_offset()+node_Tl_offset);
        y.push_back( H*fvm_node->volume() );
      }
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



void SemiconductorSimulationRegion::EBM3_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // find the node variable offset
  unsigned int node_n_offset   = ebm_variable_offset(ELECTRON);
  unsigned int node_p_offset   = ebm_variable_offset(HOLE);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);
  unsigned int node_Tn_offset  = ebm_variable_offset(E_TEMP);
  unsigned int node_Tp_offset  = ebm_variable_offset(H_TEMP);

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

  const double r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);

  // process node related terms
  // including \rho of poisson's equation and recombination term of continuation equation
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();

    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

    // process \partial t

    //second order
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false)
    {
      // electron density
      PetscScalar n         =  x[fvm_node->local_offset()+node_n_offset];
      PetscScalar dndt_BDF2 = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                              / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
      iy.push_back(fvm_node->global_offset()+node_n_offset);                     // save index in the buffer
      y.push_back( dndt_BDF2 );


      // hole density
      PetscScalar p         =  x[fvm_node->local_offset()+node_p_offset];

      PetscScalar dpdt_BDF2 = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*node_data->p() + (1-r)/r*node_data->p_last())
                              / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
      iy.push_back(fvm_node->global_offset()+node_p_offset);
      y.push_back( dpdt_BDF2 );


      // lattice temperature if required
      if(get_advanced_model()->enable_Tl())
      {
        PetscScalar Tl           =  x[fvm_node->local_offset()+node_Tl_offset];
        PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Tl);

        PetscScalar dTldt_BDF2   = -((2-r)/(1-r)*Tl - 1.0/(r*(1-r))*node_data->T() + (1-r)/r*node_data->T_last())*node_data->density()*HeatCapacity
                                   / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tl_offset);
        y.push_back( dTldt_BDF2 );
      }

      // electron temperature if required
      if(get_advanced_model()->enable_Tn())
      {
        PetscScalar n          =  x[fvm_node->local_offset()+node_n_offset];
        PetscScalar Tn         =  x[fvm_node->local_offset()+node_Tn_offset]/n;

        PetscScalar dWndt_BDF2 = -((2-r)/(1-r)*1.5*n*kb*Tn - 1.0/(r*(1-r))*1.5*node_data->n()*kb*node_data->Tn() + (1-r)/r*1.5*node_data->n_last()*kb*node_data->Tn_last())
                                 / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tn_offset);
        y.push_back( dWndt_BDF2 );
      }

      // hole temperature if required
      if(get_advanced_model()->enable_Tp())
      {
        PetscScalar p          =  x[fvm_node->local_offset()+node_p_offset];
        PetscScalar Tp         =  x[fvm_node->local_offset()+node_Tp_offset]/p;

        PetscScalar dWpdt_BDF2 = -((2-r)/(1-r)*1.5*p*kb*Tp - 1.0/(r*(1-r))*1.5*node_data->p()*kb*node_data->Tp() + (1-r)/r*1.5*node_data->p_last()*kb*node_data->Tp_last())
                                 / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tp_offset);
        y.push_back( dWpdt_BDF2 );
      }
    }
    else //first order
    {
      // electron density
      PetscScalar n         =  x[fvm_node->local_offset()+node_n_offset];
      PetscScalar dndt_BDF1 = -(n - node_data->n())/SolverSpecify::dt*fvm_node->volume();
      iy.push_back(fvm_node->global_offset()+node_n_offset);                     // save index in the buffer
      y.push_back( dndt_BDF1 );


      // hole density
      PetscScalar p         =  x[fvm_node->local_offset()+node_p_offset];
      PetscScalar dpdt_BDF1 = -(p - node_data->p())/SolverSpecify::dt*fvm_node->volume();
      iy.push_back(fvm_node->global_offset()+node_p_offset);
      y.push_back( dpdt_BDF1 );


      // lattice temperature if required
      if(get_advanced_model()->enable_Tl())
      {
        PetscScalar Tl           =  x[fvm_node->local_offset()+node_Tl_offset];
        PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Tl);
        PetscScalar dTldt_BDF1   = -(Tl - node_data->T())*node_data->density()*HeatCapacity/SolverSpecify::dt*fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tl_offset);
        y.push_back( dTldt_BDF1 );
      }

      // electron temperature if required
      if(get_advanced_model()->enable_Tn())
      {
        PetscScalar n          =  x[fvm_node->local_offset()+node_n_offset];
        PetscScalar Tn         =  x[fvm_node->local_offset()+node_Tn_offset]/n;

        PetscScalar dWndt_BDF1 = -(1.5*n*kb*Tn - 1.5*node_data->n()*kb*node_data->Tn())/SolverSpecify::dt*fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tn_offset);
        y.push_back( dWndt_BDF1 );
      }

      // hole temperature if required
      if(get_advanced_model()->enable_Tp())
      {
        PetscScalar p          =  x[fvm_node->local_offset()+node_p_offset];
        PetscScalar Tp         =  x[fvm_node->local_offset()+node_Tp_offset]/p;

        PetscScalar dWpdt_BDF1 = -(1.5*p*kb*Tp - 1.5*node_data->p()*kb*node_data->Tp())/SolverSpecify::dt*fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tp_offset);
        y.push_back( dWpdt_BDF1 );
      }
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

