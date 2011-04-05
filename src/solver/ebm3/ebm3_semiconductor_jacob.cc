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


/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 * AD is fully used here
 */
void SemiconductorSimulationRegion::EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // find the node variable offset
  unsigned int n_node_var      = ebm_n_variables();
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

    //the indepedent variable number, this->ebm_n_variables()*n_nodes
    adtl::AutoDScalar::numdir = n_node_var*elem->n_nodes();

    //synchronize with material database
    mt->set_ad_num(adtl::AutoDScalar::numdir);

    // indicate the column position of the variables in the matrix
    std::vector<PetscInt> cell_col;
    for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
    {
      const FVM_Node * fvm_node = elem->get_fvm_node(nd);
      for(unsigned int nv=0; nv<n_node_var; nv++)
        cell_col.push_back( fvm_node->global_offset() + nv );
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
          // use values in the current iteration
          V  =  x[fvm_node->local_offset()+node_psi_offset];   V.setADValue(n_node_var*nd + node_psi_offset, 1.0);
          n  =  x[fvm_node->local_offset()+node_n_offset];     n.setADValue(n_node_var*nd + node_n_offset, 1.0);
          p  =  x[fvm_node->local_offset()+node_p_offset];     p.setADValue(n_node_var*nd + node_p_offset, 1.0);
          n  +=  fvm_node_data->ni()*1e-2;
          p  +=  fvm_node_data->ni()*1e-2;
        }
        else
        {
          // n and p use previous solution value
          V  =  x[fvm_node->local_offset()+node_psi_offset];   V.setADValue(n_node_var*nd + node_psi_offset, 1.0);
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
            V_neighbor.setADValue(n_node_var*elem->n_nodes()+nd, 1.0);
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

      // partial area associated with this edge
      double partial_area = elem->partial_area_with_edge(ne);
      double partial_volume = elem->partial_volume_with_edge(ne); // partial volume associated with this edge

      double truncated_partial_area =  partial_area;
      double truncated_partial_volume =  partial_volume;
      if(truncation)
      {
        truncated_partial_area =  elem->partial_area_with_edge_truncated(ne);
        truncated_partial_volume =  elem->partial_volume_with_edge_truncated(ne);
      }

      // fvm_node of node1
      const FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes.first);
      // fvm_node of node2
      const FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes.second);

      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data() ;   genius_assert(n1_data);
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data() ;   genius_assert(n2_data);

      unsigned int n1_local_offset = fvm_n1->local_offset();
      unsigned int n2_local_offset = fvm_n2->local_offset();

      // the row position of variables in the matrix
      std::vector<PetscInt> row1, row2;
      for(unsigned int nv=0; nv<n_node_var; ++nv)  row1.push_back( fvm_n1->global_offset()+nv );
      for(unsigned int nv=0; nv<n_node_var; ++nv)  row2.push_back( fvm_n2->global_offset()+nv );


      // here we use AD again. Can we hand write it for more efficient?
      {

        //for node 1 of the edge
        mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

        AutoDScalar V1 = x[n1_local_offset + node_psi_offset];
        V1.setADValue(n_node_var*edge_nodes.first + node_psi_offset, 1.0);   // electrostatic potential

        AutoDScalar n1 = x[n1_local_offset + node_n_offset];
        n1.setADValue(n_node_var*edge_nodes.first + node_n_offset, 1.0);     // electron density

        AutoDScalar p1 = x[n1_local_offset + node_p_offset];
        p1.setADValue(n_node_var*edge_nodes.first + node_p_offset, 1.0);     // hole density

        AutoDScalar T1  =  T_external();
        AutoDScalar Tn1 =  T_external();
        AutoDScalar Tp1 =  T_external();

        // lattice temperature if required
        if(get_advanced_model()->enable_Tl())
        {
          T1 =  x[n1_local_offset + node_Tl_offset];
          T1.setADValue(n_node_var*edge_nodes.first + node_Tl_offset, 1.0);
        }

        // electron temperature if required
        if(get_advanced_model()->enable_Tn())
        {
          AutoDScalar n1Tn1 = x[n1_local_offset + node_Tn_offset];
          n1Tn1.setADValue(n_node_var*edge_nodes.first + node_Tn_offset, 1.0);
          Tn1 = n1Tn1/n1;
        }

        // hole temperature if required
        if(get_advanced_model()->enable_Tp())
        {
          AutoDScalar p1Tp1 = x[n1_local_offset + node_Tp_offset];
          p1Tp1.setADValue(n_node_var*edge_nodes.first + node_Tp_offset, 1.0);
          Tp1 = p1Tp1/p1;
        }

        AutoDScalar Ec1 =  -(e*V1 + n1_data->affinity() + kb*T1*log(mt->band->nie(p1, n1, T1)));//conduct band energy level
        AutoDScalar Ev1 =  -(e*V1 + n1_data->affinity() - kb*T1*log(mt->band->nie(p1, n1, T1)));//valence band energy level
        if(get_advanced_model()->Fermi)
        {
          Ec1 = Ec1 - kb*T1*log(gamma_f(fabs(n1)/n1_data->Nc()));
          Ev1 = Ev1 + kb*T1*log(gamma_f(fabs(p1)/n1_data->Nv()));
        }
        PetscScalar eps1 =  n1_data->eps();                        // eps
        AutoDScalar kap1 =  mt->thermal->HeatConduction(T1);
        AutoDScalar Eg1= mt->band->Eg(T1);


        //for node 2 of the edge
        mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

        AutoDScalar V2   =  x[n2_local_offset + node_psi_offset];
        V2.setADValue(n_node_var*edge_nodes.second + node_psi_offset, 1.0);   // electrostatic potential

        AutoDScalar n2   =  x[n2_local_offset + node_n_offset];
        n2.setADValue(n_node_var*edge_nodes.second + node_n_offset, 1.0);     // electron density

        AutoDScalar p2   =  x[n2_local_offset + node_p_offset];
        p2.setADValue(n_node_var*edge_nodes.second + node_p_offset, 1.0);     // hole density

        AutoDScalar T2  =  T_external();
        AutoDScalar Tn2 =  T_external();
        AutoDScalar Tp2 =  T_external();

        // lattice temperature if required
        if(get_advanced_model()->enable_Tl())
        {
          T2 =  x[n2_local_offset + node_Tl_offset];
          T2.setADValue(n_node_var*edge_nodes.second+node_Tl_offset, 1.0);
        }

        // electron temperature if required
        if(get_advanced_model()->enable_Tn())
        {
          AutoDScalar n2Tn2 = x[n2_local_offset + node_Tn_offset];
          n2Tn2.setADValue(n_node_var*edge_nodes.second+node_Tn_offset, 1.0);
          Tn2 = n2Tn2/n2;
        }

        // hole temperature if required
        if(get_advanced_model()->enable_Tp())
        {
          AutoDScalar p2Tp2 = x[n2_local_offset + node_Tp_offset];
          p2Tp2.setADValue(n_node_var*edge_nodes.second+node_Tp_offset, 1.0);
          Tp2 = p2Tp2/p2;
        }

        AutoDScalar Ec2 =  -(e*V2 + n2_data->affinity() + kb*T2*log(mt->band->nie(p2, n2, T2)));//conduct band energy level
        AutoDScalar Ev2 =  -(e*V2 + n2_data->affinity() - kb*T2*log(mt->band->nie(p2, n2, T2)));//valence band energy level
        if(get_advanced_model()->Fermi)
        {
          Ec2 = Ec2 - kb*T2*log(gamma_f(fabs(n2)/n2_data->Nc()));
          Ev2 = Ev2 + kb*T2*log(gamma_f(fabs(p2)/n2_data->Nv()));
        }
        PetscScalar eps2 =  n2_data->eps();                         // eps
        AutoDScalar kap2 =  mt->thermal->HeatConduction(T2);
        AutoDScalar Eg2= mt->band->Eg(T2);


        AutoDScalar mun1;   // electron mobility for node 1 of the edge
        AutoDScalar mup1;   // hole mobility for node 1 of the edge
        AutoDScalar mun2;   // electron mobility  for node 2 of the edge
        AutoDScalar mup2;   // hole mobility for node 2 of the edge

        if(highfield_mob)
        {
          if (get_advanced_model()->Mob_Force == ModelSpecify::ESimple && !insulator_interface_elem  )
          {
            Point _dir = (*fvm_n1->root_node() - *fvm_n2->root_node()).unit();
            VectorValue<AutoDScalar> dir(_dir(0), _dir(1), _dir(2));
            AutoDScalar Ep = fabs((V2-V1)/length);
            AutoDScalar Et = 0;//
            if(mos_channel_elem)
              Et = (E - (E*dir)*dir).size();

            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T1, Ep, Et, Tn1);
            mup1 = mt->mob->HoleMob(p1, n1, T1, Ep, Et, Tp1);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T2, Ep, Et, Tn2);
            mup2 = mt->mob->HoleMob(p2, n2, T2, Ep, Et, Tp2);
          }
          else // ModelSpecify::EJ || ModelSpecify::EQF
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

        AutoDScalar mun = 0.5*(mun1+mun2);  // the electron mobility at the mid point of the edge, use linear interpolation
        AutoDScalar mup = 0.5*(mup1+mup2);  // the hole mobility at the mid point of the edge, use linear interpolation


        PetscScalar eps = 0.5*(eps1+eps2); // eps at mid point of the edge
        AutoDScalar kap = 0.5*(kap1+kap2); // kapa at mid point of the edge


        // S-G current along the edge, call different SG scheme selected by EBM level
        AutoDScalar Jn, Jp, Sn=0, Sp=0;

        switch(Jn_level)
        {
        case 1:
          Jn =  mun*In_dd(kb*T_external()/e, (Ec2-Ec1)/e, n1, n2, length);
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
          Jp =  mup*Ip_dd(kb*T_external()/e, (Ev2-Ev1)/e, p1, p2, length);
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
        AutoDScalar H=0, Hn=0, Hp=0;

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



        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {

          AutoDScalar poisson = ( eps*(V2 - V1)/length*partial_area );
          MatSetValues(*jac, 1, &row1[node_psi_offset], cell_col.size(), &cell_col[0], poisson.getADValue(), ADD_VALUES);

          AutoDScalar electron_continuation = ( Jn*truncated_partial_area );
          MatSetValues(*jac, 1, &row1[node_n_offset], cell_col.size(), &cell_col[0], electron_continuation.getADValue(), ADD_VALUES);

          AutoDScalar hole_continuation = ( - Jp*truncated_partial_area );
          MatSetValues(*jac, 1, &row1[node_p_offset], cell_col.size(), &cell_col[0], hole_continuation.getADValue(), ADD_VALUES);

          // heat transport equation if required
          if(get_advanced_model()->enable_Tl())
          {
            AutoDScalar heating_equ = ( kap*(T2 - T1)/length*partial_area + H*truncated_partial_area);
            MatSetValues(*jac, 1, &row1[node_Tl_offset], cell_col.size(), &cell_col[0], heating_equ.getADValue(), ADD_VALUES);
          }


          // energy balance equation for electron if required
          if(get_advanced_model()->enable_Tn())
          {
            AutoDScalar electron_energy = -Sn*truncated_partial_area + Hn*truncated_partial_area;
            MatSetValues(*jac, 1, &row1[node_Tn_offset], cell_col.size(), &cell_col[0], electron_energy.getADValue(), ADD_VALUES);
          }


          // energy balance equation for hole if required
          if(get_advanced_model()->enable_Tp())
          {
            AutoDScalar hole_energy = -Sp*truncated_partial_area + Hp*truncated_partial_area;
            MatSetValues(*jac, 1, &row1[node_Tp_offset], cell_col.size(), &cell_col[0], hole_energy.getADValue(), ADD_VALUES);
          }

        }

        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {

          AutoDScalar poisson = ( eps*(V1 - V2)/length*partial_area );
          MatSetValues(*jac, 1, &row2[node_psi_offset], cell_col.size(), &cell_col[0], poisson.getADValue(), ADD_VALUES);

          AutoDScalar electron_continuation = ( - Jn*truncated_partial_area );
          MatSetValues(*jac, 1, &row2[node_n_offset], cell_col.size(), &cell_col[0], electron_continuation.getADValue(), ADD_VALUES);

          AutoDScalar hole_continuation = ( Jp*truncated_partial_area );
          MatSetValues(*jac, 1, &row2[node_p_offset], cell_col.size(), &cell_col[0], hole_continuation.getADValue(), ADD_VALUES);

          // heat transport equation if required
          if(get_advanced_model()->enable_Tl())
          {
            AutoDScalar heating_equ = ( kap*(T1 - T2)/length*partial_area + H*truncated_partial_area);
            MatSetValues(*jac, 1, &row2[node_Tl_offset], cell_col.size(), &cell_col[0], heating_equ.getADValue(), ADD_VALUES);
          }

          // energy balance equation for electron if required
          if(get_advanced_model()->enable_Tn())
          {
            AutoDScalar electron_energy = Sn*truncated_partial_area + Hn*truncated_partial_area;
            MatSetValues(*jac, 1, &row2[node_Tn_offset], cell_col.size(), &cell_col[0], electron_energy.getADValue(), ADD_VALUES);
          }

          // energy balance equation for hole if required
          if(get_advanced_model()->enable_Tp())
          {
            AutoDScalar hole_energy = Sp*truncated_partial_area + Hp*truncated_partial_area;
            MatSetValues(*jac, 1, &row2[node_Tp_offset], cell_col.size(), &cell_col[0], hole_energy.getADValue(), ADD_VALUES);
          }

        }

        if (get_advanced_model()->BandBandTunneling && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {

          AutoDScalar GBTBT1 = mt->band->BB_Tunneling(T1, E.size());
          AutoDScalar GBTBT2 = mt->band->BB_Tunneling(T2, E.size());

          if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            AutoDScalar continuity = 0.5*GBTBT1*truncated_partial_volume;
            MatSetValues(*jac, 1, &row1[node_n_offset], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row1[node_p_offset], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
          }

          if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation
            AutoDScalar continuity = 0.5*GBTBT2*truncated_partial_volume;
            MatSetValues(*jac, 1, &row2[node_n_offset], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row2[node_p_offset], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
          }
        }


        if (get_advanced_model()->ImpactIonization && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          // consider impact-ionization

          AutoDScalar IIn,IIp,GIIn,GIIp;
          AutoDScalar T,Tn,Tp,Eg;

          // FIXME should use weighted carrier temperature.
          T  = std::min(T1,T2);
          Eg = 0.5* ( Eg1 + Eg2 );
          Tn = std::min(Tn1,Tn2);
          Tp = std::min(Tp1,Tp2);

          VectorValue<PetscScalar> ev0 = (elem->point(edge_nodes.second) - elem->point(edge_nodes.first));
          VectorValue<AutoDScalar> ev;
          ev(0)=ev0(0); ev(1)=ev0(1); ev(2)=ev0(2);
          AutoDScalar riin1 = 0.5 + 0.5*(ev.unit()).dot((Jnv.unit(true)));
          AutoDScalar riin2 = 1.0 - riin1;
          AutoDScalar riip2 = 0.5 + 0.5*(ev.unit()).dot((Jpv.unit(true)));
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
            AutoDScalar electron_continuity = (riin1*GIIn+riip1*GIIp)*truncated_partial_volume ;
            AutoDScalar hole_continuity     = (riin1*GIIn+riip1*GIIp)*truncated_partial_volume ;
            MatSetValues(*jac, 1, &row1[node_n_offset], cell_col.size(), &cell_col[0], electron_continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row1[node_p_offset], cell_col.size(), &cell_col[0], hole_continuity.getADValue(), ADD_VALUES);

            if (get_advanced_model()->enable_Tn())
            {
              Hn = - (Eg+1.5*kb*Tp) * riin1*GIIn + 1.5*kb*Tn * riip1*GIIp;
              AutoDScalar electron_energy = Hn*truncated_partial_volume;
              MatSetValues(*jac, 1, &row1[node_Tn_offset], cell_col.size(), &cell_col[0], electron_energy.getADValue(), ADD_VALUES);
            }
            if (get_advanced_model()->enable_Tp())
            {
              Hp = - (Eg+1.5*kb*Tn) * riip1*GIIp + 1.5*kb*Tp * riin1*GIIn;
              AutoDScalar hole_energy = Hp*truncated_partial_volume;
              MatSetValues(*jac, 1, &row1[node_Tp_offset], cell_col.size(), &cell_col[0], hole_energy.getADValue(), ADD_VALUES);
            }
          }

          if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
          {
            // continuity equation of electron
            AutoDScalar electron_continuity = (riin2*GIIn+riip2*GIIp)*truncated_partial_volume ;
            AutoDScalar hole_continuity     = (riin2*GIIn+riip2*GIIp)*truncated_partial_volume ;
            MatSetValues(*jac, 1, &row2[node_n_offset], cell_col.size(), &cell_col[0], electron_continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row2[node_p_offset], cell_col.size(), &cell_col[0], hole_continuity.getADValue(), ADD_VALUES);

            if (get_advanced_model()->enable_Tn())
            {
              Hn = - (Eg+1.5*kb*Tp) * riin2*GIIn + 1.5*kb*Tn * riip2*GIIp;
              AutoDScalar electron_energy = Hn*truncated_partial_volume;
              MatSetValues(*jac, 1, &row2[node_Tn_offset], cell_col.size(), &cell_col[0], electron_energy.getADValue(), ADD_VALUES);
            }
            if (get_advanced_model()->enable_Tp())
            {
              Hp = - (Eg+1.5*kb*Tn) * riip2*GIIp + 1.5*kb*Tp * riin2*GIIn;
              AutoDScalar hole_energy = Hp*truncated_partial_volume;
              MatSetValues(*jac, 1, &row2[node_Tp_offset], cell_col.size(), &cell_col[0], hole_energy.getADValue(), ADD_VALUES);
            }
          }
        } // end of II

      }
    }// end of scan all edges of the cell

  }// end of scan all the cell

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // process node related terms
  // including \rho of poisson's equation, recombination term of continuation equation and heat consume due to R/G and collision

  //the indepedent variable number, n_node_var for each node
  adtl::AutoDScalar::numdir = n_node_var;

  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  // process node related terms
  // including \rho of poisson's equation and recombination term of continuation equation
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();


    std::vector<PetscInt> index;
    for(unsigned int nv=0; nv<n_node_var; ++nv)  index.push_back( fvm_node->global_offset()+nv );

    AutoDScalar V = x[fvm_node->local_offset() + node_psi_offset];
    V.setADValue(node_psi_offset, 1.0);   // electrostatic potential

    AutoDScalar n = x[fvm_node->local_offset() + node_n_offset];
    n.setADValue(node_n_offset, 1.0);     // electron density

    AutoDScalar p = x[fvm_node->local_offset() + node_p_offset];
    p.setADValue(node_p_offset, 1.0);     // hole density

    AutoDScalar T  =  T_external();
    AutoDScalar Tn =  T_external();
    AutoDScalar Tp =  T_external();

    // lattice temperature if required
    if(get_advanced_model()->enable_Tl())
    {
      T =  x[fvm_node->local_offset() + node_Tl_offset];
      T.setADValue(node_Tl_offset, 1.0);
    }

    // electron temperature if required
    if(get_advanced_model()->enable_Tn())
    {
      AutoDScalar nTn = x[fvm_node->local_offset() + node_Tn_offset];
      nTn.setADValue(node_Tn_offset, 1.0);
      Tn = nTn/n;
    }

    // hole temperature if required
    if(get_advanced_model()->enable_Tp())
    {
      AutoDScalar pTp = x[fvm_node->local_offset() + node_Tp_offset];
      pTp.setADValue(node_Tp_offset, 1.0);
      Tp = pTp/p;
    }

    // map this node and its data to material database
    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

    // the charge density for poisson's equation
    AutoDScalar rho = e*(node_data->Net_doping() + p - n)*fvm_node->volume();
    MatSetValues(*jac, 1, &index[node_psi_offset], n_node_var, &index[0], rho.getADValue(), ADD_VALUES);

    // the recombination term
    AutoDScalar R_SHR  = mt->band->R_SHR(p,n,T);
    AutoDScalar R_AUG_N  = mt->band->R_Auger_N(p,n,T);
    AutoDScalar R_AUG_P  = mt->band->R_Auger_P(p,n,T);
    AutoDScalar R_DIR  = mt->band->R_Direct(p,n,T);
    AutoDScalar R   = (R_SHR + R_AUG_N + R_AUG_P + R_DIR)*fvm_node->volume();
    AutoDScalar G   = 0;
    MatSetValues(*jac, 1, &index[node_n_offset], n_node_var, &index[0], (G-R).getADValue(),   ADD_VALUES);
    MatSetValues(*jac, 1, &index[node_p_offset], n_node_var, &index[0], (G-R).getADValue(),   ADD_VALUES);


    // process heat consume due to R/G and collision
    AutoDScalar H=0, Hn=0, Hp=0;
    AutoDScalar Eg = mt->band->Eg(T);
    AutoDScalar tao_en = mt->band->ElecEnergyRelaxTime(Tn, T);
    AutoDScalar tao_ep = mt->band->HoleEnergyRelaxTime(Tp, T);

    // lattice temperature if required
    if(get_advanced_model()->enable_Tl())
      H  = R_SHR*(Eg+1.5*kb*Tn+1.5*kb*Tp);

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
      MatSetValues(*jac, 1, &index[node_Tl_offset], n_node_var, &index[0], (H*fvm_node->volume()).getADValue(),   ADD_VALUES);
    }

    if(get_advanced_model()->enable_Tn())
    {
      MatSetValues(*jac, 1, &index[node_Tn_offset], n_node_var, &index[0], (Hn*fvm_node->volume()).getADValue(),   ADD_VALUES);
    }

    if(get_advanced_model()->enable_Tp())
    {
      MatSetValues(*jac, 1, &index[node_Tp_offset], n_node_var, &index[0], (Hp*fvm_node->volume()).getADValue(),   ADD_VALUES);
    }

    if (get_advanced_model()->Trap)
    {
      // consider charge trapping in semiconductor bulk (bulk_flag=true)

      // call the Trap MPI to calculate trap occupancy using the local carrier densities and lattice temperature
      AutoDScalar ni = mt->band->nie(p, n, T);
      mt->trap->Calculate(true,p,n,ni,T);

      // calculate the contribution of trapped charge to Poisson's equation
      AutoDScalar TrappedC = mt->trap->ChargeAD(true);
      if (TrappedC !=0)
      {
        MatSetValues(*jac, 1, &index[node_psi_offset], n_node_var, &index[0], (TrappedC*fvm_node->volume()).getADValue(), ADD_VALUES);
      }

      // calculate the rates of electron and hole capture
      AutoDScalar TrapElec = mt->trap->ElectronTrapRate(true,n,ni,T);
      AutoDScalar TrapHole = mt->trap->HoleTrapRate    (true,p,ni,T);

      MatSetValues(*jac, 1, &index[node_n_offset], n_node_var, &index[0], (-TrapElec*fvm_node->volume()).getADValue(),   ADD_VALUES);
      MatSetValues(*jac, 1, &index[node_p_offset], n_node_var, &index[0], (-TrapHole*fvm_node->volume()).getADValue(),   ADD_VALUES);

      if(get_advanced_model()->enable_Tn())
      {
        Hn = - 1.5 * kb*Tn * TrapElec;
        MatSetValues(*jac, 1, &index[node_Tn_offset], n_node_var, &index[0], (Hn*fvm_node->volume()).getADValue(),   ADD_VALUES);
      }

      if(get_advanced_model()->enable_Tp())
      {
        Hp = - 1.5 * kb*Tp * TrapHole;
        MatSetValues(*jac, 1, &index[node_Tp_offset], n_node_var, &index[0], (Hp*fvm_node->volume()).getADValue(),   ADD_VALUES);
      }

      if(get_advanced_model()->enable_Tl())
      {
        AutoDScalar EcEi = 0.5*Eg - kb*T*log(node_data->Nc()/node_data->Nv());
        AutoDScalar EiEv = 0.5*Eg + kb*T*log(node_data->Nc()/node_data->Nv());
        H = mt->trap->TrapHeat(true,p,n,ni,Tp,Tn,T,EcEi,EiEv);
        MatSetValues(*jac, 1, &index[node_Tl_offset], n_node_var, &index[0], (H*fvm_node->volume()).getADValue(),   ADD_VALUES);
      }
    }

  }


  // boundary condition should be processed later!

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
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false)
    {
      PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);

      // electron density
      PetscScalar n    =  x[fvm_node->local_offset()+node_n_offset];
      PetscScalar dndt = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                         / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
      iy.push_back(fvm_node->global_offset()+node_n_offset);                     // save index in the buffer
      y.push_back( dndt );


      // hole density
      PetscScalar p    =  x[fvm_node->local_offset()+node_p_offset];
      PetscScalar dpdt = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*node_data->p() + (1-r)/r*node_data->p_last())
                         / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
      iy.push_back(fvm_node->global_offset()+node_p_offset);
      y.push_back( dpdt );


      // lattice temperature if required
      if(get_advanced_model()->enable_Tl())
      {
        PetscScalar Tl    =  x[fvm_node->local_offset()+node_Tl_offset];
        PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Tl);
        PetscScalar dTldt = -((2-r)/(1-r)*Tl - 1.0/(r*(1-r))*node_data->T() + (1-r)/r*node_data->T_last())*node_data->density()*HeatCapacity
                            / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tl_offset);
        y.push_back( dTldt );
      }

      // electron temperature if required
      if(get_advanced_model()->enable_Tn())
      {
        PetscScalar n     =  x[fvm_node->local_offset()+node_n_offset];
        PetscScalar Tn    =  x[fvm_node->local_offset()+node_Tn_offset]/n;
        PetscScalar dWndt = -((2-r)/(1-r)*1.5*n*kb*Tn - 1.0/(r*(1-r))*1.5*node_data->n()*kb*node_data->Tn() + (1-r)/r*1.5*node_data->n_last()*kb*node_data->Tn_last())
                            / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tn_offset);
        y.push_back( dWndt );
      }

      // hole temperature if required
      if(get_advanced_model()->enable_Tp())
      {
        PetscScalar p     =  x[fvm_node->local_offset()+node_p_offset];
        PetscScalar Tp    =  x[fvm_node->local_offset()+node_Tp_offset]/p;
        PetscScalar dWpdt = -((2-r)/(1-r)*1.5*p*kb*Tp - 1.0/(r*(1-r))*1.5*node_data->p()*kb*node_data->Tp() + (1-r)/r*1.5*node_data->p_last()*kb*node_data->Tp_last())
                            / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tp_offset);
        y.push_back( dWpdt );
      }
    }
    else //first order
    {
      // electron density
      PetscScalar n    =  x[fvm_node->local_offset()+node_n_offset];
      PetscScalar dndt = -(n - node_data->n())/SolverSpecify::dt*fvm_node->volume();
      iy.push_back(fvm_node->global_offset()+node_n_offset);                     // save index in the buffer
      y.push_back( dndt );


      // hole density
      PetscScalar p    =  x[fvm_node->local_offset()+node_p_offset];
      PetscScalar dpdt = -(p - node_data->p())/SolverSpecify::dt*fvm_node->volume();
      iy.push_back(fvm_node->global_offset()+node_p_offset);
      y.push_back( dpdt );


      // lattice temperature if required
      if(get_advanced_model()->enable_Tl())
      {
        PetscScalar Tl    =  x[fvm_node->local_offset()+node_Tl_offset];
        PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Tl);
        PetscScalar dTldt = -(Tl - node_data->T())*node_data->density()*HeatCapacity/SolverSpecify::dt*fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tl_offset);
        y.push_back( dTldt );
      }

      // electron temperature if required
      if(get_advanced_model()->enable_Tn())
      {
        PetscScalar n     =  x[fvm_node->local_offset()+node_n_offset];
        PetscScalar Tn    =  x[fvm_node->local_offset()+node_Tn_offset]/n;

        PetscScalar dWndt = -(1.5*n*kb*Tn - 1.5*node_data->n()*kb*node_data->Tn())/SolverSpecify::dt*fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tn_offset);
        y.push_back( dWndt );
      }

      // hole temperature if required
      if(get_advanced_model()->enable_Tp())
      {
        PetscScalar p     =  x[fvm_node->local_offset()+node_p_offset];
        PetscScalar Tp    =  x[fvm_node->local_offset()+node_Tp_offset]/p;

        PetscScalar dWpdt = -(1.5*p*kb*Tp - 1.5*node_data->p()*kb*node_data->Tp())/SolverSpecify::dt*fvm_node->volume();
        iy.push_back(fvm_node->global_offset()+node_Tp_offset);
        y.push_back( dWpdt );
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





void SemiconductorSimulationRegion::EBM3_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // find the node variable offset
  unsigned int node_n_offset   = ebm_variable_offset(ELECTRON);
  unsigned int node_p_offset   = ebm_variable_offset(HOLE);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);
  unsigned int node_Tn_offset  = ebm_variable_offset(E_TEMP);
  unsigned int node_Tp_offset  = ebm_variable_offset(H_TEMP);

  // note, we will use ADD_VALUES to set values of matrix J
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  //the indepedent variable number, max 2 for each node
  adtl::AutoDScalar::numdir = 2;

  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

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
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false)
    {
      PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);

      // electron density
      AutoDScalar n    =  x[fvm_node->local_offset()+node_n_offset];   n.setADValue(0, 1.0);
      AutoDScalar dndt = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*node_data->n() + (1-r)/r*node_data->n_last())
                         / (SolverSpecify::dt_last+SolverSpecify::dt)*fvm_node->volume();
      MatSetValue(*jac, fvm_node->global_offset()+node_n_offset,  fvm_node->global_offset()+node_n_offset, dndt.getADValue(0), ADD_VALUES);

      // hole density
      AutoDScalar p    =  x[fvm_node->local_offset()+node_p_offset];   p.setADValue(0, 1.0);
      AutoDScalar dpdt = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*node_data->p() + (1-r)/r*node_data->p_last())
                         / (SolverSpecify::dt_last+SolverSpecify::dt)*fvm_node->volume();
      MatSetValue(*jac, fvm_node->global_offset()+node_p_offset,  fvm_node->global_offset()+node_p_offset, dpdt.getADValue(0), ADD_VALUES);

      // lattice temperature if required
      if(get_advanced_model()->enable_Tl())
      {
        AutoDScalar Tl    =  x[fvm_node->local_offset()+3];   Tl.setADValue(0, 1.0);              // lattice temperature
        AutoDScalar HeatCapacity =  mt->thermal->HeatCapacity(Tl);
        AutoDScalar dTldt = -((2-r)/(1-r)*Tl - 1.0/(r*(1-r))*node_data->T() + (1-r)/r*node_data->T_last())*node_data->density()*HeatCapacity
                            / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
        MatSetValue(*jac, fvm_node->global_offset()+node_Tl_offset,  fvm_node->global_offset()+node_Tl_offset, dTldt.getADValue(0), ADD_VALUES);
      }

      // electron temperature if required
      if(get_advanced_model()->enable_Tn())
      {
        AutoDScalar n     =  x[fvm_node->local_offset()+node_n_offset];  n.setADValue(0, 1.0);
        AutoDScalar nTn   =  x[fvm_node->local_offset()+node_Tn_offset]; nTn.setADValue(1, 1.0);
        AutoDScalar Tn    =  nTn/n;
        AutoDScalar dWndt = -((2-r)/(1-r)*1.5*n*kb*Tn - 1.0/(r*(1-r))*1.5*node_data->n()*kb*node_data->Tn() + (1-r)/r*1.5*node_data->n_last()*kb*node_data->Tn_last())
                            / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
        MatSetValue(*jac, fvm_node->global_offset()+node_Tn_offset,  fvm_node->global_offset()+node_n_offset,  dWndt.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, fvm_node->global_offset()+node_Tn_offset,  fvm_node->global_offset()+node_Tn_offset, dWndt.getADValue(1), ADD_VALUES);
      }

      // hole temperature if required
      if(get_advanced_model()->enable_Tn())
      {
        AutoDScalar p     =  x[fvm_node->local_offset()+node_p_offset];  p.setADValue(0, 1.0);
        AutoDScalar pTp   =  x[fvm_node->local_offset()+node_Tp_offset]; pTp.setADValue(1, 1.0);
        AutoDScalar Tp    =  pTp/p;
        AutoDScalar dWpdt = -((2-r)/(1-r)*1.5*p*kb*Tp - 1.0/(r*(1-r))*1.5*node_data->p()*kb*node_data->Tp() + (1-r)/r*1.5*node_data->p_last()*kb*node_data->Tp_last())
                            / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
        MatSetValue(*jac, fvm_node->global_offset()+node_Tp_offset,  fvm_node->global_offset()+node_p_offset,  dWpdt.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, fvm_node->global_offset()+node_Tp_offset,  fvm_node->global_offset()+node_Tp_offset, dWpdt.getADValue(1), ADD_VALUES);
      }

    }
    else //first order
    {
      // electron density
      AutoDScalar n    =  x[fvm_node->local_offset()+node_n_offset];   n.setADValue(0, 1.0);
      AutoDScalar dndt = -(n - node_data->n())/SolverSpecify::dt*fvm_node->volume();
      MatSetValue(*jac, fvm_node->global_offset()+node_n_offset,  fvm_node->global_offset()+node_n_offset, dndt.getADValue(0), ADD_VALUES);

      // hole density
      AutoDScalar p    =  x[fvm_node->local_offset()+node_p_offset];   p.setADValue(0, 1.0);
      AutoDScalar dpdt = -(p - node_data->p())/SolverSpecify::dt*fvm_node->volume();
      MatSetValue(*jac, fvm_node->global_offset()+node_p_offset,  fvm_node->global_offset()+node_p_offset, dpdt.getADValue(0), ADD_VALUES);

      // lattice temperature if required
      if(get_advanced_model()->enable_Tl())
      {
        AutoDScalar Tl    =  x[fvm_node->local_offset()+3];   Tl.setADValue(0, 1.0);              // lattice temperature
        AutoDScalar HeatCapacity =  mt->thermal->HeatCapacity(Tl);
        AutoDScalar dTldt = -(Tl - node_data->T())*node_data->density()*HeatCapacity/SolverSpecify::dt*fvm_node->volume();
        MatSetValue(*jac, fvm_node->global_offset()+node_Tl_offset,  fvm_node->global_offset()+node_Tl_offset, dTldt.getADValue(0), ADD_VALUES);
      }

      // electron temperature if required
      if(get_advanced_model()->enable_Tn())
      {
        AutoDScalar n     =  x[fvm_node->local_offset()+node_n_offset];  n.setADValue(0, 1.0);
        AutoDScalar nTn   =  x[fvm_node->local_offset()+node_Tn_offset]; nTn.setADValue(1, 1.0);
        AutoDScalar Tn    =  nTn/n;
        AutoDScalar dWndt = -(1.5*n*kb*Tn - 1.5*node_data->n()*kb*node_data->Tn())/SolverSpecify::dt*fvm_node->volume();
        MatSetValue(*jac, fvm_node->global_offset()+node_Tn_offset,  fvm_node->global_offset()+node_n_offset,  dWndt.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, fvm_node->global_offset()+node_Tn_offset,  fvm_node->global_offset()+node_Tn_offset, dWndt.getADValue(1), ADD_VALUES);
      }

      // hole temperature if required
      if(get_advanced_model()->enable_Tn())
      {
        AutoDScalar p     =  x[fvm_node->local_offset()+node_p_offset];  p.setADValue(0, 1.0);
        AutoDScalar pTp   =  x[fvm_node->local_offset()+node_Tp_offset]; pTp.setADValue(1, 1.0);
        AutoDScalar Tp    =  pTp/p;
        AutoDScalar dWpdt = -(1.5*p*kb*Tp - 1.5*node_data->p()*kb*node_data->Tp())/SolverSpecify::dt*fvm_node->volume();
        MatSetValue(*jac, fvm_node->global_offset()+node_Tp_offset,  fvm_node->global_offset()+node_p_offset,  dWpdt.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, fvm_node->global_offset()+node_Tp_offset,  fvm_node->global_offset()+node_Tp_offset, dWpdt.getADValue(1), ADD_VALUES);
      }

    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}




void SemiconductorSimulationRegion::EBM3_Update_Solution(PetscScalar *lxx)
{
  unsigned int node_psi_offset = ebm_variable_offset(POTENTIAL);
  unsigned int node_n_offset   = ebm_variable_offset(ELECTRON);
  unsigned int node_p_offset   = ebm_variable_offset(HOLE);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);
  unsigned int node_Tn_offset  = ebm_variable_offset(E_TEMP);
  unsigned int node_Tp_offset  = ebm_variable_offset(H_TEMP);


  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;

    PetscScalar V          = lxx[fvm_node->local_offset()+node_psi_offset];
    PetscScalar n          = lxx[fvm_node->local_offset()+node_n_offset];
    PetscScalar p          = lxx[fvm_node->local_offset()+node_p_offset];
    PetscScalar T          = T_external();

    FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data!=NULL);

    //update psi
    node_data->psi_last() =  node_data->psi();
    node_data->psi()      =  V;

    // electron density
    node_data->n_last()   = node_data->n();
    node_data->n()        =  n;

    // hole density
    node_data->p_last()   = node_data->p();
    node_data->p()        =  p;

    // lattice temperature if required
    if(get_advanced_model()->enable_Tl())
    {
      T = lxx[fvm_node->local_offset()+node_Tl_offset];
      node_data->T_last()   = node_data->T();
      node_data->T()        =  T;
    }


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


    // electron temperature if required
    if(get_advanced_model()->enable_Tn())
    {
      node_data->Tn_last()   = node_data->Tn();
      node_data->Tn()        = lxx[fvm_node->local_offset()+node_Tn_offset]/n;
    }

    // hole temperature if required
    if(get_advanced_model()->enable_Tp())
    {
      node_data->Tp_last()   = node_data->Tp();
      node_data->Tp()        = lxx[fvm_node->local_offset()+node_Tp_offset]/p;
    }

    // Update traps
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
