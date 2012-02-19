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
#include "solver_specify.h"
#include "hall/hall.h"
#include "log.h"

#include "jflux1.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::um;
using PhysicalUnit::cm;



/*-------------------------------------------------------------------------------------
 * filling solution data from FVM_NodeData into petsc vector of L3 EBM.
 */
void SemiconductorSimulationRegion::HALL_Fill_Value(Vec x, Vec L)
{
  this->DDM1_Fill_Value(x, L);
}



static void build_magnetic_field_matrix(PetscScalar mu, const VectorValue<double> & B, TensorValue<PetscScalar> & M)
{
  // note: mu*B is dimension less!
  PetscScalar a = mu*B(0), b = mu*B(1), c = mu*B(2);
  PetscScalar scale = 1.0/(1+a*a+b*b+c*c);
  M(0,0) = (1+a*a)/scale;  M(0,1)=(a*b-c)/scale;   M(0,2)=(c*a+b)/scale;
  M(1,0) = (c+a*b)/scale;  M(1,1)=(1+b*b)/scale;   M(1,2)=(b*c-a)/scale;
  M(2,0) = (c*a-b)/scale;  M(2,1)=(a+b*c)/scale;   M(2,2)=(1+c*c)/scale;
}

static void build_magnetic_field_matrix(const AutoDScalar &mu, const VectorValue<double> & B, TensorValue<AutoDScalar> & M)
{
  // note: mu*B is dimension less!
  AutoDScalar a = mu*B(0), b = mu*B(1), c = mu*B(2);
  AutoDScalar scale = 1.0/(1+a*a+b*b+c*c);
  M(0,0) = (1+a*a)/scale;  M(0,1)=(a*b-c)/scale;   M(0,2)=(c*a+b)/scale;
  M(1,0) = (c+a*b)/scale;  M(1,1)=(1+b*b)/scale;   M(1,2)=(b*c-a)/scale;
  M(2,0) = (c*a-b)/scale;  M(2,1)=(a+b*c)/scale;   M(2,2)=(1+c*c)/scale;
}

static VectorValue<AutoDScalar> matmult(const TensorValue<double> & M, const VectorValue<AutoDScalar> & v)
{
  VectorValue<AutoDScalar> rst;
  for(unsigned int m=0; m<3; m++)
    for(unsigned int n=0; n<3; n++)
      rst(m) += M(m,n)*v(n);
  return rst;
}


static AutoDScalar dot(const VectorValue<PetscScalar> &a, const VectorValue<AutoDScalar> & b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


static AutoDScalar dot(const VectorValue<AutoDScalar> &a, const VectorValue<PetscScalar> & b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static AutoDScalar dot(const VectorValue<AutoDScalar> &a, const VectorValue<AutoDScalar> & b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static VectorValue<AutoDScalar> cross(const VectorValue<PetscScalar> &a, const VectorValue<AutoDScalar> & b)
{
  return VectorValue<AutoDScalar>(  a[1]*b[2] - a[2]*b[1],
                                    -a[0]*b[2] + a[2]*b[0],
                                    a[0]*b[1] - a[1]*b[0]);
}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void SemiconductorSimulationRegion::HALL_Function(const VectorValue<double> & B, PetscScalar * x, Vec f, InsertMode &add_value_flag)
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

  iy.reserve(3*n_node());
  y.reserve(3*n_node());

  //common used variable
  const PetscScalar T   = T_external();
  const PetscScalar Vt  = kb*T/e;
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
    bool boundary_elem = elem->on_boundary() || elem->on_interface();
    bool truncation =  SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationAlways ||
        (SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationBoundary && boundary_elem) ;

    // build the gradient of psi and fermi potential in this cell.
    // which are the vector of electric field and current density.

    std::vector<PetscScalar> psi_vertex;
    std::vector<PetscScalar> phin_vertex;
    std::vector<PetscScalar> phip_vertex;

    std::vector<PetscScalar> mun_edge; //store all the edge mun
    std::vector<PetscScalar> mup_edge; //store all the edge mup
    std::vector<PetscScalar> Jp_edge; //store all the edge Jp
    std::vector<PetscScalar> Jn_edge; //store all the edge Jn
    std::vector<PetscScalar> Vp_edge; //store all the edge Vp=Jp/p
    std::vector<PetscScalar> Vn_edge; //store all the edge Vn=Jn/n

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
      for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
      {
        const FVM_Node * fvm_node = elem->get_fvm_node(nd);
        const FVM_NodeData * fvm_node_data = fvm_node->node_data();

        PetscScalar V;  // electrostatic potential
        PetscScalar n;  // electron density
        PetscScalar p;  // hole density

        if(get_advanced_model()->HighFieldMobilitySelfConsistently)
        {
          // use values in the current iteration
          V  =  x[fvm_node->local_offset()+0];
          n  =  fabs(x[fvm_node->local_offset()+1]) + fvm_node_data->ni()*1e-2;
          p  =  fabs(x[fvm_node->local_offset()+2]) + fvm_node_data->ni()*1e-2;
        }
        else
        {
          // n and p will use previous solution value
          V  =  x[fvm_node->local_offset()+0];
          n  =  fvm_node_data->n() + 1.0*std::pow(cm, -3);
          p  =  fvm_node_data->p() + 1.0*std::pow(cm, -3);
        }

        psi_vertex.push_back  ( V );
        //fermi potential
        phin_vertex.push_back ( V - Vt*log(std::abs(n)/fvm_node_data->ni()) );
        phip_vertex.push_back ( V + Vt*log(std::abs(p)/fvm_node_data->ni()) );
      }


      // compute the gradient
      E   = - elem->gradient(psi_vertex);  // E = - grad(psi)
      Jnv = - elem->gradient(phin_vertex); // we only need the direction of Jnv, here Jnv = - gradient of Fn
      Jpv = - elem->gradient(phip_vertex); // Jpv = - gradient of Fp


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
          if(E_neighbor.size()>E_insul.size())
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
      AutoPtr<Elem> edge = elem->build_edge (ne);
      genius_assert(edge->type()==EDGE2);

      std::vector<unsigned int > edge_nodes;
      elem->nodes_on_edge(ne, edge_nodes);

      double length = edge->volume();                           // the length of this edge
      VectorValue<double> dir = (edge->point(1) - edge->point(0)).unit(); // unit direction of the edge

      const FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes[0]);   assert(fvm_n1);  // fvm_node of node1
      const FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes[1]);   assert(fvm_n2);  // fvm_node of node2

      const FVM_NodeData * n1_data = fvm_n1->node_data();  assert(n1_data);            // fvm_node_data of node1
      const FVM_NodeData * n2_data = fvm_n2->node_data();  assert(n2_data);            // fvm_node_data of node2

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

      // build governing equation of DDML1
      {

        //for node 1 of the edge
        mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

        PetscScalar V1   =  x[n1_local_offset+0];                  // electrostatic potential
        PetscScalar n1   =  x[n1_local_offset+1];                  // electron density
        PetscScalar p1   =  x[n1_local_offset+2];                  // hole density

        // NOTE: Here Ec1, Ev1 are not the conduction/valence band energy.
        // They are here for the calculation of effective driving field for electrons and holes
        // They differ from the conduction/valence band energy by the term with log(Nc), which
        // takes care of the change effective DOS.
        // Ec/Ev should not be used except when its difference between two nodes.
        // The same comment applies to Ec2/Ev2.
        PetscScalar Ec1 =  -(e*V1 + n1_data->affinity() + kb*T*log(mt->band->nie(p1, n1, T)));
        PetscScalar Ev1 =  -(e*V1 + n1_data->affinity() - kb*T*log(mt->band->nie(p1, n1, T)));
        if(get_advanced_model()->Fermi)
        {
          Ec1 = Ec1 - e*Vt*log(gamma_f(fabs(n1)/n1_data->Nc()));
          Ev1 = Ev1 + e*Vt*log(gamma_f(fabs(p1)/n1_data->Nv()));
        }

        PetscScalar eps1 =  n1_data->eps();                        // eps
        PetscScalar Eg1  = mt->band->Eg(T);


        //for node 2 of the edge
        mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

        PetscScalar V2   =  x[n2_local_offset+0];                   // electrostatic potential
        PetscScalar n2   =  x[n2_local_offset+1];                   // electron density
        PetscScalar p2   =  x[n2_local_offset+2];                   // hole density

        PetscScalar Ec2 =  -(e*V2 + n2_data->affinity() + kb*T*log(mt->band->nie(p2, n2, T)));
        PetscScalar Ev2 =  -(e*V2 + n2_data->affinity() - kb*T*log(mt->band->nie(p2, n2, T)));
        if(get_advanced_model()->Fermi)
        {
          Ec2 = Ec2 - e*Vt*log(gamma_f(fabs(n2)/n2_data->Nc()));
          Ev2 = Ev2 + e*Vt*log(gamma_f(fabs(p2)/n2_data->Nv()));
        }

        PetscScalar eps2 =  n2_data->eps();                         // eps
        PetscScalar Eg2  = mt->band->Eg(T);

        PetscScalar mun1;  // electron mobility
        PetscScalar mup1;  // hole mobility
        PetscScalar mun2;   // electron mobility
        PetscScalar mup2;   // hole mobility

        if(highfield_mob)
        {
          // high field mobility
          if (get_advanced_model()->Mob_Force == ModelSpecify::ESimple && !insulator_interface_elem)
          {
            Point dir = (*fvm_n1->root_node() - *fvm_n2->root_node()).unit();
            PetscScalar Ep = fabs((V2-V1)/length);
            PetscScalar Et = 0;
            if(mos_channel_elem)
              Et = (E - dir*(E*dir)).size();

            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T, Ep, Et, T);
            mup1 = mt->mob->HoleMob(p1, n1, T, Ep, Et, T);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T, Ep, Et, T);
            mup2 = mt->mob->HoleMob(p2, n2, T, Ep, Et, T);
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
        else // low field mobility
        {
          mun1 = n1_data->mun();
          mup1 = n1_data->mup();

          mun2 = n2_data->mun();
          mup2 = n2_data->mup();
        }


        PetscScalar mun = 0.5*(mun1+mun2); // the electron mobility at the mid point of the edge, use linear interpolation
        PetscScalar mup = 0.5*(mup1+mup2); // the hole mobility at the mid point of the edge, use linear interpolation
        PetscScalar eps = 0.5*(eps1+eps2); // eps at mid point of the edge

        // S-G current along the edge
        PetscScalar Jn =  mun*In_dd(Vt,(Ec2-Ec1)/e,n1,n2,length);
        PetscScalar Jp =  mup*Ip_dd(Vt,(Ev2-Ev1)/e,p1,p2,length);

        PetscScalar nmid =  nmid_dd(Vt, Ec1/e, Ec2/e, n1, n2);
        PetscScalar pmid =  pmid_dd(Vt, Ev1/e, Ev2/e, p1, p2);

        mun_edge.push_back(mun);
        mup_edge.push_back(mup);
        Jn_edge.push_back(Jn);
        Jp_edge.push_back(Jp);
        Vn_edge.push_back(Jn/nmid);
        Vp_edge.push_back(Jp/pmid);

        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->on_processor() )
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
        }

        // for node 2.
        if( fvm_n2->on_processor() )
        {
          // poisson's equation
          iy.push_back( fvm_n2->global_offset()+0 );
          y.push_back ( eps*(V1 - V2)/length*partial_area );

          // continuity equation of electron
          iy.push_back( fvm_n2->global_offset()+1);
          y.push_back ( - Jn*truncated_partial_area );

          // continuity equation of hole
          iy.push_back( fvm_n2->global_offset()+2);
          y.push_back ( Jp*truncated_partial_area );
        }


        if (get_advanced_model()->BandBandTunneling && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          VectorValue<PetscScalar> E   = - elem->gradient(psi_vertex);  // E = - grad(psi)
          PetscScalar GBTBT1 = mt->band->BB_Tunneling(T, E.size());
          PetscScalar GBTBT2 = mt->band->BB_Tunneling(T, E.size());

          if( fvm_n1->on_processor() )
          {
            // continuity equation
            iy.push_back( fvm_n1->global_offset() + 1);
            y.push_back ( 0.5*GBTBT1*truncated_partial_volume );

            iy.push_back( fvm_n1->global_offset() + 1);
            y.push_back ( 0.5*GBTBT1*truncated_partial_volume );
          }

          if( fvm_n2->on_processor() )
          {
            // continuity equation
            iy.push_back( fvm_n2->global_offset() + 1);
            y.push_back ( 0.5*GBTBT2*truncated_partial_volume );

            iy.push_back( fvm_n2->global_offset() + 1);
            y.push_back ( 0.5*GBTBT2*truncated_partial_volume );
          }
        }

        if (get_advanced_model()->ImpactIonization && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          // consider impact-ionization

          PetscScalar IIn,IIp,GIIn,GIIp;
          PetscScalar Eg;

          Eg = 0.5* ( Eg1 + Eg2 );

          VectorValue<PetscScalar> E   = - elem->gradient(psi_vertex);  // E = - grad(psi)
          VectorValue<PetscScalar> Jnv = - elem->gradient(phin_vertex); // we only need the direction of Jnv, here Jnv = - gradient of Fn
          VectorValue<PetscScalar> Jpv = - elem->gradient(phip_vertex); // Jpv = - gradient of Fpa
          Jnv.add_scaled(VectorValue<PetscScalar>(0, 1e-20, 0), 1.0);
          Jpv.add_scaled(VectorValue<PetscScalar>(0, 1e-20, 0), 1.0);

          VectorValue<PetscScalar> ev = (elem->point(edge_nodes[1]) - elem->point(edge_nodes[0]));
          PetscScalar riin1 = 0.5 + 0.5* (ev.unit()).dot(Jnv.unit());
          PetscScalar riin2 = 1.0 - riin1;
          PetscScalar riip2 = 0.5 + 0.5* (ev.unit()).dot(Jpv.unit());
          PetscScalar riip1 = 1.0 - riip2;

          switch (get_advanced_model()->II_Force)
          {
          case ModelSpecify::IIForce_EdotJ:
            Epn = std::max(E.dot(Jnv.unit()), 0.0);
            Epp = std::max(E.dot(Jpv.unit()), 0.0);
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

          if( fvm_n1->on_processor() )
          {
            // continuity equation
            iy.push_back( fvm_n1->global_offset() + 1);
            y.push_back ( (riin1*GIIn+riip1*GIIp)*truncated_partial_volume );

            iy.push_back( fvm_n1->global_offset() + 2);
            y.push_back ( (riin1*GIIn+riip1*GIIp)*truncated_partial_volume );
          }

          if( fvm_n2->on_processor() )
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

    /*
     * process hall item
     */

    // the average cell electron/hole current density vector
    VectorValue<PetscScalar> vn_cell =  elem->reconstruct_vector(Vn_edge);
    VectorValue<PetscScalar> vp_cell =  elem->reconstruct_vector(Vp_edge);


    PetscScalar mun_cell = std::accumulate(mun_edge.begin(), mun_edge.end(), 0.0 ) / elem->n_edges();
    PetscScalar mup_cell = std::accumulate(mup_edge.begin(), mup_edge.end(), 0.0 ) / elem->n_edges();

    // magnetic field deflection matrix
    TensorValue<double> Mn,Mp;
    build_magnetic_field_matrix(mun_cell*mt->mob->RH_ELEC(), B, Mn);
    build_magnetic_field_matrix(mup_cell*mt->mob->RH_HOLE(), B, Mp);

    // current vector under magnetic field
    VectorValue<PetscScalar> Mvn = Mn*vn_cell;
    VectorValue<PetscScalar> Mvp = Mp*vp_cell;



    for(unsigned int ne=0; ne<elem->n_edges(); ++ne )
    {
      AutoPtr<Elem> edge = elem->build_edge (ne);
      VectorValue<double> dir = (edge->point(1) - edge->point(0)).unit(); // unit direction of the edge

      double length = edge->volume();
      double partial_area = elem->partial_area_with_edge(ne);  // partial area associated with this edge
      double truncated_partial_area =  partial_area;
      if(truncation)
      {
        truncated_partial_area =  this->truncated_partial_area(elem, ne);
      }

      std::vector<unsigned int > edge_nodes;
      elem->nodes_on_edge(ne, edge_nodes);

      FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes[0]);   // fvm_node of node1
      FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes[1]);   // fvm_node of node2

      FVM_NodeData * n1_data = fvm_n1->node_data();  assert(n1_data);            // fvm_node_data of node1
      FVM_NodeData * n2_data = fvm_n2->node_data();  assert(n2_data);            // fvm_node_data of node2

      mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
      PetscScalar V1   =  x[fvm_n1->local_offset()+0];                  // electrostatic potential
      PetscScalar n1   =  x[fvm_n1->local_offset()+1];                  // electron density
      PetscScalar p1   =  x[fvm_n1->local_offset()+2];                  // hole density
      PetscScalar Ec1 =  -(e*V1 + n1_data->affinity() + kb*T*log(mt->band->nie(p1, n1, T)));
      PetscScalar Ev1 =  -(e*V1 + n1_data->affinity() - kb*T*log(mt->band->nie(p1, n1, T)));
      if(get_advanced_model()->Fermi)
      {
        Ec1 = Ec1 - e*Vt*log(gamma_f(fabs(n1)/n1_data->Nc()));
        Ev1 = Ev1 + e*Vt*log(gamma_f(fabs(p1)/n1_data->Nv()));
      }

      mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
      PetscScalar V2   =  x[fvm_n2->local_offset()+0];                   // electrostatic potential
      PetscScalar n2   =  x[fvm_n2->local_offset()+1];                   // electron density
      PetscScalar p2   =  x[fvm_n2->local_offset()+2];                   // hole density
      PetscScalar Ec2 =  -(e*V2 + n2_data->affinity() + kb*T*log(mt->band->nie(p2, n2, T)));
      PetscScalar Ev2 =  -(e*V2 + n2_data->affinity() - kb*T*log(mt->band->nie(p2, n2, T)));
      if(get_advanced_model()->Fermi)
      {
        Ec2 = Ec2 - e*Vt*log(gamma_f(fabs(n2)/n2_data->Nc()));
        Ev2 = Ev2 + e*Vt*log(gamma_f(fabs(p2)/n2_data->Nv()));
      }

      PetscScalar nmid =  nmid_dd(Vt, Ec1/e, Ec2/e, n1, n2);
      PetscScalar pmid =  pmid_dd(Vt, Ev1/e, Ev2/e, p1, p2);

      // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
      if( fvm_n1->on_processor() )
      {
        // continuity equation of electron
        iy.push_back( fvm_n1->global_offset()+1 );
        y.push_back ( - (mun_cell*mt->mob->RH_ELEC())*B*dir.cross(nmid*Mvn)*truncated_partial_area );

        // continuity equation of hole
        iy.push_back( fvm_n1->global_offset()+2 );
        y.push_back ( - (mup_cell*mt->mob->RH_HOLE())*B*dir.cross(pmid*Mvp)*truncated_partial_area );
      }

      // for node 2.
      if( fvm_n2->on_processor() )
      {
        // continuity equation of electron
        iy.push_back( fvm_n2->global_offset()+1);
        y.push_back ( (mun_cell*mt->mob->RH_ELEC())*B*dir.cross(nmid*Mvn)*truncated_partial_area );

        // continuity equation of hole
        iy.push_back( fvm_n2->global_offset()+2);
        y.push_back ( (mup_cell*mt->mob->RH_HOLE())*B*dir.cross(pmid*Mvp)*truncated_partial_area );
      }
    }

    // the average cell electron/hole current density vector
    elem_data->Jn() = - elem->reconstruct_vector(Jn_edge);
    elem_data->Jp() =   elem->reconstruct_vector(Jp_edge);

  }// end of scan all the cell



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

    //PetscScalar V   =  x[fvm_node->local_offset()+0];                       // electrostatic potential
    PetscScalar n   =  x[fvm_node->local_offset()+1];                         // electron density
    PetscScalar p   =  x[fvm_node->local_offset()+2];                         // hole density

    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);      // map this node and its data to material database

    PetscScalar R   = - mt->band->Recomb(p, n, T)*fvm_node->volume();         // the recombination term
    PetscScalar rho = e*(node_data->Net_doping() + p - n)*fvm_node->volume(); // the charge density

    // consider carrier generation
    PetscScalar Field_G = node_data->Field_G()*fvm_node->volume();

    iy.push_back(fvm_node->global_offset()+0);                                // save index in the buffer
    iy.push_back(fvm_node->global_offset()+1);
    iy.push_back(fvm_node->global_offset()+2);
    y.push_back( rho );                                                       // save value in the buffer
    y.push_back( R + Field_G );
    y.push_back( R + Field_G );
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
 * build function and its jacobian for DDML1 solver
 * AD is fully used here
 */
void SemiconductorSimulationRegion::HALL_Jacobian(const VectorValue<double> & B, PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of matrix J
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  //common used variable
  const PetscScalar T   = T_external();
  const PetscScalar Vt  = kb*T/e;
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
    bool boundary_elem = elem->on_boundary() || elem->on_interface();
    bool truncation =  SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationAlways ||
        (SolverSpecify::VoronoiTruncation == SolverSpecify::VoronoiTruncationBoundary && boundary_elem) ;

    //the indepedent variable number, 3*n_nodes
    adtl::AutoDScalar::numdir = 3*elem->n_nodes();

    //synchronize with material database
    mt->set_ad_num(adtl::AutoDScalar::numdir);

    // indicate the column position of the variables in the matrix
    std::vector<PetscInt> cell_col;
    for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
    {
      const FVM_Node * fvm_node = elem->get_fvm_node(nd);
      cell_col.push_back( fvm_node->global_offset()+0 );
      cell_col.push_back( fvm_node->global_offset()+1 );
      cell_col.push_back( fvm_node->global_offset()+2 );
    }

    // first, we build the gradient of psi and fermi potential in this cell.
    // which are the vector of electric field and current density.
    // here use type AutoDScalar, we should make sure the order of independent variable keeps the
    // same all the time
    std::vector<AutoDScalar> psi_vertex;
    std::vector<AutoDScalar> phin_vertex;
    std::vector<AutoDScalar> phip_vertex;

    std::vector<AutoDScalar> mun_edge; //store all the edge mun
    std::vector<AutoDScalar> mup_edge; //store all the edge mup
    std::vector<AutoDScalar> Vp_edge; //store all the edge Jp/p
    std::vector<AutoDScalar> Vn_edge; //store all the edge Jn/n

    // E field parallel to current flow
    AutoDScalar Epn=0;
    AutoDScalar Epp=0;

    // E field vertical to current flow
    AutoDScalar Etn=0;
    AutoDScalar Etp=0;

    // first, we build the gradient of psi and fermi potential in this cell.
    VectorValue<AutoDScalar> E;
    VectorValue<AutoDScalar> Jnv;
    VectorValue<AutoDScalar> Jpv;

    // evaluate E field parallel and vertical to current flow
    if(highfield_mob)
    {
      for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
      {
        const FVM_Node * fvm_node = elem->get_fvm_node(nd);
        const FVM_NodeData * fvm_node_data = fvm_node->node_data();

        AutoDScalar V;               // electrostatic potential
        AutoDScalar n;               // electron density
        AutoDScalar p;               // hole density

        if(get_advanced_model()->HighFieldMobilitySelfConsistently)
        {
          // use values in the current iteration
          V  =  x[fvm_node->local_offset()+0];   V.setADValue(3*nd+0, 1.0);
          n  =  fabs(x[fvm_node->local_offset()+1]) + fvm_node_data->ni()*1e-2;
          p  =  fabs(x[fvm_node->local_offset()+2]) + fvm_node_data->ni()*1e-2;
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
        phin_vertex.push_back ( V - Vt*log(fabs(n)/fvm_node_data->ni()) );
        phip_vertex.push_back ( V + Vt*log(fabs(p)/fvm_node_data->ni()) );
      }

      // compute the gradient
      E   = - elem->gradient(psi_vertex);  // E = - grad(psi)
      Jnv = - elem->gradient(phin_vertex); // we only need the direction of Jnv, here Jnv = - gradient of Fn
      Jpv = - elem->gradient(phip_vertex); // the same as Jnv


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
            V_neighbor.setADValue(3*elem->n_nodes()+nd, 1.0);
            psi_vertex_neighbor.push_back(V_neighbor);
          }
          VectorValue<AutoDScalar> E_neighbor = - elem_neighbor->gradient(psi_vertex_neighbor);
          if(E_neighbor.size()>E_insul.size())
          {
            E_insul = E_neighbor;
            elem_insul = elem_neighbor;
            side_insul = sides[ne];
            region_insul = regions[ne];
          }
        }

        // we need more AD variable
        adtl::AutoDScalar::numdir += 1*elem_insul->n_nodes();
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
      AutoPtr<Elem> edge = elem->build_edge (ne);
      genius_assert(edge->type()==EDGE2);

      std::vector<unsigned int > edge_nodes;
      elem->nodes_on_edge(ne, edge_nodes);

      // the length of this edge
      double length = edge->volume();
      VectorValue<double> dir = (edge->point(1) - edge->point(0)).unit(); // unit direction of the edge

      // fvm_node of node1
      const FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes[0]);
      // fvm_node of node2
      const FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes[1]);

      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data() ;
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data() ;

      // partial area associated with this edge
      double partial_area = elem->partial_area_with_edge(ne);        // partial area associated with this edge
      double partial_volume = elem->partial_volume_with_edge(ne);    // partial volume associated with this edge
      double truncated_partial_area =  partial_area;
      double truncated_partial_volume =  partial_volume;
      if(truncation)
      {
        truncated_partial_area =  this->truncated_partial_area(elem, ne);
        truncated_partial_volume =  elem->partial_volume_with_edge_truncated(ne);
      }

      unsigned int n1_local_offset = fvm_n1->local_offset();
      unsigned int n2_local_offset = fvm_n2->local_offset();

      // the row position of variables in the matrix
      std::vector<PetscInt> row;
      for(unsigned int i=0; i<3; ++i) row.push_back( fvm_n1->global_offset()+i );
      for(unsigned int i=0; i<3; ++i) row.push_back( fvm_n2->global_offset()+i );


      // here we use AD again. Can we hand write it for more efficient?
      {
        //for node 1 of the edge
        mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

        AutoDScalar V1   =  x[n1_local_offset+0];       V1.setADValue(3*edge_nodes[0]+0, 1.0);           // electrostatic potential
        AutoDScalar n1   =  x[n1_local_offset+1];       n1.setADValue(3*edge_nodes[0]+1, 1.0);           // electron density
        AutoDScalar p1   =  x[n1_local_offset+2];       p1.setADValue(3*edge_nodes[0]+2, 1.0);           // hole density

        AutoDScalar Ec1 =  -(e*V1 + n1_data->affinity() + kb*T*log(mt->band->nie(p1, n1, T)) );
        AutoDScalar Ev1 =  -(e*V1 + n1_data->affinity() - kb*T*log(mt->band->nie(p1, n1, T)) );
        if(get_advanced_model()->Fermi)
        {
          Ec1 = Ec1 - e*Vt*log(gamma_f(fabs(n1)/n1_data->Nc()));
          Ev1 = Ev1 + e*Vt*log(gamma_f(fabs(p1)/n1_data->Nv()));
        }
        PetscScalar eps1 =  n1_data->eps();                        // eps
        PetscScalar Eg1  = mt->band->Eg(T);


        //for node 2 of the edge
        mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

        AutoDScalar V2   =  x[n2_local_offset+0];       V2.setADValue(3*edge_nodes[1]+0, 1.0);             // electrostatic potential
        AutoDScalar n2   =  x[n2_local_offset+1];       n2.setADValue(3*edge_nodes[1]+1, 1.0);             // electron density
        AutoDScalar p2   =  x[n2_local_offset+2];       p2.setADValue(3*edge_nodes[1]+2, 1.0);             // hole density

        AutoDScalar Ec2 =  -(e*V2 + n2_data->affinity() + kb*T*log(mt->band->nie(p2, n2, T)) );
        AutoDScalar Ev2 =  -(e*V2 + n2_data->affinity() - kb*T*log(mt->band->nie(p2, n2, T)) );
        if(get_advanced_model()->Fermi)
        {
          Ec2 = Ec2 - e*Vt*log(gamma_f(fabs(n2)/n2_data->Nc()));
          Ev2 = Ev2 + e*Vt*log(gamma_f(fabs(p2)/n2_data->Nv()));
        }
        PetscScalar eps2 =  n2_data->eps();                         // eps

        AutoDScalar mun1;   // electron mobility for node 1 of the edge
        AutoDScalar mup1;   // hole mobility for node 1 of the edge
        AutoDScalar mun2;   // electron mobility  for node 2 of the edge
        AutoDScalar mup2;   // hole mobility for node 2 of the edge

        if(highfield_mob)
        {
          // high field mobility
          if (get_advanced_model()->Mob_Force == ModelSpecify::ESimple && !insulator_interface_elem )
          {
            Point _dir = (*fvm_n1->root_node() - *fvm_n2->root_node()).unit();
            VectorValue<AutoDScalar> dir(_dir(0), _dir(1), _dir(2));
            AutoDScalar Ep = fabs((V2-V1)/length);
            AutoDScalar Et = 0;//
            if(mos_channel_elem)
              Et = (E - (E*dir)*dir).size();

            mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
            mun1 = mt->mob->ElecMob(p1, n1, T, Ep, Et, T);
            mup1 = mt->mob->HoleMob(p1, n1, T, Ep, Et, T);

            mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
            mun2 = mt->mob->ElecMob(p2, n2, T, Ep, Et, T);
            mup2 = mt->mob->HoleMob(p2, n2, T, Ep, Et, T);
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
        else // low field mobility
        {
          mun1 = n1_data->mun();
          mup1 = n1_data->mup();

          mun2 = n2_data->mun();
          mup2 = n2_data->mup();
        }

        AutoDScalar mun = 0.5*(mun1+mun2);  // the electron mobility at the mid point of the edge, use linear interpolation
        AutoDScalar mup = 0.5*(mup1+mup2);  // the hole mobility at the mid point of the edge, use linear interpolation


        PetscScalar eps = 0.5*(eps1+eps2);
        PetscScalar Eg2  = mt->band->Eg(T);

        // S-G current along the edge
        AutoDScalar Jn =  mun*In_dd(Vt,(Ec2-Ec1)/e,n1,n2,length);
        AutoDScalar Jp =  mup*Ip_dd(Vt,(Ev2-Ev1)/e,p1,p2,length);

        AutoDScalar nmid =  nmid_dd(Vt, Ec1/e, Ec2/e, n1, n2);
        AutoDScalar pmid =  pmid_dd(Vt, Ev1/e, Ev2/e, p1, p2);

        mun_edge.push_back(mun);
        mup_edge.push_back(mup);
        Vn_edge.push_back(Jn/nmid);
        Vp_edge.push_back(Jp/pmid);

        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->on_processor() )
        {

          AutoDScalar ff1 = ( eps*(V2 - V1)/length*partial_area );
          AutoDScalar ff2 = (   Jn*truncated_partial_area );
          AutoDScalar ff3 = ( - Jp*truncated_partial_area );

          // general coding always has some overkill... bypass it.
          MatSetValues(*jac, 1, &row[0], cell_col.size(), &cell_col[0], ff1.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[1], cell_col.size(), &cell_col[0], ff2.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[2], cell_col.size(), &cell_col[0], ff3.getADValue(), ADD_VALUES);
        }

        if( fvm_n2->on_processor() )
        {

          AutoDScalar ff1 = ( eps*(V1 - V2)/length*partial_area );
          AutoDScalar ff2 = ( - Jn*truncated_partial_area );
          AutoDScalar ff3 = (   Jp*truncated_partial_area );

          MatSetValues(*jac, 1, &row[3], cell_col.size(), &cell_col[0], ff1.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[4], cell_col.size(), &cell_col[0], ff2.getADValue(), ADD_VALUES);
          MatSetValues(*jac, 1, &row[5], cell_col.size(), &cell_col[0], ff3.getADValue(), ADD_VALUES);
        }

        if (get_advanced_model()->BandBandTunneling && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          VectorValue<AutoDScalar> E   = - elem->gradient(psi_vertex);  // E = - grad(psi)
          AutoDScalar GBTBT1 = mt->band->BB_Tunneling(T, E.size());
          AutoDScalar GBTBT2 = mt->band->BB_Tunneling(T, E.size());

          if( fvm_n1->on_processor() )
          {
            // continuity equation
            AutoDScalar continuity = 0.5*GBTBT1*truncated_partial_volume;
            MatSetValues(*jac, 1, &row[1], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row[2], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
          }

          if( fvm_n2->on_processor() )
          {
            // continuity equation
            AutoDScalar continuity = 0.5*GBTBT2*truncated_partial_volume;
            MatSetValues(*jac, 1, &row[4], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row[5], cell_col.size(), &cell_col[0], continuity.getADValue(), ADD_VALUES);
          }
        }

        if (get_advanced_model()->ImpactIonization && SolverSpecify::Type!=SolverSpecify::EQUILIBRIUM)
        {
          // consider impact-ionization

          AutoDScalar IIn,IIp,GIIn,GIIp;
          AutoDScalar Eg;

          // FIXME should use weighted carrier temperature.
          Eg = 0.5* ( Eg1 + Eg2 );

          VectorValue<AutoDScalar> E   = - elem->gradient(psi_vertex);  // E = - grad(psi)
          VectorValue<AutoDScalar> Jnv = - elem->gradient(phin_vertex); // we only need the direction of Jnv, here Jnv = - gradient of Fn
          VectorValue<AutoDScalar> Jpv = - elem->gradient(phip_vertex); // Jpv = - gradient of Fpa
          Jnv.add_scaled(VectorValue<AutoDScalar>(0, 1e-20, 0), 1.0);
          Jpv.add_scaled(VectorValue<AutoDScalar>(0, 1e-20, 0), 1.0);

          VectorValue<PetscScalar> ev0 = (elem->point(edge_nodes[1]) - elem->point(edge_nodes[0]));
          VectorValue<AutoDScalar> ev;
          ev(0)=ev0(0); ev(1)=ev0(1); ev(2)=ev0(2);
          AutoDScalar riin1 = 0.5 + 0.5* (ev.unit()).dot(Jnv.unit());
          AutoDScalar riin2 = 1.0 - riin1;
          AutoDScalar riip2 = 0.5 + 0.5* (ev.unit()).dot(Jpv.unit());
          AutoDScalar riip1 = 1.0 - riip2;

          switch (get_advanced_model()->II_Force)
          {
          case ModelSpecify::IIForce_EdotJ:
            Epn = adtl::fmax(E.dot(Jnv.unit()), 0.0);
            Epp = adtl::fmax(E.dot(Jpv.unit()), 0.0);
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

          if( fvm_n1->on_processor() )
          {
            // continuity equation
            AutoDScalar electron_continuity = (riin1*GIIn+riip1*GIIp)*truncated_partial_volume ;
            AutoDScalar hole_continuity     = (riin1*GIIn+riip1*GIIp)*truncated_partial_volume ;
            MatSetValues(*jac, 1, &row[1], cell_col.size(), &cell_col[0], electron_continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row[2], cell_col.size(), &cell_col[0], hole_continuity.getADValue(), ADD_VALUES);
          }

          if( fvm_n2->on_processor() )
          {
            // continuity equation
            AutoDScalar electron_continuity = (riin2*GIIn+riip2*GIIp)*truncated_partial_volume ;
            AutoDScalar hole_continuity     = (riin2*GIIn+riip2*GIIp)*truncated_partial_volume ;
            MatSetValues(*jac, 1, &row[4], cell_col.size(), &cell_col[0], electron_continuity.getADValue(), ADD_VALUES);
            MatSetValues(*jac, 1, &row[5], cell_col.size(), &cell_col[0], hole_continuity.getADValue(), ADD_VALUES);
          }
        }

      }
    }// end of scan all edges of the cell

    /*
     * process hall item
     */

    // the average cell electron/hole current density vector
    VectorValue<AutoDScalar> vn_cell =  elem->reconstruct_vector(Vn_edge);
    VectorValue<AutoDScalar> vp_cell =  elem->reconstruct_vector(Vp_edge);

    AutoDScalar mun_cell = std::accumulate(mun_edge.begin(), mun_edge.end(), AutoDScalar(0.0) ) / elem->n_edges();
    AutoDScalar mup_cell = std::accumulate(mup_edge.begin(), mup_edge.end(), AutoDScalar(0.0) ) / elem->n_edges();

    // magnetic field deflection matrix
    TensorValue<AutoDScalar> Mn,Mp;
    build_magnetic_field_matrix(mun_cell*mt->mob->RH_ELEC(), B, Mn);
    build_magnetic_field_matrix(mup_cell*mt->mob->RH_HOLE(), B, Mp);

    // current vector under magnetic field
    VectorValue<AutoDScalar> Mvn = Mn*vn_cell;
    VectorValue<AutoDScalar> Mvp = Mp*vp_cell;

    for(unsigned int ne=0; ne<elem->n_edges(); ++ne )
    {
      AutoPtr<Elem> edge = elem->build_edge (ne);
      VectorValue<double> dir = (edge->point(1) - edge->point(0)).unit(); // unit direction of the edge

      double length = edge->volume();
      double partial_area = elem->partial_area_with_edge(ne);  // partial area associated with this edge
      double truncated_partial_area =  partial_area;
      if(truncation)
      {
        truncated_partial_area =  this->truncated_partial_area(elem, ne);
      }

      std::vector<unsigned int > edge_nodes;
      elem->nodes_on_edge(ne, edge_nodes);

      // fvm_node of node1
      const FVM_Node * fvm_n1 = elem->get_fvm_node(edge_nodes[0]);
      // fvm_node of node2
      const FVM_Node * fvm_n2 = elem->get_fvm_node(edge_nodes[1]);

      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data() ;
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data() ;

      unsigned int n1_local_offset = fvm_n1->local_offset();
      unsigned int n2_local_offset = fvm_n2->local_offset();

      // the row position of variables in the matrix
      std::vector<PetscInt> row;
      for(unsigned int i=0; i<3; ++i) row.push_back( fvm_n1->global_offset()+i );
      for(unsigned int i=0; i<3; ++i) row.push_back( fvm_n2->global_offset()+i );

      //for node 1 of the edge
      mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);

      AutoDScalar V1   =  x[n1_local_offset+0];       V1.setADValue(3*edge_nodes[0]+0, 1.0);           // electrostatic potential
      AutoDScalar n1   =  x[n1_local_offset+1];       n1.setADValue(3*edge_nodes[0]+1, 1.0);           // electron density
      AutoDScalar p1   =  x[n1_local_offset+2];       p1.setADValue(3*edge_nodes[0]+2, 1.0);           // hole density

      AutoDScalar Ec1 =  -(e*V1 + n1_data->affinity() + kb*T*log(mt->band->nie(p1, n1, T)) );
      AutoDScalar Ev1 =  -(e*V1 + n1_data->affinity() - kb*T*log(mt->band->nie(p1, n1, T)) );
      if(get_advanced_model()->Fermi)
      {
        Ec1 = Ec1 - e*Vt*log(gamma_f(fabs(n1)/n1_data->Nc()));
        Ev1 = Ev1 + e*Vt*log(gamma_f(fabs(p1)/n1_data->Nv()));
      }

      //for node 2 of the edge
      mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);

      AutoDScalar V2   =  x[n2_local_offset+0];       V2.setADValue(3*edge_nodes[1]+0, 1.0);             // electrostatic potential
      AutoDScalar n2   =  x[n2_local_offset+1];       n2.setADValue(3*edge_nodes[1]+1, 1.0);             // electron density
      AutoDScalar p2   =  x[n2_local_offset+2];       p2.setADValue(3*edge_nodes[1]+2, 1.0);             // hole density

      AutoDScalar Ec2 =  -(e*V2 + n2_data->affinity() + kb*T*log(mt->band->nie(p2, n2, T)) );
      AutoDScalar Ev2 =  -(e*V2 + n2_data->affinity() - kb*T*log(mt->band->nie(p2, n2, T)) );
      if(get_advanced_model()->Fermi)
      {
        Ec2 = Ec2 - e*Vt*log(gamma_f(fabs(n2)/n2_data->Nc()));
        Ev2 = Ev2 + e*Vt*log(gamma_f(fabs(p2)/n2_data->Nv()));
      }

      AutoDScalar nmid =  nmid_dd(Vt, Ec1/e, Ec2/e, n1, n2);
      AutoDScalar pmid =  pmid_dd(Vt, Ev1/e, Ev2/e, p1, p2);

      // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
      if( fvm_n1->on_processor() )
      {
        AutoDScalar ff2 = - mun_cell*mt->mob->RH_ELEC()*dot(B, cross(dir, nmid*Mvn))*truncated_partial_area;
        AutoDScalar ff3 = - mup_cell*mt->mob->RH_HOLE()*dot(B, cross(dir, pmid*Mvp))*truncated_partial_area;
        MatSetValues(*jac, 1, &row[1], cell_col.size(), &cell_col[0], ff2.getADValue(), ADD_VALUES);
        MatSetValues(*jac, 1, &row[2], cell_col.size(), &cell_col[0], ff3.getADValue(), ADD_VALUES);
      }

      // for node 2.
      if( fvm_n2->on_processor() )
      {
        AutoDScalar ff2 = mun_cell*mt->mob->RH_ELEC()*dot(B, cross(dir, nmid*Mvn))*truncated_partial_area;
        AutoDScalar ff3 = mup_cell*mt->mob->RH_HOLE()*dot(B, cross(dir, pmid*Mvp))*truncated_partial_area;
        MatSetValues(*jac, 1, &row[4], cell_col.size(), &cell_col[0], ff2.getADValue(), ADD_VALUES);
        MatSetValues(*jac, 1, &row[5], cell_col.size(), &cell_col[0], ff3.getADValue(), ADD_VALUES);
      }
    }


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

    const FVM_NodeData * node_data = fvm_node->node_data();

    PetscInt index[3] = {fvm_node->global_offset()+0, fvm_node->global_offset()+1, fvm_node->global_offset()+2};

    AutoDScalar V   =  x[fvm_node->local_offset()+0];   V.setADValue(0, 1.0);              // psi
    AutoDScalar n   =  x[fvm_node->local_offset()+1];   n.setADValue(1, 1.0);              // electron density
    AutoDScalar p   =  x[fvm_node->local_offset()+2];   p.setADValue(2, 1.0);              // hole density

    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);                   // map this node and its data to material database

    AutoDScalar R   = - mt->band->Recomb(p, n, T)*fvm_node->volume();                      // the recombination term
    AutoDScalar rho = e*(node_data->Net_doping() + p - n)*fvm_node->volume();              // the charge density

    // ADD to Jacobian matrix,
    MatSetValues(*jac, 1, &index[0], 3, &index[0], rho.getADValue(), ADD_VALUES);
    MatSetValues(*jac, 1, &index[1], 3, &index[0], R.getADValue(), ADD_VALUES);
    MatSetValues(*jac, 1, &index[2], 3, &index[0], R.getADValue(), ADD_VALUES);
  }


  // boundary condition should be processed later!

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif


}




void SemiconductorSimulationRegion::HALL_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  this->DDM1_Time_Dependent_Function(x, f, add_value_flag);
}


void SemiconductorSimulationRegion::HALL_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  this->DDM1_Time_Dependent_Jacobian(x, jac, add_value_flag);
}


void SemiconductorSimulationRegion::HALL_Function_Hanging_Node(PetscScalar *x, Vec f, InsertMode &add_value_flag)
{
  this->DDM1_Function_Hanging_Node(x, f, add_value_flag);
}


void SemiconductorSimulationRegion::HALL_Jacobian_Hanging_Node(PetscScalar *x, Mat *jac, InsertMode &add_value_flag)
{
  this->DDM1_Jacobian_Hanging_Node(x, jac, add_value_flag);
}


void SemiconductorSimulationRegion::HALL_Update_Solution(PetscScalar *lxx)
{
  this->DDM1_Update_Solution(lxx);
}
