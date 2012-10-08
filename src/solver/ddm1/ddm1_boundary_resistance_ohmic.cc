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


// C++ includes
#include <algorithm>
#include <numeric>

#include "adsmt.h"

// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "resistance_region.h"
#include "insulator_region.h"
#include "boundary_condition_resistance_ohmic.h"
#include "parallel.h"
#include "mathfunc.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


/*---------------------------------------------------------------------
 * set scaling constant
 */
void IF_Metal_OhmicBC::DDM1_Fill_Value(Vec , Vec L)
{

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      switch ( region->type() )
      {
      case SemiconductorRegion:
        {
          const FVM_Node * fvm_node = ( *rnode_it ).second.second;
          VecSetValue(L, fvm_node->global_offset()+0, 1.0, INSERT_VALUES);
          VecSetValue(L, fvm_node->global_offset()+1, 1.0, INSERT_VALUES);
          VecSetValue(L, fvm_node->global_offset()+2, 1.0, INSERT_VALUES);
          break;
        }

      case InsulatorRegion:
        {
          const FVM_Node * fvm_node = ( *rnode_it ).second.second;
          VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
          break;
        }
      default: break;
      }
    }
  }

}

///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


void IF_Metal_OhmicBC::DDM1_Function( PetscScalar * x, Vec f, InsertMode &add_value_flag )
{
  if( this->_infinity_recombination() )
    return _DDM1_Function_Infinite_Recombination(x, f, add_value_flag);
  else
    return _DDM1_Function_Limited_Recombination(x, f, add_value_flag);
}


void IF_Metal_OhmicBC::DDM1_Jacobian( PetscScalar * x, Mat *jac, InsertMode &add_value_flag )
{
  if( this->_infinity_recombination() )
    return _DDM1_Jacobian_Infinite_Recombination(x, jac, add_value_flag);
  else
    return _DDM1_Jacobian_Limited_Recombination(x, jac, add_value_flag);
}



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void IF_Metal_OhmicBC::DDM1_Function_Preprocess(PetscScalar *x, Vec f, std::vector<PetscInt> &src_row,  std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  if( !this->_infinity_recombination() )
  {
    BoundaryCondition::const_node_iterator node_it = nodes_begin();
    BoundaryCondition::const_node_iterator end_it = nodes_end();
    for ( ; node_it!=end_it; ++node_it )
    {
      // skip node not belongs to this processor
      if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

      const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );

      clear_row.push_back ( semiconductor_node->global_offset() );

      if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
      {
        BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
        BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
        for ( ; rnode_it!=end_rnode_it; ++rnode_it )
        {
          const SimulationRegion * region = ( *rnode_it ).second.first;
          if ( region->type() != InsulatorRegion ) continue;
          const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
          clear_row.push_back ( insulator_node->global_offset() );
        }
      }
    }
  }



  if( this->_infinity_recombination() )
  {
    this->_current_buffer.clear();

    BoundaryCondition::const_node_iterator node_it = nodes_begin();
    BoundaryCondition::const_node_iterator end_it = nodes_end();
    for ( ; node_it!=end_it; ++node_it )
    {
      // skip node not belongs to this processor
      if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

      const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
      const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

      clear_row.push_back ( semiconductor_node->global_offset()+0 );
      clear_row.push_back ( semiconductor_node->global_offset()+1 );
      clear_row.push_back ( semiconductor_node->global_offset()+2 );

      // for conduction current
      PetscScalar inject_current = 0;
      {
        PetscInt    ix[2] = {semiconductor_node->global_offset()+1, semiconductor_node->global_offset()+2};
        // I={In, Ip} the electron and hole current flow into this boundary cell.
        // NOTE: although In has dn/dt and R items, they are zero since n is const and n=n0 holds
        // so does Ip
        PetscScalar I[2];

        VecGetValues(f, 2, ix, I);

        // the current = In - Ip;
        // inject_current is the current from semiconductor region to metal region, thus has a negative sign
        inject_current = (-(I[0] - I[1]));
      }

      // displacement current in semiconductor region
      if(SolverSpecify::TimeDependent == true)
      {
        const PetscScalar V_semiconductor  = x[semiconductor_node->local_offset() ];

        PetscScalar I_displacement = 0.0;
        FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
        for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
        {
          const FVM_Node *nb_node = (*nb_it).first;
          const FVM_NodeData * nb_node_data = nb_node->node_data();
          // the psi of neighbor node
          PetscScalar V_nb = x[nb_node->local_offset() +0];
          // distance from nb node to this node
          PetscScalar distance = semiconductor_node->distance ( nb_node );
          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = semiconductor_node->cv_surface_area ( nb_node );
          PetscScalar dEdt;
          if ( SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false ) //second order
          {
            PetscScalar r = SolverSpecify::dt_last/ ( SolverSpecify::dt_last + SolverSpecify::dt );
            dEdt = ( ( 2-r ) / ( 1-r ) * ( V_semiconductor-V_nb )
                     - 1.0/ ( r* ( 1-r ) ) * ( semiconductor_node_data->psi()-nb_node_data->psi() )
                     + ( 1-r ) /r* ( semiconductor_node_data->psi_last()-nb_node_data->psi_last() ) ) /distance/ ( SolverSpecify::dt_last+SolverSpecify::dt );
          }
          else//first order
          {
            dEdt = ( ( V_semiconductor-V_nb )- ( semiconductor_node_data->psi()-nb_node_data->psi() ) ) /distance/SolverSpecify::dt;
          }

          I_displacement += cv_boundary*semiconductor_node_data->eps() *dEdt;
        }
        inject_current -= I_displacement;
      }

      this->_current_buffer.push_back(inject_current);

      if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
      {
        BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
        BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
        for ( ; rnode_it!=end_rnode_it; ++rnode_it )
        {
          const SimulationRegion * region = ( *rnode_it ).second.first;
          if ( region->type() != InsulatorRegion ) continue;
          const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
          clear_row.push_back ( insulator_node->global_offset() );
        }
      }
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void IF_Metal_OhmicBC::_DDM1_Function_Limited_Recombination ( PetscScalar * x, Vec f, InsertMode &add_value_flag )
{
  // Ohmic boundary condition is processed here.

  const PetscScalar T = T_external();
  const PetscScalar eRecombVelocity = scalar("elec.recomb.velocity");
  const PetscScalar hRecombVelocity = scalar("hole.recomb.velocity");

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // data buffer for mesh nodes
  std::vector<PetscInt> iy;
  iy.reserve(3*n_nodes());
  std::vector<PetscScalar> y;
  y.reserve(3*n_nodes());

  std::vector<PetscScalar> current_buffer;
  current_buffer.reserve(n_nodes());

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  const PetscScalar current_scale = this->z_width();

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *> ( _r2 );

  genius_assert(semiconductor_region);
  genius_assert(resistance_region);

  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data = resistance_node->node_data();


    const PetscScalar V_resistance = x[resistance_node->local_offset() ];
    const PetscScalar V_semiconductor  = x[semiconductor_node->local_offset() ];
    const PetscScalar n = x[semiconductor_node->local_offset() +1];
    const PetscScalar p = x[semiconductor_node->local_offset() +2];

    // process semiconductor region

    // mapping this node to material library
    semiconductor_region->material()->mapping ( semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock );

    const PetscScalar nie = semiconductor_region->material()->band->nie ( p, n, T );
    const PetscScalar Nc  = semiconductor_region->material()->band->Nc ( T );
    const PetscScalar Nv  = semiconductor_region->material()->band->Nv ( T );
    const PetscScalar Eg  = semiconductor_region->material()->band->Eg ( T );

    PetscScalar  electron_density;
    PetscScalar  hole_density;
    const PetscScalar  net_dpoing = semiconductor_node_data->Net_doping();
    if ( net_dpoing <0 )                  //p-type
    {
      hole_density = ( -net_dpoing + sqrt ( net_dpoing*net_dpoing + 4*nie*nie ) ) /2.0;
      electron_density = nie*nie/hole_density;
    }
    else                                  //n-type
    {
      electron_density = ( net_dpoing + sqrt ( net_dpoing*net_dpoing + 4*nie*nie ) ) /2.0;
      hole_density = nie*nie/electron_density;
    }

    //governing equation for psi
    PetscScalar f_psi =  V_semiconductor - kb*T/e*boost::math::asinh(semiconductor_node_data->Net_doping()/(2*nie))
                         + Eg/(2*e)
                         + kb*T*log(Nc/Nv)/(2*e)
                         + semiconductor_node_data->affinity()/e
                         - (V_resistance + resistance_node_data->affinity()/e) ;
    VecSetValue ( f, semiconductor_node->global_offset(), f_psi, ADD_VALUES );


    // conservation equation of electron/hole
    PetscScalar S  = semiconductor_node->outside_boundary_surface_area();
    PetscScalar In = - eRecombVelocity * ( n-electron_density ) *S; // electron emit to resistance region
    PetscScalar Ip =   hRecombVelocity * ( p-hole_density ) *S; // hole emit to resistance region
    iy.push_back ( semiconductor_node->global_offset() +1 );
    y.push_back ( In );
    iy.push_back ( semiconductor_node->global_offset() +2 );
    y.push_back ( -Ip );

    // very confusing code, I will explain it here
    // current inject to resistance region can be expressed as I2c + I2d,
    // here I2c=In+Ip is the conduct current from semiconductor region to resistance region
    // and I2d is the displacement current from semiconductor region to resistance region
    // But, I don't know the value of I2d here!
    // As a result, I get inject current by current continuation law, that is I1c + I1d + I2c + I2d = 0
    // I1c is the current from neighbor semiconductor node flow into this cell
    // I1d is the displacement current from from neighbor semiconductor node flow into this cell
    // from electron continuation equation: dn/dt + I1c +I2c - R =0
    // we can got I1c here
    PetscScalar inject_current = In+Ip;

    // compute I1c, see above
    if(SolverSpecify::TimeDependent == true)
    {
      //second order
      if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false)
      {
        PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
        PetscScalar Tn = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*semiconductor_node_data->n() + (1-r)/r*semiconductor_node_data->n_last())
                         / (SolverSpecify::dt_last+SolverSpecify::dt) * semiconductor_node->volume();
        PetscScalar Tp = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*semiconductor_node_data->p() + (1-r)/r*semiconductor_node_data->p_last())
                         / (SolverSpecify::dt_last+SolverSpecify::dt) * semiconductor_node->volume();
        inject_current += Tn;
        inject_current += Tp;
      }
      else //first order
      {
        PetscScalar Tn = -(n - semiconductor_node_data->n())/SolverSpecify::dt*semiconductor_node->volume();
        PetscScalar Tp = -(p - semiconductor_node_data->p())/SolverSpecify::dt*semiconductor_node->volume();
        inject_current += Tn;
        inject_current += Tp;
      }
    }

    // displacement current in semiconductor region
    if(SolverSpecify::TimeDependent == true)
    {
      PetscScalar I_displacement = 0.0;
      FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
      for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
      {
        const FVM_Node *nb_node = (*nb_it).first;
        const FVM_NodeData * nb_node_data = nb_node->node_data();
        // the psi of neighbor node
        PetscScalar V_nb = x[nb_node->local_offset() +0];
        // distance from nb node to this node
        PetscScalar distance = semiconductor_node->distance ( nb_node );
        // area of out surface of control volume related with neighbor node
        PetscScalar cv_boundary = semiconductor_node->cv_surface_area ( nb_node );
        PetscScalar dEdt;
        if ( SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false ) //second order
        {
          PetscScalar r = SolverSpecify::dt_last/ ( SolverSpecify::dt_last + SolverSpecify::dt );
          dEdt = ( ( 2-r ) / ( 1-r ) * ( V_semiconductor-V_nb )
                   - 1.0/ ( r* ( 1-r ) ) * ( semiconductor_node_data->psi()-nb_node_data->psi() )
                   + ( 1-r ) /r* ( semiconductor_node_data->psi_last()-nb_node_data->psi_last() ) ) /distance/ ( SolverSpecify::dt_last+SolverSpecify::dt );
        }
        else//first order
        {
          dEdt = ( ( V_semiconductor-V_nb )- ( semiconductor_node_data->psi()-nb_node_data->psi() ) ) /distance/SolverSpecify::dt;
        }

        I_displacement += cv_boundary*semiconductor_node_data->eps() *dEdt;
      }
      inject_current -= I_displacement;
    }

    // process resistance region, the equation is \sigma J = 0
    iy.push_back ( resistance_node->global_offset() );
    y.push_back ( inject_current );//current flow into resistance region

    current_buffer.push_back ( inject_current );//current flow into resistance region

    // if we have insulator node?
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;

        const FVM_Node * insulator_node = ( *rnode_it ).second.second;
        const PetscScalar V_insulator  = x[insulator_node->local_offset() ];
        PetscScalar f_phi =  V_insulator - V_resistance;
        VecSetValue ( f, insulator_node->global_offset(), f_phi, ADD_VALUES );
      }
    }
  }

  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y[0]), ADD_VALUES);

  // for get the current, we must sum all the terms in current_buffer
  this->current() = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



void IF_Metal_OhmicBC::_DDM1_Function_Infinite_Recombination ( PetscScalar * x, Vec f, InsertMode &add_value_flag )
{
  // Ohmic boundary condition is processed here.

  const PetscScalar T = T_external();

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  std::vector<PetscScalar> current_buffer;
  current_buffer.reserve(n_nodes());

  // data buffer for mesh nodes
  std::vector<PetscInt> iy;
  iy.reserve(3*n_nodes());
  std::vector<PetscScalar> y;
  y.reserve(3*n_nodes());

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  const PetscScalar current_scale = this->z_width();

  const PetscScalar distributed_res = this->scalar("contact.resistance");

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *> ( _r2 );

  genius_assert(semiconductor_region);
  genius_assert(resistance_region);

  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for (unsigned int i=0 ; node_it!=end_it; ++node_it,++i )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data = resistance_node->node_data();


    const PetscScalar V_resistance = x[resistance_node->local_offset() ];
    const PetscScalar V_semiconductor  = x[semiconductor_node->local_offset() ];
    const PetscScalar n = x[semiconductor_node->local_offset() +1];
    const PetscScalar p = x[semiconductor_node->local_offset() +2];

    // process semiconductor region

    // mapping this node to material library
    semiconductor_region->material()->mapping ( semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock );

    const PetscScalar nie = semiconductor_region->material()->band->nie ( p, n, T );
    const PetscScalar Nc  = semiconductor_region->material()->band->Nc ( T );
    const PetscScalar Nv  = semiconductor_region->material()->band->Nv ( T );
    const PetscScalar Eg  = semiconductor_region->material()->band->Eg ( T );

    // here we should calculate current flow into this cell
    const PetscScalar inject_current = this->_current_buffer[i];//pre-computed

    //potential drop in distributed contact resistance
    const PetscScalar cv_boundary = std::abs(semiconductor_node->outside_boundary_surface_area())+1e-10;
    const PetscScalar dV = -inject_current*distributed_res/(cv_boundary*current_scale);// dV = V_resistance - V_semiconductor

    //governing equation for Ohmic contact boundary
    if(semiconductor_region->get_advanced_model()->Fermi) //Fermi
    {
      PetscScalar Ec =  -(e*V_semiconductor + semiconductor_node_data->affinity());
      PetscScalar Ev =  -(e*V_semiconductor + semiconductor_node_data->affinity()+ Eg);

      // the quasi-fermi potential equals to electrode Vapp
      PetscScalar phin = V_resistance - dV + resistance_node_data->affinity()/e;
      PetscScalar phip = V_resistance - dV + resistance_node_data->affinity()/e;

      PetscScalar etan = (-e*phin-Ec)/kb/T;
      PetscScalar etap = (Ev+e*phip)/kb/T;

      y.push_back( Nc*fermi_half(etan) - Nv*fermi_half(etap) -semiconductor_node_data->Net_doping() );
      y.push_back( n - Nc*fermi_half(etan) );
      y.push_back( p - Nv*fermi_half(etap) );

    }
    else     //Boltzmann
    {
      //governing equation for psi
      PetscScalar f_psi =  V_semiconductor - kb*T/e*boost::math::asinh(semiconductor_node_data->Net_doping()/(2*nie))
                           + Eg/(2*e)
                           + kb*T*log(Nc/Nv)/(2*e)
                           + semiconductor_node_data->affinity()/e
                           - (V_resistance - dV + resistance_node_data->affinity()/e) ;
      y.push_back(f_psi);

      PetscScalar  electron_density;
      PetscScalar  hole_density;
      PetscScalar  net_dpoing = semiconductor_node_data->Net_doping();
      if( net_dpoing <0 )                   //p-type
      {
        hole_density = (-net_dpoing + sqrt(net_dpoing*net_dpoing + 4*nie*nie))/2.0;
        electron_density = nie*nie/hole_density;
      }
      else                                  //n-type
      {
        electron_density = (net_dpoing + sqrt(net_dpoing*net_dpoing + 4*nie*nie))/2.0;
        hole_density = nie*nie/electron_density;
      }
      y.push_back( n - electron_density );  //governing equation for electron density
      y.push_back( p - hole_density );      //governing equation for hole density
    }

    // save insert position
    iy.push_back(semiconductor_node->global_offset()+0);
    iy.push_back(semiconductor_node->global_offset()+1);
    iy.push_back(semiconductor_node->global_offset()+2);


    // process resistance region, the equation is \sigma J = 0
    VecSetValue(f, resistance_node->global_offset(),  inject_current , ADD_VALUES);
    current_buffer.push_back( inject_current );

    // if we have insulator node?
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;

        const FVM_Node * insulator_node = ( *rnode_it ).second.second;
        const PetscScalar V_insulator  = x[insulator_node->local_offset() ];
        PetscScalar f_phi =  V_insulator - V_resistance;
        y.push_back( f_phi );
        iy.push_back(insulator_node->global_offset());
      }
    }
  }

  // for get the current, we must sum all the terms in current_buffer
  this->current() = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  if(iy.size())  VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);
  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void IF_Metal_OhmicBC::DDM1_Jacobian_Preprocess(PetscScalar *x, Mat *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  if( !this->_infinity_recombination() )
  {
    BoundaryCondition::const_node_iterator node_it = nodes_begin();
    BoundaryCondition::const_node_iterator end_it = nodes_end();
    for ( ; node_it!=end_it; ++node_it )
    {
      // skip node not belongs to this processor
      if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

      const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );

      clear_row.push_back ( semiconductor_node->global_offset() );

      if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
      {
        BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
        BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
        for ( ; rnode_it!=end_rnode_it; ++rnode_it )
        {
          const SimulationRegion * region = ( *rnode_it ).second.first;
          if ( region->type() != InsulatorRegion ) continue;
          const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
          clear_row.push_back ( insulator_node->global_offset() );
        }
      }
    }
  }



  if( this->_infinity_recombination() )
  {
    _buffer_rows.clear();
    _buffer_cols.clear();
    _buffer_jacobian_entries.clear();

    // search and process all the boundary nodes
    BoundaryCondition::const_node_iterator node_it;
    BoundaryCondition::const_node_iterator end_it = nodes_end();
    for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // only process nodes belong to this processor
      if( (*node_it)->processor_id() != Genius::processor_id() ) continue;

      std::map<PetscInt, PetscScalar>    J_col;

      // get the derivative of electrode current to ohmic node
      const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
      const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();
      const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );

      PetscScalar A1[3], A2[3];
      PetscInt row[3];
      row[0] = semiconductor_node->global_offset()+0;
      row[1] = semiconductor_node->global_offset()+1;
      row[2] = semiconductor_node->global_offset()+2;

      //NOTE MatGetValues only get value from local block!
      MatGetValues(*jac, 1, &row[1], 3, &row[0], &A1[0]);
      MatGetValues(*jac, 1, &row[2], 3, &row[0], &A2[0]);

      J_col[row[0]] += (-(A1[0]-A2[0]));
      J_col[row[1]] += (-(A1[1]-A2[1]));
      J_col[row[2]] += (-(A1[2]-A2[2]));

      //col_vec.insert(col_vec.end(), row.begin(), row.end());
      //J_vec.insert(J_vec.begin(), JM.begin(), JM.end());

      // get the derivative of electrode current to neighbors of ohmic node
      // NOTE neighbors and ohmic bc node may on different processor!
      FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
      FVM_Node::fvm_neighbor_node_iterator nb_it_end = semiconductor_node->neighbor_node_end();
      for(; nb_it != nb_it_end; ++nb_it)
      {
        const FVM_Node *  semiconductor_nb_node = (*nb_it).first;

        PetscInt  col[3];
        col[0] = semiconductor_nb_node->global_offset()+0;
        col[1] = semiconductor_nb_node->global_offset()+1;
        col[2] = semiconductor_nb_node->global_offset()+2;

        MatGetValues(*jac, 1, &row[1], 3, &col[0], &A1[0]);
        MatGetValues(*jac, 1, &row[2], 3, &col[0], &A2[0]);

        J_col[col[0]] += (-(A1[0]-A2[0]));
        J_col[col[1]] += (-(A1[1]-A2[1]));
        J_col[col[2]] += (-(A1[2]-A2[2]));
      }

      // displacement current
      if(SolverSpecify::TimeDependent == true)
      {
        adtl::AutoDScalar::numdir=2;
        AutoDScalar V_semiconductor  = x[semiconductor_node->local_offset() ]; V_semiconductor.setADValue(0, 1.0);
        FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
        for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
        {
          const FVM_Node *nb_node = (*nb_it).first;
          const FVM_NodeData * nb_node_data = nb_node->node_data();
          // the psi of neighbor node
          AutoDScalar V_nb = x[nb_node->local_offset() +0]; V_nb.setADValue(1, 1.0);
          // distance from nb node to this node
          PetscScalar distance = semiconductor_node->distance ( nb_node );
          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = semiconductor_node->cv_surface_area ( nb_node );
          AutoDScalar dEdt;
          if ( SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false ) //second order
          {
            PetscScalar r = SolverSpecify::dt_last/ ( SolverSpecify::dt_last + SolverSpecify::dt );
            dEdt = ( ( 2-r ) / ( 1-r ) * ( V_semiconductor-V_nb )
                     - 1.0/ ( r* ( 1-r ) ) * ( semiconductor_node_data->psi()-nb_node_data->psi() )
                     + ( 1-r ) /r* ( semiconductor_node_data->psi_last()-nb_node_data->psi_last() ) ) / distance/ ( SolverSpecify::dt_last+SolverSpecify::dt );
          }
          else//first order
          {
            dEdt = ( ( V_semiconductor-V_nb )- ( semiconductor_node_data->psi()-nb_node_data->psi() ) ) /distance/SolverSpecify::dt;
          }

          AutoDScalar I_displacement = cv_boundary*semiconductor_node_data->eps() *dEdt;
          J_col[semiconductor_node->global_offset()] += -I_displacement.getADValue(0);
          J_col[nb_node->global_offset()] += -I_displacement.getADValue(1);
        }
      }

      _buffer_rows.push_back(resistance_node->global_offset());

      std::vector<PetscInt> col_vec;
      std::vector<PetscScalar> J_vec;
      std::map<PetscInt, PetscScalar>::const_iterator it = J_col.begin();
      for(; it != J_col.end(); ++it)
      {
        col_vec.push_back(it->first);
        J_vec.push_back(it->second);
      }
      _buffer_cols.push_back(col_vec);
      _buffer_jacobian_entries.push_back(J_vec);

    }


    // then, we can zero all the rows corresponding to Ohmic bc since they have been filled previously.
    for ( node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // skip node not belongs to this processor
      if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

      const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );

      clear_row.push_back ( semiconductor_node->global_offset()+0);
      clear_row.push_back ( semiconductor_node->global_offset()+1);
      clear_row.push_back ( semiconductor_node->global_offset()+2);

      if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
      {
        BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
        BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
        for ( ; rnode_it!=end_rnode_it; ++rnode_it )
        {
          const SimulationRegion * region = ( *rnode_it ).second.first;
          if ( region->type() != InsulatorRegion ) continue;
          const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
          clear_row.push_back ( insulator_node->global_offset() );
        }
      }
    }

  }
}



/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void IF_Metal_OhmicBC::_DDM1_Jacobian_Limited_Recombination ( PetscScalar * x, Mat *jac, InsertMode &add_value_flag )
{

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  const PetscScalar T = T_external();
  const PetscScalar eRecombVelocity = scalar("elec.recomb.velocity");
  const PetscScalar hRecombVelocity = scalar("hole.recomb.velocity");

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *> ( _r2 );

  adtl::AutoDScalar::numdir=4;
  //synchronize with material database
  semiconductor_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;


    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data = resistance_node->node_data();

    AutoDScalar V_resistance = x[resistance_node->local_offset() ];  V_resistance.setADValue ( 0,1.0 );
    AutoDScalar V_semiconductor  = x[semiconductor_node->local_offset() ]; V_semiconductor.setADValue ( 1,1.0 );
    AutoDScalar n = x[semiconductor_node->local_offset() +1]; n.setADValue ( 2,1.0 );
    AutoDScalar p = x[semiconductor_node->local_offset() +2]; p.setADValue ( 3,1.0 );

    // process semiconductor region

    // mapping this node to material library
    semiconductor_region->material()->mapping ( semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock );

    AutoDScalar nie = semiconductor_region->material()->band->nie ( p, n, T );
    const PetscScalar Nc  = semiconductor_region->material()->band->Nc ( T );
    const PetscScalar Nv  = semiconductor_region->material()->band->Nv ( T );
    const PetscScalar Eg  = semiconductor_region->material()->band->Eg ( T );

    AutoDScalar  electron_density;
    AutoDScalar  hole_density;
    const PetscScalar  net_dpoing = semiconductor_node_data->Net_doping();
    if ( net_dpoing <0 )                  //p-type
    {
      hole_density = ( -net_dpoing + sqrt ( net_dpoing*net_dpoing + 4*nie*nie ) ) /2.0;
      electron_density = nie*nie/hole_density;
    }
    else                                  //n-type
    {
      electron_density = ( net_dpoing + sqrt ( net_dpoing*net_dpoing + 4*nie*nie ) ) /2.0;
      hole_density = nie*nie/electron_density;
    }

    //governing equation for psi
    AutoDScalar f_phi =  V_semiconductor - kb*T/e*adtl::asinh(semiconductor_node_data->Net_doping()/(2*nie))
                         + Eg/(2*e)
                         + kb*T*log(Nc/Nv)/(2*e)
                         + semiconductor_node_data->affinity()/e
                         - (V_resistance + resistance_node_data->affinity()/e) ;
    MatSetValue ( *jac, semiconductor_node->global_offset(), resistance_node->global_offset(), f_phi.getADValue ( 0 ), ADD_VALUES );
    MatSetValue ( *jac, semiconductor_node->global_offset(), semiconductor_node->global_offset(), f_phi.getADValue ( 1 ), ADD_VALUES );


    // conservation equation of electron/hole
    PetscScalar S  = semiconductor_node->outside_boundary_surface_area();
    AutoDScalar In = -eRecombVelocity * ( n-electron_density ) *S; // electron emit to resistance region
    AutoDScalar Ip =  hRecombVelocity * ( p-hole_density ) *S; // hole emit to resistance region
    MatSetValue ( *jac, semiconductor_node->global_offset() +1, semiconductor_node->global_offset() +1, In.getADValue ( 2 ), ADD_VALUES );
    MatSetValue ( *jac, semiconductor_node->global_offset() +2, semiconductor_node->global_offset() +2, -Ip.getADValue ( 3 ), ADD_VALUES );

    // process resistance region

    // electron/hole emit current
    MatSetValue ( *jac, resistance_node->global_offset(), semiconductor_node->global_offset() +1, In.getADValue ( 2 ), ADD_VALUES );
    MatSetValue ( *jac, resistance_node->global_offset(), semiconductor_node->global_offset() +2, Ip.getADValue ( 3 ), ADD_VALUES );

    if(SolverSpecify::TimeDependent == true)
    {
      //second order
      if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false)
      {
        PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
        AutoDScalar Tn = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*semiconductor_node_data->n() + (1-r)/r*semiconductor_node_data->n_last())
                         / (SolverSpecify::dt_last+SolverSpecify::dt)*semiconductor_node->volume();
        AutoDScalar Tp = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*semiconductor_node_data->p() + (1-r)/r*semiconductor_node_data->p_last())
                         / (SolverSpecify::dt_last+SolverSpecify::dt)*semiconductor_node->volume();
        // ADD to Jacobian matrix
        MatSetValue(*jac, resistance_node->global_offset(), semiconductor_node->global_offset() +1, Tn.getADValue(2), ADD_VALUES);
        MatSetValue(*jac, resistance_node->global_offset(), semiconductor_node->global_offset() +2, Tp.getADValue(3), ADD_VALUES);
      }
      else //first order
      {
        AutoDScalar Tn = -(n - semiconductor_node_data->n())/SolverSpecify::dt*semiconductor_node->volume();
        AutoDScalar Tp = -(p - semiconductor_node_data->p())/SolverSpecify::dt*semiconductor_node->volume();
        // ADD to Jacobian matrix
        MatSetValue(*jac, resistance_node->global_offset(), semiconductor_node->global_offset() +1, Tn.getADValue(2), ADD_VALUES);
        MatSetValue(*jac, resistance_node->global_offset(), semiconductor_node->global_offset() +2, Tp.getADValue(3), ADD_VALUES);
      }
    }

    // displacement current
    if(SolverSpecify::TimeDependent == true)
    {
      FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
      for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
      {
        const FVM_Node *nb_node = (*nb_it).first;
        const FVM_NodeData * nb_node_data = nb_node->node_data();
        // the psi of neighbor node
        AutoDScalar V_nb = x[nb_node->local_offset() +0]; V_nb.setADValue ( 2, 1.0 );
        // distance from nb node to this node
        PetscScalar distance = semiconductor_node->distance ( nb_node );
        // area of out surface of control volume related with neighbor node
        PetscScalar cv_boundary = semiconductor_node->cv_surface_area ( nb_node );
        AutoDScalar dEdt;
        if ( SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false ) //second order
        {
          PetscScalar r = SolverSpecify::dt_last/ ( SolverSpecify::dt_last + SolverSpecify::dt );
          dEdt = ( ( 2-r ) / ( 1-r ) * ( V_semiconductor-V_nb )
                   - 1.0/ ( r* ( 1-r ) ) * ( semiconductor_node_data->psi()-nb_node_data->psi() )
                   + ( 1-r ) /r* ( semiconductor_node_data->psi_last()-nb_node_data->psi_last() ) ) /distance/ ( SolverSpecify::dt_last+SolverSpecify::dt );
        }
        else//first order
        {
          dEdt = ( ( V_semiconductor-V_nb )- ( semiconductor_node_data->psi()-nb_node_data->psi() ) ) /distance/SolverSpecify::dt;
        }

        AutoDScalar I_displacement = cv_boundary*semiconductor_node_data->eps() *dEdt;
        MatSetValue ( *jac, resistance_node->global_offset(), semiconductor_node->global_offset(), -I_displacement.getADValue ( 1 ), ADD_VALUES );
        MatSetValue ( *jac, resistance_node->global_offset(), nb_node->global_offset(), -I_displacement.getADValue ( 2 ), ADD_VALUES );
      }
    }


    // if we have insulator node?
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;
        const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
        AutoDScalar V_insulator  = x[insulator_node->local_offset() ]; V_insulator.setADValue ( 1, 1.0 );
        AutoDScalar f_phi =  V_insulator - V_resistance;
        MatSetValue ( *jac, insulator_node->global_offset(), resistance_node->global_offset(), f_phi.getADValue ( 0 ), ADD_VALUES );
        MatSetValue ( *jac, insulator_node->global_offset(), insulator_node->global_offset(), f_phi.getADValue ( 1 ), ADD_VALUES );
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




void IF_Metal_OhmicBC::_DDM1_Jacobian_Infinite_Recombination ( PetscScalar * x, Mat *jac, InsertMode &add_value_flag )
{

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  const PetscScalar T = T_external();
  const PetscScalar distributed_res = this->scalar("contact.resistance");

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  const PetscScalar current_scale = this->z_width();

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *> ( _r2 );


  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for (unsigned int i=0 ; node_it!=end_it; ++node_it,++i)
  {
    assert ( ( *node_it )->processor_id() == Genius::processor_id() );

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data = resistance_node->node_data();

    // the insert position
    PetscInt row[4];
    row[0] = resistance_node->global_offset();
    row[1] = semiconductor_node->global_offset()+0;
    row[2] = semiconductor_node->global_offset()+1;
    row[3] = semiconductor_node->global_offset()+2;

    std::vector<PetscInt> col_pattern;
    col_pattern.push_back(resistance_node->global_offset());
    col_pattern.push_back(semiconductor_node->global_offset()+0);
    col_pattern.push_back(semiconductor_node->global_offset()+1);
    col_pattern.push_back(semiconductor_node->global_offset()+2);
    std::sort(col_pattern.begin(), col_pattern.end());

    AsmtDScalar V_resistance( x[resistance_node->local_offset()] );  V_resistance.setADValue(col_pattern, resistance_node->global_offset(), 1.0 );
    AsmtDScalar V_semiconductor( x[semiconductor_node->local_offset()] ); V_semiconductor.setADValue (col_pattern, semiconductor_node->global_offset()+0,1.0 );
    AsmtDScalar n(x[semiconductor_node->local_offset()+1]); n.setADValue(col_pattern, semiconductor_node->global_offset()+1,1.0 );
    AsmtDScalar p(x[semiconductor_node->local_offset()+2]); p.setADValue(col_pattern, semiconductor_node->global_offset()+2,1.0 );

    // process semiconductor region

    // mapping this node to material library
    semiconductor_region->material()->mapping ( semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock );

    const PetscScalar nie = semiconductor_region->material()->band->nie ( p.getValue(), n.getValue(), T );
    const PetscScalar Nc  = semiconductor_region->material()->band->Nc ( T );
    const PetscScalar Nv  = semiconductor_region->material()->band->Nv ( T );
    const PetscScalar Eg  = semiconductor_region->material()->band->Eg ( T );

    // d(current)/d(independent variables of bd node and its neighbors)
    const std::vector<PetscInt> & inject_current_col = _buffer_cols[i];
    const std::vector<PetscScalar> & inject_current_J = _buffer_jacobian_entries[i];

    // here we should calculate current flow into this cell
    AsmtDScalar inject_current(this->_current_buffer[i], inject_current_col, inject_current_J);//pre-computed

    //potential drop in distributed contact resistance
    const PetscScalar cv_boundary = std::abs(semiconductor_node->outside_boundary_surface_area())+1e-10;
    AsmtDScalar dV = -inject_current*distributed_res/(cv_boundary*current_scale);// dV = V_resistance - V_semiconductor

    //governing equation for Ohmic contact boundary
    if(semiconductor_region->get_advanced_model()->Fermi) //Fermi
    {
      AsmtDScalar Ec =  -(e*V_semiconductor + semiconductor_node_data->affinity());
      AsmtDScalar Ev =  -(e*V_semiconductor + semiconductor_node_data->affinity() + Eg);

      // the quasi-fermi potential equals to electrode Vapp
      AsmtDScalar phin = V_resistance - dV + resistance_node_data->affinity()/e;
      AsmtDScalar phip = V_resistance - dV + resistance_node_data->affinity()/e;

      AsmtDScalar etan = (-e*phin-Ec)/kb/T;
      AsmtDScalar etap = (Ev+e*phip)/kb/T;

      AsmtDScalar f_phi =  Nc*fermi_half(etan) - Nv*fermi_half(etap) - semiconductor_node_data->Net_doping();
      AsmtDScalar f_elec =  n - Nc*fermi_half(etan);
      AsmtDScalar f_hole =  p - Nv*fermi_half(etap);

      // set Jacobian of governing equations
      MatSetValues(*jac, 1, &row[1], f_phi.getADSize(), f_phi.getADIndex(), f_phi.getADValue(), ADD_VALUES);
      MatSetValues(*jac, 1, &row[2], f_elec.getADSize(), f_elec.getADIndex(), f_elec.getADValue(), ADD_VALUES);
      MatSetValues(*jac, 1, &row[3], f_hole.getADSize(), f_hole.getADIndex(), f_hole.getADValue(), ADD_VALUES);
    }
    else //Boltzmann
    {
      AsmtDScalar f_phi = V_semiconductor - kb*T/e*boost::math::asinh(semiconductor_node_data->Net_doping()/(2*nie))
                         + Eg/(2*e)
                         + kb*T*log(Nc/Nv)/(2*e)
                         + semiconductor_node_data->affinity()/e
                         - (V_resistance -dV + resistance_node_data->affinity()/e);

      PetscScalar  electron_density;
      PetscScalar  hole_density;
      PetscScalar  net_dpoing = semiconductor_node_data->Net_doping();
      if( net_dpoing <0 )   //p-type
      {
        hole_density = (-net_dpoing + sqrt(net_dpoing*net_dpoing + 4*nie*nie))/2.0;
        electron_density = nie*nie/hole_density;
      }
      else                               //n-type
      {
        electron_density = (net_dpoing + sqrt(net_dpoing*net_dpoing + 4*nie*nie))/2.0;
        hole_density = nie*nie/electron_density;
      }

      AsmtDScalar f_elec =  n - electron_density;  //governing equation for electron density
      AsmtDScalar f_hole =  p - hole_density;      //governing equation for hole density

      // set Jacobian of governing equations
      MatSetValues(*jac, 1, &row[1], f_phi.getADSize(),  f_phi.getADIndex(),  f_phi.getADValue(),  ADD_VALUES);
      MatSetValues(*jac, 1, &row[2], f_elec.getADSize(), f_elec.getADIndex(), f_elec.getADValue(), ADD_VALUES);
      MatSetValues(*jac, 1, &row[3], f_hole.getADSize(), f_hole.getADIndex(), f_hole.getADValue(), ADD_VALUES);
    }


    MatSetValues(*jac, 1, &row[0], inject_current.getADSize(), inject_current.getADIndex(), inject_current.getADValue(), ADD_VALUES);
    //MatSetValue(*jac, semiconductor_node->global_offset(), semiconductor_node->global_offset(), 1.0e-12, ADD_VALUES);

    // if we have insulator node?
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;
        const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
        MatSetValue ( *jac, insulator_node->global_offset(), resistance_node->global_offset(), -1.0, ADD_VALUES );
        MatSetValue ( *jac, insulator_node->global_offset(), insulator_node->global_offset(), 1.0, ADD_VALUES );
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML1 solver
 */
void IF_Metal_OhmicBC::DDM1_Jacobian_Reserve ( Mat *jac, InsertMode &add_value_flag )
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.


  // since we will use ADD_VALUES operat, check the matrix state.
  if ( ( add_value_flag != ADD_VALUES ) && ( add_value_flag != NOT_SET_VALUES ) )
  {
    MatAssemblyBegin ( *jac, MAT_FLUSH_ASSEMBLY );
    MatAssemblyEnd ( *jac, MAT_FLUSH_ASSEMBLY );
  }

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );

    MatSetValue ( *jac, semiconductor_node->global_offset(), resistance_node->global_offset(), 0, ADD_VALUES );

    // process resistance region
    MatSetValue ( *jac, resistance_node->global_offset(), semiconductor_node->global_offset() +0, 0, ADD_VALUES );
    MatSetValue ( *jac, resistance_node->global_offset(), semiconductor_node->global_offset() +1, 0, ADD_VALUES );
    MatSetValue ( *jac, resistance_node->global_offset(), semiconductor_node->global_offset() +2, 0, ADD_VALUES );

    FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
    for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
    {
      const FVM_Node *nb_node = (*nb_it).first;
      MatSetValue ( *jac, resistance_node->global_offset(), nb_node->global_offset()+0, 0, ADD_VALUES );
      MatSetValue ( *jac, resistance_node->global_offset(), nb_node->global_offset()+1, 0, ADD_VALUES );
      MatSetValue ( *jac, resistance_node->global_offset(), nb_node->global_offset()+2, 0, ADD_VALUES );
    }

    // process insulator region when necessary
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;
        const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
        MatSetValue ( *jac, insulator_node->global_offset(), resistance_node->global_offset(), 0, ADD_VALUES );
      }
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


/*---------------------------------------------------------------------
 * update electrode IV
 */
void IF_Metal_OhmicBC::DDM1_Update_Solution(PetscScalar *x)
{
  Parallel::sum(this->current());

  std::vector<PetscScalar> psi_buffer;
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;
  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;
    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const PetscScalar V_resistance = x[resistance_node->local_offset() ];
    psi_buffer.push_back ( V_resistance );
  }

  // we can only get average psi
  Parallel::allgather(psi_buffer);
  this->psi() = psi_buffer.size() ? std::accumulate(psi_buffer.begin(), psi_buffer.end(), 0.0 )/psi_buffer.size() : 0.0;
}
