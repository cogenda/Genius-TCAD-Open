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
#include <numeric>

#include "asinh.hpp" // for asinh

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




///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


void IF_Metal_OhmicBC::EBM3_Function( PetscScalar * x, Vec f, InsertMode &add_value_flag )
{
  return _EBM3_Function_Infinite_Recombination(x, f, add_value_flag);
}


void IF_Metal_OhmicBC::EBM3_Jacobian( PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag )
{
  return _EBM3_Jacobian_Infinite_Recombination(x, jac, add_value_flag);
}



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for EBM3 solver
 */
void IF_Metal_OhmicBC::EBM3_Function_Preprocess(PetscScalar *,Vec f, std::vector<PetscInt> &src_row,  std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  const SimulationRegion * semiconductor_region = bc_regions().first;
  const SimulationRegion * metal_region = bc_regions().second;

  this->_current_buffer.clear();

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), semiconductor_region );
    const FVM_Node * metal_node  = get_region_fvm_node ( ( *node_it ), metal_region );

    clear_row.push_back ( semiconductor_node->global_offset()+0 );
    clear_row.push_back ( semiconductor_node->global_offset()+1 );
    clear_row.push_back ( semiconductor_node->global_offset()+2 );

    if(semiconductor_region->get_advanced_model()->enable_Tl())
    {
      src_row.push_back(semiconductor_node->global_offset()+semiconductor_region->ebm_variable_offset(TEMPERATURE));
      dst_row.push_back(metal_node->global_offset()+metal_region->ebm_variable_offset(TEMPERATURE));
      clear_row.push_back ( semiconductor_node->global_offset()+semiconductor_region->ebm_variable_offset(TEMPERATURE) );
    }

    // electron temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tn())
    {
      clear_row.push_back(semiconductor_node->global_offset() + semiconductor_region->ebm_variable_offset(E_TEMP));
    }

    // hole temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tp())
    {
      clear_row.push_back(semiconductor_node->global_offset() + semiconductor_region->ebm_variable_offset(H_TEMP));
    }

    // for conduction current
    {
      PetscInt    ix[2] = {semiconductor_node->global_offset()+1, semiconductor_node->global_offset()+2};
      // I={In, Ip} the electron and hole current flow into this boundary cell.
      // NOTE: although In has dn/dt and R items, they are zero since n is const and n=n0 holds
      // so does Ip
      PetscScalar I[2];

      VecGetValues(f, 2, ix, I);

      // the current = In - Ip;
      this->_current_buffer.push_back(I[1] - I[0]);
    }

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

        if(region->get_advanced_model()->enable_Tl())
        {
          src_row.push_back(insulator_node->global_offset()+1);
          dst_row.push_back(metal_node->global_offset()+1);
          clear_row.push_back ( insulator_node->global_offset()+1);
        }
      }
    }
  }



}



void IF_Metal_OhmicBC::_EBM3_Function_Infinite_Recombination ( PetscScalar * x, Vec f, InsertMode &add_value_flag )
{
  // Ohmic boundary condition is processed here.

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
  iy.reserve(6*n_nodes());
  std::vector<PetscScalar> y;
  y.reserve(6*n_nodes());


  std::vector<PetscScalar> psi_buffer;
  psi_buffer.reserve(n_nodes());

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
  for (unsigned int i=0 ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data = resistance_node->node_data();


    PetscScalar V_resistance = x[resistance_node->local_offset() +0];
    PetscScalar T_resistance = T_external();;
    PetscScalar V_semiconductor  = x[semiconductor_node->local_offset() +0];
    PetscScalar n                = x[semiconductor_node->local_offset() +1];
    PetscScalar p                = x[semiconductor_node->local_offset() +2];
    PetscScalar T_semiconductor  = T_external();

    // lattice temperature
    if(resistance_region->get_advanced_model()->enable_Tl())
    {
      unsigned int node_Tl_offset  = resistance_region->ebm_variable_offset(TEMPERATURE);
      T_resistance = x[resistance_node->local_offset() + node_Tl_offset];
    }

    if(semiconductor_region->get_advanced_model()->enable_Tl())
    {
      unsigned int node_Tl_offset  = semiconductor_region->ebm_variable_offset(TEMPERATURE);
      T_semiconductor = x[semiconductor_node->local_offset() + node_Tl_offset];
    }

    psi_buffer.push_back ( V_resistance );

    // process semiconductor region

    // mapping this node to material library
    semiconductor_region->material()->mapping ( semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock );

    const PetscScalar nie = semiconductor_region->material()->band->nie ( p, n, T_semiconductor );
    const PetscScalar Nc  = semiconductor_region->material()->band->Nc ( T_semiconductor );
    const PetscScalar Nv  = semiconductor_region->material()->band->Nv ( T_semiconductor );
    const PetscScalar Eg  = semiconductor_region->material()->band->Eg ( T_semiconductor );


    //governing equation for Ohmic contact boundary
    if(semiconductor_region->get_advanced_model()->Fermi) //Fermi
    {
      PetscScalar Ec =  -(e*V_semiconductor + semiconductor_node_data->affinity() );
      PetscScalar Ev =  -(e*V_semiconductor + semiconductor_node_data->affinity() + Eg);

      // the quasi-fermi potential equals to electrode Vapp
      PetscScalar phin = V_resistance + resistance_node_data->affinity()/e;
      PetscScalar phip = V_resistance + resistance_node_data->affinity()/e;

      PetscScalar etan = (-e*phin-Ec)/kb/T_semiconductor;
      PetscScalar etap = (Ev+e*phip)/kb/T_semiconductor;

      y.push_back( Nc*fermi_half(etan) - Nv*fermi_half(etap) -semiconductor_node_data->Net_doping() );
      y.push_back( n - Nc*fermi_half(etan) );
      y.push_back( p - Nv*fermi_half(etap) );

    }
    else     //Boltzmann
    {
      //governing equation for psi
      PetscScalar f_psi =  V_semiconductor - kb*T_semiconductor/e*asinh(semiconductor_node_data->Net_doping()/(2*nie))
                           + Eg/(2*e)
                           + kb*T_semiconductor*log(Nc/Nv)/(2*e)
                           + semiconductor_node_data->affinity()/e
                           - (V_resistance + resistance_node_data->affinity()/e) ;
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

    if(semiconductor_region->get_advanced_model()->enable_Tl())
    {
      unsigned int node_Tl_offset  = semiconductor_region->ebm_variable_offset(TEMPERATURE);
      PetscScalar f_T =  T_semiconductor - T_resistance;
      iy.push_back ( semiconductor_node->global_offset() +node_Tl_offset );
      y.push_back(f_T);
    }

    // electron temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tn())
    {
      unsigned int node_Tn_offset  = semiconductor_region->ebm_variable_offset(E_TEMP);
      PetscScalar Tn = x[semiconductor_node->local_offset() + node_Tn_offset]/n;
      y.push_back(n*(Tn - T_semiconductor));
      iy.push_back(semiconductor_node->global_offset() + node_Tn_offset);
    }

    // hole temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tp())
    {
      unsigned int node_Tp_offset  = semiconductor_region->ebm_variable_offset(H_TEMP);
      PetscScalar Tp = x[semiconductor_node->local_offset() + node_Tp_offset]/p;
      y.push_back(p*(Tp - T_semiconductor));
      iy.push_back(semiconductor_node->global_offset() + node_Tp_offset);
    }


    // here we should calculate current flow into this cell
    PetscScalar inject_current=this->_current_buffer[i++];//pre-computed


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
        const PetscScalar V_insulator  = x[insulator_node->local_offset() +0];
        PetscScalar f_phi =  V_insulator - V_resistance;
        y.push_back( f_phi );
        iy.push_back(insulator_node->global_offset()+0);

        if(region->get_advanced_model()->enable_Tl())
        {
          const PetscScalar T_insulator  = x[insulator_node->local_offset() +1];
          PetscScalar f_T =  T_insulator - T_resistance;
          y.push_back( f_T );
          iy.push_back(insulator_node->global_offset()+1);
        }
      }
    }
  }


  // we first gather the electrode current
  Parallel::allgather(current_buffer);

  // for get the current, we must sum all the terms in current_buffer
  this->current() = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  // we can only get average psi
  Parallel::allgather(psi_buffer);
  this->psi() = psi_buffer.size() ? std::accumulate(psi_buffer.begin(), psi_buffer.end(), 0.0 )/psi_buffer.size() : 0.0;


  if(iy.size())  VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);
  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for EBM3 solver
 */
void IF_Metal_OhmicBC::EBM3_Jacobian_Preprocess(PetscScalar * ,SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,  std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  const SimulationRegion * semiconductor_region = bc_regions().first;
  const SimulationRegion * metal_region = bc_regions().second;


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

    // get the derivative of electrode current to ohmic node
    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), semiconductor_region );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();
    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), metal_region );

    unsigned int n_node_var  = semiconductor_region->ebm_n_variables();
    unsigned int node_n_offset = semiconductor_region->ebm_variable_offset(ELECTRON);
    unsigned int node_p_offset = semiconductor_region->ebm_variable_offset(HOLE);

    std::vector<PetscScalar> A1(n_node_var), A2(n_node_var), JM(n_node_var), JN(n_node_var);
    std::vector<PetscInt>    row(n_node_var);
    for(unsigned int nv=0; nv<n_node_var; ++nv)
      row[nv] = semiconductor_node->global_offset()+nv;

    jac->get_row(  row[node_n_offset],  n_node_var,  &row[0],  &A1[0]);
    jac->get_row(  row[node_p_offset],  n_node_var,  &row[0],  &A2[0]);

    for(unsigned int nv=0; nv<n_node_var; ++nv)
      JM[nv] = -(A1[nv]-A2[nv]);

    // get the derivative of electrode current to neighbors of ohmic node
    // NOTE neighbors and ohmic bc node may on different processor!
    FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
    FVM_Node::fvm_neighbor_node_iterator nb_it_end = semiconductor_node->neighbor_node_end();
    for(; nb_it != nb_it_end; ++nb_it)
    {
      const FVM_Node *  semiconductor_nb_node = (*nb_it).first;

      // distance from nb node to this node
      PetscScalar distance = semiconductor_node->distance( semiconductor_nb_node );
      // area of out surface of control volume related with neighbor node
      PetscScalar cv_boundary = semiconductor_node->cv_surface_area(semiconductor_nb_node);

      std::vector<PetscInt>    col(n_node_var);
      for(unsigned int i=0; i<n_node_var; ++i)
        col[i] = semiconductor_nb_node->global_offset()+i;

      jac->get_row(  row[node_n_offset],  n_node_var,  &col[0],  &A1[0]);
      jac->get_row(  row[node_p_offset],  n_node_var,  &col[0],  &A2[0]);

      for(unsigned int nv=0; nv<n_node_var; ++nv)
        JN[nv] = -(A1[nv]-A2[nv]);

      _buffer_rows.push_back(resistance_node->global_offset());
      _buffer_cols.push_back(col);
      _buffer_jacobian_entries.push_back(JN);
    }
    _buffer_rows.push_back(resistance_node->global_offset());
    _buffer_cols.push_back(row);
    _buffer_jacobian_entries.push_back(JM);
  }

  // then, we can zero all the rows corresponding to Ohmic bc since they have been filled previously.
  for ( node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), semiconductor_region );
    const FVM_Node * metal_node  = get_region_fvm_node ( ( *node_it ), metal_region );

    clear_row.push_back ( semiconductor_node->global_offset()+0 );
    clear_row.push_back ( semiconductor_node->global_offset()+1 );
    clear_row.push_back ( semiconductor_node->global_offset()+2 );

    if(semiconductor_region->get_advanced_model()->enable_Tl())
    {
      src_row.push_back(semiconductor_node->global_offset()+semiconductor_region->ebm_variable_offset(TEMPERATURE));
      dst_row.push_back(metal_node->global_offset()+metal_region->ebm_variable_offset(TEMPERATURE));
      clear_row.push_back ( semiconductor_node->global_offset()+semiconductor_region->ebm_variable_offset(TEMPERATURE) );
    }

    // electron temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tn())
    {
      clear_row.push_back(semiconductor_node->global_offset() + semiconductor_region->ebm_variable_offset(E_TEMP));
    }

    // hole temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tp())
    {
      clear_row.push_back(semiconductor_node->global_offset() + semiconductor_region->ebm_variable_offset(H_TEMP));
    }

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

        if(region->get_advanced_model()->enable_Tl())
        {
          src_row.push_back(insulator_node->global_offset()+1);
          dst_row.push_back(metal_node->global_offset()+1);
          clear_row.push_back ( insulator_node->global_offset()+1);
        }
      }
    }
  }

}



void IF_Metal_OhmicBC::_EBM3_Jacobian_Infinite_Recombination ( PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag )
{

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *> ( _r2 );

  // d(current)/d(independent variables of bd node and its neighbors)
  for(unsigned int n=0; n<_buffer_rows.size(); ++n)
  {
    jac->add_row(  _buffer_rows[n],  _buffer_cols[n].size(),  &(_buffer_cols[n])[0],  &(_buffer_jacobian_entries[n])[0] );
  }

  adtl::AutoDScalar::numdir=8;
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

    PetscInt ad_index = 0;
    std::vector<PetscInt> col;
    AutoDScalar V_resistance = x[resistance_node->local_offset() +0];        V_resistance.setADValue ( ad_index++,1.0 );
    AutoDScalar T_resistance = T_external();
    AutoDScalar V_semiconductor  = x[semiconductor_node->local_offset() +0]; V_semiconductor.setADValue ( ad_index++,1.0 );
    AutoDScalar n                = x[semiconductor_node->local_offset() +1]; n.setADValue ( ad_index++,1.0 );
    AutoDScalar p                = x[semiconductor_node->local_offset() +2]; p.setADValue ( ad_index++,1.0 );
    AutoDScalar T_semiconductor  = T_external();
    AutoDScalar Tn =  T_external();
    AutoDScalar Tp =  T_external();

    col.push_back(resistance_node->global_offset() +0);
    col.push_back(semiconductor_node->global_offset() +0);
    col.push_back(semiconductor_node->global_offset() +1);
    col.push_back(semiconductor_node->global_offset() +2);

    // lattice temperature
    if(resistance_region->get_advanced_model()->enable_Tl())
    {
      unsigned int node_Tl_offset  = resistance_region->ebm_variable_offset(TEMPERATURE);
      T_resistance = x[resistance_node->local_offset() + node_Tl_offset];
      T_resistance.setADValue ( ad_index++, 1.0 );
      col.push_back(resistance_node->global_offset() +node_Tl_offset);
    }

    if(semiconductor_region->get_advanced_model()->enable_Tl())
    {
      unsigned int node_Tl_offset  = semiconductor_region->ebm_variable_offset(TEMPERATURE);
      T_semiconductor = x[semiconductor_node->local_offset() + node_Tl_offset];
      T_semiconductor.setADValue ( ad_index++, 1.0 );
      col.push_back(semiconductor_node->global_offset() +node_Tl_offset);
    }

    // electron temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tn())
    {
      unsigned int node_Tn_offset  = semiconductor_region->ebm_variable_offset(E_TEMP);
      AutoDScalar nTn = x[semiconductor_node->local_offset() + node_Tn_offset];
      nTn.setADValue(ad_index++, 1.0);
      col.push_back(semiconductor_node->global_offset() +node_Tn_offset);
      Tn = nTn/n;
    }

    // hole temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tp())
    {
      unsigned int node_Tp_offset  = semiconductor_region->ebm_variable_offset(H_TEMP);
      AutoDScalar pTp = x[semiconductor_node->local_offset() + node_Tp_offset];
      pTp.setADValue(ad_index++, 1.0);
      col.push_back(semiconductor_node->global_offset() +node_Tp_offset);
      Tp = pTp/p;
    }

    // process semiconductor region

    // mapping this node to material library
    semiconductor_region->material()->mapping ( semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock );

    AutoDScalar nie = semiconductor_region->material()->band->nie ( p, n, T_semiconductor );
    const AutoDScalar Nc  = semiconductor_region->material()->band->Nc ( T_semiconductor );
    const AutoDScalar Nv  = semiconductor_region->material()->band->Nv ( T_semiconductor );
    const AutoDScalar Eg  = semiconductor_region->material()->band->Eg ( T_semiconductor );

    //governing equation for Ohmic contact boundary
    AutoDScalar f_phi,f_elec,f_hole;
    if(semiconductor_region->get_advanced_model()->Fermi) //Fermi
    {
      AutoDScalar Ec =  -(e*V_semiconductor + semiconductor_node_data->affinity() );
      AutoDScalar Ev =  -(e*V_semiconductor + semiconductor_node_data->affinity() + Eg);

      // the quasi-fermi potential equals to electrode Vapp
      AutoDScalar phin = V_resistance + resistance_node_data->affinity()/e;
      AutoDScalar phip = V_resistance + resistance_node_data->affinity()/e;

      AutoDScalar etan = (-e*phin-Ec)/kb/T_semiconductor;
      AutoDScalar etap = (Ev+e*phip)/kb/T_semiconductor;

      f_phi =  Nc*fermi_half(etan) - Nv*fermi_half(etap) - semiconductor_node_data->Net_doping();
      f_elec =  n - Nc*fermi_half(etan);
      f_hole =  p - Nv*fermi_half(etap);
    }
    else //Boltzmann
    {

      f_phi = V_semiconductor - kb*T_semiconductor/e*adtl::asinh(semiconductor_node_data->Net_doping()/(2*nie))
              + Eg/(2*e)
              + kb*T_semiconductor*log(Nc/Nv)/(2*e)
              + semiconductor_node_data->affinity()/e
              - (V_resistance + resistance_node_data->affinity()/e);

      AutoDScalar  electron_density;
      AutoDScalar  hole_density;
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

      f_elec =  n - electron_density;  //governing equation for electron density
      f_hole =  p - hole_density;      //governing equation for hole density
    }



    // the insert position
    PetscInt gloabl_node_phi_offset = semiconductor_node->global_offset()+0;
    PetscInt gloabl_node_n_offset = semiconductor_node->global_offset()+1;
    PetscInt gloabl_node_p_offset = semiconductor_node->global_offset()+2;

    // set Jacobian for governing equations
    jac->add_row(  gloabl_node_phi_offset,  col.size(),  &col[0],  f_phi.getADValue() );
    jac->add_row(  gloabl_node_n_offset,  col.size(),  &col[0],  f_elec.getADValue() );
    jac->add_row(  gloabl_node_p_offset,  col.size(),  &col[0],  f_hole.getADValue() );

    if(semiconductor_region->get_advanced_model()->enable_Tl())
    {
      AutoDScalar f_T =  T_semiconductor - T_resistance;
      PetscInt gloabl_node_Tl_offset  = semiconductor_node->global_offset() + semiconductor_region->ebm_variable_offset(TEMPERATURE);
      jac->add_row(  gloabl_node_Tl_offset,  col.size(),  &col[0],  f_T.getADValue() );
    }

    // electron temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tn())
    {
      AutoDScalar fTn = (n*(Tn - T_semiconductor));
      PetscInt gloabl_node_Tn_offset  = semiconductor_node->global_offset()+semiconductor_region->ebm_variable_offset(E_TEMP);
      jac->add_row(  gloabl_node_Tn_offset,  col.size(),  &col[0],  fTn.getADValue() );
    }

    // hole temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tp())
    {
      AutoDScalar fTp = (p*(Tp - T_semiconductor));
      PetscInt gloabl_node_Tp_offset  = semiconductor_node->global_offset()+semiconductor_region->ebm_variable_offset(H_TEMP);
      jac->add_row(  gloabl_node_Tp_offset,  col.size(),  &col[0],  fTp.getADValue() );
    }

    // displacement current
    if(SolverSpecify::TimeDependent == true)
    {
      AutoDScalar V_semiconductor  = x[semiconductor_node->local_offset() +0]; V_semiconductor.setADValue ( 0,1.0 );

      FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
      for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
      {
        const FVM_Node *nb_node = (*nb_it).first;
        const FVM_NodeData * nb_node_data = nb_node->node_data();
        // the psi of neighbor node
        AutoDScalar V_nb = x[nb_node->local_offset() +0]; V_nb.setADValue ( 1, 1.0 );
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
          jac->add( resistance_node->global_offset(),  semiconductor_node->global_offset(),  -I_displacement.getADValue ( 0 ) );
          jac->add( resistance_node->global_offset(),  nb_node->global_offset(),  -I_displacement.getADValue ( 1 ) );
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
        AutoDScalar V_resistance = x[resistance_node->local_offset() +0]; V_resistance.setADValue ( 0, 1.0 );
        AutoDScalar V_insulator  = x[insulator_node->local_offset() +0];  V_insulator.setADValue ( 1, 1.0 );
        AutoDScalar f_phi =  V_insulator - V_resistance;
          jac->add( insulator_node->global_offset(),  resistance_node->global_offset(),  f_phi.getADValue ( 0 ) );
          jac->add( insulator_node->global_offset(),  insulator_node->global_offset(),  f_phi.getADValue ( 1 ) );

        if(region->get_advanced_model()->enable_Tl())
        {
          AutoDScalar T_resistance = x[resistance_node->local_offset() + 1];   T_resistance.setADValue ( 0, 1.0 );
          AutoDScalar T_insulator  = x[insulator_node->local_offset() +1]; T_insulator.setADValue ( 1, 1.0 );
          AutoDScalar f_T =  T_insulator - T_resistance;
            jac->add( insulator_node->global_offset()+1,  resistance_node->global_offset() +1,  f_T.getADValue ( 0 ) );
            jac->add( insulator_node->global_offset()+1,  insulator_node->global_offset() +1,  f_T.getADValue ( 1 ) );
        }
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




