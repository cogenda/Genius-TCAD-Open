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

// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "resistance_region.h"
#include "insulator_region.h"
#include "boundary_condition_resistance_schottky.h"
#include "parallel.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;



///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for EBM3 solver
 */
void IF_Metal_SchottkyBC::EBM3_Function_Preprocess(Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  const SimulationRegion * semiconductor_region = bc_regions().first;
  const SimulationRegion * metal_region = bc_regions().second;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
      // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), semiconductor_region );
    const FVM_Node * metal_node  = get_region_fvm_node ( ( *node_it ), metal_region );

    clear_row.push_back ( semiconductor_node->global_offset() );

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


/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void IF_Metal_SchottkyBC::EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Schottky boundary condition is processed here

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // data buffer for mesh nodes
  std::vector<PetscInt> iy;
  std::vector<PetscScalar> y;

  // the electrode current, since the electrode may be partitioned into several processor,
  // we should collect it.
  std::vector<PetscScalar> current_buffer;

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
  for ( ; node_it!=end_it; ++node_it )
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


    //Schotty Barrier Lowerring
    PetscScalar deltaVB = semiconductor_region->material()->band->SchottyBarrierLowerring(semiconductor_node_data->eps(), semiconductor_node_data->E().size());
    //Schottky current
    PetscScalar S  = semiconductor_node->outside_boundary_surface_area();
    PetscScalar In = semiconductor_region->material()->band->SchottyJsn(n, T_semiconductor, resistance_node_data->affinity() - semiconductor_node_data->affinity() - deltaVB) * S;
    PetscScalar Ip = semiconductor_region->material()->band->SchottyJsp(p, T_semiconductor, resistance_node_data->affinity() - semiconductor_node_data->affinity() + deltaVB) * S;

    PetscScalar f_phi = V_semiconductor - V_resistance - deltaVB;

    // set governing equation to boundary condition of poisson's equation
    VecSetValue(f, semiconductor_node->global_offset()+0, f_phi, ADD_VALUES);

    iy.push_back ( semiconductor_node->global_offset() +1 );
    y.push_back ( In );
    iy.push_back ( semiconductor_node->global_offset() +2 );
    y.push_back ( -Ip );

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

    PetscScalar inject_current = In+Ip;
    // compute Ic
    if(SolverSpecify::TimeDependent == true)
    {
      //second order
      if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false)
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
    PetscScalar I_displacement = 0.0;
    if(SolverSpecify::TimeDependent == true)
    {
      FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
      for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
      {
        const FVM_Node *nb_node = ( *nb_it ).second;
        const FVM_NodeData * nb_node_data = nb_node->node_data();
        // the psi of neighbor node
        PetscScalar V_nb = x[nb_node->local_offset() +0];
        // distance from nb node to this node
        PetscScalar distance = semiconductor_node->distance ( nb_node );
        // area of out surface of control volume related with neighbor node
        PetscScalar cv_boundary = semiconductor_node->cv_surface_area ( nb_node->root_node() );
        PetscScalar dEdt;
        if ( SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false ) //second order
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
    }

    // process resistance region, the equation is \sigma J = 0
    iy.push_back ( resistance_node->global_offset() );
    y.push_back ( inject_current - I_displacement );


    current_buffer.push_back( inject_current - I_displacement );

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

        const PetscScalar T_insulator  = x[insulator_node->local_offset() +1];
        PetscScalar f_T =  T_insulator - T_resistance;
        y.push_back( f_T );
        iy.push_back(insulator_node->global_offset()+1);
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


  // prevent insert zero length vector
  if ( iy.size() )  VecSetValues ( f, iy.size(), &iy[0], &y[0], ADD_VALUES );


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for EBM3 solver
 */
void IF_Metal_SchottkyBC::EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.


  // since we will use ADD_VALUES operat, check the matrix state.
  if ( ( add_value_flag != ADD_VALUES ) && ( add_value_flag != NOT_SET_VALUES ) )
  {
    MatAssemblyBegin ( *jac, MAT_FLUSH_ASSEMBLY );
    MatAssemblyEnd ( *jac, MAT_FLUSH_ASSEMBLY );
  }

  const SimulationRegion * semiconductor_region = bc_regions().first;
  const SimulationRegion * metal_region = bc_regions().second;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), semiconductor_region );
    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), metal_region );
    unsigned int n_node_var_semiconductor  = semiconductor_region->ebm_n_variables();

    MatSetValue ( *jac, semiconductor_node->global_offset()+0, resistance_node->global_offset()+0, 0, ADD_VALUES );
    for(unsigned int n=0; n<n_node_var_semiconductor; ++n)
      MatSetValue ( *jac, resistance_node->global_offset(), semiconductor_node->global_offset() +n, 0, ADD_VALUES );
    FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
    for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
    {
      const FVM_Node *nb_node = ( *nb_it ).second;
      for(unsigned int n=0; n<n_node_var_semiconductor; ++n)
        MatSetValue ( *jac, resistance_node->global_offset(), nb_node->global_offset()+n, 0, ADD_VALUES );
    }

    if(semiconductor_region->get_advanced_model()->enable_Tl())
    {
      unsigned int node_Tl_offset  = semiconductor_region->ebm_variable_offset(TEMPERATURE);
      MatSetValue ( *jac, semiconductor_node->global_offset()+node_Tl_offset, resistance_node->global_offset()+1, 0, ADD_VALUES );
      for(unsigned int n=0; n<n_node_var_semiconductor; ++n)
        MatSetValue ( *jac, resistance_node->global_offset()+1, semiconductor_node->global_offset() +n, 0, ADD_VALUES );
      FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
      for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
      {
        const FVM_Node *nb_node = ( *nb_it ).second;
        for(unsigned int n=0; n<n_node_var_semiconductor; ++n)
          MatSetValue ( *jac, resistance_node->global_offset()+1, nb_node->global_offset()+n, 0, ADD_VALUES );
      }
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
        MatSetValue ( *jac, insulator_node->global_offset()+0, resistance_node->global_offset()+0, 0, ADD_VALUES );
        MatSetValue ( *jac, resistance_node->global_offset()+0, insulator_node->global_offset() +0, 0, ADD_VALUES );
        FVM_Node::fvm_neighbor_node_iterator nb_it = insulator_node->neighbor_node_begin();
        for ( ; nb_it != insulator_node->neighbor_node_end(); ++nb_it )
        {
          const FVM_Node *nb_node = ( *nb_it ).second;
          MatSetValue ( *jac, resistance_node->global_offset()+0, nb_node->global_offset()+0, 0, ADD_VALUES );
        }

        if(region->get_advanced_model()->enable_Tl())
        {
          MatSetValue ( *jac, insulator_node->global_offset()+1, resistance_node->global_offset()+1, 0, ADD_VALUES );
          MatSetValue ( *jac, resistance_node->global_offset()+1, insulator_node->global_offset() +1, 0, ADD_VALUES );
          FVM_Node::fvm_neighbor_node_iterator nb_it = insulator_node->neighbor_node_begin();
          for ( ; nb_it != insulator_node->neighbor_node_end(); ++nb_it )
          {
            const FVM_Node *nb_node = ( *nb_it ).second;
            MatSetValue ( *jac, resistance_node->global_offset()+1, nb_node->global_offset()+1, 0, ADD_VALUES );
          }
        }
      }
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for EBM3 solver
 */
void IF_Metal_SchottkyBC::EBM3_Jacobian_Preprocess(Mat *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  const SimulationRegion * semiconductor_region = bc_regions().first;
  const SimulationRegion * metal_region = bc_regions().second;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
      // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), semiconductor_region );
    const FVM_Node * metal_node  = get_region_fvm_node ( ( *node_it ), metal_region );

    clear_row.push_back ( semiconductor_node->global_offset() );

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


/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void IF_Metal_SchottkyBC::EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *> ( _r2 );

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
    AutoDScalar Nc  = semiconductor_region->material()->band->Nc ( T_semiconductor );
    AutoDScalar Nv  = semiconductor_region->material()->band->Nv ( T_semiconductor );
    AutoDScalar Eg  = semiconductor_region->material()->band->Eg ( T_semiconductor );


    //Schotty Barrier Lowerring
    PetscScalar deltaVB = semiconductor_region->material()->band->SchottyBarrierLowerring(semiconductor_node_data->eps(), semiconductor_node_data->E().size());
    //Schottky current
    PetscScalar S  = semiconductor_node->outside_boundary_surface_area();
    AutoDScalar In = semiconductor_region->material()->band->SchottyJsn(n, T_semiconductor, resistance_node_data->affinity() - semiconductor_node_data->affinity() - deltaVB) * S;
    AutoDScalar Ip = semiconductor_region->material()->band->SchottyJsp(p, T_semiconductor, resistance_node_data->affinity() - semiconductor_node_data->affinity() + deltaVB) * S;

    // the insert position
    PetscInt gloabl_node_phi_offset = semiconductor_node->global_offset()+0;
    PetscInt gloabl_node_n_offset = semiconductor_node->global_offset()+1;
    PetscInt gloabl_node_p_offset = semiconductor_node->global_offset()+2;

    // schottky boundary condition of poisson's equation
    AutoDScalar f_phi = V_semiconductor - V_resistance - deltaVB;
    MatSetValues(*jac, 1, &gloabl_node_phi_offset, col.size(), &col[0], f_phi.getADValue(), ADD_VALUES);

    // process the Jacobian of Schottky current
    MatSetValues(*jac, 1, &gloabl_node_n_offset, col.size(), &col[0], In.getADValue(), ADD_VALUES);
    MatSetValues(*jac, 1, &gloabl_node_p_offset, col.size(), &col[0], ( -Ip ).getADValue(), ADD_VALUES);

    if(semiconductor_region->get_advanced_model()->enable_Tl())
    {
      AutoDScalar f_T =  T_semiconductor - T_resistance;
      PetscInt gloabl_node_Tl_offset  = semiconductor_node->global_offset() + semiconductor_region->ebm_variable_offset(TEMPERATURE);
      MatSetValues(*jac, 1, &gloabl_node_Tl_offset, col.size(), &col[0], f_T.getADValue(), ADD_VALUES);
    }

    // electron temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tn())
    {
      AutoDScalar fTn = (n*(Tn - T_semiconductor));
      PetscInt gloabl_node_Tn_offset  = semiconductor_node->global_offset()+semiconductor_region->ebm_variable_offset(E_TEMP);
      MatSetValues(*jac, 1, &gloabl_node_Tn_offset, col.size(), &col[0], fTn.getADValue(),  ADD_VALUES);
    }

    // hole temperature if required
    if(semiconductor_region->get_advanced_model()->enable_Tp())
    {
      AutoDScalar fTp = (p*(Tp - T_semiconductor));
      PetscInt gloabl_node_Tp_offset  = semiconductor_node->global_offset()+semiconductor_region->ebm_variable_offset(H_TEMP);
      MatSetValues(*jac, 1, &gloabl_node_Tp_offset, col.size(), &col[0], fTp.getADValue(),  ADD_VALUES);
    }

    if(SolverSpecify::TimeDependent == true)
    {
      //second order
      if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false)
      {
        PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
        AutoDScalar Tn = -((2-r)/(1-r)*n - 1.0/(r*(1-r))*semiconductor_node_data->n() + (1-r)/r*semiconductor_node_data->n_last())
                         / (SolverSpecify::dt_last+SolverSpecify::dt)*semiconductor_node->volume();
        AutoDScalar Tp = -((2-r)/(1-r)*p - 1.0/(r*(1-r))*semiconductor_node_data->p() + (1-r)/r*semiconductor_node_data->p_last())
                         / (SolverSpecify::dt_last+SolverSpecify::dt)*semiconductor_node->volume();
        // ADD to Jacobian matrix
        MatSetValue(*jac, resistance_node->global_offset(), semiconductor_node->global_offset() +1, Tn.getADValue(3), ADD_VALUES);
        MatSetValue(*jac, resistance_node->global_offset(), semiconductor_node->global_offset() +2, Tp.getADValue(4), ADD_VALUES);
      }
      else //first order
      {
        AutoDScalar Tn = -(n - semiconductor_node_data->n())/SolverSpecify::dt*semiconductor_node->volume();
        AutoDScalar Tp = -(p - semiconductor_node_data->p())/SolverSpecify::dt*semiconductor_node->volume();
        // ADD to Jacobian matrix
        MatSetValue(*jac, resistance_node->global_offset(), semiconductor_node->global_offset() +1, Tn.getADValue(3), ADD_VALUES);
        MatSetValue(*jac, resistance_node->global_offset(), semiconductor_node->global_offset() +2, Tp.getADValue(4), ADD_VALUES);
      }
    }

    // process resistance region

    // displacement current
    if(SolverSpecify::TimeDependent == true)
    {
      AutoDScalar V_semiconductor  = x[semiconductor_node->local_offset() +0]; V_semiconductor.setADValue ( 0,1.0 );

      FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
      for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
      {
        const FVM_Node *nb_node = ( *nb_it ).second;
        const FVM_NodeData * nb_node_data = nb_node->node_data();
        // the psi of neighbor node
        AutoDScalar V_nb = x[nb_node->local_offset() +0]; V_nb.setADValue ( 1, 1.0 );
        // distance from nb node to this node
        PetscScalar distance = semiconductor_node->distance ( nb_node );
        // area of out surface of control volume related with neighbor node
        PetscScalar cv_boundary = semiconductor_node->cv_surface_area ( nb_node->root_node() );
        AutoDScalar dEdt;
        if ( SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false ) //second order
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
        MatSetValue ( *jac, resistance_node->global_offset(), semiconductor_node->global_offset(), -I_displacement.getADValue ( 0 ), ADD_VALUES );
        MatSetValue ( *jac, resistance_node->global_offset(), nb_node->global_offset(), -I_displacement.getADValue ( 1 ), ADD_VALUES );
      }
    }

    // electron/hole emit current
    {
      PetscInt gloabl_node_phi_offset = resistance_node->global_offset();
      MatSetValues ( *jac, 1, &gloabl_node_phi_offset, col.size(), &col[0], In.getADValue(), ADD_VALUES );
      MatSetValues ( *jac, 1, &gloabl_node_phi_offset, col.size(), &col[0], Ip.getADValue(), ADD_VALUES );
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
        MatSetValue ( *jac, insulator_node->global_offset(), resistance_node->global_offset(), f_phi.getADValue ( 0 ), ADD_VALUES );
        MatSetValue ( *jac, insulator_node->global_offset(), insulator_node->global_offset(), f_phi.getADValue ( 1 ), ADD_VALUES );

        if(region->get_advanced_model()->enable_Tl())
        {
          AutoDScalar T_resistance = x[resistance_node->local_offset() + 1];   T_resistance.setADValue ( 0, 1.0 );
          AutoDScalar T_insulator  = x[insulator_node->local_offset() +1]; T_insulator.setADValue ( 1, 1.0 );
          AutoDScalar f_T =  T_insulator - T_resistance;
          MatSetValue ( *jac, insulator_node->global_offset()+1, resistance_node->global_offset() +1, f_T.getADValue ( 0 ), ADD_VALUES );
          MatSetValue ( *jac, insulator_node->global_offset()+1, insulator_node->global_offset() +1, f_T.getADValue ( 1 ), ADD_VALUES );
        }
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


}

