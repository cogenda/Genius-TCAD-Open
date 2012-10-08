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

/*---------------------------------------------------------------------
 * set scaling constant
 */
void IF_Metal_SchottkyBC::DDM1R_Fill_Value(Vec , Vec L)
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

/*---------------------------------------------------------------------
 * do pre-process to function for DDML1 solver
 */
void IF_Metal_SchottkyBC::DDM1R_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;


    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * resistance_node     = get_region_fvm_node ( ( *node_it ), _r2 );

    src_row.push_back ( semiconductor_node->global_offset() );
    dst_row.push_back ( resistance_node->global_offset() );
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


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void IF_Metal_SchottkyBC::DDM1R_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Schottky boundary condition is processed here

  const PetscScalar T = T_external();

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


    //Schotty Barrier Lowerring
    PetscScalar deltaVB = semiconductor_region->material()->band->SchottyBarrierLowerring(semiconductor_node_data->eps(), semiconductor_node_data->E().size());
    //Schottky current
    PetscScalar S  = semiconductor_node->outside_boundary_surface_area();
    PetscScalar In = semiconductor_region->material()->band->SchottyJsn(n, T, resistance_node_data->affinity() - semiconductor_node_data->affinity() - deltaVB) * S;
    PetscScalar Ip = semiconductor_region->material()->band->SchottyJsp(p, T, resistance_node_data->affinity() - semiconductor_node_data->affinity() + deltaVB) * S;

    PetscScalar f_phi = V_semiconductor - V_resistance - deltaVB;

    // set governing equation to boundary condition of poisson's equation
    VecSetValue(f, semiconductor_node->global_offset()+0, f_phi, ADD_VALUES);

    iy.push_back ( semiconductor_node->global_offset() +1 );
    y.push_back ( In );
    iy.push_back ( semiconductor_node->global_offset() +2 );
    y.push_back ( -Ip );

    PetscScalar inject_current = -(In+Ip);
    // process resistance region, the equation is \sigma J = 0
    iy.push_back ( resistance_node->global_offset()+1 );
    y.push_back ( inject_current );


    current_buffer.push_back( -inject_current );

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

  // for get the current, we must sum all the terms in current_buffer
  this->current() = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  // prevent insert zero length vector
  if ( iy.size() )  VecSetValues ( f, iy.size(), &iy[0], &y[0], ADD_VALUES );


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML1 solver
 */
void IF_Metal_SchottkyBC::DDM1R_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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
    MatSetValue ( *jac, resistance_node->global_offset()+0, semiconductor_node->global_offset() +0, 0, ADD_VALUES );
    MatSetValue ( *jac, resistance_node->global_offset()+0, semiconductor_node->global_offset() +1, 0, ADD_VALUES );
    MatSetValue ( *jac, resistance_node->global_offset()+0, semiconductor_node->global_offset() +2, 0, ADD_VALUES );
    MatSetValue ( *jac, resistance_node->global_offset()+1, semiconductor_node->global_offset() +0, 0, ADD_VALUES );
    MatSetValue ( *jac, resistance_node->global_offset()+1, semiconductor_node->global_offset() +1, 0, ADD_VALUES );
    MatSetValue ( *jac, resistance_node->global_offset()+1, semiconductor_node->global_offset() +2, 0, ADD_VALUES );

    FVM_Node::fvm_neighbor_node_iterator nb_it = semiconductor_node->neighbor_node_begin();
    for ( ; nb_it != semiconductor_node->neighbor_node_end(); ++nb_it )
    {
      const FVM_Node *nb_node = ( *nb_it ).first;
      MatSetValue ( *jac, resistance_node->global_offset()+0, nb_node->global_offset()+0, 0, ADD_VALUES );
      MatSetValue ( *jac, resistance_node->global_offset()+0, nb_node->global_offset()+1, 0, ADD_VALUES );
      MatSetValue ( *jac, resistance_node->global_offset()+0, nb_node->global_offset()+2, 0, ADD_VALUES );
      MatSetValue ( *jac, resistance_node->global_offset()+1, nb_node->global_offset()+0, 0, ADD_VALUES );
      MatSetValue ( *jac, resistance_node->global_offset()+1, nb_node->global_offset()+1, 0, ADD_VALUES );
      MatSetValue ( *jac, resistance_node->global_offset()+1, nb_node->global_offset()+2, 0, ADD_VALUES );
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
 * do pre-process to jacobian matrix for DDML1 solver
 */
void IF_Metal_SchottkyBC::DDM1R_Jacobian_Preprocess(PetscScalar *, Mat *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * resistance_node     = get_region_fvm_node ( ( *node_it ), _r2 );

    src_row.push_back ( semiconductor_node->global_offset() );
    dst_row.push_back ( resistance_node->global_offset() );
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


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void IF_Metal_SchottkyBC::DDM1R_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  const PetscScalar T = T_external();

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


    //Schotty Barrier Lowerring
    PetscScalar deltaVB = semiconductor_region->material()->band->SchottyBarrierLowerring(semiconductor_node_data->eps(), semiconductor_node_data->E().size());
    //Schottky current
    PetscScalar S  = semiconductor_node->outside_boundary_surface_area();
    AutoDScalar In = semiconductor_region->material()->band->SchottyJsn(n, T, resistance_node_data->affinity() - semiconductor_node_data->affinity() - deltaVB) * S;
    AutoDScalar Ip = semiconductor_region->material()->band->SchottyJsp(p, T, resistance_node_data->affinity() - semiconductor_node_data->affinity() + deltaVB) * S;

    // schottky boundary condition of poisson's equation
    AutoDScalar f_phi = V_semiconductor - V_resistance - deltaVB;

    MatSetValue ( *jac, semiconductor_node->global_offset(), resistance_node->global_offset(), f_phi.getADValue ( 0 ), ADD_VALUES );
    MatSetValue ( *jac, semiconductor_node->global_offset(), semiconductor_node->global_offset(), f_phi.getADValue ( 1 ), ADD_VALUES );

    // process the Jacobian of Schottky current
    MatSetValue ( *jac, semiconductor_node->global_offset() +1, semiconductor_node->global_offset() +1, In.getADValue ( 2 ), ADD_VALUES );
    MatSetValue ( *jac, semiconductor_node->global_offset() +2, semiconductor_node->global_offset() +2, ( -Ip ).getADValue ( 3 ), ADD_VALUES );


    // electron/hole emit current
    MatSetValue ( *jac, resistance_node->global_offset(), semiconductor_node->global_offset() +1, (-In).getADValue ( 2 ), ADD_VALUES );
    MatSetValue ( *jac, resistance_node->global_offset(), semiconductor_node->global_offset() +2, (-Ip).getADValue ( 3 ), ADD_VALUES );



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


/*---------------------------------------------------------------------
 * update electrode IV
 */
void IF_Metal_SchottkyBC::DDM1R_Update_Solution(PetscScalar *x)
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

