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



#include "simulation_system.h"
#include "resistance_region.h"
#include "insulator_region.h"
#include "boundary_condition_rr.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


/*---------------------------------------------------------------------
 * set scaling constant
 */
void ResistanceResistanceBC::Poissin_Fill_Value(Vec , Vec L)
{

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * resistance_fvm_node_1  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * resistance_fvm_node_2 = get_region_fvm_node ( ( *node_it ), _r2 );

    VecSetValue(L, resistance_fvm_node_2->global_offset(), 1.0, INSERT_VALUES);

    // if we have insulator node?
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;

        const FVM_Node * insulator_fvm_node = (*rnode_it).second.second;
        VecSetValue(L, insulator_fvm_node->global_offset(), 1.0, INSERT_VALUES);
      }
    }
  }

}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * do pre-process to function for poisson solver
 */
void ResistanceResistanceBC::Poissin_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * resistance_fvm_node_1  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * resistance_fvm_node_2 = get_region_fvm_node ( ( *node_it ), _r2 );

    src_row.push_back(resistance_fvm_node_2->global_offset());
    dst_row.push_back(resistance_fvm_node_1->global_offset());

    clear_row.push_back(resistance_fvm_node_2->global_offset());

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      if( region->type() == InsulatorRegion )
      {
        clear_row.push_back(fvm_node->global_offset());
      }
    }
  }
}



/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void ResistanceResistanceBC::Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // buffer for Vec value
  std::vector<PetscInt> iy;
  std::vector<PetscScalar> y_new;

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * resistance_fvm_node_1  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * resistance_node_data_1 = resistance_fvm_node_1->node_data();

    const FVM_Node * resistance_fvm_node_2 = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data_2 = resistance_fvm_node_2->node_data();

    // the governing equation of this fvm node --

    // psi of this node
    PetscScalar V_resistance_1 = x[resistance_fvm_node_1->local_offset()];
    // psi for ghost node
    PetscScalar V_resistance_2 = x[resistance_fvm_node_2->local_offset()];


    // affinity of this node
    PetscScalar affinite_resistance_1 = resistance_node_data_1->affinity();
    // affinity for ghost node
    PetscScalar affinite_resistance_2 = resistance_node_data_2->affinity();

    // the psi of this node is equal to corresponding psi of  node in the first resistance region
    // since psi should be continuous for the interface
    PetscScalar f_psi_metal = (e*V_resistance_1 + affinite_resistance_1) - (e*V_resistance_2 + affinite_resistance_2);
    iy.push_back(resistance_fvm_node_2->global_offset());
    y_new.push_back(f_psi_metal);

    // if we have insulator node?
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;

        const FVM_Node * insulator_fvm_node = (*rnode_it).second.second;
        const FVM_NodeData * insulator_fvm_node_data = insulator_fvm_node->node_data();
        const PetscScalar V_insulator  = x[insulator_fvm_node->local_offset()];

        PetscScalar f_psi_insulator =  V_insulator - 0.5*(V_resistance_1+V_resistance_2);
        iy.push_back(insulator_fvm_node->global_offset());
        y_new.push_back(f_psi_insulator);

      }
    }
  }


  // set new value
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y_new[0]), ADD_VALUES);
  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for poisson solver
 */
void ResistanceResistanceBC::Poissin_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * resistance_fvm_node_1  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * resistance_node_data_1 = resistance_fvm_node_1->node_data();

    const FVM_Node * resistance_fvm_node_2 = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data_2 = resistance_fvm_node_2->node_data();

    // reserve items for all the ghost nodes
    FVM_Node::fvm_ghost_node_iterator gn_it = resistance_fvm_node_1->ghost_node_begin();
    for(; gn_it != resistance_fvm_node_1->ghost_node_end(); ++gn_it)
    {
      const FVM_Node * ghost_fvm_node = (*gn_it).first;
      MatSetValue(*jac, resistance_fvm_node_1->global_offset(), ghost_fvm_node->global_offset(), 0, ADD_VALUES);

      FVM_Node::fvm_neighbor_node_iterator  gnb_it = ghost_fvm_node->neighbor_node_begin();
      for(; gnb_it != ghost_fvm_node->neighbor_node_end(); ++gnb_it)
        MatSetValue(*jac, resistance_fvm_node_1->global_offset(), (*gnb_it).first->global_offset(), 0, ADD_VALUES);
    }

    MatSetValue(*jac, resistance_fvm_node_2->global_offset(), resistance_fvm_node_1->global_offset(), 0, ADD_VALUES);

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      if ( region->type() != InsulatorRegion ) continue;

      const FVM_Node * insulator_fvm_node = (*rnode_it).second.second;

      // reserve for later operator
      MatSetValue(*jac, insulator_fvm_node->global_offset(), resistance_fvm_node_1->global_offset(), 0, ADD_VALUES);
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for poisson solver
 */
void ResistanceResistanceBC::Poissin_Jacobian_Preprocess(PetscScalar *, Mat *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * resistance_fvm_node_1  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * resistance_fvm_node_2 = get_region_fvm_node ( ( *node_it ), _r2 );

    src_row.push_back(resistance_fvm_node_2->global_offset());
    dst_row.push_back(resistance_fvm_node_1->global_offset());

    clear_row.push_back(resistance_fvm_node_2->global_offset());

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      if( region->type() == InsulatorRegion )
      {
        clear_row.push_back(fvm_node->global_offset());
      }
    }
  }
}




/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void ResistanceResistanceBC::Poissin_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  //the indepedent variable number, we need 3 here.
  adtl::AutoDScalar::numdir=3;

  // after that, set values to source rows
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * resistance_fvm_node_1  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * resistance_node_data_1 = resistance_fvm_node_1->node_data();

    const FVM_Node * resistance_fvm_node_2 = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data_2 = resistance_fvm_node_2->node_data();


    // the governing equation of this fvm node --

    // psi of this node
    AutoDScalar V_resistance_1 = x[resistance_fvm_node_1->local_offset()]; V_resistance_1.setADValue(0,1.0);
    // psi for ghost node
    AutoDScalar V_resistance_2 = x[resistance_fvm_node_2->local_offset()]; V_resistance_2.setADValue(1,1.0);


    // affinity of this node
    PetscScalar affinite_resistance_1 = resistance_node_data_1->affinity();
    // affinity for ghost node
    PetscScalar affinite_resistance_2 = resistance_node_data_2->affinity();

    // the psi of this node is equal to corresponding psi of  node in the first resistance region
    // since psi should be continuous for the interface
    AutoDScalar f_psi_metal = (e*V_resistance_1 + affinite_resistance_1) - (e*V_resistance_2 + affinite_resistance_2);

    // set Jacobian of governing equation f_psi
    MatSetValue(*jac, resistance_fvm_node_2->global_offset(), resistance_fvm_node_1->global_offset(), f_psi_metal.getADValue(0), ADD_VALUES);
    MatSetValue(*jac, resistance_fvm_node_2->global_offset(), resistance_fvm_node_2->global_offset(), f_psi_metal.getADValue(1), ADD_VALUES);


    // if we have insulator node?
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;

        const FVM_Node * insulator_fvm_node = (*rnode_it).second.second;
        const FVM_NodeData * insulator_fvm_node_data = insulator_fvm_node->node_data();
        AutoDScalar V_insulator  = x[insulator_fvm_node->local_offset()]; V_insulator.setADValue(2, 1.0);

        // the psi of this node is equal to corresponding psi of  node in the first resistance region
        // since psi should be continuous for the interface
        AutoDScalar f_psi_insulator = V_insulator - - 0.5*(V_resistance_1+V_resistance_2);

        // set Jacobian of governing equation f_psi_insulator
        MatSetValue(*jac, insulator_fvm_node->global_offset(), resistance_fvm_node_1->global_offset(), f_psi_insulator.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, insulator_fvm_node->global_offset(), resistance_fvm_node_2->global_offset(), f_psi_insulator.getADValue(1), ADD_VALUES);
        MatSetValue(*jac, insulator_fvm_node->global_offset(), insulator_fvm_node->global_offset(), f_psi_insulator.getADValue(2), ADD_VALUES);
      }
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
