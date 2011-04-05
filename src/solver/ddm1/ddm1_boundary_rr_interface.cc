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
void ResistanceResistanceBC::DDM1_Fill_Value(Vec , Vec L)
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
 * do pre-process to function for DDML1 solver
 */
void ResistanceResistanceBC::DDM1_Function_Preprocess(Vec f, std::vector<PetscInt> &src_row,
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
 * build function and its jacobian for DDML1 solver
 */
void ResistanceResistanceBC::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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
    PetscScalar f_psi_metal = (V_resistance_1 + affinite_resistance_1) - (V_resistance_2 + affinite_resistance_2);
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

        PetscScalar f_psi_insulator =  V_insulator - V_resistance_1;
        iy.push_back(insulator_fvm_node->global_offset());
        y_new.push_back(f_psi_insulator);


        //for transiet simulation, we should consider displacement current from insulator region
        if(SolverSpecify::TimeDependent == true)
        {
          PetscScalar I_displacement = 0.0;
          FVM_Node::fvm_neighbor_node_iterator nb_it = insulator_fvm_node->neighbor_node_begin();
          for(; nb_it != insulator_fvm_node->neighbor_node_end(); ++nb_it)
          {
            const FVM_Node *nb_node = (*nb_it).second;
            const FVM_NodeData * nb_node_data = nb_node->node_data();
            // the psi of neighbor node
            PetscScalar V_nb = x[nb_node->local_offset()];
            // distance from nb node to this node
            PetscScalar distance = insulator_fvm_node->distance(nb_node);
            // area of out surface of control volume related with neighbor node
            PetscScalar cv_boundary = insulator_fvm_node->cv_surface_area(nb_node->root_node());
            PetscScalar dEdt;
            if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false) //second order
            {
              PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
              dEdt = ( (2-r)/(1-r)*(V_insulator-V_nb)
                       - 1.0/(r*(1-r))*(insulator_fvm_node_data->psi()-nb_node_data->psi())
                       + (1-r)/r*(insulator_fvm_node_data->psi_last()-nb_node_data->psi_last()))/distance/(SolverSpecify::dt_last+SolverSpecify::dt);
            }
            else//first order
            {
              dEdt = ((V_insulator-V_nb)-(insulator_fvm_node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;
            }

            I_displacement += cv_boundary*insulator_fvm_node_data->eps()*dEdt;
          }
          VecSetValue ( f, resistance_fvm_node_1->global_offset(), -I_displacement, ADD_VALUES );
        }
      }
    }
  }

  // set new value
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y_new[0]), ADD_VALUES);
  add_value_flag = ADD_VALUES;
}



/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML1 solver
 */
void ResistanceResistanceBC::DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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
        MatSetValue(*jac, resistance_fvm_node_1->global_offset(), (*gnb_it).second->global_offset(), 0, ADD_VALUES);
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

      //reserve for displacement current
      MatSetValue(*jac, resistance_fvm_node_1->global_offset(), insulator_fvm_node->global_offset(), 0, ADD_VALUES);
      FVM_Node::fvm_neighbor_node_iterator nb_it = insulator_fvm_node->neighbor_node_begin();
      for(; nb_it != insulator_fvm_node->neighbor_node_end(); ++nb_it)
      {
        const FVM_Node *nb_node = (*nb_it).second;
        MatSetValue(*jac, resistance_fvm_node_1->global_offset(), nb_node->global_offset(), 0, ADD_VALUES);
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void ResistanceResistanceBC::DDM1_Jacobian_Preprocess(Mat *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  // here we do several things:
  // add some row to other, clear some row, insert some value to row
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
 * build function and its jacobian for DDML1 solver
 */
void ResistanceResistanceBC::DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  //the indepedent variable number, we need 2 here.
  adtl::AutoDScalar::numdir=2;

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
    AutoDScalar f_psi_metal = (V_resistance_1 + affinite_resistance_1) - (V_resistance_2 + affinite_resistance_2);

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
        AutoDScalar V_insulator  = x[insulator_fvm_node->local_offset()]; V_insulator.setADValue(1, 1.0);

        // the psi of this node is equal to corresponding psi of  node in the first resistance region
        // since psi should be continuous for the interface
        AutoDScalar f_psi_insulator = V_insulator - V_resistance_1;

        // set Jacobian of governing equation ff
        MatSetValue(*jac, insulator_fvm_node->global_offset(), resistance_fvm_node_1->global_offset(), f_psi_insulator.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, insulator_fvm_node->global_offset(), insulator_fvm_node->global_offset(), f_psi_insulator.getADValue(1), ADD_VALUES);

        // displacement current
        if(SolverSpecify::TimeDependent == true)
        {
          FVM_Node::fvm_neighbor_node_iterator nb_it = insulator_fvm_node->neighbor_node_begin();
          for(; nb_it != insulator_fvm_node->neighbor_node_end(); ++nb_it)
          {
            const FVM_Node *nb_node = (*nb_it).second;
            const FVM_NodeData * nb_node_data = nb_node->node_data();

            // the psi of neighbor node
            AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(0, 1.0);

            // distance from nb node to this node
            PetscScalar distance = insulator_fvm_node->distance(nb_node);

            // area of out surface of control volume related with neighbor node
            PetscScalar cv_boundary = insulator_fvm_node->cv_surface_area(nb_node->root_node());
            AutoDScalar dEdt;
            if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false) //second order
            {
              PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
              dEdt = ( (2-r)/(1-r)*(V_insulator-V_nb)
                       - 1.0/(r*(1-r))*(insulator_fvm_node_data->psi()-nb_node_data->psi())
                       + (1-r)/r*(insulator_fvm_node_data->psi_last()-nb_node_data->psi_last()))/distance/(SolverSpecify::dt_last+SolverSpecify::dt);
            }
            else//first order
            {
              dEdt = ((V_insulator-V_nb)-(insulator_fvm_node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;
            }

            AutoDScalar current_disp = cv_boundary*insulator_fvm_node_data->eps()*dEdt;

            MatSetValue(*jac, resistance_fvm_node_1->global_offset(), insulator_fvm_node->global_offset(), -current_disp.getADValue(1), ADD_VALUES);
            MatSetValue(*jac, resistance_fvm_node_1->global_offset(), nb_node->global_offset(), -current_disp.getADValue(0), ADD_VALUES);
          }
        }
      }
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
