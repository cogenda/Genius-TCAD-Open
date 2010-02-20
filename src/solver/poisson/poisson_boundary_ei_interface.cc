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

//  $Id: poisson_boundary_ei_interface.cc,v 1.7 2008/07/09 05:58:16 gdiso Exp $


#include "simulation_system.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void ElectrodeInsulatorInterfaceBC::Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Electrode-Insulator interface is processed here

  // note, we will use INSERT_VALUES to set values of vec f
  // if the previous operator is not insert_VALUES, we should assembly the vec
  if( (add_value_flag != INSERT_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // buffer for Vec location
  std::vector<PetscInt> src_row;

  // buffer for Vec value
  std::vector<PetscScalar> y_new;

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
        // Electrode-Insulator interface at Insulator side
      case InsulatorRegion:
        {
          // Insulator region should be the first region
          genius_assert(i==0);

          // a node on the Interface of Insulator-Semiconductor should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);

          // record the source row
          src_row.push_back(fvm_nodes[i]->global_offset());

          // find the position of ghost node
          // since we know only one ghost node exit, there is ghost_node_begin()
          FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
          const FVM_Node * ghost_fvm_node = (*gn_it).first;

          genius_assert(fvm_nodes[i]->root_node()->processor_id() == ghost_fvm_node->root_node()->processor_id() );

          // the governing equation of this fvm node

          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()];

          // since the region is sorted, we know region[0] is Insulator region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for Insulator region
          // and x[ghost_fvm_node->local_offset()] is psi for corresponding conductor region
          PetscScalar V_elec = x[ghost_fvm_node->local_offset()];

          // the psi of this node is equal to corresponding psi of conductor node
          // since psi should be continuous for the interface
          PetscScalar ff = V - V_elec;

          y_new.push_back(ff);

          genius_assert(src_row.size()==y_new.size());
          break;
        }
        // Electrode-Insulator interface at Conductor side
      case ConductorRegion:
        {
          // Conductor region should be the second region
          genius_assert(i==1);
          // a node on the Interface of Electrode-Insulator should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);

          // do nothing
          break;
        }
      case VacuumRegion:
        break;

      default: genius_error(); //we should never reach here
      }
    }

  }

  // insert new value to src row
  if( src_row.size() )
    VecSetValues(f, src_row.size(), &(src_row[0]), &(y_new[0]), INSERT_VALUES);

  add_value_flag = INSERT_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for poisson solver
 */
void ElectrodeInsulatorInterfaceBC::Poissin_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {

        // Electrode-Insulator interface at Insulator side
      case InsulatorRegion:
        {
          // a node on the Interface of Insulator-Semiconductor should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);

          // find the position of ghost node
          // since we know only one ghost node exit, there is ghost_node_begin()
          FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
          const FVM_Node * ghost_fvm_node = (*gn_it).first;

          // reserve for later operator
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), ghost_fvm_node->global_offset(), 0, ADD_VALUES);

          break;
        }
        // Electrode-Insulator interface at Conductor side
      case ConductorRegion:
        {
          // Conductor region should be the second region
          genius_assert(i==1);
          // a node on the Interface of Electrode-Insulator should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);
          break;
        }
      case VacuumRegion:
        break;

      default: genius_error(); //we should never reach here
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}






/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void ElectrodeInsulatorInterfaceBC::Poissin_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // Jacobian of Electrode-Insulator interface is processed here

  //note! MatZeroRows should be excuted on all the processor
  //no matter whether it owns this row!
  MatAssemblyBegin(*jac, MAT_FINAL_ASSEMBLY);

  // buffer for mat rows which should be removed
  std::vector<PetscInt> row_for_clear;

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
        // Electrode-Insulator interface at Insulator side
      case InsulatorRegion:
        {
          // a node on the Interface of Insulator-Semiconductor should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);

          // record the source row
          row_for_clear.push_back(fvm_nodes[i]->global_offset());

          break;
        }
        // Electrode-Insulator interface at Conductor side
      case ConductorRegion:
        {
          //do nothing
          break;
        }
      case VacuumRegion:
        break;

      default: genius_error(); //we should never reach here
      }
    }

  }

  // for efficient resion, we separate MatAssemblyBegin and MatAssemblyEnd
  MatAssemblyEnd(*jac, MAT_FINAL_ASSEMBLY);

  // clear required rows
  MatZeroRows(*jac, row_for_clear.size(), row_for_clear.empty() ? NULL : &row_for_clear[0], 0.0);


  // after that, set new Jacobian entrance to source rows
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
        // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
      case InsulatorRegion:
        {

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          // find the position of ghost node
          // since we know only one ghost node exit, there is ghost_node_begin()
          FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
          const FVM_Node * ghost_fvm_node = (*gn_it).first;

          genius_assert(fvm_nodes[i]->root_node()->processor_id() == ghost_fvm_node->root_node()->processor_id() );

          // psi of this node
          AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

          // since the region is sorted, we know region[0] is Insulator region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for Insulator region
          // and x[ghost_fvm_node->local_offset()] is psi for corresponding conductor region
          AutoDScalar V_elec = x[ghost_fvm_node->local_offset()]; V_elec.setADValue(1,1.0);

          // the psi of this node is equal to corresponding psi of conductor node
          // since psi should be continuous for the interface
          AutoDScalar ff = V - V_elec;

          // set Jacobian of governing equation ff
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), ghost_fvm_node->global_offset(), ff.getADValue(1), ADD_VALUES);

          break;

        }
        // Electrode-Insulator interface at Conductor side
      case ConductorRegion:
        {
          //do nothing
          break;
        }

      case VacuumRegion:
        break;

      default: genius_error(); //we should never reach here
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
