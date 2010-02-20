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

//  $Id: poisson_boundary_gate.cc,v 1.8 2008/07/09 07:53:36 gdiso Exp $


#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

// this file is nearly the same as poisson_boundary_ohmic.

/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void GateContactBC::Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Gate boundary condition is processed here

  // note, we will use INSERT_VALUES to set values of vec f
  // if the previous operator is not insert_VALUES, we should assembly the vec
  if( (add_value_flag != INSERT_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

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
        // insulator region.
      case InsulatorRegion:
        {

          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()];

          // the governing equation
          PetscScalar ff = V + Work_Function() - ext_circuit()->Vapp();
          // set governing equation to function vector
          VecSetValue(f, fvm_nodes[i]->global_offset(), ff, INSERT_VALUES);

          break;
        }


        // conductor region
      case ConductorRegion:
        {

          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()];

          // since the region is sorted, we know region[0] is Insulator region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding Insulator region

          PetscScalar V_in = x[fvm_nodes[0]->local_offset()];

          // the psi of this node is equal to corresponding psi of Insulator node
          PetscScalar ff = V - V_in;

          // set governing equation to function vector
          VecSetValue(f, fvm_nodes[i]->global_offset(), ff, INSERT_VALUES);

          break;
        }

      case VacuumRegion:
        break;

      default: genius_error(); //we should never reach here
      }
    }

  }

  // the last operator is INSERT_VALUES
  add_value_flag = INSERT_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for poisson solver
 */
void GateContactBC::Poissin_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

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
      case InsulatorRegion:
        {
          // do nothing
          break;
        }
      case ConductorRegion:
        {
          // insert none zero pattern
          //MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), 0, ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), 0, ADD_VALUES);

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
void GateContactBC::Poissin_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of GateContact boundary condition is processed here
  // we use AD again. no matter it is overkill here.

  // first, we zero all the rows corresponding to GateContact bc.

  // the indicator which rows should be set to zero
  std::vector<PetscInt> id;
  id.reserve(n_nodes());

  //note! MatZeroRows should be excuted on all the processor
  //no matter whether it owns this row!
  MatAssemblyBegin(*jac, MAT_FINAL_ASSEMBLY);

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    // should clear all the rows related with this boundary condition
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;

      PetscInt row = fvm_node->global_offset();
      id.push_back(row);
    }
  }

  // for efficient resion, we separate MatAssemblyBegin and MatAssemblyEnd
  MatAssemblyEnd(*jac, MAT_FINAL_ASSEMBLY);

  // zero required rows
  MatZeroRows(*jac, id.size(), id.empty() ? NULL : &id[0], 0.0);


  // after that, we should insert some value to it
  // loop again
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
      case InsulatorRegion:
        {

          //the indepedent variable number, we need 1 here.
          adtl::AutoDScalar::numdir=1;

          // psi of this node
          AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

          // the governing equation
          AutoDScalar ff = V + Work_Function() - ext_circuit()->Vapp();

          // set Jacobian of governing equation ff
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0), INSERT_VALUES);

          break;
        }
        // conductor region which has an interface with Schottky Contact boundary to semiconductor region
      case ConductorRegion:
        {

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          // psi of this node
          AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

          // since the region is sorted, we know region[0] is Insulator region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding Insulator region
          AutoDScalar  V_in = x[fvm_nodes[0]->local_offset()]; V_in.setADValue(1,1.0);

          // the psi of this node is equal to corresponding psi of Insulator node
          AutoDScalar  ff = V - V_in;

          // set Jacobian of governing equation ff
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0), INSERT_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), ff.getADValue(1), INSERT_VALUES);

          break;
        }

      case VacuumRegion:
        break;

      default: genius_error(); //we should never reach here
      }
    }

  }

  // the last operator is INSERT_VALUES
  add_value_flag = INSERT_VALUES;

}



/*---------------------------------------------------------------------
 * update electrode potential
 */
void GateContactBC::Poissin_Update_Solution(PetscScalar *)
{
  this->ext_circuit()->potential() = this->ext_circuit()->Vapp();
}
