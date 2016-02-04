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
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_poly_gate.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;



/*---------------------------------------------------------------------
 * set scaling constant
 */
void PolyGateContactBC::Poissin_Fill_Value(Vec , Vec L)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      switch ( region->type() )
      {
          case SemiconductorRegion :
          {
            const FVM_Node * fvm_node = ( *rnode_it ).second.second;
            VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
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

// this file is nearly the same as poisson_boundary_ohmic.


/*---------------------------------------------------------------------
 * do pre-process to function for poisson solver
 */
void PolyGateContactBC::Poissin_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
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
      clear_row.push_back(row);
    }
  }
}



/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void PolyGateContactBC::Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Gate boundary condition is processed here

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  const PetscScalar Work_Function = this->scalar("workfunction");

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
          // semiconductor region.
          case SemiconductorRegion:
          {
            // psi of this node
            PetscScalar V = x[fvm_nodes[i]->local_offset()];

            // the governing equation
            PetscScalar ff = V + Work_Function - ext_circuit()->Vapp();
            // set governing equation to function vector
            VecSetValue(f, fvm_nodes[i]->global_offset(), ff, ADD_VALUES);

            break;
          }

          // insulator region
          case InsulatorRegion:
          {

            // psi of this node
            PetscScalar V = x[fvm_nodes[i]->local_offset()];

            PetscScalar V_0 = x[fvm_nodes[0]->local_offset()];

            // the psi of this node is equal to corresponding psi of semiconductor node
            PetscScalar ff = V - V_0;

            // set governing equation to function vector
            VecSetValue(f, fvm_nodes[i]->global_offset(), ff, ADD_VALUES);

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
 * do pre-process to jacobian matrix for poisson solver
 */
void PolyGateContactBC::Poissin_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
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
      clear_row.push_back(row);
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void PolyGateContactBC::Poissin_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // the Jacobian of GateContact boundary condition is processed here
  // we use AD again. no matter it is overkill here.

  const PetscScalar Work_Function = this->scalar("workfunction");

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
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
          case SemiconductorRegion:
          {

            //the indepedent variable number, we need 1 here.
            adtl::AutoDScalar::numdir=1;

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

            // the governing equation
            AutoDScalar ff = V + Work_Function - ext_circuit()->Vapp();

            // set Jacobian of governing equation ff
            jac->add(fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0));

            break;
          }
          case InsulatorRegion:
          {

            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

            AutoDScalar  V_0 = x[fvm_nodes[0]->local_offset()]; V_0.setADValue(1,1.0);

            // the psi of this node is equal to corresponding psi of Insulator node
            AutoDScalar  ff = V - V_0;

            // set Jacobian of governing equation ff
            jac->add(fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0));
            jac->add(fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), ff.getADValue(1));

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
 * update electrode potential
 */
void PolyGateContactBC::Poissin_Update_Solution(PetscScalar *)
{
  this->ext_circuit()->potential() = this->ext_circuit()->Vapp();
}
