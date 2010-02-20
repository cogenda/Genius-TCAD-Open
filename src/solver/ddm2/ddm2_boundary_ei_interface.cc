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
#include "boundary_condition.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

// for DDML2 solver, poisson's equation and lattice temperature should be considered in Electrode-Insulator interface


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void ElectrodeInsulatorInterfaceBC::DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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
  std::vector<PetscInt> dst_row;

  // buffer for Vec location
  std::vector<PetscInt> iy;

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

          // find the position of ghost node
          // since we know only one ghost node exit, there is ghost_node_begin()
          FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
          const FVM_Node * ghost_fvm_node = (*gn_it).first;

          genius_assert(fvm_nodes[i]->root_node()->processor_id() == ghost_fvm_node->root_node()->processor_id() );


          // the governing equation of this fvm node

          PetscScalar V = x[fvm_nodes[i]->local_offset()+0]; // psi of this node
          PetscScalar T = x[fvm_nodes[i]->local_offset()+1]; // T of this node

          // since the region is sorted, we know region[0] is Insulator region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for Insulator region
          // and x[ghost_fvm_node->local_offset()] is psi for corresponding conductor region
          PetscScalar V_elec = x[ghost_fvm_node->local_offset()+0];
          PetscScalar T_elec = x[ghost_fvm_node->local_offset()+1];

          // the psi of this node is equal to corresponding psi of conductor node
          // since psi should be continuous for the interface
          PetscScalar ff1 = V - V_elec;
          // record the source row
          iy.push_back(fvm_nodes[i]->global_offset()+0);
          y_new.push_back(ff1);

          // the T of this node is equal to corresponding T of conductor node
          // by assuming no heat resistance between 2 region
          PetscScalar ff2 = T - T_elec;
          // record the source row
          iy.push_back(fvm_nodes[i]->global_offset()+1);
          y_new.push_back(ff2);

          break;
        }
        // Electrode-Insulator interface at Conductor side
      case ConductorRegion:
        {
          // Conductor region should be the second region
          genius_assert(i==1);
          // a node on the Interface of Electrode-Insulator should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);

          // record the source row and dst row
          src_row.push_back(fvm_nodes[0]->global_offset()+0);
          src_row.push_back(fvm_nodes[0]->global_offset()+1);

          dst_row.push_back(fvm_nodes[i]->global_offset()+0);
          dst_row.push_back(fvm_nodes[i]->global_offset()+1);

          break;
        }
      case VacuumRegion:
        break;
      default: genius_error(); //we should never reach here
      }
    }

  }

  // add src row to dst row, it will assemble vec automatically
  PetscUtils::VecAddRowToRow(f, src_row, dst_row);

  // do insert here
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y_new[0]), INSERT_VALUES);

  add_value_flag = INSERT_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML2 solver
 */
void ElectrodeInsulatorInterfaceBC::DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, ghost_fvm_node->global_offset()+0, 0, ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, ghost_fvm_node->global_offset()+1, 0, ADD_VALUES);

          break;
        }
        // Electrode-Insulator interface at Conductor side
      case ConductorRegion:
        {
          // Conductor region should be the second region
          genius_assert(i==1);
          // a node on the Interface of Electrode-Insulator should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);

          // reserve items for all the ghost nodes
          FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
          for(; gn_it != fvm_nodes[i]->ghost_node_end(); ++gn_it)
          {
            const FVM_Node * ghost_fvm_node = (*gn_it).first;
            MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, ghost_fvm_node->global_offset()+0, 0, ADD_VALUES);
            MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, ghost_fvm_node->global_offset()+1, 0, ADD_VALUES);

            FVM_Node::fvm_neighbor_node_iterator  gnb_it = ghost_fvm_node->neighbor_node_begin();
            for(; gnb_it != ghost_fvm_node->neighbor_node_end(); ++gnb_it)
            {
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, (*gnb_it).second->global_offset()+0, 0, ADD_VALUES);
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, (*gnb_it).second->global_offset()+1, 0, ADD_VALUES);
            }
          }

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
 * build function and its jacobian for DDML2 solver
 */
void ElectrodeInsulatorInterfaceBC::DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{


  // Jacobian of Electrode-Insulator interface is processed here

  {
    // buffer for mat rows which should be added to other row
    std::vector<PetscInt> src_row;
    std::vector<PetscInt> dst_row;

    // buffer for mat rows which should be removed
    std::vector<PetscInt> rm_row;

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
            rm_row.push_back(fvm_nodes[i]->global_offset()+0);
            rm_row.push_back(fvm_nodes[i]->global_offset()+1);

            break;
          }
          // Electrode-Insulator interface at Conductor side
        case ConductorRegion:
          {
            // record the source row and dst row
            src_row.push_back(fvm_nodes[0]->global_offset()+0);
            src_row.push_back(fvm_nodes[0]->global_offset()+1);

            dst_row.push_back(fvm_nodes[i]->global_offset()+0);
            dst_row.push_back(fvm_nodes[i]->global_offset()+1);

            break;
          }
        case VacuumRegion:
          break;
        default: genius_error(); //we should never reach here
        }
      }

    }

    //ok, we add source rows to destination rows
    PetscUtils::MatAddRowToRow(*jac, src_row, dst_row);

    // clear source rows
    MatZeroRows(*jac, rm_row.size(), rm_row.empty() ? NULL : &rm_row[0], 0.0);

  }



  // after that, set new Jacobian entrance to source rows
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
          AutoDScalar ff1 = V - V_elec;

          // set Jacobian of governing equation ff
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff1.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), ghost_fvm_node->global_offset(), ff1.getADValue(1), ADD_VALUES);


          // T of this node
          AutoDScalar  T = x[fvm_nodes[i]->local_offset()+1]; T.setADValue(0,1.0);

          // T for corresponding conductor region
          AutoDScalar  T_elec = x[ghost_fvm_node->local_offset()+1]; T_elec.setADValue(1,1.0);

          // the T of this node is equal to corresponding T of conductor node
          // we assuming T is continuous for the interface
          AutoDScalar ff2 = T - T_elec;

          // set Jacobian of governing equation ff2
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, fvm_nodes[i]->global_offset()+1, ff2.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, ghost_fvm_node->global_offset()+1, ff2.getADValue(1), ADD_VALUES);

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
