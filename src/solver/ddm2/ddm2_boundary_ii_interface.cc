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


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void InsulatorInsulatorInterfaceBC::DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // buffer for Vec location
  std::vector<PetscInt> src_row;
  std::vector<PetscInt> dst_row;

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
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      if(!fvm_node->is_valid()) continue;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // we may have several regions with insulator material

      // the first insulator region
      if(i==0)
      {
        // do nothing.
        // however, we will add fvm integral of other regions to it.
      }
      // other insulator regions
      else
      {
        // record the source row and dst row
        src_row.push_back(fvm_nodes[i]->global_offset()+0);
        src_row.push_back(fvm_nodes[i]->global_offset()+1);

        // find the position ff will be add to
        // ff will be added to fvm_nodes[0]

        const FVM_Node * ghost_fvm_node = fvm_nodes[0];
        dst_row.push_back(ghost_fvm_node->global_offset()+0);
        dst_row.push_back(ghost_fvm_node->global_offset()+1);


        // the ghost node should have the same processor_id with me
        genius_assert(fvm_nodes[i]->root_node()->processor_id() == ghost_fvm_node->root_node()->processor_id() );

        // the governing equation of this fvm node --


        PetscScalar V = x[fvm_nodes[i]->local_offset()+0]; // psi of this node
        PetscScalar T = x[fvm_nodes[i]->local_offset()+1]; // T of this node


        PetscScalar V_in = x[fvm_nodes[0]->local_offset()+0]; // psi for ghost node
        PetscScalar T_in = x[fvm_nodes[0]->local_offset()+1]; // T for ghost node

        // the psi of this node is equal to corresponding psi of node in the first insulator region
        // since psi should be continuous for the interface
        PetscScalar ff1 = V - V_in;
        y_new.push_back(ff1);

        // the T of this node is equal to corresponding T of node in the first insulator region
        // by assuming no heat resistance between 2 region
        PetscScalar ff2 = T - T_in;
        y_new.push_back(ff2);

        genius_assert(src_row.size()==y_new.size());

      }

    }

  }

  // add src row to dst row, it will assemble vec automatically
  PetscUtils::VecAddRowToRow(f, src_row, dst_row);

  // insert new value to src row
  if( src_row.size() )
    VecSetValues(f, src_row.size(), &(src_row[0]), &(y_new[0]), INSERT_VALUES);

  add_value_flag = INSERT_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML2 solver
 */
void InsulatorInsulatorInterfaceBC::DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      if(!fvm_node->is_valid()) continue;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first insulator region
      if(i==0)
      {
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
      }

      // other insulator region
      else
      {
        // reserve for later operator
        MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, fvm_nodes[0]->global_offset()+0, 0, ADD_VALUES);
        MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, fvm_nodes[0]->global_offset()+1, 0, ADD_VALUES);
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}






/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void InsulatorInsulatorInterfaceBC::DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // here we do several things:
  // add some row to other, clear some row, insert some value to row
  // I wonder if there are some more efficient way to do these.

  {
    // buffer for mat rows which should be added to other row
    std::vector<PetscInt> src_row;
    std::vector<PetscInt> dst_row;

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
        const SimulationRegion * region = (*rnode_it).second.first;
        const FVM_Node * fvm_node = (*rnode_it).second.second;
        if(!fvm_node->is_valid()) continue;

        regions.push_back( region );
        fvm_nodes.push_back( fvm_node );

        // the first insulator region
        if(i==0) continue;

        // other insulator region
        else
        {
          // record the source row and dst row

          src_row.push_back(fvm_nodes[i]->global_offset()+0);
          src_row.push_back(fvm_nodes[i]->global_offset()+1);

          dst_row.push_back(fvm_nodes[0]->global_offset()+0);
          dst_row.push_back(fvm_nodes[0]->global_offset()+1);
        }
      }
    }

    //ok, we add source rows to destination rows
    PetscUtils::MatAddRowToRow(*jac, src_row, dst_row);

    // clear source rows
    MatZeroRows(*jac, src_row.size(), src_row.empty() ? NULL : &src_row[0], 0.0);

  }



  // after that, set values to source rows
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
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      if(!fvm_node->is_valid()) continue;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );


      // the first insulator region
      if(i==0) continue;

      // other insulator region
      else
      {

        //the indepedent variable number, we need 2 here.
        adtl::AutoDScalar::numdir=2;

        // psi of this node
        AutoDScalar  V    = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

        // psi for ghost node
        AutoDScalar  V_in = x[fvm_nodes[0]->local_offset()]; V_in.setADValue(1,1.0);

        // the psi of this node is equal to corresponding psi of insulator node int the other region
        AutoDScalar  ff1 = V - V_in;

        // set Jacobian of governing equation ff
        MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff1.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), ff1.getADValue(1), ADD_VALUES);


        // T of this node
        AutoDScalar  T = x[fvm_nodes[i]->local_offset()+1]; T.setADValue(0,1.0);

        // T for corresponding conductor region
        AutoDScalar  T_elec = x[fvm_nodes[0]->local_offset()+1]; T_elec.setADValue(1,1.0);

        // the T of this node is equal to corresponding T of conductor node
        // we assuming T is continuous for the interface
        AutoDScalar ff2 = T - T_elec;

        // set Jacobian of governing equation ff2
        MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, fvm_nodes[i]->global_offset()+1, ff2.getADValue(0), ADD_VALUES);
        MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, fvm_nodes[0]->global_offset()+1, ff2.getADValue(1), ADD_VALUES);

      }

    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
