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
#include "boundary_condition_hetero.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


/*---------------------------------------------------------------------
 * set scaling constant
 */
void HeteroInterfaceBC::Poissin_Fill_Value(Vec , Vec L)
{
  // search for all the node with this boundary type
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
      // the first semiconductor region
      if(i==0)
      {
        // do nothing.
      }
      // other semiconductor regions with same material
      else
      {
        const FVM_Node * fvm_node = ( *rnode_it ).second.second;
        VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
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
void HeteroInterfaceBC::Poissin_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
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

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0) continue;

      // other semiconductor regions
      else
      {
        // record the source row and dst row
        src_row.push_back(fvm_nodes[i]->global_offset());

        const FVM_Node * ghost_fvm_node = fvm_nodes[0];
        dst_row.push_back(ghost_fvm_node->global_offset());

        clear_row.push_back(fvm_nodes[i]->global_offset());
      }
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void HeteroInterfaceBC::Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // buffer for Vec index
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
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      if(!fvm_node->is_valid()) continue;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // we may have several regions with same semiconductor material

      // the first semiconductor region
      if(i==0)
      {
        // do nothing.
        // however, we will add fvm integral of other regions to it.
      }
      // other semiconductor regions with different material
      else
      {
        iy.push_back(fvm_nodes[i]->global_offset());

        // the governing equation of this fvm node --

        // psi of this node
        PetscScalar V = x[fvm_nodes[i]->local_offset()];

        PetscScalar V_semi = x[fvm_nodes[0]->local_offset()];

        // the psi of this node is equal to corresponding psi of  node in the first semiconductor region
        // since psi should be continuous for the interface
        PetscScalar ff = V - V_semi;

        y_new.push_back(ff);

      }

    }

  }

  // set new value to src row
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y_new[0]), ADD_VALUES);

  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for poisson solver
 */
void HeteroInterfaceBC::Poissin_Jacobian_Preprocess(PetscScalar *, SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
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

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0) continue;

      // other semiconductor regions
      else
      {
        // record the source row and dst row
        src_row.push_back(fvm_nodes[i]->global_offset());

        const FVM_Node * ghost_fvm_node = fvm_nodes[0];
        dst_row.push_back(ghost_fvm_node->global_offset());

        clear_row.push_back(fvm_nodes[i]->global_offset());
      }
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void HeteroInterfaceBC::Poissin_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{

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

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );


      // the first semiconductor region
      if(i==0) continue;

      // other semiconductor regions
      else
      {

        //the indepedent variable number, we need 2 here.
        adtl::AutoDScalar::numdir=2;

        // psi of this node
        AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

        AutoDScalar  V_semi = x[fvm_nodes[0]->local_offset()]; V_semi.setADValue(1,1.0);

        // the psi of this node is equal to corresponding psi of semiconductor node int the other region
        AutoDScalar  ff = V - V_semi;

        // set Jacobian of governing equation ff
        jac->add(fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0));
        jac->add(fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), ff.getADValue(1));

      }

    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
