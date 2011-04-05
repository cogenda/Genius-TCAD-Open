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

//  $Id: poisson_boundary_is_interface.cc,v 1.14 2008/07/09 05:58:16 gdiso Exp $


#include "simulation_system.h"
#include "semiconductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_is.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;



/*---------------------------------------------------------------------
 * set scaling constant
 */
void InsulatorSemiconductorInterfaceBC::Poissin_Fill_Value(Vec , Vec L)
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
      const SimulationRegion * region = (*rnode_it).second.first;
      if( region->type() != InsulatorRegion ) continue;

      const FVM_Node * fvm_node = (*rnode_it).second.second;
      VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
    }
  }
}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for poisson solver
 */
void InsulatorSemiconductorInterfaceBC::Poissin_Function_Preprocess(Vec f, std::vector<PetscInt> &src_row,
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
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
          // Insulator-Semiconductor interface at Semiconductor side, do nothing
          case SemiconductorRegion:  break;

          // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
          case InsulatorRegion:
          {
            // record the source row and dst row
            src_row.push_back(fvm_nodes[i]->global_offset());
            dst_row.push_back(fvm_nodes[0]->global_offset());
            clear_row.push_back(fvm_nodes[i]->global_offset());
            break;
          }
          case VacuumRegion:
          break;

          default: genius_error(); //we should never reach here
      }
    }

  }

}



/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void InsulatorSemiconductorInterfaceBC::Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Insulator-Semiconductor interface is processed here

  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // buffer for Vec location
  std::vector<PetscInt> iy;
  iy.reserve(n_nodes());
  // buffer for Vec value
  std::vector<PetscScalar> y;
  y.reserve(n_nodes());

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
          // Insulator-Semiconductor interface at Semiconductor side, do nothing
          case SemiconductorRegion:
          {
            // semiconductor region should be the first region
            genius_assert(i==0);

            // process interface fixed charge density
            PetscScalar boundary_area = fvm_nodes[i]->outside_boundary_surface_area();
            VecSetValue(f, fvm_nodes[i]->global_offset(), this->Qf()*boundary_area, ADD_VALUES);

            break;
          }

          // Insulator-Semiconductor interface at Insulator side, we get the previous f for insulator region,
          // as well as surface charge density, plus to corresponding function of semiconductor. and then
          // force the potential equal to corresponding point in semiconductor
          case InsulatorRegion:
          {
            // insulator region should be the second region
            genius_assert(i==1);

            // record the source row and dst row
            iy.push_back(fvm_nodes[i]->global_offset());

            // the governing equation of this fvm node

            // psi of this node
            PetscScalar V = x[fvm_nodes[i]->local_offset()];

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            //genius_assert( regions[0]->type()==SemiconductorRegion );
            PetscScalar V_semi = x[fvm_nodes[0]->local_offset()];

            // the psi of this node is equal to corresponding psi of semiconductor node
            // since psi should be continuous for the interface
            PetscScalar ff = V - V_semi;

            y.push_back(ff);

            genius_assert(iy.size()==y.size());
            break;
          }

          case VacuumRegion:
          break;

          default: genius_error(); //we should never reach here
      }
    }

  }

  // insert new value to src row
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y[0]), ADD_VALUES);

  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for poisson solver
 */
void InsulatorSemiconductorInterfaceBC::Poissin_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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
          // Insulator-Semiconductor interface at Semiconductor side, we should reserve entrance for later add operator
          case SemiconductorRegion:
          {
            // semiconductor region should be the first region
            genius_assert(i==0);

            // since we know only one ghost node exit, there is ghost_node_begin()
            FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
            const FVM_Node * ghost_fvm_node = (*gn_it).first;
            MatSetValue(*jac, fvm_nodes[i]->global_offset(), ghost_fvm_node->global_offset(), 0, ADD_VALUES);

            FVM_Node::fvm_neighbor_node_iterator  gnb_it = ghost_fvm_node->neighbor_node_begin();
            for(; gnb_it != ghost_fvm_node->neighbor_node_end(); ++gnb_it)
              MatSetValue(*jac, fvm_nodes[i]->global_offset(), (*gnb_it).second->global_offset(), 0, ADD_VALUES);

            break;
          }

          // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
          case InsulatorRegion:
          {
            genius_assert(i==1);
            // reserve for later operator
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
 * do pre-process to jacobian matrix for poisson solver
 */
void InsulatorSemiconductorInterfaceBC::Poissin_Jacobian_Preprocess(Mat *jac, std::vector<PetscInt> &src_row,
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
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
          // Insulator-Semiconductor interface at Semiconductor side, do nothing
          case SemiconductorRegion:  break;

          // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
          case InsulatorRegion:
          {
            // record the source row and dst row
            src_row.push_back(fvm_nodes[i]->global_offset());
            dst_row.push_back(fvm_nodes[0]->global_offset());
            clear_row.push_back(fvm_nodes[i]->global_offset());
            break;
          }
          case VacuumRegion:
          break;

          default: genius_error(); //we should never reach here
      }
    }

  }

}




/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void InsulatorSemiconductorInterfaceBC::Poissin_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  //the indepedent variable number, we need 2 here.
  adtl::AutoDScalar::numdir=2;

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
          // Insulator-Semiconductor interface at Semiconductor side, do nothing
          case SemiconductorRegion:  break;

          // Insulator-Semiconductor interface at Insulator side, we should add the rows to semiconductor region
          case InsulatorRegion:
          {

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

            // since the region is sorted, we know region[0] is semiconductor region
            // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
            genius_assert( regions[0]->type()==SemiconductorRegion );
            AutoDScalar  V_semi = x[fvm_nodes[0]->local_offset()]; V_semi.setADValue(1,1.0);

            // the psi of this node is equal to corresponding psi of semiconductor node
            AutoDScalar  ff = V - V_semi;

            // set Jacobian of governing equation ff
            MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0), ADD_VALUES);
            MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[0]->global_offset(), ff.getADValue(1), ADD_VALUES);

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
