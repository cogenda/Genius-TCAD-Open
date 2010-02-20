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

#include "simulation_system.h"
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition.h"
#include "parallel.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * build function and its jacobian for DDM1 solver
 */
void ChargedContactBC::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // charged float metal is processed here

  // note, we will use INSERT_VALUES to set values of vec f
  // if the previous operator is not insert_VALUES, we should assembly the vec
  if( (add_value_flag != INSERT_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar z_width = this->z_width();

  // buffer for float metal surface potential equation \integral {\eps \frac{\partial P}{\partial n} } + \sigma = 0
  std::vector<PetscScalar> surface_equ;

  // the float metal potential in current iteration
  genius_assert( local_offset()!=invalid_uint );
  PetscScalar Vm = x[this->local_offset()];

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
        // FloatMetal-Insulator interface at Insulator side
      case InsulatorRegion:
        {
          // Insulator region should be the first region
          genius_assert(i==0);
          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

          // the governing equation of this fvm node

          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()];

          // the psi of this node is equal to psi of float metal
          PetscScalar ff = V - Vm;

          VecSetValue(f,  fvm_nodes[i]->global_offset(), ff, INSERT_VALUES);

          // process float metal surface equation
          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
          FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_nodes[i]->neighbor_node_end();
          for(; nb_it != nb_it_end; ++nb_it)
          {
            const FVM_Node *nb_node = (*nb_it).second;
            // the psi of neighbor node
            PetscScalar V_nb = x[nb_node->local_offset()];
            // distance from nb node to this node
            PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();
            // area of out surface of control volume related with neighbor node,
            // here we should consider the difference of 2D/3D by multiply z_width
            PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node->root_node())*z_width;
            // surface electric field
            PetscScalar E = (V_nb-V)/distance;
            PetscScalar D = node_data->eps()*E;
            // insert surface integral of electric displacement of this node into surface_equ buffer
            surface_equ.push_back(cv_boundary*D);
          }

          break;
        }
        // FloatMetal-Insulator interface at Conductor side
      case ConductorRegion:
        {
          // Conductor region should be the second region
          genius_assert(i==1);

          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()];

          // the psi of this node is equal to psi of float metal
          PetscScalar ff = V - Vm;

          VecSetValue(f,  fvm_nodes[i]->global_offset(), ff, INSERT_VALUES);

          break;
        }
      default: genius_error(); //we should never reach here
      }
    }

  }


  // we first gather the nodal surface integral of electric displacement from all the processors
  Parallel::allgather(surface_equ);

  // we sum all the terms in surface_equ buffer, that is surface integral of electric displacement for the whole float metal
  PetscScalar surface_integral_electric_displacement = std::accumulate(surface_equ.begin(), surface_equ.end(), 0.0 );

  if(Genius::processor_id() == Genius::n_processors() -1)
  {
    // the governing equation of float metal surface, process it only on last processor
    PetscScalar f_ext = surface_integral_electric_displacement + this->Qf();

    VecSetValue(f, this->global_offset(), f_ext, INSERT_VALUES);
  }

  add_value_flag = INSERT_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDM1 solver
 */
void ChargedContactBC::DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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

        // FloatMetal-Insulator interface at Insulator side
      case InsulatorRegion:
        {
          // a node on the Interface of Insulator-Semiconductor should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);

          // reserve for later operator
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), this->global_offset(), 0, ADD_VALUES);

          break;
        }
        // FloatMetal-Insulator interface at Conductor side
      case ConductorRegion:
        {
          // Conductor region should be the second region
          genius_assert(i==1);

          // reserve for later operator
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), this->global_offset(), 0, ADD_VALUES);

          break;
        }
      default: genius_error(); //we should never reach here
      }
    }

  }

  // reserve jacobian entries for the surface potential equation of float metal
  if(Genius::processor_id() == Genius::n_processors() -1)
  {
    MatSetValue(*jac, this->global_offset(), this->global_offset(), 0, ADD_VALUES);

    for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // skip node not belongs to this processor
      if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

      // get the derivative of surface potential equation to insulator
      const FVM_Node *  fvm_node = get_region_fvm_node(*node_it, InsulatorRegion);
      MatSetValue(*jac, this->global_offset(), fvm_node->global_offset(), 0, ADD_VALUES);

      // get the derivative of surface potential equation to neighbors of insulator node

      FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
      FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
      for(; nb_it != nb_it_end; ++nb_it)
      {
        const FVM_Node *  fvm_nb_node = (*nb_it).second;
        MatSetValue(*jac, this->global_offset(), fvm_nb_node->global_offset(), 0, ADD_VALUES);
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}






/*---------------------------------------------------------------------
 * build function and its jacobian for DDM1 solver
 */
void ChargedContactBC::DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // Jacobian of Electrode-Insulator interface is processed here

  {

    MatAssemblyBegin(*jac, MAT_FINAL_ASSEMBLY);

    // buffer for mat rows which should be removed
    std::vector<PetscInt> src_row;

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
          // FloatMetal-Insulator interface at Insulator side
        case InsulatorRegion:
          {
            // record the source row
            src_row.push_back(fvm_nodes[i]->global_offset());

            break;
          }
          // FloatMetal-Insulator interface at Conductor side
        case ConductorRegion:
          {
            // record the source row
            src_row.push_back(fvm_nodes[i]->global_offset());

            break;
          }
        default: genius_error(); //we should never reach here
        }
      }

    }

    // for efficient resion, we separate MatAssemblyBegin and MatAssemblyEnd
    MatAssemblyEnd(*jac, MAT_FINAL_ASSEMBLY);

    // clear source rows
    MatZeroRows(*jac, src_row.size(), src_row.empty() ? NULL : &src_row[0], 0.0);

  }


  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar z_width = this->z_width();

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
          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          // psi of this node
          AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

          // psi of float metal
          AutoDScalar  Vm = x[this->local_offset()];        Vm.setADValue(1,1.0);

          // the psi of this node is equal to corresponding psi of float metal
          AutoDScalar ff = V - Vm;

          // set Jacobian of governing equation ff
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), this->global_offset(), ff.getADValue(1), ADD_VALUES);


          // process float metal surface equation

          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
          FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_nodes[i]->neighbor_node_end();
          for(; nb_it != nb_it_end; ++nb_it)
          {
            const FVM_Node *nb_node = (*nb_it).second;
            // the psi of neighbor node
            AutoDScalar V_nb = x[nb_node->local_offset()];  V_nb.setADValue(1,1.0);
            // distance from nb node to this node
            PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();
            // area of out surface of control volume related with neighbor node,
            // here we should consider the difference of 2D/3D by multiply z_width
            PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node->root_node())*z_width;
            // surface electric field and electrc displacement
            AutoDScalar E = (V_nb-V)/distance;
            AutoDScalar D = node_data->eps()*E;

            // since the governing equation is \sum \integral D + \sigma = 0
            // and \integral D = \sum cv_boundary*D
            // here we use MatSetValue ADD_VALUES to process \sum \integral D one by one
            MatSetValue(*jac, this->global_offset(), fvm_nodes[i]->global_offset(), (cv_boundary*D).getADValue(0), ADD_VALUES);
            MatSetValue(*jac, this->global_offset(), nb_node->global_offset(), (cv_boundary*D).getADValue(1), ADD_VALUES);
          }

          break;

        }
        // Electrode-Insulator interface at Conductor side
      case ConductorRegion:
        {

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          // psi of this node
          AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0,1.0);

          // psi of float metal
          AutoDScalar  Vm = x[this->local_offset()];        Vm.setADValue(1,1.0);

          // the psi of this node is equal to corresponding psi of float metal
          AutoDScalar ff = V - Vm;

          // set Jacobian of governing equation ff
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), fvm_nodes[i]->global_offset(), ff.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_nodes[i]->global_offset(), this->global_offset(), ff.getADValue(1), ADD_VALUES);

          break;
        }

      default: genius_error(); //we should never reach here
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}


/*---------------------------------------------------------------------
 * update potential of float metal
 */
void ChargedContactBC::DDM1_Update_Solution(PetscScalar *lxx)
{
  this->psi() = lxx[this->local_offset()];
}

