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

// Local includes
#include "simulation_system.h"
#include "resistance_region.h"
#include "boundary_condition_solderpad.h"
#include "parallel.h"

using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::kb;
using PhysicalUnit::e;


/*---------------------------------------------------------------------
 * fill scaling value
 */
void SolderPadBC::Mix_DDM1_Fill_Value(Vec x, Vec L)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;
      VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
    }
  }
}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for Mixed DDML1 solver
 */
void SolderPadBC::Mix_DDM1_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;

      PetscInt row = fvm_node->global_offset();
      clear_row.push_back(row);
    }
  }
}

/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML1 solver
 */
void SolderPadBC::Mix_DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Gate boundary condition is processed here


  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // the electrode current, since the electrode may be partitioned into several processor,
  // we should collect it.
  std::vector<double> current_buffer;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width()/A;

  // the electrode potential in current iteration
  PetscScalar Ve = x[this->local_offset()];

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const MetalSimulationRegion * resistance_region = 0;
  if( _r1 && _r1->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r1);
  if( _r2 && _r2->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r2);
  genius_assert(resistance_region);

  const double workfunction = resistance_region->material()->basic->Affinity(T_external());
  const double sigma = resistance_region->material()->basic->Conductance();

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;
      const FVM_NodeData * node_data = fvm_node->node_data();

      switch ( region->type() )
      {
          case MetalRegion :
          {
            // psi of this node
            PetscScalar V = x[fvm_node->local_offset()+0];
            PetscScalar ff = V + node_data->affinity()/e - Ve;

            // set governing equation to function vector
            VecSetValue(f, fvm_node->global_offset(), ff, ADD_VALUES);

            // conductance current
            FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
            for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
            {
              const FVM_Node *nb_node = (*nb_it).first;
              const FVM_NodeData * nb_node_data = nb_node->node_data();
              // the psi of neighbor node
              PetscScalar V_nb = x[nb_node->local_offset()+0];
              // distance from nb node to this node
              PetscScalar distance = fvm_node->distance(nb_node);
              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = std::abs(fvm_node->cv_surface_area(nb_node));

              // current flow
              current_buffer.push_back( cv_boundary*sigma*(V-V_nb)/distance*current_scale );
            }
            break;
          }

          case InsulatorRegion:
          {
            // psi of this node
            PetscScalar V = x[fvm_node->local_offset()];
            PetscScalar f_psi = (V + workfunction - Ve);

            // set governing equation to function vector
            VecSetValue(f, fvm_node->global_offset(), f_psi, ADD_VALUES);
            break;
          }
          default: genius_error();

      }
    }
  }


  // for get the current, we must sum all the terms in current_buffer
  PetscScalar current = std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  // Add current to spice node
  VecSetValue(f, this->global_offset(), current, ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for Mixed DDML1 solver
 */
void SolderPadBC::Mix_DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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

    BoundaryCondition::const_node_iterator node_it = nodes_begin();
    BoundaryCondition::const_node_iterator end_it = nodes_end();
    for(; node_it!=end_it; ++node_it )
    {
      // skip node not belongs to this processor
      if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

      const FVM_Node * fvm_node = (*region_node_begin(*node_it)).second.second;
      MatSetValue(*jac, fvm_node->global_offset()+0, this->global_offset(), 0, ADD_VALUES);
    }


    // reserve jacobian entries for the circuit equation
    {

      std::vector<PetscInt> bc_node_reserve;
      for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
      {
        // get the derivative of electrode current to ohmic node
        const FVM_Node *  fvm_node = get_region_fvm_node(*node_it, MetalRegion);
        if(fvm_node->on_processor())
        {
          bc_node_reserve.push_back(fvm_node->global_offset()+0);

          // get the derivative of electrode current to neighbors of bc node
          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
          FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
          for(; nb_it != nb_it_end; ++nb_it)
          {
            const FVM_Node *  fvm_nb_node = (*nb_it).first;
            bc_node_reserve.push_back(fvm_nb_node->global_offset()+0);
          }
        }
      }
      Parallel::allgather(bc_node_reserve);

      if(Genius::processor_id() == Genius::n_processors()-1)
      {
        PetscInt bc_global_offset = this->global_offset();

        MatSetValue(*jac, bc_global_offset, this->global_offset(), 0, ADD_VALUES);

        if(bc_node_reserve.size())
        {
          std::vector<PetscScalar> bc_node_reserve_zero(bc_node_reserve.size(), 0.0);
          MatSetValues(*jac, 1, &bc_global_offset, bc_node_reserve.size(), &bc_node_reserve[0], &bc_node_reserve_zero[0], ADD_VALUES);
        }
      }
    }

  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void SolderPadBC::Mix_DDM1_Jacobian_Preprocess(PetscScalar *, Mat *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;

      PetscInt row = fvm_node->global_offset();
      clear_row.push_back(row);
    }
  }
}




/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML1 solver
 */
void SolderPadBC::Mix_DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of GateContact boundary condition is processed here
  // we use AD again. no matter it is overkill here.

  PetscInt bc_global_offset = this->global_offset();

  // after that, we should do gate boundary process here

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width()/A;

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const MetalSimulationRegion * resistance_region = 0;
  if( _r1 && _r1->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r1);
  if( _r2 && _r2->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r2);
  genius_assert(resistance_region);

  const double workfunction = resistance_region->material()->basic->Affinity(T_external());
  const double sigma = resistance_region->material()->basic->Conductance();


  //the indepedent variable number, we only need 2 here.
  adtl::AutoDScalar::numdir=2;

  // loop again
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;
      const FVM_NodeData * node_data = fvm_node->node_data();

      switch ( region->type() )
      {
          case MetalRegion :
          {
            // psi of this node
            AutoDScalar V = x[fvm_node->local_offset()];  V.setADValue(0, 1.0);

            // the electrode potential in current iteration
            AutoDScalar Ve = x[this->local_offset()];     Ve.setADValue(1, 1.0);

            AutoDScalar ff = V + node_data->affinity()/e - Ve;

            //governing equation
            MatSetValue(*jac, fvm_node->global_offset(), fvm_node->global_offset(), ff.getADValue(0), ADD_VALUES);
            MatSetValue(*jac, fvm_node->global_offset(), bc_global_offset, ff.getADValue(1), ADD_VALUES);


            // conductance current
            FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
            for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
            {
              const FVM_Node *nb_node = (*nb_it).first;
              const FVM_NodeData * nb_node_data = nb_node->node_data();

              // the psi of neighbor node
              AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);

              // distance from nb node to this node
              PetscScalar distance = fvm_node->distance(nb_node);

              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = std::abs(fvm_node->cv_surface_area(nb_node));


              AutoDScalar current = cv_boundary*sigma*(V-V_nb)/distance*current_scale;

              MatSetValue(*jac, bc_global_offset, fvm_node->global_offset(), current.getADValue(0), ADD_VALUES);
              MatSetValue(*jac, bc_global_offset, nb_node->global_offset(), current.getADValue(1), ADD_VALUES);
            }
            break;
          }

          case InsulatorRegion :
          {
            // psi of this node
            AutoDScalar V = x[fvm_node->local_offset()];  V.setADValue(0, 1.0);

            // the electrode potential in current iteration
            AutoDScalar Ve = x[this->local_offset()];     Ve.setADValue(1, 1.0);

            AutoDScalar ff = (V + workfunction - Ve);

            //governing equation
            MatSetValue(*jac, fvm_node->global_offset(), fvm_node->global_offset(), ff.getADValue(0), ADD_VALUES);
            MatSetValue(*jac, fvm_node->global_offset(), bc_global_offset, ff.getADValue(1), ADD_VALUES);
            break;
          }
          default: genius_error();

      }
    }
  }
  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


