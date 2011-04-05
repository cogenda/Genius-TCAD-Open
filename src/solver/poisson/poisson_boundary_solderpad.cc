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
#include "resistance_region.h"
#include "boundary_condition_solderpad.h"
#include "parallel.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;


/*---------------------------------------------------------------------
 * set scaling constant
 */
void SolderPadBC::Poissin_Fill_Value(Vec , Vec L)
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
      switch ( region->type() )
      {
        case MetalRegion :
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


/*---------------------------------------------------------------------
 * do pre-process to function evaluation for poisson solver
 */
void SolderPadBC::Poissin_Function_Preprocess(Vec f, std::vector<PetscInt> &src_row,
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
      switch ( region->type() )
      {
        case MetalRegion :
        {
          const FVM_Node * fvm_node = ( *rnode_it ).second.second;
          PetscInt row = fvm_node->global_offset();
          clear_row.push_back(row);
          break;
        }

        case InsulatorRegion:
        {
          const FVM_Node * fvm_node = ( *rnode_it ).second.second;
          PetscInt row = fvm_node->global_offset();
          clear_row.push_back(row);
          break;
        }
        default: break;
      }
    }
  }
}


/*---------------------------------------------------------------------
 * build function for poisson solver
 */
void SolderPadBC::Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // get the workfunction and sigma
  const SimulationRegion * region = bc_regions().first; genius_assert(region);
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *>(region); genius_assert(resistance_region);
  const PetscScalar workfunction = resistance_region->material()->basic->Affinity(T_external());

  // the electrode potential should be zero
  PetscScalar Ve = ext_circuit()->Vapp();


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
      switch ( region->type() )
      {
          case MetalRegion :
          {

            const FVM_Node * fvm_node = ( *rnode_it ).second.second;
            const FVM_NodeData * node_data = fvm_node->node_data();
            // phi of this node
            PetscScalar V = x[fvm_node->local_offset()+0];
            PetscScalar f_phi = V + workfunction - Ve;

            // set governing equation to function vector
            VecSetValue(f, fvm_node->global_offset(), f_phi, ADD_VALUES);
            break;
          }

          case InsulatorRegion:
          {
            const FVM_Node * fvm_node = ( *rnode_it ).second.second;
            // phi of this node
            PetscScalar V = x[fvm_node->local_offset()];
            PetscScalar f_phi = (V + workfunction - Ve);

            // set governing equation to function vector
            VecSetValue(f, fvm_node->global_offset(), f_phi, ADD_VALUES);
            break;
          }
          default: genius_error();
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for poisson solver
 */
void SolderPadBC::Poissin_Jacobian_Preprocess(Mat *jac, std::vector<PetscInt> &src_row,
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
 * build function and its jacobian for poisson solver
 */
void SolderPadBC::Poissin_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // the Jacobian of SolderPad boundary condition is processed here

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  const SimulationRegion * region = bc_regions().first; genius_assert(region);
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *>(region);
  const PetscScalar workfunction = resistance_region->material()->basic->Affinity(T_external());

  PetscScalar Ve = ext_circuit()->Vapp();

  // we use AD again. no matter it is overkill here.
  //the indepedent variable number, we only need 1 here.
  adtl::AutoDScalar::numdir=1;

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

      switch ( region->type() )
      {
          case MetalRegion :
          {
            const FVM_Node * fvm_node = ( *rnode_it ).second.second;
            const FVM_NodeData * node_data = fvm_node->node_data();

            // phi of this node
            AutoDScalar V = x[fvm_node->local_offset()];  V.setADValue(0, 1.0);
            AutoDScalar f_phi = V + workfunction - Ve;

            //governing equation
            MatSetValue(*jac, fvm_node->global_offset(), fvm_node->global_offset(), f_phi.getADValue(0), ADD_VALUES);
            break;
          }
          case InsulatorRegion :
          {
            const FVM_Node * fvm_node = ( *rnode_it ).second.second;

            // phi of this node
            AutoDScalar V = x[fvm_node->local_offset()];  V.setADValue(0, 1.0);
            AutoDScalar f_phi = (V + workfunction - Ve);

            //governing equation
            MatSetValue(*jac, fvm_node->global_offset(), fvm_node->global_offset(), f_phi.getADValue(0), ADD_VALUES);
            break;
          }
          default: genius_error();

      }
    }
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * update electrode potential
 */
void SolderPadBC::Poissin_Update_Solution(PetscScalar *)
{
  this->ext_circuit()->potential() = this->ext_circuit()->Vapp();
}

