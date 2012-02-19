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
#include "boundary_condition_ohmic.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


void OhmicContactBC::LinearPoissin_Reserve(Mat A, InsertMode &add_value_flag)
{}

void OhmicContactBC::LinearPoissin_Matrix(Mat A, InsertMode &add_value_flag)
{
 // the indicator which rows should be set to zero
  std::vector<PetscInt> id;
  id.reserve(n_nodes());

  //note! PetscUtils::MatZeroRows should be excuted on all the processor
  //no matter whether it owns this row!
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    // should clear all the rows related with ohmic boundary condition
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node *  fvm_node = (*rnode_it).second.second ;
      switch ( region->type() )
      {
        case SemiconductorRegion:
        {
          PetscInt row = fvm_node->global_offset();
          id.push_back(row);
          break;
        }
        case ElectrodeRegion:
        case InsulatorRegion:
        {
          PetscInt row = fvm_node->global_offset();
          id.push_back(row);
          break;
        }
        case VacuumRegion:
          break;
        default: genius_error(); //we should never reach here
      }
    }
  }

  // for efficient resion, we separate MatAssemblyBegin and MatAssemblyEnd
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  // zero required rows
  PetscUtils::MatZeroRows(A, id.size(), id.empty() ? NULL : &id[0], 1.0);
}


void OhmicContactBC::LinearPoissin_RHS(Vec b, InsertMode &add_value_flag)
{
 // Ohmic boundary condition is processed here

  // note, we will use INSERT_VALUES to set values of vec f
  // if the previous operator is not insert_VALUES, we should assembly the vec
  if( (add_value_flag != INSERT_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
  }

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
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      const FVM_NodeData * node_data = fvm_node->node_data();

      switch ( region->type() )
      {
        // Semiconductor Region of course owns OhmicContactBC
        case SemiconductorRegion:
        {
          // semiconductor region should be the first region
          genius_assert(i==0);

          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);

          PetscScalar T = T_external();
          semi_region->material()->mapping(fvm_node->root_node(), node_data, 0.0);
          PetscScalar ni  = semi_region->material()->band->ni(T);

          PetscScalar Nc  = node_data->Nc();
          PetscScalar Nv  = node_data->Nv();
          PetscScalar Eg  = node_data->Eg();

          //governing equation for Ohmic contact boundary
          PetscScalar V =  kb*T/e*boost::math::asinh(node_data->Net_doping()/(2*ni))
                         - Eg/(2*e)
                         - kb*T*log(Nc/Nv)/(2*e)
                         - node_data->affinity()/e
                         + ext_circuit()->Vapp();

          // set governing equation to function vector
          VecSetValue(b, fvm_node->global_offset(), V, INSERT_VALUES);
          break;
        }
        // conductor region which has an interface with OhmicContact boundary to semiconductor region
        case ElectrodeRegion:
        {
          // set governing equation to function vector
          VecSetValue(b, fvm_node->global_offset(), 0.0, INSERT_VALUES);

          break;
        }
        // insulator region. if a corner where semiconductor region, insulator region and  conductor region meet.
        // the boundary for the corner point may be Ohmic. (not a nice behavier)
        case InsulatorRegion:
        {
          // set governing equation to function vector
          VecSetValue(b, fvm_node->global_offset(), 0.0, INSERT_VALUES);

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


