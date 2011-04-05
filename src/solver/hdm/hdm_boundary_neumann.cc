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
#include "boundary_condition_neumann.h"

#include "hdm_flux.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


void NeumannBC::HDM_Ghostcell_Volume( Vec vol )
{
  // search and process all the boundary nodes
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
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      const SimulationRegion * region = (*rnode_it).second.first;

      switch ( region->type() )
      {
        // Semiconductor Region of course owns OhmicContactBC
        case SemiconductorRegion:
        {
          PetscScalar vol_init[8]; //(1.0/(2.0*fvm_node->volume()), 8);
          PetscInt index[8];
          for(unsigned int i=0; i<8; ++i)
          {
            vol_init[i] = 1.0/(2.0*fvm_node->volume());
            index[i] = fvm_node->global_offset() + i;
          }
          VecSetValues(vol, 8, index, vol_init, INSERT_VALUES);

          break;
        }
        default: break;
      }
    }
  }

}



void NeumannBC::HDM_Boundary( const PetscScalar * x, Vec flux, InsertMode &  )
{

  // search and process all the boundary nodes
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
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      const SimulationRegion * region = (*rnode_it).second.first;

      switch ( region->type() )
      {
        // Semiconductor Region of course owns OhmicContactBC
      case SemiconductorRegion:
        {
          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);

          const PetscScalar mn = semi_region->material()->band->EffecElecMass(T_external());
          const PetscScalar mp = semi_region->material()->band->EffecHoleMass(T_external());

          // first, find the norm vector
          const VectorValue<Real> & norm = fvm_node->norm();
          //for all the neighbor FVM_Node, build the ghost node with same density but mirrow'd velocity
          // that is Vg=Vi-2(Vi*n)n

          PetscInt loc[8];
          for(int i=0; i<8; ++i)
            loc[i] = fvm_node->global_offset() + i;

          HDMVector Un1(&x[fvm_node->local_offset()]);
          HDMVector Up1(&x[fvm_node->local_offset()+4]);


          FVM_Node::fvm_neighbor_node_iterator  neighbor_begin = fvm_node->neighbor_node_begin();
          FVM_Node::fvm_neighbor_node_iterator  neighbor_end = fvm_node->neighbor_node_end();
          for(; neighbor_begin!=neighbor_end; ++neighbor_begin )
          {
            const FVM_Node * neigbor_node = (*neighbor_begin).second;
            genius_assert(neigbor_node->on_local());

            Real d = (*(fvm_node->root_node()) - *(neigbor_node->root_node())).size();
            Real S = fvm_node->cv_surface_area(neigbor_node->root_node());
            VectorValue<double> dir = (*(neigbor_node->root_node()) - *(fvm_node->root_node())).unit();
            dir -= 2*(dir*norm)*norm;

            PetscScalar local_dt1, local_dt2;

            HDMVector Un2(&x[neigbor_node->local_offset()]);   Un2.mirror(norm);
            HDMVector Up2(&x[neigbor_node->local_offset()+4]); Up2.mirror(norm);

            HDMVector fn = AUSM_if_flux(Un1, Un2, mn, kb, d, S, dir, local_dt1);
            HDMVector fp = AUSM_if_flux(Up1, Up2, mp, kb, d, S, dir, local_dt2);

            VecSetValues(flux, 4, &loc[0], &(fn[0]), ADD_VALUES);
            VecSetValues(flux, 4, &loc[4], &(fp[0]), ADD_VALUES);
          }

          break;
        }
      default: break;
      }
    }
  }

}



