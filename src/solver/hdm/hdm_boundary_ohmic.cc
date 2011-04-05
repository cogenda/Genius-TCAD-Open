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
#include "boundary_condition_ohmic.h"

void OhmicContactBC::HDM_Boundary( const PetscScalar * , Vec x, InsertMode &add_value_flag  )
{

  PetscScalar T = T_external();

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
      const FVM_NodeData * node_data = fvm_node->node_data();

      switch ( region->type() )
      {
        // Semiconductor Region of course owns OhmicContactBC
      case SemiconductorRegion:
        {
          // semiconductor region should be the first region
          genius_assert(i==0);

          const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(region);
          // mapping this node to material library
          semi_region->material()->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);
          PetscScalar ni = semi_region->material()->band->ni(T);

          PetscInt loc[8];
          for(int i=0; i<8; ++i)
            loc[i] = fvm_node->global_offset() + i;

          std::vector<PetscScalar> value(8, 0);
          PetscScalar  & electron_density = value[0];
          PetscScalar  & hole_density = value[4];

          PetscScalar  net_dpoing = node_data->Net_doping();
          if( net_dpoing < 0 )                   //p-type
          {
            hole_density = (-net_dpoing + sqrt(net_dpoing*net_dpoing + 4*ni*ni))/2.0;
            electron_density = ni*ni/hole_density;
          }
          else                                  //n-type
          {
            electron_density = (net_dpoing + sqrt(net_dpoing*net_dpoing + 4*ni*ni))/2.0;
            hole_density = ni*ni/electron_density;
          }

          VecSetValues(x, 8, &loc[0], &value[0], INSERT_VALUES);

          break;
        }
      default: break;
      }
    }
  }

  add_value_flag = INSERT_VALUES;
}


