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
#include "insulator_region.h"
#include "boundary_condition_ii.h"
#include "petsc_utils.h"



void InsulatorInsulatorInterfaceBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag )
{

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first insulator region
      if(i==0)
      {
        //load Jacobian entry of this node from J and fill into AC matrix A
        regions[0]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[0], POTENTIAL, A, b, J, omega, add_value_flag);

        if(regions[i]->get_advanced_model()->enable_Tl())
          regions[0]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[0], TEMPERATURE, A, b, J, omega, add_value_flag);
      }

      // other insulator region
      else
      {
        // load Jacobian entry of this node from J and fill into AC matrix A with the location of fvm_node[0]
        regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], POTENTIAL, A, b, J, omega, add_value_flag, regions[0], fvm_nodes[0]);
        regions[i]->DDMAC_Force_equal(fvm_nodes[i], POTENTIAL, A, add_value_flag, regions[0], fvm_nodes[0]);

        if(regions[i]->get_advanced_model()->enable_Tl())
        {
          regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], TEMPERATURE, A, b, J, omega, add_value_flag, regions[0], fvm_nodes[0]);
          regions[i]->DDMAC_Force_equal(fvm_nodes[i], TEMPERATURE, A, add_value_flag, regions[0], fvm_nodes[0]);
        }

      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
