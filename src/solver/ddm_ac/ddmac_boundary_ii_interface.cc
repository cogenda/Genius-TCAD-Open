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
#include "boundary_condition.h"
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
      if(!fvm_node->is_valid()) continue;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first insulator region
      if(i==0)
      {

        /*
        * load Jacobian entry of this node from J and fill into AC matrix A
         */

        regions[0]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[0], A, b, J, omega, add_value_flag);

      }

      // other insulator region
      else
      {

        /*
        * load Jacobian entry of this node from J and fill into AC matrix A with the location of fvm_node[0]
         */

        regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], A, b, J, omega, add_value_flag, regions[0], fvm_nodes[0]);

        /*
         *  let psi and T of node in this region equal to psi(T) of node in the first insulator region
         */

        unsigned int n_variables     = regions[i]->ebm_n_variables();
        unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
        unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

        //the indepedent variable number, we need 2 here.
        adtl::AutoDScalar::numdir=2;

        {
          // psi of this node
          AutoDScalar  V    = fvm_nodes[i]->node_data()->psi(); V.setADValue(0,1.0);

          // psi for ghost node
          AutoDScalar  V_in = fvm_nodes[0]->node_data()->psi(); V_in.setADValue(1,1.0);

          // the psi of this node is equal to corresponding psi of insulator node in the other region
          AutoDScalar  ff1 = V - V_in;

          // set Jacobian of governing equation ff1
          PetscInt real_row = fvm_nodes[i]->global_offset() + node_psi_offset;
          PetscInt imag_row = fvm_nodes[i]->global_offset() + n_variables + node_psi_offset;
          PetscInt real_col[2]={real_row, fvm_nodes[0]->global_offset() + node_psi_offset};
          PetscInt imag_col[2]={imag_row, fvm_nodes[0]->global_offset() + n_variables + node_psi_offset};

          MatSetValues(A, 1, &real_row, 2, real_col, ff1.getADValue(), ADD_VALUES);
          MatSetValues(A, 1, &imag_row, 2, imag_col, ff1.getADValue(), ADD_VALUES);
        }

        if(regions[i]->get_advanced_model()->enable_Tl())
        {
          // T of this node
          AutoDScalar  T = fvm_nodes[i]->node_data()->T(); T.setADValue(0,1.0);

          // T of insulator node in the other region
          AutoDScalar  T_in = fvm_nodes[0]->node_data()->T(); T_in.setADValue(1,1.0);

          // the T of this node is equal to corresponding T of insulator node in the other region
          AutoDScalar ff2 = T - T_in;

          // set Jacobian of governing equation ff1
          PetscInt real_row = fvm_nodes[i]->global_offset() + node_Tl_offset;
          PetscInt imag_row = fvm_nodes[i]->global_offset() + n_variables + node_Tl_offset;
          PetscInt real_col[2]={fvm_nodes[i]->global_offset() + node_Tl_offset, fvm_nodes[0]->global_offset() + node_Tl_offset};
          PetscInt imag_col[2]={fvm_nodes[i]->global_offset() + n_variables + node_Tl_offset, fvm_nodes[0]->global_offset() + n_variables + node_Tl_offset};

          MatSetValues(A, 1, &real_row, 2, real_col, ff2.getADValue(), ADD_VALUES);
          MatSetValues(A, 1, &imag_row, 2, imag_col, ff2.getADValue(), ADD_VALUES);
        }

      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
