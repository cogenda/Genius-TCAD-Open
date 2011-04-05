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
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_charge.h"
#include "petsc_utils.h"



void ChargedContactBC::DDMAC_Fill_Matrix_Vector( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag )
{

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar z_width = this->z_width();

  // Jacobian entrance
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node
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
          unsigned int n_variables     = regions[i]->ebm_n_variables();
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          // psi of this node
          AutoDScalar  V = node_data->psi();        V.setADValue(0,1.0);

          // psi of float metal
          AutoDScalar  Vm = this->psi();            Vm.setADValue(1,1.0);

          // the psi of this node is equal to corresponding psi of float metal
          AutoDScalar ff1 = V - Vm;

          // set Jacobian of governing equation ff
          PetscInt real_row = fvm_nodes[i]->global_offset() + node_psi_offset;
          PetscInt imag_row = fvm_nodes[i]->global_offset() + node_psi_offset + n_variables;
          PetscInt real_col[2]={real_row, this->global_offset()};
          PetscInt imag_col[2]={imag_row, this->global_offset()+1};

          MatSetValues(A, 1, &real_row, 2, real_col, ff1.getADValue(), ADD_VALUES);
          MatSetValues(A, 1, &real_row, 2, imag_col, ff1.getADValue(), ADD_VALUES);

          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            // find the ghost region and ghost node
            // since we know only one ghost node exit, there is ghost_node_begin()
            const SimulationRegion *  ghost_region = get_fvm_node_region(fvm_nodes[i]->root_node(), ElectrodeRegion);
            FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
            const FVM_Node * ghost_fvm_node = (*gn_it).first;

            regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], TEMPERATURE, A, b, J, omega, add_value_flag, ghost_region, ghost_fvm_node);

            unsigned int ghost_node_Tl_offset  = ghost_region->ebm_variable_offset(TEMPERATURE);

            // T of this node
            AutoDScalar  T = node_data->T(); T.setADValue(0,1.0);

            // T for corresponding float metal region
            AutoDScalar  T_elec = ghost_fvm_node->node_data()->T(); T_elec.setADValue(1,1.0);

            // the T of this node is equal to corresponding T of conductor node
            // we assuming T is continuous for the interface
            AutoDScalar ff2 = T - T_elec;

            // set Jacobian of governing equation ff2
            PetscInt real_row = fvm_nodes[i]->global_offset() + node_Tl_offset;
            PetscInt imag_row = fvm_nodes[i]->global_offset() + node_Tl_offset + n_variables;
            PetscInt real_col[2]={real_row, ghost_fvm_node->global_offset()+ghost_node_Tl_offset};
            PetscInt imag_col[2]={imag_row, ghost_fvm_node->global_offset()+ghost_node_Tl_offset+ghost_region->ebm_n_variables()};

            MatSetValues(A, 1, &real_row, 2, real_col, ff2.getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &imag_row, 2, imag_col, ff2.getADValue(), ADD_VALUES);
          }

          // process float metal surface equation

          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
          FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_nodes[i]->neighbor_node_end();
          for(; nb_it != nb_it_end; ++nb_it)
          {
            const FVM_Node *nb_node = (*nb_it).second;
            // the psi of neighbor node
            AutoDScalar V_nb = nb_node->node_data()->psi();  V_nb.setADValue(1,1.0);
            // distance from nb node to this node
            PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();
            // area of out surface of control volume related with neighbor node,
            // here we should consider the difference of 2D/3D by multiply z_width
            PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node->root_node())*z_width;
            // surface electric field and electrc displacement
            AutoDScalar D = node_data->eps()*(V_nb-V)/distance;
            // since the governing equation is \sum \integral D + \sigma = 0
            // and \integral D = \sum cv_boundary*D
            // here we use MatSetValue ADD_VALUES to process \sum \integral D one by one
            PetscInt real_row = this->global_offset();
            PetscInt imag_row = this->global_offset()+1;
            PetscInt real_col[2]={fvm_nodes[i]->global_offset()+node_psi_offset, nb_node->global_offset()+node_psi_offset};
            PetscInt imag_col[2]={fvm_nodes[i]->global_offset()+node_psi_offset+n_variables, nb_node->global_offset()+node_psi_offset+n_variables};

            MatSetValues(A, 1, &real_row, 2, real_col, (cv_boundary*D).getADValue(), ADD_VALUES);
            MatSetValues(A, 1, &imag_row, 2, imag_col, (cv_boundary*D).getADValue(), ADD_VALUES);
          }

          break;

        }
        // Electrode-Insulator interface at Conductor side
      case ElectrodeRegion:
        {
	  unsigned int n_variables     = regions[i]->ebm_n_variables();
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          // psi of this node
          AutoDScalar  V = fvm_nodes[i]->node_data()->psi(); V.setADValue(0,1.0);

          // psi of float metal
          AutoDScalar  Vm = this->psi();        Vm.setADValue(1,1.0);

          // the psi of this node is equal to corresponding psi of float metal
          AutoDScalar ff = V - Vm;

          // set Jacobian of governing equation ff
	  PetscInt real_row = fvm_nodes[i]->global_offset() + node_psi_offset;
          PetscInt imag_row = fvm_nodes[i]->global_offset() + node_psi_offset + n_variables;
          PetscInt real_col[2]={real_row, this->global_offset()};
          PetscInt imag_col[2]={imag_row, this->global_offset()+1};

          MatSetValues(A, 1, &real_row, 2, real_col, ff.getADValue(), ADD_VALUES);
          MatSetValues(A, 1, &imag_row, 2, imag_col, ff.getADValue(), ADD_VALUES);

	  if(regions[i]->get_advanced_model()->enable_Tl())
          {
            regions[i]->DDMAC_Fill_Nodal_Matrix_Vector(fvm_nodes[i], TEMPERATURE, A, b, J, omega, add_value_flag);
	  }

          break;
        }

      default: genius_error(); //we should never reach here
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
