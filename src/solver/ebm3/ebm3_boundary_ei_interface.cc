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
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////



/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void ElectrodeInsulatorInterfaceBC::EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Electrode-Insulator interface is processed here

  // note, we will use INSERT_VALUES to set values of vec f
  // if the previous operator is not insert_VALUES, we should assembly the vec
  if( (add_value_flag != INSERT_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // buffer for Vec location
  std::vector<PetscInt> src_row;
  std::vector<PetscInt> dst_row;

  // buffer for Vec location
  std::vector<PetscInt> iy;

  // buffer for Vec value
  std::vector<PetscScalar> y_new;

  // search for all the node with this boundary type
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
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
        // Electrode-Insulator interface at Insulator side
      case InsulatorRegion:
        {
          // Insulator region should be the first region
          genius_assert(i==0);

          // a node on the Interface of Insulator-Semiconductor should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);

          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          // record the source row and dst row
          dst_row.push_back(fvm_nodes[i]->global_offset()+node_psi_offset);

          if(regions[i]->get_advanced_model()->enable_Tl())
            dst_row.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);

          break;
        }
        // Electrode-Insulator interface at Conductor side
      case ConductorRegion:
        {
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          src_row.push_back(fvm_nodes[i]->global_offset()+node_psi_offset);

          // psi of this node
          PetscScalar V = x[fvm_nodes[i]->local_offset()+node_psi_offset];
          // since the region is sorted, we know region[0] is Insulator region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding Insulator region
          PetscScalar V_in = x[fvm_nodes[0]->local_offset()+regions[0]->ebm_variable_offset(POTENTIAL)];
          // the psi of this node is equal to corresponding psi of Insulator node
          PetscScalar ff1 = V - V_in;

          // set governing equation to function vector
          iy.push_back(fvm_nodes[i]->global_offset()+node_psi_offset);
          y_new.push_back(ff1);

          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            src_row.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);

            // lattice temperature of this node
            PetscScalar T = x[fvm_nodes[i]->local_offset()+node_Tl_offset];
            PetscScalar T_in = x[fvm_nodes[0]->local_offset()+regions[0]->ebm_variable_offset(TEMPERATURE)];
            // the T of this node is equal to corresponding T of Insulator node
            PetscScalar ff2 = T - T_in;

            iy.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);
            y_new.push_back(ff2);
          }

          break;
        }
      case VacuumRegion:
        break;
      default: genius_error(); //we should never reach here
      }
    }

  }

  // add src row to dst row, it will assemble vec automatically
  PetscUtils::VecAddRowToRow(f, src_row, dst_row);

  // insert new value to src row
  if( iy.size() )
    VecSetValues(f, iy.size(), &(iy[0]), &(y_new[0]), INSERT_VALUES);

  add_value_flag = INSERT_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for EBM3 solver
 */
void ElectrodeInsulatorInterfaceBC::EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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


    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {

        // Electrode-Insulator interface at Insulator side
      case InsulatorRegion:
        {
          // a node on the Interface of Insulator-Semiconductor should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);

          unsigned int global_offset   = fvm_nodes[i]->global_offset();
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          // reserve items for all the ghost nodes
          FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
          for(; gn_it != fvm_nodes[i]->ghost_node_end(); ++gn_it)
          {
            const FVM_Node * ghost_fvm_node = (*gn_it).first;

            const SimulationRegion * ghost_region = this->system().region((*gn_it).second.first);
            unsigned int ghostregion_node_psi_offset = ghost_region->ebm_variable_offset(POTENTIAL);
            unsigned int ghostregion_node_Tl_offset  = ghost_region->ebm_variable_offset(TEMPERATURE);

            MatSetValue(*jac, global_offset+node_psi_offset, ghost_fvm_node->global_offset()+ghostregion_node_psi_offset, 0, ADD_VALUES);

            if(regions[i]->get_advanced_model()->enable_Tl())
              MatSetValue(*jac, global_offset+node_Tl_offset, ghost_fvm_node->global_offset()+ghostregion_node_Tl_offset, 0, ADD_VALUES);

            FVM_Node::fvm_neighbor_node_iterator  gnb_it = ghost_fvm_node->neighbor_node_begin();
            for(; gnb_it != ghost_fvm_node->neighbor_node_end(); ++gnb_it)
            {
              MatSetValue(*jac, global_offset+node_psi_offset, (*gnb_it).second->global_offset()+ghostregion_node_psi_offset, 0, ADD_VALUES);

              if(regions[i]->get_advanced_model()->enable_Tl())
                MatSetValue(*jac, global_offset+node_Tl_offset, (*gnb_it).second->global_offset()+ghostregion_node_Tl_offset, 0, ADD_VALUES);
            }
          }

          break;
        }
        // Electrode-Insulator interface at Conductor side
      case ConductorRegion:
        {
          // Conductor region should be the second region
          genius_assert(i==1);
          // a node on the Interface of Electrode-Insulator should only have one ghost node.
          genius_assert(fvm_nodes[i]->n_ghost_node()==1);

          unsigned int global_offset   = fvm_nodes[i]->global_offset();
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);
          // insert none zero pattern
          MatSetValue(*jac, global_offset+node_psi_offset, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(POTENTIAL), 0, ADD_VALUES);

          if(regions[i]->get_advanced_model()->enable_Tl())
            MatSetValue(*jac, global_offset+node_Tl_offset,  fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE), 0, ADD_VALUES);

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
 * build function and its jacobian for EBM3 solver
 */
void ElectrodeInsulatorInterfaceBC::EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // Jacobian of Electrode-Insulator interface is processed here

  {
    // buffer for mat rows which should be added to other row
    std::vector<PetscInt> src_row;
    std::vector<PetscInt> dst_row;

    // buffer for mat rows which should be removed
    std::vector<PetscInt> rm_row;

    // search for all the node with this boundary type
    BoundaryCondition::const_node_iterator node_it = nodes_begin();
    BoundaryCondition::const_node_iterator end_it = nodes_end();
    for(; node_it!=end_it; ++node_it )
    {
      // skip node not belongs to this processor
      if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

      // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
      // but in different region.
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
      for(; rnode_it!=end_rnode_it;  ++rnode_it  )
      {
        const FVM_Node *  fvm_node = (*rnode_it).second.second ;
        const SimulationRegion * region = (*rnode_it).second.first;

        switch ( region->type() )
        {
          // Electrode-Insulator interface at Insulator side
        case InsulatorRegion:
          {
            // a node on the Interface of Insulator-Semiconductor should only have one ghost node.
            genius_assert(fvm_node->n_ghost_node()==1);

            // record the source row and dst row

            dst_row.push_back(fvm_node->global_offset() + region->ebm_variable_offset(POTENTIAL));

            if(region->get_advanced_model()->enable_Tl())
              dst_row.push_back(fvm_node->global_offset() + region->ebm_variable_offset(TEMPERATURE));


            break;
          }
          // Electrode-Insulator interface at Conductor side
        case ConductorRegion:
          {

            src_row.push_back(fvm_node->global_offset() + region->ebm_variable_offset(POTENTIAL));
            if(region->get_advanced_model()->enable_Tl())
              src_row.push_back(fvm_node->global_offset()+ region->ebm_variable_offset(TEMPERATURE));


            rm_row.push_back(fvm_node->global_offset() + region->ebm_variable_offset(POTENTIAL));
            if(region->get_advanced_model()->enable_Tl())
              rm_row.push_back(fvm_node->global_offset() + region->ebm_variable_offset(TEMPERATURE));

            break;
          }
        case VacuumRegion:
          break;
        default: genius_error(); //we should never reach here
        }
      }

    }

    //we add source rows to destination rows
    PetscUtils::MatAddRowToRow(*jac, src_row, dst_row);

    // clear rows
    MatZeroRows(*jac, rm_row.size(), rm_row.empty() ? NULL : &rm_row[0], 0.0);

  }


  // after that, set new Jacobian entrance to source rows
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
        // Electrode-Insulator interface at Insulator side
      case InsulatorRegion:
        {
          break;
        }
        // Electrode-Insulator interface at Conductor side
      case ConductorRegion:
        {
          unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
          unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

          //the indepedent variable number, we need 2 here.
          adtl::AutoDScalar::numdir=2;

          AutoDScalar  V = x[fvm_nodes[i]->local_offset()+node_psi_offset]; V.setADValue(0,1.0); // psi of this node
          // since the region is sorted, we know region[0] is insulator region
          // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding insulator region
          AutoDScalar  V_in = x[fvm_nodes[0]->local_offset()+regions[0]->ebm_variable_offset(POTENTIAL)]; V_in.setADValue(1,1.0);
          // the psi of this node is equal to corresponding psi of insulator node
          AutoDScalar  ff1 = V - V_in;
          PetscInt row     = fvm_nodes[i]->global_offset()+node_psi_offset;
          PetscInt cols[2] = {row, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(POTENTIAL)};
          MatSetValues(*jac, 1, &row, 2, &cols[0], ff1.getADValue(), ADD_VALUES);

          // the T of this node is equal to corresponding T of insulator node
          if(regions[i]->get_advanced_model()->enable_Tl())
          {
            // lattice temperature of this node
            AutoDScalar  T = x[fvm_nodes[i]->local_offset()+node_Tl_offset]; T.setADValue(0,1.0);
            AutoDScalar  T_in = x[fvm_nodes[0]->local_offset()+regions[0]->ebm_variable_offset(TEMPERATURE)]; T_in.setADValue(1,1.0);
            AutoDScalar  ff2 = T - T_in;
            PetscInt row     = fvm_nodes[i]->global_offset()+node_Tl_offset;
            PetscInt cols[2] = {row, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE)};
            MatSetValues(*jac, 1, &row, 2, &cols[0], ff2.getADValue(), ADD_VALUES);
          }

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
