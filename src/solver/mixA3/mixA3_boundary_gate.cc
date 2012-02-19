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
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_gate.h"
#include "spice_ckt.h"
#include "parallel.h"
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::V;
using PhysicalUnit::A;




///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * do pre-process to function for EBM3 solver
 */
void GateContactBC::MixA_EBM3_Function_Preprocess(PetscScalar *,Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    // should clear all the rows related with this boundary condition
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
          // insulator region.
          case InsulatorRegion:
          {
            clear_row.push_back(fvm_nodes[i]->global_offset()+0);
            break;
          }

          case ElectrodeRegion:
          {
            clear_row.push_back(fvm_nodes[i]->global_offset()+0);

            if(regions[i]->get_advanced_model()->enable_Tl())
            {
              src_row.push_back(fvm_nodes[i]->global_offset()+1);
              dst_row.push_back(fvm_nodes[0]->global_offset()+1);
              clear_row.push_back(fvm_nodes[i]->global_offset()+1);
            }
            break;
          }

      }
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void GateContactBC::MixA_EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Gate boundary condition is processed here
  SPICE_CKT * ckt = this->system().get_circuit();
  unsigned int spice_node_index = ckt->get_spice_node_by_bc(this);

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }


  // data buffer for INSERT_VALUES
  std::vector<int> iy;
  std::vector<PetscScalar> y;


  // the electrode current, since the electrode may be partitioned into several processor,
  // we should collect it.
  std::vector<double> current_buffer;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width()/A;

  const PetscScalar Work_Function = this->scalar("workfunction");
  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");

  // the electrode potential in current iteration
  PetscScalar Ve = x[ckt->local_offset_x(spice_node_index)];

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
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
          // insulator region.
          case InsulatorRegion:
          {
            unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
            unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

            const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

            // psi of this node
            PetscScalar V = x[fvm_nodes[i]->local_offset()+node_psi_offset];
            // the governing equation of potential
            PetscScalar ff1 = V + Work_Function - Ve;
            // set governing equation to function vector
            iy.push_back(fvm_nodes[i]->global_offset() + node_psi_offset);
            y.push_back( ff1 );

            // add heat flux out of gate boundary to lattice temperature equatiuon
            if( regions[i]->get_advanced_model()->enable_Tl() && (node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion)) )
            {
              // lattice temperature
              PetscScalar T = x[fvm_nodes[i]->local_offset() + node_Tl_offset];
              PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
              PetscScalar fT = Heat_Transfer*(T_external()-T)*S;
              VecSetValue(f, fvm_nodes[i]->global_offset()+node_Tl_offset, fT, ADD_VALUES);
            }

            // MOS gate can have displacement current and tunneling current

            // displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                FVM_Node *nb_node = (*nb_it).second;
                FVM_NodeData * nb_node_data = nb_node->node_data();
                // the psi of neighbor node
                PetscScalar V_nb = x[nb_node->local_offset()+0];
                // distance from nb node to this node
                PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();
                // area of out surface of control volume related with neighbor node
                PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node->root_node());
                PetscScalar dEdt;
                if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false) //second order
                {
                  PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
                  dEdt = ( (2-r)/(1-r)*(V-V_nb)
                           - 1.0/(r*(1-r))*(node_data->psi()-nb_node_data->psi())
                           + (1-r)/r*(node_data->psi_last()-nb_node_data->psi_last()))/distance/(SolverSpecify::dt_last+SolverSpecify::dt);
                }
                else//first order
                {
                  dEdt = ((V-V_nb)-(node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;
                }

                current_buffer.push_back( cv_boundary*node_data->eps()*dEdt*current_scale );
              }
            }
            // FIXME: tunneling current.

            break;
          }


          // conductor region
          case ElectrodeRegion:
          {
            unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
            unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

            unsigned int insuregion_node_psi_offset = regions[0]->ebm_variable_offset(POTENTIAL);
            unsigned int insuregion_node_Tl_offset  = regions[0]->ebm_variable_offset(TEMPERATURE);

            {
              // psi of this node
              PetscScalar V = x[fvm_nodes[i]->local_offset()+node_psi_offset];
              // since the region is sorted, we know region[0] is Insulator region
              // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding Insulator region
              PetscScalar V_in = x[fvm_nodes[0]->local_offset()+insuregion_node_psi_offset];
              // the psi of this node is equal to corresponding psi of Insulator node
              PetscScalar ff1 = V - V_in;
              // set governing equation to function vector
              y.push_back( ff1 );
              iy.push_back(fvm_nodes[i]->global_offset()+node_psi_offset);
            }

            if(regions[i]->get_advanced_model()->enable_Tl())
            {
              // lattice temperature of this node
              PetscScalar T = x[fvm_nodes[i]->local_offset()+node_Tl_offset];
              PetscScalar T_in = x[fvm_nodes[0]->local_offset()+insuregion_node_Tl_offset];

              // the T of this node is equal to corresponding T of Insulator node
              PetscScalar ff2 = T - T_in;
              y.push_back( ff2 );
              iy.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);

            }

            break;
          }
          case VacuumRegion:
          break;
          default: genius_error(); //we should never reach here
      }
    }

  }


  if(iy.size()) VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES) ;

  // for get the current, we must sum all the terms in current_buffer
  PetscScalar current = std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  // Add current to spice node
  VecSetValue(f, this->global_offset(), current, ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for EBM3 solver
 */
void GateContactBC::MixA_EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.
  SPICE_CKT * ckt = this->system().get_circuit();
  unsigned int spice_node_index = ckt->get_spice_node_by_bc(this);

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
          case InsulatorRegion:
          {
            genius_assert(i==0);
            // insert none zero pattern
            // none zero pattern includes bd node and their neighbors!
            unsigned int n_node_var  = regions[i]->ebm_n_variables();
            unsigned int global_offset   = fvm_nodes[i]->global_offset();
            unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

            for(unsigned int nv=0; nv<n_node_var; ++nv)
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+nv, this->global_offset(), 0, ADD_VALUES);

            // reserve for heat transport equation
            if(regions[i]->get_advanced_model()->enable_Tl())
            {
              FVM_Node::fvm_ghost_node_iterator gn_it = fvm_nodes[i]->ghost_node_begin();
              FVM_Node::fvm_ghost_node_iterator gn_it_end = fvm_nodes[i]->ghost_node_end();
              for(; gn_it != gn_it_end; ++gn_it)
              {
                const FVM_Node * ghost_fvm_node = (*gn_it).first;
                // skip NULL neighbor which means the node is on Neumann boundary
                if(ghost_fvm_node==NULL) continue;

                const SimulationRegion * ghost_region = this->system().region((*gn_it).second.first);
                genius_assert(ghost_region!=NULL);
                unsigned int ghostregion_node_Tl_offset  = ghost_region->ebm_variable_offset(TEMPERATURE);

                MatSetValue(*jac, global_offset+node_Tl_offset, ghost_fvm_node->global_offset()+ghostregion_node_Tl_offset, 0,ADD_VALUES);

                FVM_Node::fvm_neighbor_node_iterator  gnb_it = ghost_fvm_node->neighbor_node_begin();
                for(; gnb_it != ghost_fvm_node->neighbor_node_end(); ++gnb_it)
                  MatSetValue(*jac, global_offset+node_Tl_offset, (*gnb_it).second->global_offset()+ghostregion_node_Tl_offset, 0, ADD_VALUES);
              }
            }

            break;
          }
          case ElectrodeRegion:
          {
            unsigned int global_offset   = fvm_nodes[i]->global_offset();
            unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
            unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

            unsigned int insuregion_node_psi_offset = regions[0]->ebm_variable_offset(POTENTIAL);
            unsigned int insuregion_node_Tl_offset  = regions[0]->ebm_variable_offset(TEMPERATURE);

            // insert none zero pattern
            MatSetValue(*jac, global_offset+node_psi_offset, fvm_nodes[0]->global_offset()+insuregion_node_psi_offset, 0, ADD_VALUES);

            if(regions[i]->get_advanced_model()->enable_Tl())
              MatSetValue(*jac, global_offset+node_Tl_offset,  fvm_nodes[0]->global_offset()+insuregion_node_Tl_offset, 0, ADD_VALUES);

            break;
          }
          case VacuumRegion:
          break;
          default: genius_error(); //we should never reach here
      }
    }
  }

  // reserve jacobian entries for the circuit equation of gate electrode
  {
    std::vector<PetscInt> bc_node_reserve;
    for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // get the derivative of electrode current to gate node
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
      for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
      {
        const SimulationRegion * region = ( (*rnode_it).second.first );
        if(region->type() != InsulatorRegion) continue;

        const FVM_Node * fvm_node = (*rnode_it).second.second;

        if(fvm_node->on_processor())
        {
          for(unsigned int nv=0; nv<region->ebm_n_variables(); ++nv)
            bc_node_reserve.push_back(fvm_node->global_offset()+nv);

          FVM_Node::fvm_neighbor_node_iterator nb_it     =  fvm_node->neighbor_node_begin();
          FVM_Node::fvm_neighbor_node_iterator nb_it_end =  fvm_node->neighbor_node_end();
          for(; nb_it!=nb_it_end; ++nb_it)
          {
            const FVM_Node *  fvm_nb_node = (*nb_it).second;
            for(unsigned int nv=0; nv<region->ebm_n_variables(); ++nv)
              bc_node_reserve.push_back(fvm_nb_node->global_offset()+nv);
          }
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

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for EBM3 solver
 */
void GateContactBC::MixA_EBM3_Jacobian_Preprocess(PetscScalar * ,Mat *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);

    // should clear all the rows related with this boundary condition
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      regions.push_back( (*rnode_it).second.first );
      fvm_nodes.push_back( (*rnode_it).second.second );

      switch ( regions[i]->type() )
      {
          // insulator region.
          case InsulatorRegion:
          {
            clear_row.push_back(fvm_nodes[i]->global_offset()+0);
            break;
          }

          case ElectrodeRegion:
          {
            clear_row.push_back(fvm_nodes[i]->global_offset()+0);

            if(regions[i]->get_advanced_model()->enable_Tl())
            {
              src_row.push_back(fvm_nodes[i]->global_offset()+1);
              dst_row.push_back(fvm_nodes[0]->global_offset()+1);
              clear_row.push_back(fvm_nodes[i]->global_offset()+1);
            }
            break;
          }

      }
    }
  }
}




/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void GateContactBC::MixA_EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of GateContact boundary condition is processed here
  // we use AD again. no matter it is overkill here.

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  SPICE_CKT * ckt = this->system().get_circuit();
  unsigned int spice_node_index = ckt->get_spice_node_by_bc(this);
  PetscInt bc_global_offset = this->global_offset();


  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width()/A;

  const PetscScalar Work_Function = this->scalar("workfunction");
  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");

  // loop again
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
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
          case InsulatorRegion:
          {
            unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
            unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

            const FVM_NodeData * node_data = fvm_nodes[i]->node_data();

            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            // psi of this node
            AutoDScalar  V = x[fvm_nodes[i]->local_offset()+node_psi_offset]; V.setADValue(0, 1.0);

            // the electrode potential in current iteration
            AutoDScalar Ve = x[ckt->local_offset_x(spice_node_index)];                         Ve.setADValue(1, 1.0);

            // the governing equation of potential
            AutoDScalar ff = V + Work_Function - Ve;

            // the insert position
            PetscInt row     = fvm_nodes[i]->global_offset()+node_psi_offset;
            PetscInt cols[2] = {fvm_nodes[i]->global_offset()+node_psi_offset, this->global_offset()};
            // process the Jacobian of governing equation of potential
            MatSetValues(*jac, 1, &row, 2, &cols[0], ff.getADValue(), ADD_VALUES);

            // process the Jacobian of equation of T
            if(regions[i]->get_advanced_model()->enable_Tl() && (node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion)) )
            {
              //the indepedent variable number, we need 1 here.
              adtl::AutoDScalar::numdir=1;

              AutoDScalar T = x[fvm_nodes[i]->local_offset()+node_Tl_offset]; T.setADValue(0, 1.0); // psi of this node

              PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
              AutoDScalar fT = Heat_Transfer*(T_external()-T)*S;

              PetscInt row = fvm_nodes[i]->global_offset()+node_Tl_offset;
              MatSetValue(*jac, row, row, fT.getADValue(0), ADD_VALUES);
            }

            /*
             * process the Jacobian of current flow out of gate electrode
             */

            // compute displacement current

            //the indepedent variable number, we need 2 here.
            adtl::AutoDScalar::numdir=2;

            // displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                FVM_Node *nb_node = (*nb_it).second;
                FVM_NodeData * nb_node_data = nb_node->node_data();

                // the psi of this node
                AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0, 1.0);
                // the psi of neighbor node
                AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);

                // distance from nb node to this node
                PetscScalar distance = (*(fvm_nodes[i]->root_node()) - *(nb_node->root_node())).size();

                // area of out surface of control volume related with neighbor node
                PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node->root_node());
                AutoDScalar dEdt;
                if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false) //second order
                {
                  PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
                  dEdt = ( (2-r)/(1-r)*(V-V_nb)
                           - 1.0/(r*(1-r))*(node_data->psi()-nb_node_data->psi())
                           + (1-r)/r*(node_data->psi_last()-nb_node_data->psi_last()))/distance/(SolverSpecify::dt_last+SolverSpecify::dt);
                }
                else//first order
                {
                  dEdt = ((V-V_nb)-(node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;
                }

                AutoDScalar current_disp = cv_boundary*node_data->eps()*dEdt*current_scale;

                // consider electrode connect
                MatSetValue(*jac, bc_global_offset, fvm_nodes[i]->global_offset(), current_disp.getADValue(0), ADD_VALUES);
                MatSetValue(*jac, bc_global_offset, nb_node->global_offset(), current_disp.getADValue(1), ADD_VALUES);
              }
            }
            //FIXME tunneling current should be considered here

            break;
          }
          // conductor region (gate) which has an interface with insulator region
          case ElectrodeRegion:
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
            PetscInt cols[2] = {fvm_nodes[i]->global_offset()+node_psi_offset, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(POTENTIAL)};
            MatSetValues(*jac, 1, &row, 2, &cols[0], ff1.getADValue(), ADD_VALUES);

            // the T of this node is equal to corresponding T of insulator node
            if(regions[i]->get_advanced_model()->enable_Tl())
            {
              // lattice temperature of this node
              AutoDScalar  T = x[fvm_nodes[i]->local_offset()+node_Tl_offset]; T.setADValue(0,1.0);
              AutoDScalar  T_in = x[fvm_nodes[0]->local_offset()+regions[0]->ebm_variable_offset(TEMPERATURE)]; T_in.setADValue(1,1.0);
              AutoDScalar  ff2 = T - T_in;
              PetscInt row     = fvm_nodes[i]->global_offset()+node_Tl_offset;
              PetscInt cols[2] = {fvm_nodes[i]->global_offset()+node_Tl_offset, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE)};
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




