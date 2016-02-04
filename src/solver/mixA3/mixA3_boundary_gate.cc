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
          case MetalRegion:
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
                FVM_Node *nb_node = (*nb_it).first;
                FVM_NodeData * nb_node_data = nb_node->node_data();
                // the psi of neighbor node
                PetscScalar V_nb = x[nb_node->local_offset()+0];
                // distance from nb node to this node
                PetscScalar distance = fvm_nodes[i]->distance(nb_node);
                // area of out surface of control volume related with neighbor node
                PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node);
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
          case MetalRegion:
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
  this->current() = current*A;

  // Add current to spice node
  VecSetValue(f, this->global_offset(), current, ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}





/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for EBM3 solver
 */
void GateContactBC::MixA_EBM3_Jacobian_Preprocess(PetscScalar * ,SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
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
          case MetalRegion:
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
void GateContactBC::MixA_EBM3_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // the Jacobian of GateContact boundary condition is processed here
  // we use AD again. no matter it is overkill here.


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
            jac->add_row(  row,  2,  &cols[0],  ff.getADValue() );

            // process the Jacobian of equation of T
            if(regions[i]->get_advanced_model()->enable_Tl() && (node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion)) )
            {
              //the indepedent variable number, we need 1 here.
              adtl::AutoDScalar::numdir=1;

              AutoDScalar T = x[fvm_nodes[i]->local_offset()+node_Tl_offset]; T.setADValue(0, 1.0); // psi of this node

              PetscScalar S  = fvm_nodes[i]->outside_boundary_surface_area();
              AutoDScalar fT = Heat_Transfer*(T_external()-T)*S;

              PetscInt row = fvm_nodes[i]->global_offset()+node_Tl_offset;
              jac->add( row,  row,  fT.getADValue(0) );
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
                FVM_Node *nb_node = (*nb_it).first;
                FVM_NodeData * nb_node_data = nb_node->node_data();

                // the psi of this node
                AutoDScalar  V = x[fvm_nodes[i]->local_offset()]; V.setADValue(0, 1.0);
                // the psi of neighbor node
                AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);

                // distance from nb node to this node
                PetscScalar distance = fvm_nodes[i]->distance(nb_node);

                // area of out surface of control volume related with neighbor node
                PetscScalar cv_boundary = fvm_nodes[i]->cv_surface_area(nb_node);
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
                jac->add( bc_global_offset,  fvm_nodes[i]->global_offset(),  current_disp.getADValue(0) );
                jac->add( bc_global_offset,  nb_node->global_offset(),  current_disp.getADValue(1) );
              }
            }
            //FIXME tunneling current should be considered here

            break;
          }
          // conductor region (gate) which has an interface with insulator region
          case MetalRegion:
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
            jac->add_row(  row,  2,  &cols[0],  ff1.getADValue() );

            // the T of this node is equal to corresponding T of insulator node
            if(regions[i]->get_advanced_model()->enable_Tl())
            {
              // lattice temperature of this node
              AutoDScalar  T = x[fvm_nodes[i]->local_offset()+node_Tl_offset]; T.setADValue(0,1.0);
              AutoDScalar  T_in = x[fvm_nodes[0]->local_offset()+regions[0]->ebm_variable_offset(TEMPERATURE)]; T_in.setADValue(1,1.0);
              AutoDScalar  ff2 = T - T_in;
              PetscInt row     = fvm_nodes[i]->global_offset()+node_Tl_offset;
              PetscInt cols[2] = {fvm_nodes[i]->global_offset()+node_Tl_offset, fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE)};
              jac->add_row(  row,  2,  &cols[0],  ff2.getADValue() );
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



/*---------------------------------------------------------------------
 * update electrode IV
 */
void GateContactBC::MixA_EBM3_Update_Solution(PetscScalar *)
{
  Parallel::sum(this->current());
}



