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
#include "semiconductor_region.h"
#include "conductor_region.h"
#include "insulator_region.h"
#include "boundary_condition_ir.h"
#include "parallel.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


/*---------------------------------------------------------------------
 * set scaling constant
 */
void ResistanceInsulatorBC::DDM1_Fill_Value(Vec , Vec L)
{
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
      if( region->type() != InsulatorRegion ) continue;

      const FVM_Node * fvm_node = (*rnode_it).second.second;
      VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
    }
  }
}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


/*---------------------------------------------------------------------
 * do pre-process to function for DDML1 solver
 */
void ResistanceInsulatorBC::DDM1_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      if( region->type() != InsulatorRegion ) continue;

      const FVM_Node * insulator_fvm_node = (*rnode_it).second.second;
      clear_row.push_back(insulator_fvm_node->global_offset());
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void ResistanceInsulatorBC::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  // Resistance-Insulator boundary condition is processed here

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // buffer for Vec location
  std::vector<PetscInt> index;
  index.reserve(n_nodes());
  // buffer for Vec value
  std::vector<PetscScalar> current_buffer;
  current_buffer.reserve(n_nodes());

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  const PetscScalar current_scale = this->z_width();

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // we should only have one resistance region
    const FVM_Node * resistance_fvm_node = get_region_fvm_node((*node_it), MetalRegion);
    const PetscScalar V_resistance = x[resistance_fvm_node->local_offset()];

    // however, this resistance region may contact several insulator region

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      if( region->type() != InsulatorRegion ) continue;

      const FVM_Node * insulator_fvm_node = (*rnode_it).second.second;
      const FVM_NodeData * insulator_fvm_node_data = insulator_fvm_node->node_data();
      const PetscScalar V_insulator  = x[insulator_fvm_node->local_offset()];

      // in any case, the phi in insulator region should equal to phi in resistance region.
      PetscScalar f_phi =  V_insulator - V_resistance;
      VecSetValue(f, insulator_fvm_node->global_offset(), f_phi, ADD_VALUES);

      //for transiet simulation, we should consider displacement current from insulator region
      if(SolverSpecify::TimeDependent == true)
      {
        PetscScalar I_displacement = 0.0;
        FVM_Node::fvm_neighbor_node_iterator nb_it = insulator_fvm_node->neighbor_node_begin();
        for(; nb_it != insulator_fvm_node->neighbor_node_end(); ++nb_it)
        {
          const FVM_Node *nb_node = (*nb_it).second;
          const FVM_NodeData * nb_node_data = nb_node->node_data();
          // the psi of neighbor node
          PetscScalar V_nb = x[nb_node->local_offset()];
          // distance from nb node to this node
          PetscScalar distance = insulator_fvm_node->distance(nb_node);
          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = insulator_fvm_node->cv_surface_area(nb_node->root_node());
          PetscScalar dEdt;
          if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false) //second order
          {
            PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
            dEdt = ( (2-r)/(1-r)*(V_insulator-V_nb)
                     - 1.0/(r*(1-r))*(insulator_fvm_node_data->psi()-nb_node_data->psi())
                     + (1-r)/r*(insulator_fvm_node_data->psi_last()-nb_node_data->psi_last()))/distance/(SolverSpecify::dt_last+SolverSpecify::dt);
          }
          else//first order
          {
            dEdt = ((V_insulator-V_nb)-(insulator_fvm_node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;
          }

          I_displacement += cv_boundary*insulator_fvm_node_data->eps()*dEdt;
        }
        index.push_back(resistance_fvm_node->global_offset());
        current_buffer.push_back(-I_displacement);
      }

    }
  }

  if( index.size() )
    VecSetValues(f, index.size(), &(index[0]), &(current_buffer[0]), ADD_VALUES);


  // for get the current, we must sum all the terms in current_buffer
  this->current() = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML1 solver
 */
void ResistanceInsulatorBC::DDM1_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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

    // we should only have one resistance region
    const FVM_Node * resistance_fvm_node = get_region_fvm_node((*node_it), MetalRegion);

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      if( region->type() != InsulatorRegion ) continue;
      // reserve for later operator
      {
        MatSetValue(*jac, fvm_node->global_offset(), resistance_fvm_node->global_offset(), 0, ADD_VALUES);

        //reserve for displacement current
        MatSetValue(*jac, resistance_fvm_node->global_offset(), fvm_node->global_offset(), 0, ADD_VALUES);
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
        {
          const FVM_Node *nb_node = (*nb_it).second;
          MatSetValue(*jac, resistance_fvm_node->global_offset(), nb_node->global_offset(), 0, ADD_VALUES);
        }
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML1 solver
 */
void ResistanceInsulatorBC::DDM1_Jacobian_Preprocess(PetscScalar *, Mat *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(; node_it!=end_it; ++node_it )
  {
      // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      if( region->type() != InsulatorRegion ) continue;

      const FVM_Node * insulator_fvm_node = (*rnode_it).second.second;
      clear_row.push_back(insulator_fvm_node->global_offset());
    }
  }
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void ResistanceInsulatorBC::DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of Resistance-Insulator boundary condition is processed here
  // we use AD again. no matter it is overkill here.

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  PetscScalar dt = SolverSpecify::dt;

  //the indepedent variable number, we need 2 here.
  adtl::AutoDScalar::numdir=2;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * resistance_fvm_node = get_region_fvm_node((*node_it), MetalRegion);
    AutoDScalar V_resistance = x[resistance_fvm_node->local_offset()]; V_resistance.setADValue(0, 1.0);

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      if( region->type() != InsulatorRegion ) continue;
      const FVM_Node * insulator_fvm_node = (*rnode_it).second.second;
      const FVM_NodeData * insulator_fvm_node_data = insulator_fvm_node->node_data();

      AutoDScalar V_insulator  = x[insulator_fvm_node->local_offset()];  V_insulator.setADValue(1, 1.0);

      // the psi of this node is equal to corresponding psi of Insulator node
      AutoDScalar f_phi =  V_insulator - V_resistance;

      // set Jacobian of governing equation ff. since these rows are erased, ADD_VALUES is equal to INSERT_VALUES
      MatSetValue(*jac, insulator_fvm_node->global_offset(), resistance_fvm_node->global_offset(), f_phi.getADValue(0), ADD_VALUES);
      MatSetValue(*jac, insulator_fvm_node->global_offset(), insulator_fvm_node->global_offset(), f_phi.getADValue(1), ADD_VALUES);


      // displacement current
      if(SolverSpecify::TimeDependent == true)
      {
        FVM_Node::fvm_neighbor_node_iterator nb_it = insulator_fvm_node->neighbor_node_begin();
        for(; nb_it != insulator_fvm_node->neighbor_node_end(); ++nb_it)
        {
          const FVM_Node *nb_node = (*nb_it).second;
          const FVM_NodeData * nb_node_data = nb_node->node_data();

          // the psi of neighbor node
          AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(0, 1.0);

          // distance from nb node to this node
          PetscScalar distance = insulator_fvm_node->distance(nb_node);

          // area of out surface of control volume related with neighbor node
          PetscScalar cv_boundary = insulator_fvm_node->cv_surface_area(nb_node->root_node());
          AutoDScalar dEdt;
          if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_LowerOrder==false) //second order
          {
            PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
            dEdt = ( (2-r)/(1-r)*(V_insulator-V_nb)
                     - 1.0/(r*(1-r))*(insulator_fvm_node_data->psi()-nb_node_data->psi())
                     + (1-r)/r*(insulator_fvm_node_data->psi_last()-nb_node_data->psi_last()))/distance/(SolverSpecify::dt_last+SolverSpecify::dt);
          }
          else//first order
          {
            dEdt = ((V_insulator-V_nb)-(insulator_fvm_node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;
          }

          AutoDScalar current_disp = cv_boundary*insulator_fvm_node_data->eps()*dEdt;

          MatSetValue(*jac, resistance_fvm_node->global_offset(), insulator_fvm_node->global_offset(), -current_disp.getADValue(1), ADD_VALUES);
          MatSetValue(*jac, resistance_fvm_node->global_offset(), nb_node->global_offset(), -current_disp.getADValue(0), ADD_VALUES);
        }
      }

    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * update electrode IV
 */
void ResistanceInsulatorBC::DDM1_Update_Solution(PetscScalar *)
{
  Parallel::sum(this->current());
}


