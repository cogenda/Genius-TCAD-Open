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
#include "resistance_region.h"
#include "boundary_condition_solderpad.h"
#include "spice_ckt.h"
#include "parallel.h"

using PhysicalUnit::V;
using PhysicalUnit::A;
using PhysicalUnit::kb;
using PhysicalUnit::e;


/*---------------------------------------------------------------------
 * fill initial and scaling value
 */
void SolderPadBC::MixA_DDM2_Fill_Value(Vec x, Vec L)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;
      VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
    }
  }
}




///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for Mixed DDML2 solver
 */
void SolderPadBC::MixA_DDM2_Function_Preprocess(PetscScalar * ,Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;

      PetscInt row = fvm_node->global_offset();
      clear_row.push_back(row);
    }
  }
}

/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML2 solver
 */
void SolderPadBC::MixA_DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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

  // the electrode current, since the electrode may be partitioned into several processor,
  // we should collect it.
  std::vector<double> current_buffer;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  const PetscScalar current_scale = this->z_width()/A;

  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");

  // the electrode potential in current iteration
  PetscScalar Ve = x[ckt->local_offset_x(spice_node_index)];

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const MetalSimulationRegion * resistance_region = 0;
  if( _r1 && _r1->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r1);
  if( _r2 && _r2->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r2);
  genius_assert(resistance_region);

  const double workfunction = resistance_region->material()->basic->Affinity(T_external());
  const double sigma = resistance_region->material()->basic->Conductance();

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;
      const FVM_NodeData * node_data = fvm_node->node_data();

      switch ( region->type() )
      {
      case MetalRegion :
        {
          // psi of this node
          PetscScalar V = x[fvm_node->local_offset()+0];
          // T of this node
          PetscScalar T = x[fvm_node->local_offset()+1];

          PetscScalar f_psi = V + node_data->affinity()/e - Ve;

          // add heat flux out of boundary to lattice temperature equatiuon
          PetscScalar f_q = Heat_Transfer*(T_external()-T)*fvm_node->outside_boundary_surface_area();

          // set governing equation to function vector
          VecSetValue(f, fvm_node->global_offset()+0, f_psi, ADD_VALUES);
          VecSetValue(f, fvm_node->global_offset()+1, f_q, ADD_VALUES);

          // conductance current
          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
          for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
          {
            const FVM_Node *nb_node = (*nb_it).first;
            const FVM_NodeData * nb_node_data = nb_node->node_data();
            // psi of neighbor node
            PetscScalar V_nb = x[nb_node->local_offset()];
            // T of neighbor node
            PetscScalar T_nb = x[nb_node->local_offset()+1];
            // distance from nb node to this node
            PetscScalar distance = fvm_node->distance(nb_node);
            // area of out surface of control volume related with neighbor node
            PetscScalar cv_boundary = std::abs(fvm_node->cv_surface_area(nb_node));
            // current density 
            PetscScalar current_density = resistance_region->material()->basic->CurrentDensity((V-V_nb)/distance, 0.5*(T+T_nb));

            // current flow
            current_buffer.push_back( cv_boundary*current_density*current_scale );
          }
          break;
        }

      case InsulatorRegion:
        {
          const FVM_Node * fvm_node = ( *rnode_it ).second.second;
          // psi of this node
          PetscScalar V = x[fvm_node->local_offset()];
          PetscScalar f_psi = (V + workfunction - Ve);

          // assume heat flux out of boundary is zero

          // set governing equation to function vector
          VecSetValue(f, fvm_node->global_offset(), f_psi, ADD_VALUES);

          // displacement current
          if(SolverSpecify::TimeDependent == true)
          {
            FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
            for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
            {
              const FVM_Node *nb_node = (*nb_it).first;
              const FVM_NodeData * nb_node_data = nb_node->node_data();
              // the psi of neighbor node
              PetscScalar V_nb = x[nb_node->local_offset()+0];
              // distance from nb node to this node
              PetscScalar distance = (*(fvm_node->root_node()) - *(nb_node->root_node())).size();
              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = fvm_node->cv_surface_area(nb_node);
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

          break;
        }
      default: genius_error();
      }
    }
  }

  // for get the current, we must sum all the terms in current_buffer
  PetscScalar current = std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );
  this->current() = current*A;

  VecSetValue(f, this->global_offset(), current, ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}






/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for DDML2 solver
 */
void SolderPadBC::MixA_DDM2_Jacobian_Preprocess(PetscScalar *,SparseMatrix<PetscScalar> *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;

      PetscInt row = fvm_node->global_offset();
      clear_row.push_back(row);
    }
  }
}




/*---------------------------------------------------------------------
 * build function and its jacobian for Mixed DDML2 solver
 */
void SolderPadBC::MixA_DDM2_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // the Jacobian of GateContact boundary condition is processed here
  // we use AD again. no matter it is overkill here.
  SPICE_CKT * ckt = this->system().get_circuit();
  unsigned int spice_node_index = ckt->get_spice_node_by_bc(this);

  PetscInt bc_global_offset = this->global_offset();

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  const PetscScalar current_scale = this->z_width()/A;

  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");
  
  //the indepedent variable number, we only need 4 here.
  adtl::AutoDScalar::numdir=4;

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const MetalSimulationRegion * resistance_region = 0;
  if( _r1 && _r1->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r1);
  if( _r2 && _r2->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r2);
  genius_assert(resistance_region);
  
  resistance_region->material()->set_ad_num(adtl::AutoDScalar::numdir);

  const double workfunction = resistance_region->material()->basic->Affinity(T_external());




  // loop again
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;
      const FVM_NodeData * node_data = fvm_node->node_data();

      switch ( region->type() )
      {
      case MetalRegion :
        {
          // psi of this node
          AutoDScalar V = x[fvm_node->local_offset()];  V.setADValue(0, 1.0);

          // the electrode potential in current iteration
          genius_assert( local_offset()!=invalid_uint );
          AutoDScalar Ve = x[ckt->local_offset_x(spice_node_index)];     Ve.setADValue(1, 1.0);

          AutoDScalar ff = V + node_data->affinity()/e - Ve;

          //governing equation of psi
          jac->add( fvm_node->global_offset(),  fvm_node->global_offset(),  ff.getADValue(0) );
          jac->add( fvm_node->global_offset(),  bc_global_offset,  ff.getADValue(1) );


          // T of this node
          AutoDScalar T = x[fvm_node->local_offset()+1];  T.setADValue(1, 1.0);

          // add heat flux out of boundary to lattice temperature equatiuon
          AutoDScalar f_q = Heat_Transfer*(T_external()-T)*fvm_node->outside_boundary_surface_area();
          //governing equation of T
          jac->add( fvm_node->global_offset()+1,  fvm_node->global_offset()+1,  f_q.getADValue(1) );

          // conductance current
          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
          for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
          {
            const FVM_Node *nb_node = (*nb_it).first;
            const FVM_NodeData * nb_node_data = nb_node->node_data();

            // the psi of neighbor node
            AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(2, 1.0);
            // T of neighbor node
            AutoDScalar T_nb = x[nb_node->local_offset()+1]; T_nb.setADValue(3, 1.0);

            // distance from nb node to this node
            PetscScalar distance = fvm_node->distance(nb_node);

            // area of out surface of control volume related with neighbor node
            PetscScalar cv_boundary = std::abs(fvm_node->cv_surface_area(nb_node));
            
            // current density
            AutoDScalar current_density = resistance_region->material()->basic->CurrentDensity((V-V_nb)/distance, 0.5*(T+T_nb));
            
            AutoDScalar current = cv_boundary*current_density*current_scale;

            jac->add( bc_global_offset,  fvm_node->global_offset(),   current.getADValue(0) );
            jac->add( bc_global_offset,  fvm_node->global_offset()+1, current.getADValue(1) );
            jac->add( bc_global_offset,  nb_node->global_offset(),    current.getADValue(2) );
            jac->add( bc_global_offset,  nb_node->global_offset()+1,  current.getADValue(3) );
            
          }
          break;
        }
      case InsulatorRegion :
        {
          // psi of this node
          AutoDScalar V = x[fvm_node->local_offset()];  V.setADValue(0, 1.0);

          // the electrode potential in current iteration
          genius_assert( local_offset()!=invalid_uint );
          AutoDScalar Ve = x[ckt->local_offset_x(spice_node_index)];     Ve.setADValue(1, 1.0);

          AutoDScalar ff = (V + workfunction - Ve);

          //governing equation
          jac->add( fvm_node->global_offset(),  fvm_node->global_offset(),  ff.getADValue(0) );
          jac->add( fvm_node->global_offset(),  bc_global_offset,  ff.getADValue(1) );

          // displacement current
          if(SolverSpecify::TimeDependent == true)
          {
            FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
            for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
            {
              const FVM_Node *nb_node = (*nb_it).first;
              const FVM_NodeData * nb_node_data = nb_node->node_data();

              // the psi of this node
              AutoDScalar  V = x[fvm_node->local_offset()]; V.setADValue(0, 1.0);
              // the psi of neighbor node
              AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);

              // distance from nb node to this node
              PetscScalar distance = fvm_node->distance(nb_node);

              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = fvm_node->cv_surface_area(nb_node);
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

              jac->add( bc_global_offset,  fvm_node->global_offset(),  current_disp.getADValue(0) );
              jac->add( bc_global_offset,  nb_node->global_offset(),  current_disp.getADValue(1) );
            }
          }
          break;
        }
      default: genius_error();

      }
    }
  }
  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * update electrode IV
 */
void SolderPadBC::MixA_DDM2_Update_Solution(PetscScalar *)
{
  Parallel::sum(this->current());
}


