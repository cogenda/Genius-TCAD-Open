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
#include "resistance_region.h"
#include "boundary_condition_solderpad.h"
#include "parallel.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;


/*---------------------------------------------------------------------
 * fill electrode potential into initial vector
 */
void SolderPadBC::EBM3_Fill_Value(Vec x, Vec L)
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

  const PetscScalar current_scale = this->z_width();

  if(Genius::processor_id() == Genius::n_processors() -1)
  {
    VecSetValue(x, this->global_offset(), this->ext_circuit()->potential(), INSERT_VALUES);

    if(this->is_inter_connect_bc())
    {
      VecSetValue(L, this->global_offset(), 1.0, INSERT_VALUES);
    }
    //for stand alone electrode
    else
    {
      if(ext_circuit()->is_voltage_driven())
      {
        const PetscScalar R = this->ext_circuit()->R();
        VecSetValue(L, this->global_offset(), 1.0/((1.0+R)*current_scale), INSERT_VALUES);
      }

      if(ext_circuit()->is_current_driven())
      {
        VecSetValue(L, this->global_offset(), 1.0/(current_scale), INSERT_VALUES);
      }
    }

  }

}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for EBM L3 solver
 */
void SolderPadBC::EBM3_Function_Preprocess(Vec f, std::vector<PetscInt> &src_row,
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
 * build function and its jacobian for EBM L3 solver
 */
void SolderPadBC::EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

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
  PetscScalar current_scale = this->z_width();

  // get the workfunction and sigma
  const SimulationRegion * region = bc_regions().first; genius_assert(region);
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *>(region); genius_assert(resistance_region);
  const PetscScalar workfunction = resistance_region->material()->basic->Affinity(T_external());
  const PetscScalar sigma = resistance_region->material()->basic->Conductance();

  // the electrode potential in current iteration
  genius_assert( local_offset()!=invalid_uint );
  PetscScalar Ve = x[this->local_offset()];


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
      switch ( region->type() )
      {
          case MetalRegion :
          {
            const FVM_Node * fvm_node = ( *rnode_it ).second.second;
            const FVM_NodeData * node_data = fvm_node->node_data();

            // psi of this node
            PetscScalar V = x[fvm_node->local_offset()+0];
            PetscScalar f_psi = V + workfunction - Ve;
            VecSetValue(f, fvm_node->global_offset()+0, f_psi, ADD_VALUES);

            if(region->get_advanced_model()->enable_Tl())
            {
              // T of this node
              PetscScalar T = x[fvm_node->local_offset()+1];
              // add heat flux out of boundary to lattice temperature equatiuon
              PetscScalar f_q =this->Heat_Transfer()*(T_external()-T)*fvm_node->outside_boundary_surface_area();
              // set governing equation to function vector
              VecSetValue(f, fvm_node->global_offset()+1, f_q, ADD_VALUES);
            }

            // conductance current
            FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
            for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
            {
              const FVM_Node *nb_node = (*nb_it).second;
              const FVM_NodeData * nb_node_data = nb_node->node_data();
              // the psi of neighbor node
              PetscScalar V_nb = x[nb_node->local_offset()];
              // distance from nb node to this node
              PetscScalar distance = (*(fvm_node->root_node()) - *(nb_node->root_node())).size();
              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = fvm_node->cv_surface_area(nb_node->root_node());

              // current flow
              current_buffer.push_back( cv_boundary*sigma*(V-V_nb)/distance*current_scale );
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
            break;
          }
          default: genius_error();
      }
    }
  }

  // the extra equation of gate boundary
  // For voltage driven
  //
  //          _____                Ve
  //    -----|_____|----/\/\/\/\-------> to gate electrode (Ve, I)
  //    | +     R          L       |
  //   Vapp                     C ===
  //    | -                        |
  //    |__________________________|

  //           GND
  //
  // And for current driven
  // NOTE: It is dangerous to attach current source to MOS gate!
  //                               Ve
  //    -->-----------------------------> to gate electrode (Ve, I)
  //    |                          |
  //   Iapp                     C ===
  //    |__________________________|
  //           GND

  // Or for inter connect
  //
  //          _____                Ve
  //    -----|_____|-------------------> to gate electrode (Ve, I)
  //    |       R
  //    |
  // V_inter_connect
  //
  //

  // we first gather the electrode current
  Parallel::allgather(current_buffer);

  // for get the current, we must sum all the terms in current_buffer
  PetscScalar current = std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );


  if(Genius::is_last_processor())
  {
    // here we process the external circuit

    //for inter connect electrode
    if(this->is_inter_connect_bc())
    {
      PetscScalar V_ic = x[this->inter_connect_hub()->local_offset()];  // potential at inter connect node
      PetscScalar R = ext_circuit()->R();                               // resistance
      PetscScalar f_ext = Ve - V_ic;
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }
    // for stand alone electrode
    else
    {
      if(ext_circuit()->is_voltage_driven())
      {
        PetscScalar Vapp = ext_circuit()->Vapp();       // application voltage
        PetscScalar R = ext_circuit()->R();             // resistance
        PetscScalar C = ext_circuit()->C();             // capacitance
        PetscScalar L = ext_circuit()->L();             // inductance
        PetscScalar I = ext_circuit()->current();       // the previous step current flow into electrode
        PetscScalar Ic = ext_circuit()->cap_current();  // the previous step current flow pass though cap to ground.
        PetscScalar P  = ext_circuit()->potential();    // the previous step potential of the electrode
        PetscScalar dt = SolverSpecify::dt;
        PetscScalar f_ext = (Ve-Vapp) + (L/dt+R)*C/dt*Ve - (L/dt+R)*C/dt*P - L/dt*(I+Ic);
        VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
      }

      if(ext_circuit()->is_current_driven())
      {
        PetscScalar Iapp = ext_circuit()->Iapp();         // application current
        PetscScalar Ic   = ext_circuit()->cap_current();  // the previous step current flow pass though cap to ground.
        PetscScalar f_ext = Ic - Iapp;
        VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
      }
    }

  }

  // save the IV of current iteration
  ext_circuit()->current_itering() = current;
  ext_circuit()->potential_itering() = Ve;

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for EBM L3 solver
 */
void SolderPadBC::EBM3_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.

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

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;
      // bd node, psi = Ve
      MatSetValue(*jac, fvm_node->global_offset(), this->global_offset(), 0, ADD_VALUES);
    }
  }

  // reserve jacobian entries for the circuit equation of SolderPadBC
  {

    std::vector<PetscInt> bc_node_reserve;
    for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // skip node not belongs to this processor
      if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != MetalRegion ) continue;

        const FVM_Node * fvm_node = ( *rnode_it ).second.second;
        bc_node_reserve.push_back(fvm_node->global_offset());

        // get the derivative of electrode current to neighbors of bc node
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
        for(; nb_it != nb_it_end; ++nb_it)
        {
          const FVM_Node *  fvm_nb_node = (*nb_it).second;
          bc_node_reserve.push_back(fvm_nb_node->global_offset());
        }
      }
    }
    Parallel::allgather(bc_node_reserve);

    if(Genius::processor_id() == Genius::n_processors()-1)
    {
      PetscInt bc_global_offset = this->global_offset();

      MatSetValue(*jac, bc_global_offset, bc_global_offset, 0, ADD_VALUES);

      if(this->is_inter_connect_bc())
        MatSetValue(*jac, bc_global_offset, this->inter_connect_hub()->global_offset(), 0, ADD_VALUES);

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
 * do pre-process to jacobian matrix for EBM L3 solver
 */
void SolderPadBC::EBM3_Jacobian_Preprocess(Mat *jac, std::vector<PetscInt> &src_row,
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
 * build function and its jacobian for EBM L3 solver
 */
void SolderPadBC::EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // the Jacobian of SolderPad boundary condition is processed here

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  const PetscInt bc_global_offset = this->global_offset();
  const PetscScalar R             = this->ext_circuit()->R();  // resistance
  const PetscScalar C             = this->ext_circuit()->C();  // capacitance
  const PetscScalar L             = this->ext_circuit()->L();  // inductance
  const PetscScalar dt            = SolverSpecify::dt;

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();

  const SimulationRegion * region = bc_regions().first; genius_assert(region);
  const MetalSimulationRegion * resistance_region = dynamic_cast<const MetalSimulationRegion *>(region);
  const PetscScalar workfunction = resistance_region->material()->basic->Affinity(T_external());
  const PetscScalar sigma = resistance_region->material()->basic->Conductance();


  // we use AD again. no matter it is overkill here.
  //the indepedent variable number, we only need 2 here.
  adtl::AutoDScalar::numdir=2;

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
      switch ( region->type() )
      {
          case MetalRegion :
          {

            const FVM_Node * fvm_node = ( *rnode_it ).second.second;
            const FVM_NodeData * node_data = fvm_node->node_data();

            // psi of this node
            AutoDScalar V = x[fvm_node->local_offset()+0];  V.setADValue(0, 1.0);
            // the electrode potential in current iteration
            genius_assert( local_offset()!=invalid_uint );
            AutoDScalar Ve = x[this->local_offset()];     Ve.setADValue(1, 1.0);

            AutoDScalar f_psi = V + workfunction - Ve;
            //governing equation
            MatSetValue(*jac, fvm_node->global_offset(), fvm_node->global_offset(), f_psi.getADValue(0), ADD_VALUES);
            MatSetValue(*jac, fvm_node->global_offset(), bc_global_offset, f_psi.getADValue(1), ADD_VALUES);

            if(region->get_advanced_model()->enable_Tl())
            {
              // T of this node
              AutoDScalar T = x[fvm_node->local_offset()+1];  T.setADValue(0, 1.0);
              // add heat flux out of boundary to lattice temperature equatiuon
              AutoDScalar f_q =this->Heat_Transfer()*(T_external()-T)*fvm_node->outside_boundary_surface_area();

              //governing equation
              MatSetValue(*jac, fvm_node->global_offset()+1, fvm_node->global_offset()+1, f_q.getADValue(0), ADD_VALUES);
            }

            // conductance current
            FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
            for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
            {
              const FVM_Node *nb_node = (*nb_it).second;
              const FVM_NodeData * nb_node_data = nb_node->node_data();

              // the psi of neighbor node
              AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);

              // distance from nb node to this node
              PetscScalar distance = fvm_node->distance(nb_node);

              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = fvm_node->cv_surface_area(nb_node->root_node());


              AutoDScalar current = cv_boundary*sigma*(V-V_nb)/distance*current_scale;

              // consider electrode connect

              //for inter connect electrode
              if(this->is_inter_connect_bc())
                current = R*current;
              // for stand alone electrode
              else
              {
                if(ext_circuit()->is_voltage_driven())
                  current = (L/dt+R)*current;
                else
                  current = current;
              }

              MatSetValue(*jac, bc_global_offset, fvm_node->global_offset(), current.getADValue(0), ADD_VALUES);
              MatSetValue(*jac, bc_global_offset, nb_node->global_offset(), current.getADValue(1), ADD_VALUES);

            }

            break;
          }

          case InsulatorRegion :
          {
            const FVM_Node * fvm_node = ( *rnode_it ).second.second;

            // psi of this node
            AutoDScalar V = x[fvm_node->local_offset()];  V.setADValue(0, 1.0);

            // the electrode potential in current iteration
            genius_assert( local_offset()!=invalid_uint );
            AutoDScalar Ve = x[this->local_offset()];     Ve.setADValue(1, 1.0);

            AutoDScalar f_psi = (V + workfunction - Ve);

            //governing equation
            MatSetValue(*jac, fvm_node->global_offset(), fvm_node->global_offset(), f_psi.getADValue(0), ADD_VALUES);
            MatSetValue(*jac, fvm_node->global_offset(), bc_global_offset, f_psi.getADValue(1), ADD_VALUES);
            break;
          }
          default: genius_error();
      }
    }

  }


  // the extra equation of gate boundary
  // For voltage driven
  //
  //          _____                Ve
  //    -----|_____|----/\/\/\/\-------> to gate electrode (Ve, I)
  //    |       R          L       |
  //   Vapp                     C ===
  //    |__________________________|
  //           GND
  //
  // And for current driven
  // NOTE: It is dangerous to attach current source to MOS gate!
  //                               Ve
  //    --------------------------------> to gate electrode (Ve, I)
  //    |                          |
  //   Iapp                     C ===
  //    |__________________________|
  //           GND
  //
  // Or for inter connect
  //
  //          _____                Ve
  //    -----|_____|-------------------> to gate electrode (Ve, I)
  //    |       R
  //    |
  // V_inter_connect
  //
  //

  if(Genius::is_last_processor())
  {
    // here we process the external circuit, we do not use AD here
    // NOTE current item such as (L/dt+R)*current and current has already been processed before

    //for inter connect electrode
    if(this->is_inter_connect_bc())
    {
      // the external electrode equation is:
      // f_ext = Ve - V_ic + R*current;

      // d(f_ext)/d(Ve)
      MatSetValue(*jac, bc_global_offset, bc_global_offset, 1.0, ADD_VALUES);
      // d(f_ext)/d(V_ic)
      MatSetValue(*jac, bc_global_offset, this->inter_connect_hub()->global_offset(), -1.0, ADD_VALUES);
    }
    //for stand alone electrode
    else
    {
      if(ext_circuit()->is_voltage_driven())
      {
        // the external electrode equation is:
        // f_ext = (L/dt+R)*current + (Ve-Vapp) + (L/dt+R)*C/dt*Ve - (L/dt+R)*C/dt*P - L/dt*(I+Ic);

        // d(f_ext)/d(Ve)
        MatSetValue(*jac, bc_global_offset, bc_global_offset, 1+(L/dt+R)*C/dt, ADD_VALUES);
      }

      if(ext_circuit()->is_current_driven())
      {
        // the external electrode equation is:
        // f_ext = current + Ic - Iapp;
        // so nothing to do
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * update electrode potential
 */
void SolderPadBC::EBM3_Update_Solution(PetscScalar *)
{
  this->ext_circuit()->update();
}
