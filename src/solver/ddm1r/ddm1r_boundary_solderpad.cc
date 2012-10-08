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
#include "petsc_utils.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::Ohm;

/*---------------------------------------------------------------------
 * fill electrode potential into initial vector
 */
void SolderPadBC::DDM1R_Fill_Value(Vec x, Vec L)
{
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
      const PetscScalar s = ext_circuit()->electrode_scaling(SolverSpecify::dt);
      VecSetValue(L, this->global_offset(), s, INSERT_VALUES);
    }

  }

}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for DDML1 solver
 */
void SolderPadBC::DDM1R_Function_Preprocess(PetscScalar *, Vec f, std::vector<PetscInt> &src_row,
                                              std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  this->_current_buffer.clear();

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

      if ( region->type() == MetalRegion)
      {
        PetscInt    ix = fvm_node->global_offset()+1;
        PetscScalar I;
        VecGetValues(f, 1, &ix, &I);
        this->_current_buffer.push_back(I);

        clear_row.push_back(fvm_node->global_offset()+0);
        clear_row.push_back(fvm_node->global_offset()+1);
      }

      if ( region->type() == InsulatorRegion)
      {
        clear_row.push_back(fvm_node->global_offset()+0);
      }
    }
  }
}



/*---------------------------------------------------------------------
 * build function and its jacobian for DDM L1 solver
 */
void SolderPadBC::DDM1R_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();
  std::vector<PetscScalar>   displacement_current_buffer;

  // the electrode potential in current iteration
  genius_assert( local_offset()!=invalid_uint );
  PetscScalar Ve = x[this->local_offset()];

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const MetalSimulationRegion * resistance_region = 0;
  if( _r1 && _r1->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r1);
  if( _r2 && _r2->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r2);
  genius_assert(resistance_region);

  const double workfunction = resistance_region->material()->basic->Affinity(T_external());
  const double sigma = resistance_region->material()->basic->Conductance();
  const PetscScalar ion   = resistance_region->material()->basic->IonDensity(T_external()); //ion density

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
            PetscScalar V = x[fvm_node->local_offset()];
            PetscScalar f_psi = V + node_data->affinity()/e - Ve;
            VecSetValue(f, fvm_node->global_offset(), f_psi, ADD_VALUES);

            PetscScalar n = x[fvm_node->local_offset()+1];
            PetscScalar f_n = n-ion;
            VecSetValue(f, fvm_node->global_offset()+1, f_n, ADD_VALUES);

            // displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              PetscScalar I_displacement = 0.0;
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
                PetscScalar cv_boundary = std::abs(fvm_node->cv_surface_area(nb_node));
                PetscScalar dEdt = ((V-V_nb)-(node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;
                I_displacement += cv_boundary*node_data->eps() *dEdt;
              }
              displacement_current_buffer.push_back( I_displacement );
            }
            break;
          }

          case InsulatorRegion:
          {
            // psi of this node
            PetscScalar V = x[fvm_node->local_offset()];
            PetscScalar f_psi = (V + workfunction - Ve);

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

  // for get the current, we must sum all the terms in current_buffer
  // NOTE: only statistic current flow belongs to on processor node
  const std::vector<double> & current_buffer = this->_current_buffer;
  PetscScalar current = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );
  PetscScalar current_displacement = current_scale*std::accumulate(displacement_current_buffer.begin(), displacement_current_buffer.end(), 0.0 );
  current += current_displacement;

  ext_circuit()->potential() = Ve;
  ext_circuit()->current() = current;

  PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);

  //for inter connect electrode
  if(this->is_inter_connect_bc())
  {
    PetscScalar R = std::max(1e-3*Ohm, ext_circuit()->serial_resistance());                               // resistance
    PetscScalar f_ext = R*current;
    VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
  }
  // for stand alone electrode
  else
  {
    PetscScalar f_ext = mna_scaling*current;
    VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
  }

  if(Genius::is_last_processor())
  {
    //for inter connect electrode
    if(this->is_inter_connect_bc())
    {
      PetscScalar V_ic = x[this->inter_connect_hub()->local_offset()];  // potential at inter connect node
      PetscScalar f_ext = Ve - V_ic;
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }
    // for stand alone electrode
    else
    {
      PetscScalar f_ext = ext_circuit()->mna_function(SolverSpecify::dt);
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for DDML1 solver
 */
void SolderPadBC::DDM1R_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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
      // get the derivative of electrode current to ohmic node
      // skip node not belongs to this processor
      if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

      // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
      // but belong to different regions in logic.
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
      for(; rnode_it!=end_rnode_it; ++rnode_it  )
      {
        const SimulationRegion * region = (*rnode_it).second.first;
        const FVM_Node *  fvm_node = (*rnode_it).second.second;

        if ( region->type() == MetalRegion)
        {
          bc_node_reserve.push_back(fvm_node->global_offset()+0);
          bc_node_reserve.push_back(fvm_node->global_offset()+1);

          // get the derivative of electrode current to neighbors of ohmic node

          FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
          FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
          for(; nb_it != nb_it_end; ++nb_it)
          {
            const FVM_Node *  fvm_nb_node = (*nb_it).first;
            bc_node_reserve.push_back(fvm_nb_node->global_offset()+0);
            bc_node_reserve.push_back(fvm_nb_node->global_offset()+1);
          }
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
 * do pre-process to jacobian matrix for DDML1 solver
 */
void SolderPadBC::DDM1R_Jacobian_Preprocess(PetscScalar *, Mat *jac, std::vector<PetscInt> &src_row,
                                           std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{
  _buffer_cols.clear();
  _buffer_jacobian_entries.clear();

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

      if ( region->type() == MetalRegion)
      {
        const FVM_NodeData * node_data = fvm_node->node_data();

        PetscScalar A1[2], A2[2];
        std::vector<PetscScalar> JM(2), JN(2);
        std::vector<PetscInt>    row(2);
        row[0] = fvm_node->global_offset()+0;
        row[1] = fvm_node->global_offset()+1;

        //NOTE MatGetValues only get value from local block!
        MatGetValues(*jac, 1, &row[1], 2, &row[0], &JM[0]);

        // get the derivative of electrode current to neighbors of ohmic node
        // NOTE neighbors and ohmic bc node may on different processor!
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
        for(; nb_it != nb_it_end; ++nb_it)
        {
          const FVM_Node *  fvm_nb_node = (*nb_it).first;

          std::vector<PetscInt>    col(2);
          col[0] = fvm_nb_node->global_offset()+0;
          col[1] = fvm_nb_node->global_offset()+1;

          MatGetValues(*jac, 1, &row[1], 2, &col[0], &JN[0]);

          _buffer_cols.push_back(col);
          _buffer_jacobian_entries.push_back(JN);
        }

        _buffer_cols.push_back(row);
        _buffer_jacobian_entries.push_back(JM);


        clear_row.push_back(fvm_node->global_offset()+0);
        clear_row.push_back(fvm_node->global_offset()+1);
      }

      if ( region->type() == InsulatorRegion)
      {
        clear_row.push_back(fvm_node->global_offset()+0);
      }

    }
  }
}



/*---------------------------------------------------------------------
 * build function and its jacobian for DDM L1 solver
 */
void SolderPadBC::DDM1R_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // the Jacobian of SolderPad boundary condition is processed here

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  const PetscInt bc_global_offset = this->global_offset();

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();
  // d(current)/d(independent variables of bd node and its neighbors)
  {
    PetscScalar bc_current_scale = 1.0;
    if(this->is_inter_connect_bc())
    {
      PetscScalar R = std::max(1e-3*Ohm, ext_circuit()->serial_resistance());
      bc_current_scale *= R*current_scale;
    }
    // for stand alone electrode
    else
    {
      PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);
      bc_current_scale *= mna_scaling*current_scale;
    }
    for(unsigned int n=0; n<_buffer_cols.size(); ++n)
    {
      std::vector<PetscScalar> bc_current_jacobian = _buffer_jacobian_entries[n];
      for(unsigned int k=0; k<bc_current_jacobian.size(); ++k)
        bc_current_jacobian[k] *= bc_current_scale;
      MatSetValues(*jac, 1, &bc_global_offset, _buffer_cols[n].size(), &(_buffer_cols[n])[0], &(bc_current_jacobian)[0], ADD_VALUES);
    }
  }

  // we use AD again. no matter it is overkill here.
  //the indepedent variable number, we only need 2 here.
  adtl::AutoDScalar::numdir=2;

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  const MetalSimulationRegion * resistance_region = 0;
  if( _r1 && _r1->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r1);
  if( _r2 && _r2->type() == MetalRegion ) resistance_region = dynamic_cast<const MetalSimulationRegion *>(_r2);
  genius_assert(resistance_region);

  const double workfunction = resistance_region->material()->basic->Affinity(T_external());
  const double sigma = resistance_region->material()->basic->Conductance();
  const PetscScalar ion   = resistance_region->material()->basic->IonDensity(T_external()); //ion density

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
            AutoDScalar V = x[fvm_node->local_offset()];  V.setADValue(0, 1.0);
            // the electrode potential in current iteration
            AutoDScalar Ve = x[this->local_offset()];     Ve.setADValue(1, 1.0);
            AutoDScalar f_psi = V + node_data->affinity()/e - Ve;

            //governing equation
            MatSetValue(*jac, fvm_node->global_offset(), fvm_node->global_offset(), f_psi.getADValue(0), ADD_VALUES);
            MatSetValue(*jac, fvm_node->global_offset(), bc_global_offset, f_psi.getADValue(1), ADD_VALUES);

            AutoDScalar n = x[fvm_node->local_offset()+1];  n.setADValue(0, 1.0);
            AutoDScalar f_n = n-ion;
            MatSetValue(*jac, fvm_node->global_offset()+1, fvm_node->global_offset()+1, f_n.getADValue(0), ADD_VALUES);

            // displacement current
            if(SolverSpecify::TimeDependent == true)
            {
              FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
              for(; nb_it != fvm_node->neighbor_node_end(); ++nb_it)
              {
                FVM_Node *nb_node = (*nb_it).first;
                FVM_NodeData * nb_node_data = nb_node->node_data();
                // the psi of neighbor node
                AutoDScalar V_nb = x[nb_node->local_offset()+0]; V_nb.setADValue(1, 1.0);
                // distance from nb node to this node
                PetscScalar distance = (*(fvm_node->root_node()) - *(nb_node->root_node())).size();
                // area of out surface of control volume related with neighbor node
                PetscScalar cv_boundary = std::abs(fvm_node->cv_surface_area(nb_node));

                AutoDScalar dEdt = ((V-V_nb)-(node_data->psi()-nb_node_data->psi()))/distance/SolverSpecify::dt;

                AutoDScalar current_disp = cv_boundary*node_data->eps()*dEdt*current_scale;

                //for inter connect electrode
                if(this->is_inter_connect_bc())
                {
                  PetscScalar R = std::max(1e-3*Ohm, ext_circuit()->serial_resistance());
                  current_disp=R*current_disp;
                }
                //for stand alone electrode
                if(ext_circuit()->is_voltage_driven())
                {
                  PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);
                  current_disp=mna_scaling*current_disp;
                }
                MatSetValue(*jac, bc_global_offset, fvm_node->global_offset()+0, current_disp.getADValue(0), ADD_VALUES);
                MatSetValue(*jac, bc_global_offset, nb_node->global_offset()+0, current_disp.getADValue(1), ADD_VALUES);
              }
            }

            break;
          }

          case InsulatorRegion :
          {
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
      ext_circuit()->potential() = x[this->local_offset()];
      MatSetValue(*jac, bc_global_offset, bc_global_offset, ext_circuit()->mna_jacobian(SolverSpecify::dt), ADD_VALUES);
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


void SolderPadBC::DDM1R_Electrode_Trace(Vec lx, Mat *jac, Vec pdI_pdx, Vec pdF_pdV)
{
  VecZeroEntries(pdI_pdx);
  VecZeroEntries(pdF_pdV);

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();
  for(unsigned int n=0; n<_buffer_cols.size(); ++n)
  {
    std::vector<PetscScalar> bc_current_jacobian = _buffer_jacobian_entries[n];
    for(unsigned int k=0; k<bc_current_jacobian.size(); ++k)
      bc_current_jacobian[k] *= current_scale;
    VecSetValues(pdI_pdx, _buffer_cols[n].size(), &(_buffer_cols[n])[0], &(bc_current_jacobian)[0], ADD_VALUES);
  }


  BoundaryCondition::const_node_iterator node_it;
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
      if( region->type() !=  MetalRegion) continue;

      const FVM_Node * fvm_node = ( *rnode_it ).second.second;
      VecSetValue( pdF_pdV, fvm_node->global_offset(), 1.0, ADD_VALUES);
    }
  }

  VecAssemblyBegin(pdI_pdx);
  VecAssemblyBegin(pdF_pdV);

  VecAssemblyEnd(pdI_pdx);
  VecAssemblyEnd(pdF_pdV);


  //delete electrode current equation, omit the effect of external resistance
  PetscInt bc_global_offset = this->global_offset();
  PetscUtils::MatZeroRows(*jac, 1, &bc_global_offset, 1.0);
}


/*---------------------------------------------------------------------
 * update electrode IV
 */
void SolderPadBC::DDM1R_Update_Solution(PetscScalar *)
{
  Parallel::sum(ext_circuit()->current());
  this->ext_circuit()->update();
}
