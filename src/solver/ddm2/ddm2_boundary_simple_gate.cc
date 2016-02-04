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
#include "semiconductor_region.h"
#include "boundary_condition_simplegate.h"
#include "parallel.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::Ohm;


/*---------------------------------------------------------------------
 * fill electrode potential into initial vector
 */
void SimpleGateContactBC::DDM2_Fill_Value(Vec x, Vec L)
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
 * build function and its jacobian for DDM L2 solver
 */
void SimpleGateContactBC::DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // Simple gate boundary condition is processed here

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

  // the electrode potential in current iteration
  genius_assert( this->local_offset()!=invalid_uint );
  PetscScalar Ve = x[this->local_offset()];

  const PetscScalar q = e*this->scalar("qf");            // surface change density
  const PetscScalar Thick = this->scalar("thickness");   // the thickness of gate oxide
  const PetscScalar eps_ox = this->scalar("eps");        // the permittivity of gate material
  const PetscScalar Work_Function = this->scalar("workfunction");
  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * fvm_node = (*region_node_begin(*node_it)).second.second;
    const FVM_NodeData * node_data = fvm_node->node_data();


    PetscScalar V = x[fvm_node->local_offset()+0]; // psi of this node
    PetscScalar S = fvm_node->outside_boundary_surface_area();
    PetscScalar dP = S*(eps_ox*(Ve - Work_Function-V)/Thick + q);
    // set governing equation of psi
    VecSetValue(f, fvm_node->global_offset()+0, dP, ADD_VALUES);

    // add heat flux out of gate boundary to lattice temperature equatiuon
    // when this bc is external boundary
    if( node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion))
    {
      PetscScalar T = x[fvm_node->local_offset()+3]; // T of this node
      PetscScalar S  = fvm_node->outside_boundary_surface_area();
      PetscScalar fT = Heat_Transfer*(T_external()-T)*S;
      VecSetValue(f, fvm_node->global_offset()+3, fT, ADD_VALUES);
    }

    // displacement current, only first order in time.
    if(SolverSpecify::TimeDependent == true)
    {
      PetscScalar dEdt = ( (Ve-this->ext_circuit()->potential_old()) - (V-node_data->psi()) ) /Thick/SolverSpecify::dt;
      current_buffer.push_back( S*eps_ox*dEdt );
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
  PetscScalar current = current_scale*std::accumulate(current_buffer.begin(), current_buffer.end(), 0.0 );

  ext_circuit()->potential() = Ve;
  ext_circuit()->current() = current;

  PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);

  //for inter connect electrode
  if(this->is_inter_connect_bc())
  {
    PetscScalar R = ext_circuit()->inter_connect_resistance();                               // resistance
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
 * build function and its jacobian for DDM L2 solver
 */
void SimpleGateContactBC::DDM2_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{
  // the Jacobian of simple gate boundary condition is processed here

  genius_assert( this->global_offset()!=invalid_uint );

  PetscInt bc_global_offset = this->global_offset();

  const PetscScalar q = e*this->scalar("qf");            // surface change density
  const PetscScalar Thick = this->scalar("thickness");   // the thickness of gate oxide
  const PetscScalar eps_ox = this->scalar("eps");        // the permittivity of gate material
  const PetscScalar Work_Function = this->scalar("workfunction");
  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");


  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();

  // we use AD again. no matter it is overkill here.
  //the indepedent variable number, we only need 2 here.
  adtl::AutoDScalar::numdir=2;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    const FVM_Node * fvm_node = (*region_node_begin(*node_it)).second.second;
    const FVM_NodeData * node_data = fvm_node->node_data();

    // psi of this node
    AutoDScalar V = x[fvm_node->local_offset()];  V.setADValue(0, 1.0);

    // the electrode potential in current iteration
    genius_assert( this->local_offset()!=invalid_uint );
    AutoDScalar Ve = x[this->local_offset()];     Ve.setADValue(1, 1.0);

    PetscScalar S = fvm_node->outside_boundary_surface_area();

    AutoDScalar dP = S*(eps_ox*(Ve - Work_Function-V)/Thick + q);

    //governing equation of psi
    jac->add( fvm_node->global_offset(),  fvm_node->global_offset(),  dP.getADValue(0) );
    jac->add( fvm_node->global_offset(),  bc_global_offset,  dP.getADValue(1) );

    // process the Jacobian of equation of T
    // if this gate bc is external boundary, set heat flux here
    if( node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion))
    {
      AutoDScalar T = x[fvm_node->local_offset()+3]; T.setADValue(0, 1.0); // psi of this node
      PetscScalar S  = fvm_node->outside_boundary_surface_area();
      AutoDScalar fT = Heat_Transfer*(T_external()-T)*S;
      jac->add( fvm_node->global_offset()+3,  fvm_node->global_offset()+3,  fT.getADValue(0) );
    }

    // displacement current, only first order in time.
    if(SolverSpecify::TimeDependent == true)
    {
      AutoDScalar dEdt = ( (Ve-this->ext_circuit()->potential_old()) - (V-node_data->psi()) ) /Thick/SolverSpecify::dt;
      AutoDScalar current_disp = S*eps_ox*dEdt*current_scale;
      //for inter connect electrode
      if(this->is_inter_connect_bc())
      {
        PetscScalar R = ext_circuit()->inter_connect_resistance();
        current_disp = R*current_disp;
      }
      // for stand alone electrode
      else
      {
        PetscScalar mna_scaling = ext_circuit()->mna_scaling(SolverSpecify::dt);
        current_disp = mna_scaling*current_disp;
      }
      jac->add( bc_global_offset,  fvm_node->global_offset(),  current_disp.getADValue(0) );
      jac->add( bc_global_offset,  bc_global_offset,  current_disp.getADValue(1) );
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
      jac->add( bc_global_offset,  bc_global_offset,  1.0 );
      // d(f_ext)/d(V_ic)
      jac->add( bc_global_offset,  this->inter_connect_hub()->global_offset(),  -1.0 );
    }
    //for stand alone electrode
    else
    {
      ext_circuit()->potential() = x[this->local_offset()];
      jac->add( bc_global_offset,  bc_global_offset,  ext_circuit()->mna_jacobian(SolverSpecify::dt) );
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * update electrode potential
 */
void SimpleGateContactBC::DDM2_Update_Solution(PetscScalar *)
{
  Parallel::sum(ext_circuit()->current());
  this->ext_circuit()->update();
}

