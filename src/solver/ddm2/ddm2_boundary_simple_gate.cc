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

  PetscScalar q = e*this->Qf();            // surface change density
  PetscScalar Thick = this->Thickness();   // the thickness of gate oxide
  PetscScalar eps_ox = this->eps();        // the permittivity of gate material

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
    PetscScalar dP = S*(eps_ox*(Ve - this->Work_Function()-V)/Thick + q);
    // set governing equation of psi
    VecSetValue(f, fvm_node->global_offset()+0, dP, ADD_VALUES);

    // add heat flux out of gate boundary to lattice temperature equatiuon
    // when this bc is external boundary
    if( node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion))
    {
      PetscScalar T = x[fvm_node->local_offset()+3]; // T of this node
      PetscScalar h = this->Heat_Transfer();
      PetscScalar S  = fvm_node->outside_boundary_surface_area();
      PetscScalar fT = h*(T_external()-T)*S;
      VecSetValue(f, fvm_node->global_offset()+3, fT, ADD_VALUES);
    }

    // displacement current, only first order in time.
    PetscScalar dEdt = ( (Ve-this->ext_circuit()->potential_old()) - (V-node_data->psi()) ) /Thick;
    current_buffer.push_back( S*eps_ox*dEdt );

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

  //for inter connect electrode
  if(this->is_inter_connect_bc())
  {
    PetscScalar R = ext_circuit()->R();                               // resistance
    PetscScalar f_ext = R*current;
    VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
  }
  // for stand alone electrode
  else
  {
    if(ext_circuit()->is_voltage_driven())
    {
      PetscScalar R = ext_circuit()->R();             // resistance
      PetscScalar L = ext_circuit()->L();             // inductance
      PetscScalar dt = SolverSpecify::dt;
      PetscScalar f_ext = (L/dt+R)*current;
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }

    if(ext_circuit()->is_current_driven())
    {
      PetscScalar f_ext = current;
      VecSetValue(f, this->global_offset(), f_ext, ADD_VALUES);
    }
  }


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
 * reserve non zero pattern in jacobian matrix for DDML2 solver
 */
void SimpleGateContactBC::DDM2_Jacobian_Reserve(Mat *jac, InsertMode &add_value_flag)
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

    const FVM_Node * fvm_node = (*region_node_begin(*node_it)).second.second;
    // bd node, psi = Ve
    MatSetValue(*jac, fvm_node->global_offset()+0, this->global_offset(), 0, ADD_VALUES);
  }

  // reserve jacobian entries for the circuit equation of simple gate electrode
  {

    std::vector<PetscInt> bc_node_reserve;
    for(node_it = nodes_begin(); node_it!=end_it; ++node_it )
    {
      // get the derivative of electrode current to ohmic node
      const FVM_Node *  fvm_node = get_region_fvm_node(*node_it, SemiconductorRegion);
      if(fvm_node->on_processor())
      {
        bc_node_reserve.push_back(fvm_node->global_offset()+0);

        // get the derivative of electrode current to neighbors of bc node
        FVM_Node::fvm_neighbor_node_iterator nb_it = fvm_node->neighbor_node_begin();
        FVM_Node::fvm_neighbor_node_iterator nb_it_end = fvm_node->neighbor_node_end();
        for(; nb_it != nb_it_end; ++nb_it)
        {
          const FVM_Node *  fvm_nb_node = (*nb_it).second;
          bc_node_reserve.push_back(fvm_nb_node->global_offset()+0);
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
 * build function and its jacobian for DDM L2 solver
 */
void SimpleGateContactBC::DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  // the Jacobian of simple gate boundary condition is processed here

  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  genius_assert( this->global_offset()!=invalid_uint );
  PetscInt bc_global_offset = this->global_offset();
  PetscScalar q             = e*this->Qf();              // surface change density
  PetscScalar Thick         = this->Thickness();         // the thickness of gate oxide
  PetscScalar eps_ox        = this->eps();               // the permittivity of gate material
  PetscScalar R             = this->ext_circuit()->R();  // resistance
  PetscScalar C             = this->ext_circuit()->C();  // capacitance
  PetscScalar L             = this->ext_circuit()->L();  // inductance
  PetscScalar dt            = SolverSpecify::dt;

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

    AutoDScalar dP = S*(eps_ox*(Ve - this->Work_Function()-V)/Thick + q);

    //governing equation of psi
    MatSetValue(*jac, fvm_node->global_offset(), fvm_node->global_offset(), dP.getADValue(0), ADD_VALUES);
    MatSetValue(*jac, fvm_node->global_offset(), bc_global_offset, dP.getADValue(1), ADD_VALUES);

    // process the Jacobian of equation of T
    // if this gate bc is external boundary, set heat flux here
    if( node_on_boundary(*node_it) || has_associated_region(*node_it, VacuumRegion))
    {
      AutoDScalar T = x[fvm_node->local_offset()+3]; T.setADValue(0, 1.0); // psi of this node
      PetscScalar h = this->Heat_Transfer();
      PetscScalar S  = fvm_node->outside_boundary_surface_area();
      AutoDScalar fT = h*(T_external()-T)*S;
      MatSetValue(*jac, fvm_node->global_offset()+3, fvm_node->global_offset()+3, fT.getADValue(0), ADD_VALUES);
    }

    // displacement current, only first order in time.
    {
      AutoDScalar dEdt = ( (Ve-this->ext_circuit()->potential_old()) - (V-node_data->psi()) ) /Thick;
      AutoDScalar current_disp = S*eps_ox*dEdt*current_scale;
      //for inter connect electrode
      if(this->is_inter_connect_bc())
        current_disp = R*current_disp;
      // for stand alone electrode
      else
      {
        if(ext_circuit()->is_voltage_driven())
          current_disp = (L/dt+R)*current_disp;
        else
          current_disp = current_disp;
      }
      MatSetValue(*jac, bc_global_offset, fvm_node->global_offset(), current_disp.getADValue(0), ADD_VALUES);
      MatSetValue(*jac, bc_global_offset, bc_global_offset, current_disp.getADValue(1), ADD_VALUES);
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
void SimpleGateContactBC::DDM2_Update_Solution(PetscScalar *)
{
  Parallel::sum(ext_circuit()->current_itering());
  this->ext_circuit()->update();
}
