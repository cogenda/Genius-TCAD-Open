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


// Local includes
#include "simulation_system.h"
#include "boundary_condition_electrode_interconnect.h"


using PhysicalUnit::V;
using PhysicalUnit::A;

/*---------------------------------------------------------------------
 * fill initial value to inter-connect node
 */
void ElectrodeInterConnectBC::inter_connect_fill_value(Vec x, Vec L)
{
  if(Genius::is_last_processor())
  {
    VecSetValue(x, this->global_offset(), this->ext_circuit()->potential(), INSERT_VALUES);

    const std::vector<BoundaryCondition * > & electrodes = this->inter_connect();
    PetscScalar rmin = 1e30;
    for(unsigned int n=0; n<electrodes.size(); ++n)
    {
      const BoundaryCondition * electrode = electrodes[n];
      PetscScalar r  = electrode->ext_circuit()->R();
      rmin = std::min(r, rmin);
    }

    if(this->ext_circuit()->is_voltage_driven())
    {
      PetscScalar R = this->ext_circuit()->R();
      rmin = std::min(R, rmin);
    }

    VecSetValue(L, this->global_offset(), rmin, INSERT_VALUES);
  }

}



// The inter connect node can be voltage driven, current driven and float
//
// voltage driven
//          _____                Vic
//    -----|_____|----/\/\/\/\-------> inter_connect
//    | +     R          L       |
//   Vapp                     C ===
//    | -                        |
//    |__________________________|
//           GND
//
// current driven
//                               Vic
//    -->-----------------------------> inter_connect
//    |                          |
//   Iapp                     C ===
//    |                          |
//    |__________________________|
//           GND
//
// float
//
//                              Vic
//    -----| |-----------------------> inter_connect
//    |     C
//   GND
//
//


/*---------------------------------------------------------------------
 * set governing equation for inter-connect node. use nodal analysis method
 */
void ElectrodeInterConnectBC::inter_connect_function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
{
  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }


  // potential of inter-connect node
  PetscScalar Vic = x[this->local_offset()];
  // total current flow into each electrodes
  PetscScalar I = 0;
  const std::vector<BoundaryCondition * > & electrodes = this->inter_connect();
  for(unsigned int n=0; n<electrodes.size(); ++n)
  {
    const BoundaryCondition * electrode = electrodes[n];
    PetscScalar Ve = x[electrode->local_offset()];
    PetscScalar r  = electrode->ext_circuit()->R();
    I += (Vic-Ve)/r;
  }

    // only last processor do this
  if(Genius::is_last_processor())
  {

    // the inter connect node has external voltage source
    if(this->ext_circuit()->is_voltage_driven())
    {
      PetscScalar Vapp = this->ext_circuit()->Vapp();
      PetscScalar R = this->ext_circuit()->R();
      PetscScalar C = this->ext_circuit()->C();
      PetscScalar L = this->ext_circuit()->L();
      PetscScalar dt = SolverSpecify::dt;

      // current pass through cap to ground
      PetscScalar Ic = C/dt*(Vic - ext_circuit()->potential());
      // current in R-L branch
      PetscScalar Irl = (Vic-Vapp)/(R+L/dt) - L/dt*(ext_circuit()->current() + ext_circuit()->cap_current());
      // total current should be zero
      PetscScalar ff = I + Ic + Irl;
      VecSetValue(f, this->global_offset(), ff, ADD_VALUES);
    }

    // the inter connect node is driven by external current source
    if(this->ext_circuit()->is_current_driven())
    {
      // external current source, here we should take care of the current direction
      PetscScalar Iapp = -this->ext_circuit()->Iapp();
      PetscScalar C = this->ext_circuit()->C();
      PetscScalar dt = SolverSpecify::dt;

      // current flow into cap
      PetscScalar Ic = C/dt*(Vic-ext_circuit()->potential());
      // total current should be zero
      PetscScalar ff = I + Ic + Iapp;
      VecSetValue(f, this->global_offset(), ff, ADD_VALUES);
    }

    // the inter connect node is float
    if(this->ext_circuit()->is_float())
    {
      PetscScalar C = this->ext_circuit()->C();
      PetscScalar dt = SolverSpecify::dt;

      // current flow into cap
      PetscScalar Ic = C/dt*(Vic-ext_circuit()->potential());

      // total current should be zero
      PetscScalar ff = I + Ic;
      VecSetValue(f, this->global_offset(), ff, ADD_VALUES);
    }

    PetscScalar Is = (Vic - ext_circuit()->potential())*SolverSpecify::Gmin/(V/A);
    VecSetValue(f, this->global_offset(), Is, ADD_VALUES);
  }

  // save the IV of current iteration
  ext_circuit()->current_itering() = I;
  ext_circuit()->potential_itering() = Vic;

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve matrix entries
 */
void ElectrodeInterConnectBC::inter_connect_reserve(Mat *jac, InsertMode &add_value_flag)
{
  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  // only past processor do this
  if(Genius::is_last_processor())
  {
    // current flow into each electrode
    const std::vector<BoundaryCondition * > & electrodes = this->inter_connect();
    for(unsigned int n=0; n<electrodes.size(); ++n)
    {
      const BoundaryCondition * electrode = electrodes[n];
      MatSetValue(*jac, this->global_offset(), this->global_offset(), 0, ADD_VALUES);
      MatSetValue(*jac, this->global_offset(), electrode->global_offset(), 0, ADD_VALUES);
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * set jacobian matrix entries
 */
void ElectrodeInterConnectBC::inter_connect_jacobian(PetscScalar * , Mat *jac, InsertMode &add_value_flag)
{
  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  // only last processor do this
  if(Genius::is_last_processor())
  {
    // current flow into each electrode
    const std::vector<BoundaryCondition * > & electrodes = this->inter_connect();
    for(unsigned int n=0; n<electrodes.size(); ++n)
    {
      const BoundaryCondition * electrode = electrodes[n];
      PetscScalar r  = electrode->ext_circuit()->R();
      //I += (Vic-Ve)/r;
      MatSetValue(*jac, this->global_offset(), this->global_offset(), 1.0/r, ADD_VALUES);
      MatSetValue(*jac, this->global_offset(), electrode->global_offset(), -1.0/r, ADD_VALUES);
    }

    // the inter connect node has external voltage source
    if(this->ext_circuit()->is_voltage_driven())
    {
      PetscScalar R = this->ext_circuit()->R();
      PetscScalar C = this->ext_circuit()->C();
      PetscScalar L = this->ext_circuit()->L();
      PetscScalar dt = SolverSpecify::dt;

      // current pass through cap to ground
      // Ic = C/dt*(Vic - ext_circuit()->potential());
      // current in R-L branch
      // Irl = (Vic-Vapp)/(R+L/dt) - L/dt*(ext_circuit()->current() + ext_circuit()->cap_current());
      // total current should be zero
      // ff = I + Ic + Irl;

      MatSetValue(*jac, this->global_offset(), this->global_offset(), C/dt + 1.0/(R+L/dt), ADD_VALUES);
    }

    // the inter connect node is driven by external current source
    if(this->ext_circuit()->is_current_driven())
    {
      PetscScalar C = this->ext_circuit()->C();
      PetscScalar dt = SolverSpecify::dt;

      // current flow into cap
      // Ic = C/dt*(Vic-ext_circuit()->potential());
      // total current should be zero
      // ff = I + Ic + Iapp;
      MatSetValue(*jac, this->global_offset(), this->global_offset(), C/dt, ADD_VALUES);
    }

    // the inter connect node is float
    if(this->ext_circuit()->is_float())
    {
      PetscScalar C = this->ext_circuit()->C();
      PetscScalar dt = SolverSpecify::dt;

      // current flow into cap
      // Ic = C/dt*(Vic-ext_circuit()->potential());
      // total current should be zero
      // ff = I + Ic;
      MatSetValue(*jac, this->global_offset(), this->global_offset(), C/dt, ADD_VALUES);
    }

    //PetscScalar Is = (Vic - ext_circuit()->potential())*SolverSpecify::Gmin/(V/A);
    MatSetValue(*jac, this->global_offset(), this->global_offset(), SolverSpecify::Gmin/(V/A), ADD_VALUES);
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


}



/*---------------------------------------------------------------------
 * update solution data
 */
void ElectrodeInterConnectBC::inter_connect_update_solution(PetscScalar *)
{
  this->ext_circuit()->update();
}

