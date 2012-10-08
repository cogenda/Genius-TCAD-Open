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
using PhysicalUnit::Ohm;

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
      PetscScalar r  = electrode->ext_circuit()->inter_connect_resistance();
      rmin = std::min(r, rmin);
    }

    if(this->ext_circuit()->is_voltage_driven())
    {
      PetscScalar R = this->ext_circuit()->serial_resistance();
      rmin = std::min(std::max(R, 1e-3*Ohm), rmin);
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

  if(Genius::is_last_processor())
  {

    const PetscScalar mna_scaling = this->ext_circuit()->mna_scaling(SolverSpecify::dt);

    // potential of inter-connect node
    PetscScalar Vic = x[this->local_offset()];
    ext_circuit()->potential() = Vic;

    // total current flow into each electrodes
    PetscScalar I = 0;
    const std::vector<BoundaryCondition * > & electrodes = this->inter_connect();
    for(unsigned int n=0; n<electrodes.size(); ++n)
    {
      const BoundaryCondition * electrode = electrodes[n];
      PetscScalar Ve = x[electrode->local_offset()];
      PetscScalar r  = electrode->ext_circuit()->inter_connect_resistance();
      I += mna_scaling*(Vic-Ve)/r;
    }

    // total current should be zero
    PetscScalar ff = I + this->ext_circuit()->mna_function(SolverSpecify::dt);
    VecSetValue(f, this->global_offset(), ff, ADD_VALUES);

    PetscScalar Is = mna_scaling*(Vic - ext_circuit()->potential())*SolverSpecify::Gmin/(V/A);
    VecSetValue(f, this->global_offset(), Is, ADD_VALUES);

  }

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
void ElectrodeInterConnectBC::inter_connect_jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
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
    const PetscScalar mna_scaling = this->ext_circuit()->mna_scaling(SolverSpecify::dt);

    // potential of inter-connect node
    PetscScalar Vic = x[this->local_offset()];
    ext_circuit()->potential() = Vic;

    // current flow into each electrode
    const std::vector<BoundaryCondition * > & electrodes = this->inter_connect();
    for(unsigned int n=0; n<electrodes.size(); ++n)
    {
      const BoundaryCondition * electrode = electrodes[n];
      PetscScalar r  = electrode->ext_circuit()->inter_connect_resistance();
      //I += mna_scaling*(Vic-Ve)/r;
      MatSetValue(*jac, this->global_offset(), this->global_offset(), mna_scaling/r, ADD_VALUES);
      MatSetValue(*jac, this->global_offset(), electrode->global_offset(), -mna_scaling/r, ADD_VALUES);
    }

    MatSetValue(*jac, this->global_offset(), this->global_offset(), ext_circuit()->mna_jacobian(SolverSpecify::dt), ADD_VALUES);

    //PetscScalar Is = mna_scaling*(Vic - ext_circuit()->potential())*SolverSpecify::Gmin/(V/A);
    MatSetValue(*jac, this->global_offset(), this->global_offset(), mna_scaling*SolverSpecify::Gmin/(V/A), ADD_VALUES);
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * update solution data
 */
void ElectrodeInterConnectBC::inter_connect_update_solution(PetscScalar *x)
{
  // potential of inter-connect node
  PetscScalar Vic = x[this->local_offset()];
  // total current flow into each electrodes
  PetscScalar I = 0;
  const std::vector<BoundaryCondition * > & electrodes = this->inter_connect();
  for(unsigned int n=0; n<electrodes.size(); ++n)
  {
    const BoundaryCondition * electrode = electrodes[n];
    PetscScalar Ve = x[electrode->local_offset()];
    PetscScalar r  = electrode->ext_circuit()->inter_connect_resistance();
    I += (Vic-Ve)/r;
  }

  // save the IV of current iteration
  ext_circuit()->current() = I;
  ext_circuit()->potential() = Vic;

  this->ext_circuit()->update();
}

