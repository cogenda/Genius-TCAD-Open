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
#include "spice_ckt.h"

using PhysicalUnit::V;
using PhysicalUnit::A;


void ElectrodeInterConnectBC::mix_inter_connect_fill_value(Vec x, Vec L)
{
  if(Genius::is_last_processor())
  {
    VecSetValue(L, this->global_offset(), 1.0, INSERT_VALUES);
  }
}



void ElectrodeInterConnectBC::mix_inter_connect_function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
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

  // only last processor do this
  if(Genius::is_last_processor())
  {
    // potential of inter-connect node
    PetscScalar Vic = x[ckt->local_offset_x(spice_node_index)];

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

    VecSetValue(f, this->global_offset(), I/A, ADD_VALUES);
  }

}



void ElectrodeInterConnectBC::mix_inter_connect_reserve(Mat *jac, InsertMode &add_value_flag)
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



void ElectrodeInterConnectBC::mix_inter_connect_jacobian(PetscScalar *x , Mat *jac, InsertMode &add_value_flag)
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
      PetscScalar r  = electrode->ext_circuit()->inter_connect_resistance();
      //I += (Vic-Ve)/r/A;
      MatSetValue(*jac, this->global_offset(), this->global_offset(), 1.0/r/A, ADD_VALUES);
      MatSetValue(*jac, this->global_offset(), electrode->global_offset(), -1.0/r/A, ADD_VALUES);
    }

    //PetscScalar Is = (Vic - ext_circuit()->potential())*SolverSpecify::Gmin/(V/A);
    //MatSetValue(*jac, this->global_offset(), this->global_offset(), SolverSpecify::Gmin/(V/A), ADD_VALUES);
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


