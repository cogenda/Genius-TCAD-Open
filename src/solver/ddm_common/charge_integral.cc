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
#include "boundary_condition_charge_integral.h"


/*---------------------------------------------------------------------
 * fill initial value to inter-connect node
 */
void ChargeIntegralBC::charge_integral_fill_value(Vec x, Vec L)
{
  if(Genius::is_last_processor())
  {
    VecSetValue(x, this->global_offset(), 0.0, INSERT_VALUES);
    VecSetValue(L, this->global_offset(), 1.0, INSERT_VALUES);
  }
}



/*---------------------------------------------------------------------
 * set governing equation for inter-connect node. use nodal analysis method
 */
void ChargeIntegralBC::charge_integral_function(PetscScalar *x , Vec f, InsertMode &add_value_flag)
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
    // add to ChargeIntegralBC
    VecSetValue(f, this->global_offset(), this->scalar("qf"), ADD_VALUES);
  }
}



/*---------------------------------------------------------------------
 * set jacobian matrix entries
 */
void ChargeIntegralBC::charge_integral_jacobian(PetscScalar * , SparseMatrix<PetscScalar> *jac, InsertMode &)
{
  //if(Genius::is_last_processor())
  //  MatSetValue(*jac, this->global_offset(), this->global_offset(), 0.0, ADD_VALUES);
}


/*---------------------------------------------------------------------
 * update solution data
 */
void ChargeIntegralBC::charge_integral_update_solution(PetscScalar *lxx)
{
  this->psi() = lxx[this->local_offset()];
}

