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



#include "simulation_system.h"
#include "insulator_region.h"
#include "solver_specify.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;




///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////



void InsulatorSimulationRegion::HALL_Fill_Value(Vec x, Vec L)
{
  this->DDM1_Fill_Value(x, L);
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDM1 solver
 */
void InsulatorSimulationRegion::HALL_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  this->DDM1_Function(x, f, add_value_flag);
}


/*---------------------------------------------------------------------
 * build function and its jacobian for DDM1 solver
 */
void InsulatorSimulationRegion::HALL_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  this->DDM1_Jacobian(x, jac, add_value_flag);
}


void InsulatorSimulationRegion::HALL_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{
  this->DDM1_Time_Dependent_Function(x, f, add_value_flag);
}


void InsulatorSimulationRegion::HALL_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{
  this->DDM1_Time_Dependent_Jacobian(x, jac, add_value_flag);
}


void InsulatorSimulationRegion::HALL_Function_Hanging_Node(PetscScalar *x, Vec f, InsertMode &add_value_flag)
{
  this->DDM1_Function_Hanging_Node(x, f, add_value_flag);
}


void InsulatorSimulationRegion::HALL_Jacobian_Hanging_Node(PetscScalar *x, Mat *jac, InsertMode &add_value_flag)
{
  this->DDM1_Jacobian_Hanging_Node(x, jac, add_value_flag);
}


void InsulatorSimulationRegion::HALL_Update_Solution(PetscScalar *lxx)
{
  this->DDM1_Update_Solution(lxx);
}


