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

//  $Id: fvm_data.cc,v 1.2 2008/07/09 05:58:16 gdiso Exp $


#include "fvm_node_data_semiconductor.h"
#include "fvm_node_data_insulator.h"
#include "fvm_node_data_conductor.h"
#include "fvm_node_data_resistance.h"
#include "physical_unit.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;


// the static member in simulation region
PetscScalar FVM_NodeData::_scalar_dummy_ = 0.0;

std::complex<PetscScalar> FVM_NodeData::_complex_dummy_ = std::complex<PetscScalar>(0.0, 0.0);

VectorValue<PetscScalar> FVM_NodeData::_vector_dummy_ = VectorValue<PetscScalar>();

TensorValue<PetscScalar> FVM_NodeData::_tensor_dummy_ = TensorValue<PetscScalar>();

/*----------------------------------------------------------------
 * @return the intrinsic carrier concentration.
 * @note will not consider bandgap narrowing
 */
PetscScalar  FVM_Semiconductor_NodeData::ni()           const
{ return sqrt(Nc()*Nv())*exp(-Eg()/(2*kb*T())); }


//-----------

PetscScalar  FVM_Insulator_NodeData::Ec()           const
{ return -(e*psi() + affinity()) ; }


PetscScalar  FVM_Insulator_NodeData::Ev()           const
{ return -(e*psi() + affinity() + Eg() ); }


PetscScalar  FVM_Insulator_NodeData::qFn()           const
{ return 0.5*(Ec() + Ev()); }


PetscScalar  FVM_Insulator_NodeData::qFp()           const
{ return 0.5*(Ec() + Ev()); }

//----------

PetscScalar  FVM_Conductor_NodeData::Ec()           const
{ return -(e*psi() + affinity()) ; }


PetscScalar  FVM_Conductor_NodeData::Ev()           const
{ return -(e*psi() + affinity()); }


PetscScalar  FVM_Conductor_NodeData::qFn()           const
{ return Ec(); }


PetscScalar  FVM_Conductor_NodeData::qFp()           const
{ return Ev(); }

//----------

PetscScalar  FVM_Resistance_NodeData::Ec()           const
{ return -(e*psi() + affinity()) ; }


PetscScalar  FVM_Resistance_NodeData::Ev()           const
{ return -(e*psi() + affinity()); }


PetscScalar  FVM_Resistance_NodeData::qFn()           const
{ return Ec(); }


PetscScalar  FVM_Resistance_NodeData::qFp()           const
{ return Ev(); }

