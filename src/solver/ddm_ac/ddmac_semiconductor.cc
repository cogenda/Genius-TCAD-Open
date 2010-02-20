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
#include "semiconductor_region.h"
#include "solver_specify.h"
#include "boundary_info.h"
#include "boundary_condition_collector.h"


using PhysicalUnit::kb;


void SemiconductorSimulationRegion::DDMAC_Fill_Value(Vec x, Vec L) const
{
  // find the node variable offset
  unsigned int n_variables     = this->ebm_n_variables();
  unsigned int node_psi_offset = this->ebm_variable_offset(POTENTIAL);
  unsigned int node_n_offset   = this->ebm_variable_offset(ELECTRON);
  unsigned int node_p_offset   = this->ebm_variable_offset(HOLE);
  unsigned int node_Tl_offset  = this->ebm_variable_offset(TEMPERATURE);
  unsigned int node_Tn_offset  = this->ebm_variable_offset(E_TEMP);
  unsigned int node_Tp_offset  = this->ebm_variable_offset(H_TEMP);

  // data buffer
  std::vector<int> ix;
  std::vector<PetscScalar> y;
  std::vector<PetscScalar> s;

  // reserve menory for data buffer
  PetscInt n_local_dofs;
  VecGetLocalSize(x, &n_local_dofs);
  ix.reserve(n_local_dofs);
  y.reserve(n_local_dofs);
  s.reserve(n_local_dofs);

  // for all the on processor node, insert value to petsc vector
  const_node_iterator it = nodes_begin();
  const_node_iterator it_end = nodes_end();
  for(; it!=it_end; ++it)
  {

    const FVM_Node * node = (*it).second;
    //if this node NOT belongs to this processor, continue
    if( node->root_node()->processor_id() != Genius::processor_id() ) continue;

    const FVM_NodeData * node_data = node->node_data();

    // the first variable, psi
    genius_assert(node_psi_offset!=invalid_uint);
    ix.push_back(node->global_offset() + node_psi_offset);
    y.push_back(node_data->psi());
    s.push_back(1.0/node->volume());

    ix.push_back(node->global_offset() + n_variables + node_psi_offset);
    y.push_back(node_data->psi());
    s.push_back(1.0/node->volume());

    // the second variable, n
    genius_assert(node_n_offset!=invalid_uint);
    ix.push_back(node->global_offset() + node_n_offset);
    y.push_back(node_data->n());
    s.push_back(1.0/node->volume());

    ix.push_back(node->global_offset() + n_variables + node_n_offset);
    y.push_back(node_data->n());
    s.push_back(1.0/node->volume());

    // the third variable, p
    genius_assert(node_p_offset!=invalid_uint);
    ix.push_back(node->global_offset() + node_p_offset);
    y.push_back(node_data->p());
    s.push_back(1.0/node->volume());

    ix.push_back(node->global_offset() + n_variables + node_p_offset);
    y.push_back(node_data->p());
    s.push_back(1.0/node->volume());


    // for extra Tl temperature equations
    if(get_advanced_model()->enable_Tl())
    {
      genius_assert(node_Tl_offset!=invalid_uint);
      ix.push_back(node->global_offset() + node_Tl_offset);
      y.push_back(node_data->T());
      s.push_back(1.0/node->volume());

      ix.push_back(node->global_offset() + n_variables + node_Tl_offset);
      y.push_back(node_data->T());
      s.push_back(1.0/node->volume());
    }

    // for extra Tn temperature equations, use n*Tn as indepedent variable
    if(get_advanced_model()->enable_Tn())
    {
      genius_assert(node_Tn_offset!=invalid_uint);
      ix.push_back(node->global_offset() + node_Tn_offset);
      y.push_back(node_data->n()*node_data->Tn());
      s.push_back(1.0/node->volume());

      ix.push_back(node->global_offset() + n_variables + node_Tn_offset);
      y.push_back(node_data->n()*node_data->Tn());
      s.push_back(1.0/node->volume());
    }

    // for extra Tp temperature equations, use p*Tp as indepedent variable
    if(get_advanced_model()->enable_Tp())
    {
      genius_assert(node_Tp_offset!=invalid_uint);
      ix.push_back(node->global_offset() + node_Tp_offset);
      y.push_back(node_data->p()*node_data->Tp());
      s.push_back(1.0/node->volume());

      ix.push_back(node->global_offset() + n_variables + node_Tp_offset);
      y.push_back(node_data->p()*node_data->Tp());
      s.push_back(1.0/node->volume());
    }

  }

  // call petsc VecSetValues routine to insert bufferred value
  if( ix.size() )
  {
    VecSetValues(x, ix.size(), &ix[0], &y[0], INSERT_VALUES) ;
    VecSetValues(L, ix.size(), &ix[0], &s[0], INSERT_VALUES) ;
  }

}




void SemiconductorSimulationRegion::DDMAC_Fill_Matrix_Vector(Mat A, Vec b, const Mat J, const double omega, InsertMode &add_value_flag) const
{

  // note, we will use ADD_VALUES to set values of matrix A
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( add_value_flag != ADD_VALUES && add_value_flag != NOT_SET_VALUES)
  {
    MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);
  }

  //
  const_node_iterator it = nodes_begin();
  const_node_iterator it_end = nodes_end();
  for(; it!=it_end; ++it)
  {
    const FVM_Node * fvm_node = (*it).second;

    //if this node NOT belongs to this processor, continue
    if( fvm_node->root_node()->processor_id() != Genius::processor_id() ) continue;

    // if the node on boundary, continue
    if( fvm_node->boundary_id() != BoundaryInfo::invalid_id ) continue;

    this->DDMAC_Fill_Nodal_Matrix_Vector(fvm_node, A, b, J, omega, add_value_flag);
  }

}





void SemiconductorSimulationRegion::DDMAC_Fill_Nodal_Matrix_Vector(
  const FVM_Node *fvm_node, Mat A, Vec, const Mat J, const double omega,
  InsertMode & add_value_flag,
  const SimulationRegion * adjacent_region,
  const FVM_Node * adjacent_fvm_node) const
{

  // each indepedent variable
  std::vector<SolutionVariable> indepedent_variable;

  indepedent_variable.push_back( POTENTIAL );
  indepedent_variable.push_back( ELECTRON );
  indepedent_variable.push_back( HOLE );

  if(this->get_advanced_model()->enable_Tl())
    indepedent_variable.push_back( TEMPERATURE );

  if(this->get_advanced_model()->enable_Tn())
    indepedent_variable.push_back( E_TEMP );

  if(this->get_advanced_model()->enable_Tp())
    indepedent_variable.push_back( H_TEMP );

  genius_assert( fvm_node->root_node()->processor_id() == Genius::processor_id() );


  const FVM_NodeData * node_data = fvm_node->node_data();
  this->mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

  for(unsigned int n=0; n<indepedent_variable.size(); ++n)
  {
    PetscInt jacobian_row = fvm_node->global_offset() + this->ebm_variable_offset(indepedent_variable[n]);

    // get derivative from Jacobian matrix
    PetscInt ncols;
    const PetscInt    * row_cols_pointer;
    const PetscScalar * row_vals_pointer;
    MatGetRow(J, jacobian_row, &ncols, &row_cols_pointer, &row_vals_pointer);

    PetscInt ac_real = fvm_node->global_offset() + this->ebm_variable_offset(indepedent_variable[n]);
    PetscInt ac_imag = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(indepedent_variable[n]);
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL)
    {
      ac_real = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(indepedent_variable[n]);
      ac_imag = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(indepedent_variable[n]);
    }

    std::vector<PetscInt> real_row_entry, imag_row_entry;
    for(int n=0; n<ncols; ++n)
    {
      real_row_entry.push_back(row_cols_pointer[n]);
      imag_row_entry.push_back(row_cols_pointer[n]+this->ebm_n_variables());
    }

    // fill real part to AC matrix
    MatSetValues(A, 1, &ac_real, real_row_entry.size(), &real_row_entry[0], row_vals_pointer, ADD_VALUES );

    // fill image part to AC matrix
    MatSetValues(A, 1, &ac_imag, imag_row_entry.size(), &imag_row_entry[0], row_vals_pointer, ADD_VALUES );

    // restore row pointers of Jacobian matrix
    MatRestoreRow(J, jacobian_row, &ncols, &row_cols_pointer, &row_vals_pointer);
  }

  // process omega item of electron continuity equ
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(ELECTRON);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(ELECTRON) ;
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL)
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(ELECTRON);
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(ELECTRON) ;
    }
    MatSetValue(A, real_row, imag_row,  omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue(A, imag_row, real_row, -omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of hole continuity equ
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(HOLE);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(HOLE) ;
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL)
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(HOLE);
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(HOLE) ;
    }
    MatSetValue(A, real_row, imag_row,  omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue(A, imag_row, real_row, -omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of lattice temperature equ
  if(this->get_advanced_model()->enable_Tl())
  {
    PetscScalar HeatCapacity =  this->mt->thermal->HeatCapacity(node_data->T());

    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(TEMPERATURE);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(TEMPERATURE) ;
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL)
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(TEMPERATURE);
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(TEMPERATURE) ;
    }
    MatSetValue(A, real_row, imag_row,  node_data->density()*HeatCapacity*omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue(A, imag_row, real_row, -node_data->density()*HeatCapacity*omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of electron temperature equ
  if(this->get_advanced_model()->enable_Tn())
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(E_TEMP);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(E_TEMP) ;
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      if( adjacent_region->ebm_variable_offset(E_TEMP)!=invalid_uint )
      {
        real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(E_TEMP);
        imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(E_TEMP) ;
        MatSetValue(A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
        MatSetValue(A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      }
    }
    else
    {
      MatSetValue(A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      MatSetValue(A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
    }
  }

  // process omega item of hole temperature equ
  if(this->get_advanced_model()->enable_Tn())
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(H_TEMP);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(H_TEMP) ;
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL)
    {
      if( adjacent_region->ebm_variable_offset(H_TEMP)!=invalid_uint )
      {
        real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(H_TEMP);
        imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(H_TEMP) ;
        MatSetValue(A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
        MatSetValue(A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      }
    }
    else
    {
      MatSetValue(A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      MatSetValue(A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


}





void SemiconductorSimulationRegion::DDMAC_Fill_Nodal_Matrix_Vector(
  const FVM_Node *fvm_node, const SolutionVariable var,
  Mat A, Vec , const Mat J, const double omega, InsertMode & add_value_flag,
  const SimulationRegion * adjacent_region,
  const FVM_Node * adjacent_fvm_node) const
{

  genius_assert( fvm_node->root_node()->processor_id() == Genius::processor_id() );

  const FVM_NodeData * node_data = fvm_node->node_data();
  this->mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

  PetscInt ncols;
  const PetscInt    * row_cols_pointer;
  const PetscScalar * row_vals_pointer;

  PetscInt jacobian_row = fvm_node->global_offset() + this->ebm_variable_offset(var);

  // get derivative from Jacobian matrix
  MatGetRow(J, jacobian_row, &ncols, &row_cols_pointer, &row_vals_pointer);

  PetscInt ac_real = fvm_node->global_offset() + this->ebm_variable_offset(var);
  PetscInt ac_imag = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(var);
  if(adjacent_region!=NULL && adjacent_fvm_node!=NULL)
  {
    ac_real = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(var);
    ac_imag = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(var);
  }

  std::vector<PetscInt> real_row_entry, imag_row_entry;
  for(int n=0; n<ncols; ++n)
  {
    real_row_entry.push_back(row_cols_pointer[n]);
    imag_row_entry.push_back(row_cols_pointer[n]+this->ebm_n_variables());
  }

  // fill real part to AC matrix
  MatSetValues(A, 1, &ac_real, real_row_entry.size(), &real_row_entry[0], row_vals_pointer, ADD_VALUES );

  // fill image part to AC matrix
  MatSetValues(A, 1, &ac_imag, imag_row_entry.size(), &imag_row_entry[0], row_vals_pointer, ADD_VALUES );

  // restore row pointers of Jacobian matrix
  MatRestoreRow(J, jacobian_row, &ncols, &row_cols_pointer, &row_vals_pointer);


  // process omega item of electron continuity equ
  if(var == ELECTRON)
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(ELECTRON);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(ELECTRON) ;
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL)
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(ELECTRON);
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(ELECTRON) ;
    }
    MatSetValue(A, real_row, imag_row,  omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue(A, imag_row, real_row, -omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of hole continuity equ
  if(var == HOLE)
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(HOLE);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(HOLE) ;
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL)
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(HOLE);
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(HOLE) ;
    }
    MatSetValue(A, real_row, imag_row,  omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue(A, imag_row, real_row, -omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of lattice temperature equ
  if( var == TEMPERATURE && this->get_advanced_model()->enable_Tl())
  {
    PetscScalar HeatCapacity =  this->mt->thermal->HeatCapacity(node_data->T());

    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(TEMPERATURE);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(TEMPERATURE) ;
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL)
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(TEMPERATURE);
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(TEMPERATURE) ;
    }
    MatSetValue(A, real_row, imag_row,  node_data->density()*HeatCapacity*omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue(A, imag_row, real_row, -node_data->density()*HeatCapacity*omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of electron temperature equ
  if( var == E_TEMP && this->get_advanced_model()->enable_Tn())
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(E_TEMP);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(E_TEMP) ;
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      if( adjacent_region->ebm_variable_offset(E_TEMP)!=invalid_uint )
      {
        real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(E_TEMP);
        imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(E_TEMP) ;
        MatSetValue(A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
        MatSetValue(A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      }
    }
    else
    {
      MatSetValue(A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      MatSetValue(A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
    }
  }

  // process omega item of hole temperature equ
  if( var == H_TEMP && this->get_advanced_model()->enable_Tp())
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(H_TEMP);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(H_TEMP) ;
    if(adjacent_region!=NULL && adjacent_fvm_node!=NULL)
    {
      if( adjacent_region->ebm_variable_offset(H_TEMP)!=invalid_uint )
      {
        real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(H_TEMP);
        imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(H_TEMP) ;
        MatSetValue(A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
        MatSetValue(A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      }
    }
    else
    {
      MatSetValue(A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      MatSetValue(A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


}




void SemiconductorSimulationRegion::DDMAC_Update_Solution(PetscScalar *lxx)
{
  unsigned int node_psi_offset = ebm_variable_offset(POTENTIAL);
  unsigned int node_n_offset   = ebm_variable_offset(ELECTRON);
  unsigned int node_p_offset   = ebm_variable_offset(HOLE);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);
  unsigned int node_Tn_offset  = ebm_variable_offset(E_TEMP);
  unsigned int node_Tp_offset  = ebm_variable_offset(H_TEMP);


  node_iterator it = nodes_begin();
  node_iterator it_end = nodes_end();
  for(; it!=it_end; ++it)
  {
    FVM_Node * fvm_node = (*it).second;

    // NOTE: here, solution for all the local node should be updated!
    if( !fvm_node->root_node()->on_local() ) continue;

    FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data!=NULL);


    //update psi
    {
      PetscScalar psi_real = lxx[fvm_node->local_offset() + node_psi_offset];
      PetscScalar psi_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_psi_offset];
      node_data->psi_ac()  = std::complex<PetscScalar>(psi_real, psi_imag);
    }

    // electron density
    {
      PetscScalar n_real = lxx[fvm_node->local_offset() + node_n_offset];
      PetscScalar n_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_n_offset];
      node_data->n_ac()  = std::complex<PetscScalar>(n_real, n_imag);
    }

    // hole density
    {
      PetscScalar p_real = lxx[fvm_node->local_offset() + node_p_offset];
      PetscScalar p_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_p_offset];
      node_data->p_ac()  = std::complex<PetscScalar>(p_real, p_imag);
    }

    // lattice temperature if required
    if(this->get_advanced_model()->enable_Tl())
    {
      PetscScalar T_real = lxx[fvm_node->local_offset() + node_Tl_offset];
      PetscScalar T_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_Tl_offset];
      node_data->T_ac()  = std::complex<PetscScalar>(T_real, T_imag);
    }

    // electron temperature if required
    if(this->get_advanced_model()->enable_Tn())
    {
      PetscScalar nTn_real = lxx[fvm_node->local_offset() + node_Tn_offset];
      PetscScalar nTn_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_Tn_offset];
      node_data->Tn_ac()   = std::complex<PetscScalar>(nTn_real, nTn_imag)/node_data->n_ac();
    }

    // hole temperature if required
    if(this->get_advanced_model()->enable_Tp())
    {
      PetscScalar pTp_real = lxx[fvm_node->local_offset() + node_Tp_offset];
      PetscScalar pTp_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_Tp_offset];
      node_data->Tp_ac()   = std::complex<PetscScalar>(pTp_real, pTp_imag)/node_data->p_ac();
    }
  }

}
