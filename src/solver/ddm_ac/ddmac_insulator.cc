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
#include "boundary_info.h"



void InsulatorSimulationRegion::DDMAC_Fill_Value(Vec x, Vec L) const
{

  // find the node variable offset
  unsigned int n_variables     = this->ebm_n_variables();
  unsigned int node_psi_offset = this->ebm_variable_offset(POTENTIAL);
  unsigned int node_Tl_offset  = this->ebm_variable_offset(TEMPERATURE);


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

    // psi, real part
    genius_assert(node_psi_offset!=invalid_uint);
    ix.push_back(node->global_offset() + node_psi_offset);
    y.push_back(node_data->psi());
    s.push_back(1.0/node->volume());

    // image part
    ix.push_back(node->global_offset() + n_variables + node_psi_offset);
    y.push_back(node_data->psi());
    s.push_back(1.0/node->volume());


    // for extra Tl temperature equations
    if(get_advanced_model()->enable_Tl())
    {
      // T, real part
      genius_assert(node_Tl_offset!=invalid_uint);
      ix.push_back(node->global_offset() + node_Tl_offset);
      y.push_back(node_data->T());
      s.push_back(1.0/node->volume());

      // T, image part
      ix.push_back(node->global_offset() + n_variables + node_Tl_offset);
      y.push_back(node_data->T());
      s.push_back(1.0/node->volume());
    }
  }

  if( ix.size() )
  {
    VecSetValues(x, ix.size(), &ix[0], &y[0], INSERT_VALUES) ;
    VecSetValues(L, ix.size(), &ix[0], &s[0], INSERT_VALUES) ;
  }

}






void InsulatorSimulationRegion::DDMAC_Fill_Matrix_Vector(Mat A, Vec b, const Mat J, const double omega, InsertMode &add_value_flag) const
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

    // if the node on electrode boundary, continue
    if( fvm_node->boundary_id() != BoundaryInfo::invalid_id ) continue;

    this->DDMAC_Fill_Nodal_Matrix_Vector(fvm_node, A, b, J, omega, add_value_flag);
  }

}





void InsulatorSimulationRegion::DDMAC_Fill_Nodal_Matrix_Vector(
  const FVM_Node *fvm_node, Mat A, Vec, const Mat J, const double omega,
  InsertMode & add_value_flag,
  const SimulationRegion * adjacent_region,
  const FVM_Node * adjacent_fvm_node) const
{

  // each indepedent variable
  std::vector<SolutionVariable> indepedent_variable;

  indepedent_variable.push_back( POTENTIAL );

  if(this->get_advanced_model()->enable_Tl())
    indepedent_variable.push_back( TEMPERATURE );

  genius_assert( fvm_node->root_node()->processor_id() == Genius::processor_id() );


  const FVM_NodeData * node_data = fvm_node->node_data();
  this->mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

  for(unsigned int n=0; n<indepedent_variable.size(); ++n)
  {
    PetscInt ncols;
    const PetscInt    * row_cols_pointer;
    const PetscScalar * row_vals_pointer;

    PetscInt jacobian_row = fvm_node->global_offset() + this->ebm_variable_offset(indepedent_variable[n]);

    // get derivative from Jacobian matrix
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

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


}





void InsulatorSimulationRegion::DDMAC_Fill_Nodal_Matrix_Vector(
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


  // process omega item of lattice temperature equ
  if(var == TEMPERATURE && this->get_advanced_model()->enable_Tl() )
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

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




void InsulatorSimulationRegion::DDMAC_Update_Solution(PetscScalar *lxx)
{
  unsigned int node_psi_offset = ebm_variable_offset(POTENTIAL);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);

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

    // lattice temperature if required
    if(this->get_advanced_model()->enable_Tl())
    {
      PetscScalar T_real = lxx[fvm_node->local_offset() + node_Tl_offset];
      PetscScalar T_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_Tl_offset];
      node_data->T_ac()  = std::complex<PetscScalar>(T_real, T_imag);
    }
  }

}

