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


void SemiconductorSimulationRegion::DDMAC_Fill_Value ( Vec x, Vec L ) const
{
  // find the node variable offset
  unsigned int n_variables     = this->ebm_n_variables();
  unsigned int node_psi_offset = this->ebm_variable_offset ( POTENTIAL );
  unsigned int node_n_offset   = this->ebm_variable_offset ( ELECTRON );
  unsigned int node_p_offset   = this->ebm_variable_offset ( HOLE );
  unsigned int node_Tl_offset  = this->ebm_variable_offset ( TEMPERATURE );
  unsigned int node_Tn_offset  = this->ebm_variable_offset ( E_TEMP );
  unsigned int node_Tp_offset  = this->ebm_variable_offset ( H_TEMP );

  // data buffer
  std::vector<int> ix;
  std::vector<PetscScalar> y;
  std::vector<PetscScalar> s;

  // reserve menory for data buffer
  PetscInt n_local_dofs;
  VecGetLocalSize ( x, &n_local_dofs );
  ix.reserve ( n_local_dofs );
  y.reserve ( n_local_dofs );
  s.reserve ( n_local_dofs );

  // for all the on processor node, insert value to petsc vector
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;

    const FVM_NodeData * node_data = fvm_node->node_data();

    // the first variable, psi
    genius_assert ( node_psi_offset!=invalid_uint );
    ix.push_back ( fvm_node->global_offset() + node_psi_offset );
    y.push_back ( node_data->psi() );
    s.push_back ( 1.0/fvm_node->volume() );

    ix.push_back ( fvm_node->global_offset() + n_variables + node_psi_offset );
    y.push_back ( node_data->psi() );
    s.push_back ( 1.0/fvm_node->volume() );

    // the second variable, n
    genius_assert ( node_n_offset!=invalid_uint );
    ix.push_back ( fvm_node->global_offset() + node_n_offset );
    y.push_back ( node_data->n() );
    s.push_back ( 1.0/fvm_node->volume() );

    ix.push_back ( fvm_node->global_offset() + n_variables + node_n_offset );
    y.push_back ( node_data->n() );
    s.push_back ( 1.0/fvm_node->volume() );

    // the third variable, p
    genius_assert ( node_p_offset!=invalid_uint );
    ix.push_back ( fvm_node->global_offset() + node_p_offset );
    y.push_back ( node_data->p() );
    s.push_back ( 1.0/fvm_node->volume() );

    ix.push_back ( fvm_node->global_offset() + n_variables + node_p_offset );
    y.push_back ( node_data->p() );
    s.push_back ( 1.0/fvm_node->volume() );


    // for extra Tl temperature equations
    if ( get_advanced_model()->enable_Tl() )
    {
      genius_assert ( node_Tl_offset!=invalid_uint );
      ix.push_back ( fvm_node->global_offset() + node_Tl_offset );
      y.push_back ( node_data->T() );
      s.push_back ( 1.0/fvm_node->volume() );

      ix.push_back ( fvm_node->global_offset() + n_variables + node_Tl_offset );
      y.push_back ( node_data->T() );
      s.push_back ( 1.0/fvm_node->volume() );
    }

    // for extra Tn temperature equations, use n*Tn as indepedent variable
    if ( get_advanced_model()->enable_Tn() )
    {
      genius_assert ( node_Tn_offset!=invalid_uint );
      ix.push_back ( fvm_node->global_offset() + node_Tn_offset );
      y.push_back ( node_data->n() *node_data->Tn() );
      s.push_back ( 1.0/fvm_node->volume() );

      ix.push_back ( fvm_node->global_offset() + n_variables + node_Tn_offset );
      y.push_back ( node_data->n() *node_data->Tn() );
      s.push_back ( 1.0/fvm_node->volume() );
    }

    // for extra Tp temperature equations, use p*Tp as indepedent variable
    if ( get_advanced_model()->enable_Tp() )
    {
      genius_assert ( node_Tp_offset!=invalid_uint );
      ix.push_back ( fvm_node->global_offset() + node_Tp_offset );
      y.push_back ( node_data->p() *node_data->Tp() );
      s.push_back ( 1.0/fvm_node->volume() );

      ix.push_back ( fvm_node->global_offset() + n_variables + node_Tp_offset );
      y.push_back ( node_data->p() *node_data->Tp() );
      s.push_back ( 1.0/fvm_node->volume() );
    }

  }

  // call petsc VecSetValues routine to insert bufferred value
  if ( ix.size() )
  {
    VecSetValues ( x, ix.size(), &ix[0], &y[0], INSERT_VALUES ) ;
    VecSetValues ( L, ix.size(), &ix[0], &s[0], INSERT_VALUES ) ;
  }

}




void SemiconductorSimulationRegion::DDMAC_Fill_Matrix_Vector ( Mat A, Vec b, const Mat J, const double omega, InsertMode &add_value_flag ) const
{

  // note, we will use ADD_VALUES to set values of matrix A
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if ( add_value_flag != ADD_VALUES && add_value_flag != NOT_SET_VALUES )
  {
    MatAssemblyBegin ( A, MAT_FLUSH_ASSEMBLY );
    MatAssemblyEnd ( A, MAT_FLUSH_ASSEMBLY );
  }

  //
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;

    // if the node on boundary, continue
    if ( fvm_node->boundary_id() != BoundaryInfo::invalid_id ) continue;

    this->DDMAC_Fill_Nodal_Matrix_Vector ( fvm_node, A, b, J, omega, add_value_flag );
  }

}


void SemiconductorSimulationRegion::DDMAC_Fill_Transformation_Matrix ( Mat T, const Mat J, const double omega,  InsertMode &add_value_flag ) const
{

  // note, we will use ADD_VALUES to set values of matrix A
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if ( add_value_flag != ADD_VALUES && add_value_flag != NOT_SET_VALUES )
  {
    MatAssemblyBegin ( T, MAT_FLUSH_ASSEMBLY );
    MatAssemblyEnd ( T, MAT_FLUSH_ASSEMBLY );
  }

  // each indepedent variable
  std::vector<SolutionVariable> indepedent_variable;

  indepedent_variable.push_back ( POTENTIAL );
  indepedent_variable.push_back ( ELECTRON );
  indepedent_variable.push_back ( HOLE );

  if ( this->get_advanced_model()->enable_Tl() )
    indepedent_variable.push_back ( TEMPERATURE );

  if ( this->get_advanced_model()->enable_Tn() )
    indepedent_variable.push_back ( E_TEMP );

  if ( this->get_advanced_model()->enable_Tp() )
    indepedent_variable.push_back ( H_TEMP );

  //
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;

    for ( unsigned int n=0; n<indepedent_variable.size(); ++n )
    {
      PetscInt diag_index = fvm_node->global_offset() + this->ebm_variable_offset ( indepedent_variable[n] );
      PetscScalar diag_value;

      MatGetValues(J, 1, &diag_index, 1, &diag_index, &diag_value);

      PetscInt ac_real = diag_index;
      PetscInt ac_imag = diag_index + this->ebm_n_variables();
      MatSetValue ( T, ac_real, ac_real, 1.0, ADD_VALUES );
      if( indepedent_variable[n] != POTENTIAL )
      {
        MatSetValue ( T, ac_real, ac_imag, omega/diag_value, ADD_VALUES );
        MatSetValue ( T, ac_imag, ac_real, -omega/diag_value, ADD_VALUES );
      }
      MatSetValue ( T, ac_imag, ac_imag, 1.0, ADD_VALUES );
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


void SemiconductorSimulationRegion::DDMAC_Fill_Nodal_Matrix_Vector (
  const FVM_Node *fvm_node, Mat A, Vec, const Mat J, const double omega,
  InsertMode & add_value_flag,
  const SimulationRegion * adjacent_region,
  const FVM_Node * adjacent_fvm_node ) const
{

  // each indepedent variable
  std::vector<SolutionVariable> indepedent_variable;

  indepedent_variable.push_back ( POTENTIAL );
  indepedent_variable.push_back ( ELECTRON );
  indepedent_variable.push_back ( HOLE );

  if ( this->get_advanced_model()->enable_Tl() )
    indepedent_variable.push_back ( TEMPERATURE );

  if ( this->get_advanced_model()->enable_Tn() )
    indepedent_variable.push_back ( E_TEMP );

  if ( this->get_advanced_model()->enable_Tp() )
    indepedent_variable.push_back ( H_TEMP );

  genius_assert ( fvm_node->root_node()->processor_id() == Genius::processor_id() );


  const FVM_NodeData * node_data = fvm_node->node_data();
  this->mt->mapping ( fvm_node->root_node(), node_data, SolverSpecify::clock );

  for ( unsigned int n=0; n<indepedent_variable.size(); ++n )
  {
    PetscInt jacobian_row = fvm_node->global_offset() + this->ebm_variable_offset ( indepedent_variable[n] );

    // get derivative from Jacobian matrix
    PetscInt ncols;
    const PetscInt    * row_cols_pointer;
    const PetscScalar * row_vals_pointer;
    MatGetRow ( J, jacobian_row, &ncols, &row_cols_pointer, &row_vals_pointer );

    PetscInt ac_real = fvm_node->global_offset() + this->ebm_variable_offset ( indepedent_variable[n] );
    PetscInt ac_imag = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( indepedent_variable[n] );
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      ac_real = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( indepedent_variable[n] );
      ac_imag = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( indepedent_variable[n] );
    }

    std::vector<PetscInt> real_row_entry, imag_row_entry;
    for ( int n=0; n<ncols; ++n )
    {
      real_row_entry.push_back ( row_cols_pointer[n] );
      imag_row_entry.push_back ( row_cols_pointer[n]+this->ebm_n_variables() );
    }

    // fill real part to AC matrix
    MatSetValues ( A, 1, &ac_real, real_row_entry.size(), &real_row_entry[0], row_vals_pointer, ADD_VALUES );

    // fill image part to AC matrix
    MatSetValues ( A, 1, &ac_imag, imag_row_entry.size(), &imag_row_entry[0], row_vals_pointer, ADD_VALUES );

    // restore row pointers of Jacobian matrix
    MatRestoreRow ( J, jacobian_row, &ncols, &row_cols_pointer, &row_vals_pointer );
  }

  // process omega item of electron continuity equ
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset ( ELECTRON );
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( ELECTRON ) ;
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( ELECTRON );
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( ELECTRON ) ;
    }
    MatSetValue ( A, real_row, imag_row,  omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue ( A, imag_row, real_row, -omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of hole continuity equ
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset ( HOLE );
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( HOLE ) ;
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( HOLE );
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( HOLE ) ;
    }
    MatSetValue ( A, real_row, imag_row,  omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue ( A, imag_row, real_row, -omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of lattice temperature equ
  if ( this->get_advanced_model()->enable_Tl() )
  {
    PetscScalar HeatCapacity =  this->mt->thermal->HeatCapacity ( node_data->T() );

    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset ( TEMPERATURE );
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( TEMPERATURE ) ;
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( TEMPERATURE );
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( TEMPERATURE ) ;
    }
    MatSetValue ( A, real_row, imag_row,  node_data->density() *HeatCapacity*omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue ( A, imag_row, real_row, -node_data->density() *HeatCapacity*omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of electron temperature equ
  if ( this->get_advanced_model()->enable_Tn() )
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset ( E_TEMP );
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( E_TEMP ) ;
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      if ( adjacent_region->ebm_variable_offset ( E_TEMP ) !=invalid_uint )
      {
        real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( E_TEMP );
        imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( E_TEMP ) ;
        MatSetValue ( A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
        MatSetValue ( A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      }
    }
    else
    {
      MatSetValue ( A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      MatSetValue ( A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
    }
  }

  // process omega item of hole temperature equ
  if ( this->get_advanced_model()->enable_Tn() )
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset ( H_TEMP );
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( H_TEMP ) ;
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      if ( adjacent_region->ebm_variable_offset ( H_TEMP ) !=invalid_uint )
      {
        real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( H_TEMP );
        imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( H_TEMP ) ;
        MatSetValue ( A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
        MatSetValue ( A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      }
    }
    else
    {
      MatSetValue ( A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      MatSetValue ( A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


}





void SemiconductorSimulationRegion::DDMAC_Fill_Nodal_Matrix_Vector (
  const FVM_Node *fvm_node, const SolutionVariable var,
  Mat A, Vec , const Mat J, const double omega, InsertMode & add_value_flag,
  const SimulationRegion * adjacent_region,
  const FVM_Node * adjacent_fvm_node ) const
{

  genius_assert ( fvm_node->root_node()->processor_id() == Genius::processor_id() );

  const FVM_NodeData * node_data = fvm_node->node_data();
  this->mt->mapping ( fvm_node->root_node(), node_data, SolverSpecify::clock );

  PetscInt ncols;
  const PetscInt    * row_cols_pointer;
  const PetscScalar * row_vals_pointer;

  PetscInt jacobian_row = fvm_node->global_offset() + this->ebm_variable_offset ( var );

  // get derivative from Jacobian matrix
  MatGetRow ( J, jacobian_row, &ncols, &row_cols_pointer, &row_vals_pointer );

  PetscInt ac_real = fvm_node->global_offset() + this->ebm_variable_offset ( var );
  PetscInt ac_imag = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( var );
  if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
  {
    ac_real = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( var );
    ac_imag = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( var );
  }

  std::vector<PetscInt> real_row_entry, imag_row_entry;
  for ( int n=0; n<ncols; ++n )
  {
    real_row_entry.push_back ( row_cols_pointer[n] );
    imag_row_entry.push_back ( row_cols_pointer[n]+this->ebm_n_variables() );
  }

  // fill real part to AC matrix
  MatSetValues ( A, 1, &ac_real, real_row_entry.size(), &real_row_entry[0], row_vals_pointer, ADD_VALUES );

  // fill image part to AC matrix
  MatSetValues ( A, 1, &ac_imag, imag_row_entry.size(), &imag_row_entry[0], row_vals_pointer, ADD_VALUES );

  // restore row pointers of Jacobian matrix
  MatRestoreRow ( J, jacobian_row, &ncols, &row_cols_pointer, &row_vals_pointer );


  // process omega item of electron continuity equ
  if ( var == ELECTRON )
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset ( ELECTRON );
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( ELECTRON ) ;
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( ELECTRON );
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( ELECTRON ) ;
    }
    MatSetValue ( A, real_row, imag_row,  omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue ( A, imag_row, real_row, -omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of hole continuity equ
  if ( var == HOLE )
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset ( HOLE );
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( HOLE ) ;
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( HOLE );
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( HOLE ) ;
    }
    MatSetValue ( A, real_row, imag_row,  omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue ( A, imag_row, real_row, -omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of lattice temperature equ
  if ( var == TEMPERATURE && this->get_advanced_model()->enable_Tl() )
  {
    PetscScalar HeatCapacity =  this->mt->thermal->HeatCapacity ( node_data->T() );

    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset ( TEMPERATURE );
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( TEMPERATURE ) ;
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( TEMPERATURE );
      imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( TEMPERATURE ) ;
    }
    MatSetValue ( A, real_row, imag_row,  node_data->density() *HeatCapacity*omega*fvm_node->volume(), ADD_VALUES );
    MatSetValue ( A, imag_row, real_row, -node_data->density() *HeatCapacity*omega*fvm_node->volume(), ADD_VALUES );
  }

  // process omega item of electron temperature equ
  if ( var == E_TEMP && this->get_advanced_model()->enable_Tn() )
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset ( E_TEMP );
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( E_TEMP ) ;
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      if ( adjacent_region->ebm_variable_offset ( E_TEMP ) !=invalid_uint )
      {
        real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( E_TEMP );
        imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( E_TEMP ) ;
        MatSetValue ( A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
        MatSetValue ( A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      }
    }
    else
    {
      MatSetValue ( A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      MatSetValue ( A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
    }
  }

  // process omega item of hole temperature equ
  if ( var == H_TEMP && this->get_advanced_model()->enable_Tp() )
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset ( H_TEMP );
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset ( H_TEMP ) ;
    if ( adjacent_region!=NULL && adjacent_fvm_node!=NULL )
    {
      if ( adjacent_region->ebm_variable_offset ( H_TEMP ) !=invalid_uint )
      {
        real_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset ( H_TEMP );
        imag_row = adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset ( H_TEMP ) ;
        MatSetValue ( A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
        MatSetValue ( A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      }
    }
    else
    {
      MatSetValue ( A, real_row, imag_row,  1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
      MatSetValue ( A, imag_row, real_row, -1.5*kb*omega*fvm_node->volume(), ADD_VALUES );
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


}



void SemiconductorSimulationRegion::DDMAC_Force_equal(const FVM_Node *fvm_node, Mat A, InsertMode & add_value_flag,
                                              const SimulationRegion * adjacent_region, const FVM_Node * adjacent_fvm_node) const
{
  genius_assert(add_value_flag == ADD_VALUES);

  // set entry on diag of A
  PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(POTENTIAL);
  PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(POTENTIAL);
  PetscInt real_col[2]={real_row, adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(POTENTIAL)};
  PetscInt imag_col[2]={imag_row, adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(POTENTIAL)};
  PetscScalar diag[2] = {1.0, -1.0};
  MatSetValues(A, 1, &real_row, 2, real_col, diag, ADD_VALUES);
  MatSetValues(A, 1, &imag_row, 2, imag_col, diag, ADD_VALUES);

  if( adjacent_region->ebm_variable_offset ( ELECTRON ) != invalid_uint )
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(ELECTRON);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(ELECTRON);
    PetscInt real_col[2]={real_row, adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(ELECTRON)};
    PetscInt imag_col[2]={imag_row, adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(ELECTRON)};
    PetscScalar diag[2] = {1.0, -1.0};
    MatSetValues(A, 1, &real_row, 2, real_col, diag, ADD_VALUES);
    MatSetValues(A, 1, &imag_row, 2, imag_col, diag, ADD_VALUES);
  }


  if( adjacent_region->ebm_variable_offset ( HOLE ) != invalid_uint )
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(HOLE);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(HOLE);
    PetscInt real_col[2]={real_row, adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(HOLE)};
    PetscInt imag_col[2]={imag_row, adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(HOLE)};
    PetscScalar diag[2] = {1.0, -1.0};
    MatSetValues(A, 1, &real_row, 2, real_col, diag, ADD_VALUES);
    MatSetValues(A, 1, &imag_row, 2, imag_col, diag, ADD_VALUES);
  }

  if(this->get_advanced_model()->enable_Tl())
  {
    PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(TEMPERATURE);
    PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(TEMPERATURE);
    PetscInt real_col[2]={real_row, adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(TEMPERATURE)};
    PetscInt imag_col[2]={imag_row, adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(TEMPERATURE)};
    PetscScalar diag[2] = {1.0, -1.0};
    MatSetValues(A, 1, &real_row, 2, real_col, diag, ADD_VALUES);
    MatSetValues(A, 1, &imag_row, 2, imag_col, diag, ADD_VALUES);
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}



void SemiconductorSimulationRegion::DDMAC_Force_equal(const FVM_Node *fvm_node, const SolutionVariable var, Mat A, InsertMode & add_value_flag,
                                              const SimulationRegion * adjacent_region, const FVM_Node * adjacent_fvm_node) const
{
  genius_assert(add_value_flag == ADD_VALUES);
  genius_assert(this->ebm_variable_offset(var)!=invalid_uint);
  genius_assert(adjacent_region->ebm_variable_offset(var)!=invalid_uint);

  // set entry on diag of A
  PetscInt real_row = fvm_node->global_offset() + this->ebm_variable_offset(var);
  PetscInt imag_row = fvm_node->global_offset() + this->ebm_n_variables() + this->ebm_variable_offset(var);
  PetscInt real_col[2]={real_row, adjacent_fvm_node->global_offset() + adjacent_region->ebm_variable_offset(var)};
  PetscInt imag_col[2]={imag_row, adjacent_fvm_node->global_offset() + adjacent_region->ebm_n_variables() + adjacent_region->ebm_variable_offset(var)};
  PetscScalar diag[2] = {1.0, -1.0};
  MatSetValues(A, 1, &real_row, 2, real_col, diag, ADD_VALUES);
  MatSetValues(A, 1, &imag_row, 2, imag_col, diag, ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}



void SemiconductorSimulationRegion::DDMAC_Update_Solution ( PetscScalar *lxx )
{
  unsigned int node_psi_offset = ebm_variable_offset ( POTENTIAL );
  unsigned int node_n_offset   = ebm_variable_offset ( ELECTRON );
  unsigned int node_p_offset   = ebm_variable_offset ( HOLE );
  unsigned int node_Tl_offset  = ebm_variable_offset ( TEMPERATURE );
  unsigned int node_Tn_offset  = ebm_variable_offset ( E_TEMP );
  unsigned int node_Tp_offset  = ebm_variable_offset ( H_TEMP );


  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;

    FVM_NodeData * node_data = fvm_node->node_data();
    genius_assert ( node_data!=NULL );


    //update psi
    {
      PetscScalar psi_real = lxx[fvm_node->local_offset() + node_psi_offset];
      PetscScalar psi_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_psi_offset];
      node_data->psi_ac()  = std::complex<PetscScalar> ( psi_real, psi_imag );
    }

    // electron density
    {
      PetscScalar n_real = lxx[fvm_node->local_offset() + node_n_offset];
      PetscScalar n_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_n_offset];
      node_data->n_ac()  = std::complex<PetscScalar> ( n_real, n_imag );
    }

    // hole density
    {
      PetscScalar p_real = lxx[fvm_node->local_offset() + node_p_offset];
      PetscScalar p_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_p_offset];
      node_data->p_ac()  = std::complex<PetscScalar> ( p_real, p_imag );
    }

    // lattice temperature if required
    if ( this->get_advanced_model()->enable_Tl() )
    {
      PetscScalar T_real = lxx[fvm_node->local_offset() + node_Tl_offset];
      PetscScalar T_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_Tl_offset];
      node_data->T_ac()  = std::complex<PetscScalar> ( T_real, T_imag );
    }

    // electron temperature if required
    if ( this->get_advanced_model()->enable_Tn() )
    {
      PetscScalar nTn_real = lxx[fvm_node->local_offset() + node_Tn_offset];
      PetscScalar nTn_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_Tn_offset];
      node_data->Tn_ac()   = std::complex<PetscScalar> ( nTn_real, nTn_imag ) /node_data->n_ac();
    }

    // hole temperature if required
    if ( this->get_advanced_model()->enable_Tp() )
    {
      PetscScalar pTp_real = lxx[fvm_node->local_offset() + node_Tp_offset];
      PetscScalar pTp_imag = lxx[fvm_node->local_offset() + this->ebm_n_variables() + node_Tp_offset];
      node_data->Tp_ac()   = std::complex<PetscScalar> ( pTp_real, pTp_imag ) /node_data->p_ac();
    }
  }

}
