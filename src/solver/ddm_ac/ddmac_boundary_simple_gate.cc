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


void SimpleGateContactBC::DDMAC_Fill_Matrix_Vector ( Mat A, Vec b, const Mat J, const double omega, InsertMode & add_value_flag )
{

  if ( ( add_value_flag != ADD_VALUES ) && ( add_value_flag != NOT_SET_VALUES ) )
  {
    MatAssemblyBegin ( A, MAT_FLUSH_ASSEMBLY );
    MatAssemblyEnd ( A, MAT_FLUSH_ASSEMBLY );
  }

  PetscInt bc_global_offset_re = this->global_offset();
  PetscInt bc_global_offset_im = this->global_offset() +1;

  const PetscScalar q = e*this->scalar("qf");            // surface change density
  const PetscScalar Thick = this->scalar("thickness");   // the thickness of gate oxide
  const PetscScalar eps_ox = this->scalar("eps");        // the permittivity of gate material
  const PetscScalar Work_Function = this->scalar("workfunction");
  const PetscScalar Heat_Transfer = this->scalar("heat.transfer");


  // impedance at frequency omega
  std::complex <PetscScalar> j ( 0.0, 1.0 );

  // for 2D mesh, z_width() is the device dimension in Z direction; for 3D mesh, z_width() is 1.0
  PetscScalar current_scale = this->z_width();

  // we use AD again. no matter it is overkill here.
  //the indepedent variable number, we only need 2 here.
  adtl::AutoDScalar::numdir=2;

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const SimulationRegion * region = ( *region_node_begin ( *node_it ) ).second.first;
    const FVM_Node * fvm_node = ( *region_node_begin ( *node_it ) ).second.second;
    const FVM_NodeData * node_data = fvm_node->node_data();

    // fill A with J
    region->DDMAC_Fill_Nodal_Matrix_Vector ( fvm_node, POTENTIAL, A, b, J, omega, add_value_flag );
    region->DDMAC_Fill_Nodal_Matrix_Vector ( fvm_node, ELECTRON, A, b, J, omega, add_value_flag );
    region->DDMAC_Fill_Nodal_Matrix_Vector ( fvm_node, HOLE, A, b, J, omega, add_value_flag );

    if ( region->get_advanced_model()->enable_Tn() )
      region->DDMAC_Fill_Nodal_Matrix_Vector ( fvm_node, E_TEMP, A, b, J, omega, add_value_flag );

    if ( region->get_advanced_model()->enable_Tp() )
      region->DDMAC_Fill_Nodal_Matrix_Vector ( fvm_node, H_TEMP, A, b, J, omega, add_value_flag );


    // psi of this node
    AutoDScalar V = node_data->psi();
    V.setADValue ( 0, 1.0 );

    // the electrode potential
    genius_assert ( local_offset() !=invalid_uint );
    AutoDScalar Ve = this->ext_circuit()->Vac();
    Ve.setADValue ( 1, 1.0 );

    // area of external surface
    PetscScalar S = fvm_node->outside_boundary_surface_area();

    {
      AutoDScalar dP = S* ( eps_ox* ( Ve - Work_Function-V ) /Thick + q );

      //governing equation of psi

      PetscInt index_re = fvm_node->global_offset() + region->ebm_variable_offset ( POTENTIAL );
      PetscInt col_re[2] = {index_re , bc_global_offset_re};

      PetscInt index_im = fvm_node->global_offset() + region->ebm_n_variables() + region->ebm_variable_offset ( POTENTIAL );
      PetscInt col_im[2] = {index_im , bc_global_offset_im};

      MatSetValues ( A, 1, &index_re, 2, col_re, dP.getADValue(), ADD_VALUES );
      MatSetValues ( A, 1, &index_im, 2, col_im, dP.getADValue(), ADD_VALUES );
    }

    // process the Jacobian of equation of T
    // if this gate bc is external boundary, set heat flux here
    if ( region->get_advanced_model()->enable_Tl() && ( node_on_boundary ( *node_it ) || has_associated_region ( *node_it, VacuumRegion ) ) )
    {
      region->DDMAC_Fill_Nodal_Matrix_Vector ( fvm_node, TEMPERATURE, A, b, J, omega, add_value_flag );

      AutoDScalar T = node_data->T();
      T.setADValue ( 0, 1.0 ); // T of this node
      AutoDScalar fT = Heat_Transfer* ( T_external()-T ) *S;

      PetscInt index_re = fvm_node->global_offset() + region->ebm_variable_offset ( TEMPERATURE );
      PetscInt col_re   = index_re;
      MatSetValue ( A, index_re, col_re, fT.getADValue ( 0 ), ADD_VALUES );

      PetscInt index_im = fvm_node->global_offset() + region->ebm_n_variables() + region->ebm_variable_offset ( TEMPERATURE );
      PetscInt col_im   = index_im;
      MatSetValue ( A, index_im, col_im, fT.getADValue ( 0 ), ADD_VALUES );
    }


    // displacement current

    // area of out surface of control volume related with neighbor node
    AutoDScalar D = eps_ox* ( Ve-V ) /Thick;

    // the 1/dt is replaced by j*omega.
    std::complex<PetscScalar> mna_scaling = ext_circuit()->mna_ac_scaling(omega);
    std::complex <PetscScalar> dJdisp_dV  = mna_scaling*S*D.getADValue ( 0 ) *j*omega*current_scale;
    std::complex <PetscScalar> dJdisp_dVe = mna_scaling*S*D.getADValue ( 1 ) *j*omega*current_scale;

    // V
    MatSetValue ( A, bc_global_offset_re, fvm_node->global_offset(), dJdisp_dV.real(), ADD_VALUES );
    MatSetValue ( A, bc_global_offset_re, fvm_node->global_offset() +region->ebm_n_variables(), -dJdisp_dV.imag(), ADD_VALUES );

    MatSetValue ( A, bc_global_offset_im, fvm_node->global_offset(), dJdisp_dV.imag(), ADD_VALUES );
    MatSetValue ( A, bc_global_offset_im, fvm_node->global_offset() +region->ebm_n_variables(),  dJdisp_dV.real(), ADD_VALUES );

    // Ve
    MatSetValue ( A, bc_global_offset_re, bc_global_offset_re,  dJdisp_dVe.real(), ADD_VALUES );
    MatSetValue ( A, bc_global_offset_re, bc_global_offset_im, -dJdisp_dVe.imag(), ADD_VALUES );

    MatSetValue ( A, bc_global_offset_im, bc_global_offset_re,  dJdisp_dVe.imag(), ADD_VALUES );
    MatSetValue ( A, bc_global_offset_im, bc_global_offset_im,  dJdisp_dVe.real(), ADD_VALUES );

  }


  // the extra equation of gate boundary
  // For ac scan
  //
  //          _____  (Z1)          Ve
  //    -----|_____|----/\/\/\/\-------> to gate electrode (Ve, I)
  //    |       R          L       |
  //   Vac                      C === (Y2)
  //    |__________________________|
  //           GND
  //

  if(Genius::is_last_processor())
  {
    // here we process the external circuit, we do not use AD here

    std::complex <PetscScalar> K = this->ext_circuit()->mna_ac_jacobian(omega);

    MatSetValue(A, bc_global_offset_re, bc_global_offset_re,  K.real(), ADD_VALUES);
    MatSetValue(A, bc_global_offset_re, bc_global_offset_im, -K.imag(), ADD_VALUES);
    MatSetValue(A, bc_global_offset_im, bc_global_offset_re,  K.imag(), ADD_VALUES);
    MatSetValue(A, bc_global_offset_im, bc_global_offset_im,  K.real(), ADD_VALUES);

    VecSetValue(b, bc_global_offset_re, this->ext_circuit()->Vac(), ADD_VALUES);
    VecSetValue(b, bc_global_offset_im, 0.0, ADD_VALUES);
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}





void SimpleGateContactBC::DDMAC_Update_Solution ( const PetscScalar * lxx, const Mat, const double omega )
{

  std::complex<PetscScalar> Iac ( 0.0, 0.0 );
  std::complex<PetscScalar> j ( 0.0, 1.0 );

  // for 2D mesh, system().z_width() is the device dimension in Z direction; for 3D mesh, system().z_width() is 1.0
  PetscScalar current_scale = system().z_width();

  const PetscScalar Thick = this->scalar("thickness");   // the thickness of gate oxide
  const PetscScalar eps_ox = this->scalar("eps");        // the permittivity of gate material

  std::complex<PetscScalar> Ve ( lxx[this->local_offset() ], lxx[this->local_offset() +1] );


  BoundaryCondition::const_node_iterator node_it;
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const SimulationRegion * region = ( *region_node_begin ( *node_it ) ).second.first;
    const FVM_Node * fvm_node = ( *region_node_begin ( *node_it ) ).second.second;
    const FVM_NodeData * node_data = fvm_node->node_data();

    // displacement current
    PetscScalar V_re = lxx[fvm_node->local_offset() + region->ebm_variable_offset ( POTENTIAL ) ];
    PetscScalar V_im = lxx[fvm_node->local_offset() + region->ebm_n_variables() + region->ebm_variable_offset ( POTENTIAL ) ];
    std::complex<PetscScalar> V = std::complex<PetscScalar> ( V_re, V_im );

    PetscScalar S  = fvm_node->outside_boundary_surface_area();
    std::complex<PetscScalar> D = eps_ox* ( Ve-V ) /Thick;

    Iac += D*S*j*omega*current_scale;
  }

  Parallel::sum ( Iac );

  this->ext_circuit()->current_ac()   = Iac;
  this->ext_circuit()->potential_ac() = Ve;
}
