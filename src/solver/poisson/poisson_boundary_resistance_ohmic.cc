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

#include "asinh.hpp" // for asinh

// Local includes
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "resistance_region.h"
#include "insulator_region.h"
#include "boundary_condition_resistance_ohmic.h"

#include "parallel.h"
#include "mathfunc.h"

using PhysicalUnit::cm;
using PhysicalUnit::kb;
using PhysicalUnit::e;


/*---------------------------------------------------------------------
 * set scaling constant
 */
void IF_Metal_OhmicBC::Poissin_Fill_Value(Vec , Vec L)
{

  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( ; rnode_it!=end_rnode_it; ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      switch ( region->type() )
      {
          case SemiconductorRegion:
          {
            const FVM_Node * fvm_node = ( *rnode_it ).second.second;
            VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
            break;
          }

          case MetalRegion :
          {
            const FVM_Node * fvm_node = ( *rnode_it ).second.second;
            VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
            break;
          }

          case InsulatorRegion:
          {
            const FVM_Node * fvm_node = ( *rnode_it ).second.second;
            VecSetValue(L, fvm_node->global_offset(), 1.0, INSERT_VALUES);
            break;
          }
          default: break;
      }
    }
  }

}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------
 * do pre-process to function for poisson solver
 */
void IF_Metal_OhmicBC::Poissin_Function_Preprocess(Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;


  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;


    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );

    clear_row.push_back ( semiconductor_node->global_offset() );
    clear_row.push_back ( resistance_node->global_offset() );

    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;
        const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
        clear_row.push_back ( insulator_node->global_offset() );
      }
    }
  }

}


/*---------------------------------------------------------------------
 * build function and its jacobian for Poissin solver
 */
void IF_Metal_OhmicBC::Poissin_Function ( PetscScalar * x, Vec f, InsertMode &add_value_flag )
{
  // Ohmic boundary condition is processed here.

  const PetscScalar T = T_external();

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  genius_assert(_r1->type() == SemiconductorRegion);
  genius_assert(_r2->type() == MetalRegion);

  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;


    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data = resistance_node->node_data();

    const PetscScalar V_resistance = x[resistance_node->local_offset() ];
    const PetscScalar V_semiconductor  = x[semiconductor_node->local_offset() ];

    // process semiconductor region
    {
      const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(_r1);
      semi_region->material()->mapping(semiconductor_node->root_node(), semiconductor_node_data, 0.0);
      PetscScalar ni  = semi_region->material()->band->ni(T);

      PetscScalar Nc  = semiconductor_node_data->Nc();
      PetscScalar Nv  = semiconductor_node_data->Nv();
      PetscScalar Eg  = semiconductor_node_data->Eg();

      //governing equation for Ohmic contact boundary
      PetscScalar ff = V_semiconductor - kb*T/e*boost::math::asinh(semiconductor_node_data->Net_doping()/(2*ni))
                       + Eg/(2*e)
                       + kb*T*log(Nc/Nv)/(2*e)
                       + semiconductor_node_data->affinity();

      // set governing equation to function vector
      VecSetValue(f, semiconductor_node->global_offset(), ff, ADD_VALUES);
    }

    // process metal region
    {
      //governing equation for Ohmic contact boundary
      PetscScalar ff = V_resistance + resistance_node_data->affinity();

      // set governing equation to function vector
      VecSetValue(f, resistance_node->global_offset(), ff, ADD_VALUES);
    }

    // if we have insulator node?
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;

        const FVM_Node * insulator_node = ( *rnode_it ).second.second;
        const PetscScalar V_insulator  = x[insulator_node->local_offset() ];
        PetscScalar f_phi =  V_insulator - V_resistance;
        VecSetValue ( f, insulator_node->global_offset(), f_phi, ADD_VALUES );
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for Poissin solver
 */
void IF_Metal_OhmicBC::Poissin_Jacobian_Reserve ( Mat *jac, InsertMode &add_value_flag )
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.


  // since we will use ADD_VALUES operat, check the matrix state.
  if ( ( add_value_flag != ADD_VALUES ) && ( add_value_flag != NOT_SET_VALUES ) )
  {
    MatAssemblyBegin ( *jac, MAT_FLUSH_ASSEMBLY );
    MatAssemblyEnd ( *jac, MAT_FLUSH_ASSEMBLY );
  }

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;


  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;


    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );

    // process insulator region when necessary
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;
        const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
        MatSetValue ( *jac, insulator_node->global_offset(), resistance_node->global_offset(), 0, ADD_VALUES );
      }
    }
  }


  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}



/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for poisson solver
 */
void IF_Metal_OhmicBC::Poissin_Jacobian_Preprocess(Mat *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;


  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;


    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );

    clear_row.push_back ( semiconductor_node->global_offset() );
    clear_row.push_back ( resistance_node->global_offset() );

    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;
        const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
        clear_row.push_back ( insulator_node->global_offset() );
      }
    }
  }

}


/*---------------------------------------------------------------------
 * build function and its jacobian for Poissin solver
 */
void IF_Metal_OhmicBC::Poissin_Jacobian ( PetscScalar * x, Mat *jac, InsertMode &add_value_flag )
{
  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  const PetscScalar T = T_external();

  const SimulationRegion * _r1 = bc_regions().first;
  const SimulationRegion * _r2 = bc_regions().second;

  adtl::AutoDScalar::numdir=2;

  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    const FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
    const FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

    const FVM_Node * resistance_node = get_region_fvm_node ( ( *node_it ), _r2 );
    const FVM_NodeData * resistance_node_data = resistance_node->node_data();


    // process semiconductor region
    {
      AutoDScalar V_semiconductor  = x[semiconductor_node->local_offset() ]; V_semiconductor.setADValue (0, 1.0 );
      const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *>(_r1);
      semi_region->material()->mapping(semiconductor_node->root_node(), semiconductor_node_data, 0.0);
      PetscScalar ni  = semi_region->material()->band->ni(T);

      PetscScalar Nc  = semiconductor_node_data->Nc();
      PetscScalar Nv  = semiconductor_node_data->Nv();
      PetscScalar Eg  = semiconductor_node_data->Eg();

      //governing equation for Ohmic contact boundary
      AutoDScalar ff = V_semiconductor - kb*T/e*boost::math::asinh(semiconductor_node_data->Net_doping()/(2*ni))
                       + Eg/(2*e)
                       + kb*T*log(Nc/Nv)/(2*e)
                       + semiconductor_node_data->affinity();

      // set governing equation to function vector
      MatSetValue(*jac, semiconductor_node->global_offset(), semiconductor_node->global_offset(), ff.getADValue ( 0 ), ADD_VALUES);
    }


    // process resistance region
    {
      AutoDScalar V_resistance = x[resistance_node->local_offset() ];  V_resistance.setADValue ( 0, 1.0 );
      //governing equation for Ohmic contact boundary
      AutoDScalar ff = V_resistance + resistance_node_data->affinity();
      // set governing equation to function vector
      MatSetValue(*jac, resistance_node->global_offset(), resistance_node->global_offset(), ff.getADValue ( 0 ), ADD_VALUES);
    }

    // if we have insulator node?
    if ( has_associated_region ( ( *node_it ), InsulatorRegion ) )
    {
      AutoDScalar V_resistance = x[resistance_node->local_offset() ];  V_resistance.setADValue ( 0, 1.0 );

      BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
      BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
      for ( ; rnode_it!=end_rnode_it; ++rnode_it )
      {
        const SimulationRegion * region = ( *rnode_it ).second.first;
        if ( region->type() != InsulatorRegion ) continue;
        const FVM_Node * insulator_node  = ( *rnode_it ).second.second;
        AutoDScalar V_insulator  = x[insulator_node->local_offset() ]; V_insulator.setADValue ( 1, 1.0 );
        AutoDScalar f_phi =  V_insulator - V_resistance;
        MatSetValue ( *jac, insulator_node->global_offset(), resistance_node->global_offset(), f_phi.getADValue ( 0 ), ADD_VALUES );
        MatSetValue ( *jac, insulator_node->global_offset(), insulator_node->global_offset(), f_phi.getADValue ( 1 ), ADD_VALUES );
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}


/*---------------------------------------------------------------------
 * update
 */
void IF_Metal_OhmicBC::Poissin_Update_Solution ( PetscScalar * lxx )
{
#if 0
  PetscScalar T   = T_external();

  // search and process all the boundary nodes
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( ; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );
    for ( unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it )
    {
      const SimulationRegion * _r1 = bc_regions().first;
      const SemiconductorSimulationRegion * semiconductor_region = dynamic_cast<const SemiconductorSimulationRegion *> ( _r1 );
      FVM_Node * semiconductor_node  = get_region_fvm_node ( ( *node_it ), _r1 );
      FVM_NodeData * semiconductor_node_data = semiconductor_node->node_data();

      // mapping this node to material library
      semiconductor_region->material()->mapping ( semiconductor_node->root_node(), semiconductor_node_data, SolverSpecify::clock );

      PetscScalar affinity = semiconductor_node_data->affinity();
      PetscScalar ni  = semiconductor_region->material()->band->ni ( T );
      PetscScalar Nc  = semiconductor_region->material()->band->Nc ( T );
      PetscScalar Nv  = semiconductor_region->material()->band->Nv ( T );
      PetscScalar Eg  = semiconductor_region->material()->band->Eg ( T );

      PetscScalar  electron_density;
      PetscScalar  hole_density;
      PetscScalar  net_dpoing = semiconductor_node_data->Net_doping();
      if ( net_dpoing <0 )                  //p-type
      {
        hole_density = ( -net_dpoing + sqrt ( net_dpoing*net_dpoing + 4*ni*ni ) ) /2.0;
        electron_density = ni*ni/hole_density;
      }
      else                                  //n-type
      {
        electron_density = ( net_dpoing + sqrt ( net_dpoing*net_dpoing + 4*ni*ni ) ) /2.0;
        hole_density = ni*ni/electron_density;
      }

      PetscScalar V   = lxx[semiconductor_node->local_offset()];

      // intrinsic Fermi potential.
      PetscScalar V_i = V
                        + affinity
                        + Eg / ( 2*e )
                        + kb*T*log ( Nc /Nv ) / ( 2*e );
      // electron density
      PetscScalar n    =  ni*exp ( e/ ( kb*T ) *V_i );
      semiconductor_node_data->n()   =  n;
      // hole density
      PetscScalar p    =  ni*exp ( -e/ ( kb*T ) *V_i );
      semiconductor_node_data->p()   =  p;
    }
  }
#endif
}

