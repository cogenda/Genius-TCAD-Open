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
#include "boundary_condition_hetero.h"
#include "petsc_utils.h"
#include "mathfunc.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

//FIXME how to compute current pass through hetero junction?
// Medici says there are thermal emit and tunneling current
// However, dessis only considers thermal emit current.
// what shoud I do?


/*---------------------------------------------------------------------
 * do pre-process to function for EBM3 solver
 */
void HeteroInterfaceBC::EBM3_Function_Preprocess(PetscScalar *,Vec f, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      if(!fvm_node->is_valid()) continue;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0)
      {
        // do nothing.
        // however, we will add fvm integral of other regions to it.
      }

      // other  region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(POTENTIAL));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(POTENTIAL));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(POTENTIAL));

              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(ELECTRON));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(ELECTRON));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(ELECTRON));

              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(HOLE));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(HOLE));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(HOLE));

              // when have the lattice temperature equation, we add dst to src
              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                // lattice temperature equation should be global
                genius_assert(regions[0]->get_advanced_model()->enable_Tl());
                src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
                dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              }

              // when both 2 regions have the electron/hole temperature equation, we add dst to src
              if(regions[i]->get_advanced_model()->enable_Tn())
              {
                if(regions[0]->get_advanced_model()->enable_Tn())
                {
                  src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(E_TEMP));
                  dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(E_TEMP));
                }
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(E_TEMP));
              }

              if(regions[i]->get_advanced_model()->enable_Tp())
              {
                if(regions[0]->get_advanced_model()->enable_Tp())
                {
                  src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(H_TEMP));
                  dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(H_TEMP));
                }
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(H_TEMP));
              }

              break;
            }

            case InsulatorRegion:
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+0);

              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
                dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              }
              break;
            }
            default: genius_error();
        }
      }
    }

  }

}



/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void HeteroInterfaceBC::EBM3_Function ( PetscScalar * x, Vec f, InsertMode &add_value_flag )
{

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if ( ( add_value_flag != ADD_VALUES ) && ( add_value_flag != NOT_SET_VALUES ) )
  {
    VecAssemblyBegin ( f );
    VecAssemblyEnd ( f );
  }


  // buffer for Vec value
  std::vector<PetscInt> iy;
  std::vector<PetscScalar> y_new;

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for ( ; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;


    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;
    // the variable for first region
    PetscScalar V0;
    PetscScalar n0;
    PetscScalar p0;
    PetscScalar T0;
    PetscScalar Tn0;
    PetscScalar Tp0;
    Material::MaterialSemiconductor * mt0;
    PetscScalar Ec0,Ev0;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );

    for ( unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;
      if ( !fvm_node->is_valid() ) continue;

      regions.push_back ( region );
      fvm_nodes.push_back ( fvm_node );

      // the first semiconductor region
      if ( i==0 )
      {
        const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *> ( regions[i] );
        const FVM_NodeData * n0_data = fvm_nodes[i]->node_data();

        unsigned int node_psi_offset = regions[i]->ebm_variable_offset ( POTENTIAL );
        unsigned int node_n_offset   = regions[i]->ebm_variable_offset ( ELECTRON );
        unsigned int node_p_offset   = regions[i]->ebm_variable_offset ( HOLE );
        unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset ( TEMPERATURE );
        unsigned int node_Tn_offset  = regions[i]->ebm_variable_offset ( E_TEMP );
        unsigned int node_Tp_offset  = regions[i]->ebm_variable_offset ( H_TEMP );

        V0 = x[fvm_nodes[i]->local_offset() +node_psi_offset];
        n0 = x[fvm_nodes[i]->local_offset() +node_n_offset]; // electron density
        p0 = x[fvm_nodes[i]->local_offset() +node_p_offset]; // hole density

        if ( regions[i]->get_advanced_model()->enable_Tl() )
          T0 = x[fvm_nodes[i]->local_offset() +node_Tl_offset];
        else
          T0 = T_external();

        if ( regions[i]->get_advanced_model()->enable_Tn() )
          Tn0 = x[fvm_nodes[i]->local_offset() + node_Tn_offset]/n0;
        else
          Tn0 = T0;

        if ( regions[i]->get_advanced_model()->enable_Tp() )
          Tp0 = x[fvm_nodes[i]->local_offset() + node_Tp_offset]/p0;
        else
          Tp0 = T0;

        mt0 = semi_region->material();
        mt0->mapping ( fvm_nodes[i]->root_node(), n0_data, SolverSpecify::clock );

        Ec0 =  - ( e*V0 + n0_data->affinity() + mt0->band->EgNarrowToEc ( p0, n0, T0 ) + kb*T0*log ( n0_data->Nc() ) );
        Ev0 =  - ( e*V0 + n0_data->affinity() - mt0->band->EgNarrowToEv ( p0, n0, T0 ) - kb*T0*log ( n0_data->Nv() ) + mt0->band->Eg ( T0 ) );
        if ( semi_region->get_advanced_model()->Fermi )
        {
          Ec0 = Ec0 - kb*T0*log ( gamma_f ( fabs ( n0 ) /n0_data->Nc() ) );
          Ev0 = Ev0 + kb*T0*log ( gamma_f ( fabs ( p0 ) /n0_data->Nv() ) );
        }
      }

      // other region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {

              unsigned int node_psi_offset = regions[i]->ebm_variable_offset ( POTENTIAL );
              unsigned int node_n_offset   = regions[i]->ebm_variable_offset ( ELECTRON );
              unsigned int node_p_offset   = regions[i]->ebm_variable_offset ( HOLE );
              unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset ( TEMPERATURE );
              unsigned int node_Tn_offset  = regions[i]->ebm_variable_offset ( E_TEMP );
              unsigned int node_Tp_offset  = regions[i]->ebm_variable_offset ( H_TEMP );

              const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *> ( regions[i] );
              const FVM_NodeData * n_data = fvm_nodes[i]->node_data();

              PetscScalar V = x[fvm_nodes[i]->local_offset() +node_psi_offset]; // psi of this node
              PetscScalar n = x[fvm_nodes[i]->local_offset() +node_n_offset]; // electron density
              PetscScalar p = x[fvm_nodes[i]->local_offset() +node_p_offset]; // hole density
              PetscScalar T = T_external();
              PetscScalar Tn = T_external();
              PetscScalar Tp = T_external();

              if ( regions[i]->get_advanced_model()->enable_Tl() )
                T = x[fvm_nodes[i]->local_offset() + node_Tl_offset];

              if ( regions[i]->get_advanced_model()->enable_Tn() )
                Tn = x[fvm_nodes[i]->local_offset() + node_Tn_offset]/n;

              if ( regions[i]->get_advanced_model()->enable_Tp() )
                Tp = x[fvm_nodes[i]->local_offset() + node_Tp_offset]/p;

              // mapping this node to material library
              Material::MaterialSemiconductor *mt = semi_region->material();
              mt->mapping ( fvm_nodes[i]->root_node(), n_data, SolverSpecify::clock );
              PetscScalar Ec =  - ( e*V + n_data->affinity() + mt->band->EgNarrowToEc ( p, n, T ) + kb*T*log ( n_data->Nc() ) );
              PetscScalar Ev =  - ( e*V + n_data->affinity() - mt->band->EgNarrowToEv ( p, n, T ) - kb*T*log ( n_data->Nv() ) + mt->band->Eg ( T ) );
              if ( semi_region->get_advanced_model()->Fermi )
              {
                Ec = Ec - kb*T*log ( gamma_f ( fabs ( n ) /n_data->Nc() ) );
                Ev = Ev + kb*T*log ( gamma_f ( fabs ( p ) /n_data->Nv() ) );
              }


              // the psi of this node is equal to corresponding node psi in the first semiconductor region
              PetscScalar ff1 = V - V0;
              iy.push_back ( fvm_nodes[i]->global_offset() +node_psi_offset );
              y_new.push_back ( ff1 );


              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = fvm_nodes[i]->outside_boundary_surface_area();


              if ( Ec0 > Ec )
              {
                // electrons are leaving region 0
                PetscScalar pm = mt0->band->EffecElecMass ( T ) /mt->band->EffecElecMass ( T );
                PetscScalar Jn = 2*( mt0->band->ThermalVn ( T ) *n0 - pm*mt->band->ThermalVn ( T ) *n*exp ( - ( Ec0-Ec ) / ( kb*T ) ) ) *cv_boundary;
                VecSetValue ( f, fvm_nodes[i]->global_offset() +node_n_offset,  Jn, ADD_VALUES );
                VecSetValue ( f, fvm_nodes[0]->global_offset() +node_n_offset, -Jn, ADD_VALUES );

                // electron energy flux
                PetscScalar Sn = -2*e* ( mt0->band->ThermalVn ( T ) *n0*2.5*kb*Tn0 - pm*mt->band->ThermalVn ( T ) *n*exp ( - ( Ec0-Ec ) / ( kb*T ) ) *2.5*kb*Tn ) *cv_boundary;
                if ( regions[i]->get_advanced_model()->enable_Tn() )
                  VecSetValue ( f, fvm_nodes[i]->global_offset() +node_Tn_offset,  -Sn, ADD_VALUES );
                if ( regions[0]->get_advanced_model()->enable_Tn() )
                  VecSetValue ( f, fvm_nodes[0]->global_offset() +regions[0]->ebm_variable_offset ( E_TEMP ), Sn, ADD_VALUES );
              }
              else
              {
                // electrons are leaving this region
                PetscScalar pm = mt->band->EffecElecMass ( T ) /mt0->band->EffecElecMass ( T );
                PetscScalar Jn = 2*( mt->band->ThermalVn ( T ) *n - pm*mt0->band->ThermalVn ( T ) *n0*exp ( - ( Ec-Ec0 ) / ( kb*T ) ) ) *cv_boundary;
                VecSetValue ( f, fvm_nodes[i]->global_offset() +node_n_offset, -Jn, ADD_VALUES );
                VecSetValue ( f, fvm_nodes[0]->global_offset() +node_n_offset,  Jn, ADD_VALUES );

                // electron energy flux
                PetscScalar Sn = -2*e* ( mt->band->ThermalVn ( T ) *n*2.5*kb*Tn - pm*mt0->band->ThermalVn ( T ) *n0*exp ( - ( Ec-Ec0 ) / ( kb*T ) ) *2.5*kb*Tn0 ) *cv_boundary;
                if ( regions[i]->get_advanced_model()->enable_Tn() )
                  VecSetValue ( f, fvm_nodes[i]->global_offset() +node_Tn_offset,  Sn, ADD_VALUES );
                if ( regions[0]->get_advanced_model()->enable_Tn() )
                  VecSetValue ( f, fvm_nodes[0]->global_offset() +regions[0]->ebm_variable_offset ( E_TEMP ), -Sn, ADD_VALUES );
              }

              if ( Ev0 < Ev )
              {
                // holes are leaving region 0
                PetscScalar pm = mt0->band->EffecHoleMass ( T ) /mt->band->EffecHoleMass ( T );
                PetscScalar Jp = 2*( mt0->band->ThermalVp ( T ) *p0 - pm*mt->band->ThermalVp ( T ) *p*exp ( ( Ev0-Ev ) / ( kb*T ) ) ) *cv_boundary;
                VecSetValue ( f, fvm_nodes[i]->global_offset() +node_p_offset,  Jp, ADD_VALUES );
                VecSetValue ( f, fvm_nodes[0]->global_offset() +node_p_offset, -Jp, ADD_VALUES );

                // electron energy flux
                PetscScalar Sp =  2*e* ( mt0->band->ThermalVp ( T ) *p0*2.5*kb*Tp0 - pm*mt->band->ThermalVp ( T ) *p*exp ( ( Ev0-Ev ) / ( kb*T ) ) *2.5*kb*Tp ) *cv_boundary;
                if ( regions[i]->get_advanced_model()->enable_Tp() )
                  VecSetValue ( f, fvm_nodes[i]->global_offset() +node_Tp_offset,  -Sp, ADD_VALUES );
                if ( regions[0]->get_advanced_model()->enable_Tp() )
                  VecSetValue ( f, fvm_nodes[0]->global_offset() +regions[0]->ebm_variable_offset ( H_TEMP ), Sp, ADD_VALUES );
              }
              else
              {
                // holes are leaving this region
                PetscScalar pm = mt->band->EffecHoleMass ( T ) /mt0->band->EffecHoleMass ( T );
                PetscScalar Jp = 2*( mt->band->ThermalVp ( T ) *p - pm*mt0->band->ThermalVp ( T ) *p0*exp ( ( Ev-Ev0 ) / ( kb*T ) ) ) *cv_boundary;
                VecSetValue ( f, fvm_nodes[i]->global_offset() +node_p_offset,  -Jp, ADD_VALUES );
                VecSetValue ( f, fvm_nodes[0]->global_offset() +node_p_offset,   Jp, ADD_VALUES );

                // electron energy flux
                PetscScalar Sp =  2*e* ( mt->band->ThermalVp ( T ) *p*2.5*kb*Tp - pm*mt0->band->ThermalVp ( T ) *p0*exp ( ( Ev-Ev0 ) / ( kb*T ) ) *2.5*kb*Tp0 ) *cv_boundary;
                if ( regions[i]->get_advanced_model()->enable_Tp() )
                  VecSetValue ( f, fvm_nodes[i]->global_offset() +node_Tp_offset,  Sp, ADD_VALUES );
                if ( regions[0]->get_advanced_model()->enable_Tp() )
                  VecSetValue ( f, fvm_nodes[0]->global_offset() +regions[0]->ebm_variable_offset ( H_TEMP ), -Sp, ADD_VALUES );
              }


              if ( regions[i]->get_advanced_model()->enable_Tl() )
              {
                PetscScalar ff4 = T - T0;
                iy.push_back ( fvm_nodes[i]->global_offset() +node_Tl_offset );
                y_new.push_back ( ff4 );
              }

              break;
            }
            case InsulatorRegion:
            {
              unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
              unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

              // the governing equation of this fvm node

              // psi of this node
              PetscScalar V = x[fvm_nodes[i]->local_offset()+0];

              // since the region is sorted, we know region[0] is semiconductor region
              // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
              //genius_assert( regions[0]->type()==SemiconductorRegion );
              PetscScalar V_semi = x[fvm_nodes[0]->local_offset()+0];

              // the psi of this node is equal to corresponding psi of semiconductor node
              // since psi should be continuous for the interface
              PetscScalar ff1 = V - V_semi;

              // find the position ff will be add to
              iy.push_back(fvm_nodes[i]->global_offset()+node_psi_offset);
              y_new.push_back(ff1);


              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                // T of this node
                PetscScalar T = x[fvm_nodes[i]->local_offset()+node_Tl_offset];

                // T for corresponding semiconductor region
                PetscScalar T_semi = x[fvm_nodes[0]->local_offset()+regions[0]->ebm_variable_offset(TEMPERATURE)];

                // the T of this node is equal to corresponding T of semiconductor node
                // by assuming no heat resistance between 2 region
                PetscScalar ff2 = T - T_semi;

                iy.push_back(fvm_nodes[i]->global_offset()+node_Tl_offset);
                y_new.push_back(ff2);
              }
              break;
            }
            default: genius_error();
        }
      }
    }

  }

  // insert new value
  if ( iy.size() )
    VecSetValues ( f, iy.size(), & ( iy[0] ), & ( y_new[0] ), ADD_VALUES );

  add_value_flag = ADD_VALUES;
}




/*---------------------------------------------------------------------
 * reserve non zero pattern in jacobian matrix for EBM3 solver
 */
void HeteroInterfaceBC::EBM3_Jacobian_Reserve ( Mat *jac, InsertMode &add_value_flag )
{

  // ADD 0 to some position of Jacobian matrix to prevent MatAssembly expurgation these position.

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();

  for(; node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    // buffer for saving regions and fvm_nodes this *node_it involves
    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these fvm_nodes have the same location in geometry,
    // but belong to different regions in logic.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0)
      {
        genius_assert( region->type() == SemiconductorRegion );
        // do nothing.
      }

      // other semiconductor region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              // reserve items for all the ghost nodes
              std::vector<int> rows, cols;
              rows.push_back(fvm_nodes[0]->global_offset()+0);
              rows.push_back(fvm_nodes[0]->global_offset()+1);
              rows.push_back(fvm_nodes[0]->global_offset()+2);

              cols.push_back(fvm_nodes[i]->global_offset()+0);
              cols.push_back(fvm_nodes[i]->global_offset()+1);
              cols.push_back(fvm_nodes[i]->global_offset()+2);

              if ( region->get_advanced_model()->enable_Tl() )
              {
                rows.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
                cols.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              }

              if(regions[0]->get_advanced_model()->enable_Tn())
                rows.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(E_TEMP));

              if(regions[0]->get_advanced_model()->enable_Tp())
                rows.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(H_TEMP));

              FVM_Node::fvm_neighbor_node_iterator  nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                cols.push_back((*nb_it).second->global_offset()+0);
                cols.push_back((*nb_it).second->global_offset()+1);
                cols.push_back((*nb_it).second->global_offset()+2);
                if ( regions[i]->get_advanced_model()->enable_Tl() )
                  cols.push_back((*nb_it).second->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
                if ( regions[i]->get_advanced_model()->enable_Tn() )
                  cols.push_back((*nb_it).second->global_offset()+regions[i]->ebm_variable_offset(E_TEMP));
                if ( regions[i]->get_advanced_model()->enable_Tp() )
                  cols.push_back((*nb_it).second->global_offset()+regions[i]->ebm_variable_offset(H_TEMP));
              }

              std::vector<PetscScalar> value(rows.size()*cols.size(),0);

              MatSetValues(*jac, rows.size(), &rows[0], cols.size(), &cols[0], &value[0], ADD_VALUES);

              // reserve for later operator
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, fvm_nodes[0]->global_offset()+0, 0, ADD_VALUES);
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+1, fvm_nodes[0]->global_offset()+1, 0, ADD_VALUES);
              MatSetValue(*jac, fvm_nodes[i]->global_offset()+2, fvm_nodes[0]->global_offset()+2, 0, ADD_VALUES);
              if ( regions[i]->get_advanced_model()->enable_Tl() )
                MatSetValue(*jac, fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE),
                            fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE), 0, ADD_VALUES);
              if(regions[i]->get_advanced_model()->enable_Tn() && regions[0]->get_advanced_model()->enable_Tn())
                MatSetValue(*jac, fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(E_TEMP),
                            fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(E_TEMP), 0, ADD_VALUES);
              if(regions[i]->get_advanced_model()->enable_Tn() && regions[0]->get_advanced_model()->enable_Tl())
                MatSetValue(*jac, fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(E_TEMP),
                            fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE), 0, ADD_VALUES);
              if(regions[i]->get_advanced_model()->enable_Tp() && regions[0]->get_advanced_model()->enable_Tp())
                MatSetValue(*jac, fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(H_TEMP),
                            fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(H_TEMP), 0, ADD_VALUES);
              if(regions[i]->get_advanced_model()->enable_Tp() && regions[0]->get_advanced_model()->enable_Tl())
                MatSetValue(*jac, fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(H_TEMP),
                            fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE), 0, ADD_VALUES);
              break;
            }
            case InsulatorRegion:
            {
              // reserve items for all the ghost nodes
              std::vector<int> rows, cols;
              rows.push_back(fvm_nodes[0]->global_offset()+0);
              cols.push_back(fvm_nodes[i]->global_offset()+0);

              if ( region->get_advanced_model()->enable_Tl() )
              {
                rows.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
                cols.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              }

              FVM_Node::fvm_neighbor_node_iterator  nb_it = fvm_nodes[i]->neighbor_node_begin();
              for(; nb_it != fvm_nodes[i]->neighbor_node_end(); ++nb_it)
              {
                cols.push_back((*nb_it).second->global_offset()+0);
                if ( region->get_advanced_model()->enable_Tl() )
                  cols.push_back((*nb_it).second->global_offset()+1);
              }

              std::vector<PetscScalar> value(rows.size()*cols.size(),0);

              MatSetValues(*jac, rows.size(), &rows[0], cols.size(), &cols[0], &value[0], ADD_VALUES);

              MatSetValue(*jac, fvm_nodes[i]->global_offset()+0, fvm_nodes[0]->global_offset()+0, 0, ADD_VALUES);
              if ( regions[i]->get_advanced_model()->enable_Tl() )
                MatSetValue(*jac, fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE),
                            fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE), 0, ADD_VALUES);
              break;
            }
            default: genius_error();
        }
      }
    }

  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




/*---------------------------------------------------------------------
 * do pre-process to jacobian matrix for EBM3 solver
 */
void HeteroInterfaceBC::EBM3_Jacobian_Preprocess(PetscScalar * ,Mat *jac, std::vector<PetscInt> &src_row,
    std::vector<PetscInt> &dst_row, std::vector<PetscInt> &clear_row)
{

  // search for all the node with this boundary type
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for(; node_it!=end_it; ++node_it )
  {

    // skip node not belongs to this processor
    if( (*node_it)->processor_id()!=Genius::processor_id() ) continue;

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin(*node_it);
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end(*node_it);
    for(unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it  )
    {
      const SimulationRegion * region = (*rnode_it).second.first;
      const FVM_Node * fvm_node = (*rnode_it).second.second;
      if(!fvm_node->is_valid()) continue;

      regions.push_back( region );
      fvm_nodes.push_back( fvm_node );

      // the first semiconductor region
      if(i==0)
      {
        // do nothing.
        // however, we will add fvm integral of other regions to it.
      }

      // other  region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(POTENTIAL));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(POTENTIAL));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(POTENTIAL));

              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(ELECTRON));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(ELECTRON));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(ELECTRON));

              src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(HOLE));
              dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(HOLE));
              clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(HOLE));

              // when have the lattice temperature equation, we add dst to src
              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                // lattice temperature equation should be global
                genius_assert(regions[0]->get_advanced_model()->enable_Tl());
                src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
                dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              }

              // when both 2 regions have the electron/hole temperature equation, we add dst to src
              if(regions[i]->get_advanced_model()->enable_Tn())
              {
                if(regions[0]->get_advanced_model()->enable_Tn())
                {
                  src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(E_TEMP));
                  dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(E_TEMP));
                }
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(E_TEMP));
              }

              if(regions[i]->get_advanced_model()->enable_Tp())
              {
                if(regions[0]->get_advanced_model()->enable_Tp())
                {
                  src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(H_TEMP));
                  dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(H_TEMP));
                }
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(H_TEMP));
              }

              break;
            }

            case InsulatorRegion:
            {
              // record the source row and dst row
              src_row.push_back(fvm_nodes[i]->global_offset()+0);
              dst_row.push_back(fvm_nodes[0]->global_offset()+0);
              clear_row.push_back(fvm_nodes[i]->global_offset()+0);

              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                src_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
                dst_row.push_back(fvm_nodes[0]->global_offset()+regions[0]->ebm_variable_offset(TEMPERATURE));
                clear_row.push_back(fvm_nodes[i]->global_offset()+regions[i]->ebm_variable_offset(TEMPERATURE));
              }
              break;
            }
            default: genius_error();
        }
      }
    }

  }


}



/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void HeteroInterfaceBC::EBM3_Jacobian ( PetscScalar * x, Mat *jac, InsertMode &add_value_flag )
{

  // since we will use ADD_VALUES operat, check the matrix state.
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  //the indepedent variable number, we need max 12 here.
  adtl::AutoDScalar::numdir=12;

  // after that, set values to source rows
  BoundaryCondition::const_node_iterator node_it = nodes_begin();
  BoundaryCondition::const_node_iterator end_it = nodes_end();
  for ( node_it = nodes_begin(); node_it!=end_it; ++node_it )
  {
    // skip node not belongs to this processor
    if ( ( *node_it )->processor_id() !=Genius::processor_id() ) continue;

    std::vector<const SimulationRegion *> regions;
    std::vector<const FVM_Node *> fvm_nodes;

    // the variable for first region
    AutoDScalar V0;
    AutoDScalar n0;
    AutoDScalar p0;
    AutoDScalar T0;
    AutoDScalar Tn0;
    AutoDScalar Tp0;
    Material::MaterialSemiconductor * mt0;
    AutoDScalar Ec0,Ev0;

    // search all the fvm_node which has *node_it as root node, these nodes are the same in geometry,
    // but in different region.
    BoundaryCondition::region_node_iterator  rnode_it     = region_node_begin ( *node_it );
    BoundaryCondition::region_node_iterator  end_rnode_it = region_node_end ( *node_it );

    for ( unsigned int i=0 ; rnode_it!=end_rnode_it; ++i, ++rnode_it )
    {
      const SimulationRegion * region = ( *rnode_it ).second.first;
      const FVM_Node * fvm_node = ( *rnode_it ).second.second;
      if ( !fvm_node->is_valid() ) continue;

      regions.push_back ( region );
      fvm_nodes.push_back ( fvm_node );

      // the first semiconductor region
      if ( i==0 )
      {
        const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *> ( regions[i] );

        const FVM_NodeData * n0_data = fvm_nodes[i]->node_data();

        unsigned int node_psi_offset = regions[i]->ebm_variable_offset ( POTENTIAL );
        unsigned int node_n_offset   = regions[i]->ebm_variable_offset ( ELECTRON );
        unsigned int node_p_offset   = regions[i]->ebm_variable_offset ( HOLE );
        unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset ( TEMPERATURE );
        unsigned int node_Tn_offset  = regions[i]->ebm_variable_offset ( E_TEMP );
        unsigned int node_Tp_offset  = regions[i]->ebm_variable_offset ( H_TEMP );

        V0 = x[fvm_nodes[i]->local_offset() +node_psi_offset];
        V0.setADValue ( node_psi_offset,1.0 );
        n0 = x[fvm_nodes[i]->local_offset() +node_n_offset];
        n0.setADValue ( node_n_offset,1.0 );  // electron density
        p0 = x[fvm_nodes[i]->local_offset() +node_p_offset];
        p0.setADValue ( node_p_offset,1.0 );  // hole density

        if ( regions[i]->get_advanced_model()->enable_Tl() )
        {
          T0 = x[fvm_nodes[i]->local_offset() + node_Tl_offset];
          T0.setADValue ( node_Tl_offset,1.0 );
        }
        else
          T0 = T_external();

        if ( regions[i]->get_advanced_model()->enable_Tn() )
        {
          AutoDScalar nTn0 = x[fvm_nodes[i]->local_offset() + node_Tn_offset];
          nTn0.setADValue ( node_Tn_offset, 1.0 );
          Tn0 = nTn0/n0;
        }
        else
          Tn0 = T0;

        if ( regions[i]->get_advanced_model()->enable_Tp() )
        {
          AutoDScalar pTp0 = x[fvm_nodes[i]->local_offset() + node_Tp_offset];
          pTp0.setADValue ( node_Tp_offset, 1.0 );
          Tp0 = pTp0/p0;
        }
        else
          Tp0 = T0;

        mt0 = semi_region->material();
        mt0->set_ad_num(adtl::AutoDScalar::numdir);
        mt0->mapping ( fvm_nodes[i]->root_node(), n0_data, SolverSpecify::clock );

        Ec0 =  - ( e*V0 + n0_data->affinity() + mt0->band->EgNarrowToEc ( p0, n0, T0 ) + kb*T0*log ( n0_data->Nc() ) );
        Ev0 =  - ( e*V0 + n0_data->affinity() - mt0->band->EgNarrowToEv ( p0, n0, T0 ) - kb*T0*log ( n0_data->Nv() ) + mt0->band->Eg ( T0 ) );
        if ( semi_region->get_advanced_model()->Fermi )
        {
          Ec0 = Ec0 - kb*T0*log ( gamma_f ( fabs ( n0 ) /n0_data->Nc() ) );
          Ev0 = Ev0 + kb*T0*log ( gamma_f ( fabs ( p0 ) /n0_data->Nv() ) );
        }
      }

      // other region
      else
      {
        switch( region->type() )
        {
            case SemiconductorRegion :
            {

              const SemiconductorSimulationRegion * semi_region = dynamic_cast<const SemiconductorSimulationRegion *> ( regions[i] );
              const FVM_NodeData * n_data = fvm_nodes[i]->node_data();

              unsigned int n_node_var_0    = regions[0]->ebm_n_variables();
              unsigned int n_node_var      = regions[i]->ebm_n_variables();
              unsigned int node_psi_offset = regions[i]->ebm_variable_offset ( POTENTIAL );
              unsigned int node_n_offset   = regions[i]->ebm_variable_offset ( ELECTRON );
              unsigned int node_p_offset   = regions[i]->ebm_variable_offset ( HOLE );
              unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset ( TEMPERATURE );
              unsigned int node_Tn_offset  = regions[i]->ebm_variable_offset ( E_TEMP );
              unsigned int node_Tp_offset  = regions[i]->ebm_variable_offset ( H_TEMP );

              std::vector<int> rows, cols;
              for ( unsigned int nv=0; nv<n_node_var_0; ++nv )
                rows.push_back ( fvm_nodes[0]->global_offset() +nv );
              for ( unsigned int nv=0; nv<n_node_var; ++nv )
                rows.push_back ( fvm_nodes[i]->global_offset() +nv );
              cols = rows;

              AutoDScalar V = x[fvm_nodes[i]->local_offset() +node_psi_offset];
              V.setADValue ( n_node_var_0+node_psi_offset, 1.0 ); // psi of this node
              AutoDScalar n = x[fvm_nodes[i]->local_offset() +node_n_offset];
              n.setADValue ( n_node_var_0+node_n_offset, 1.0 );   // electron density
              AutoDScalar p = x[fvm_nodes[i]->local_offset() +node_p_offset];
              p.setADValue ( n_node_var_0+node_p_offset, 1.0 );   // hole density

              AutoDScalar T  =  T_external();
              AutoDScalar Tn =  T_external();
              AutoDScalar Tp =  T_external();

              // lattice temperature if required
              if ( regions[i]->get_advanced_model()->enable_Tl() )
              {
                T =  x[fvm_nodes[i]->local_offset() + node_Tl_offset];
                T.setADValue ( n_node_var_0+node_Tl_offset, 1.0 );
              }

              // electron temperature if required
              if ( regions[i]->get_advanced_model()->enable_Tn() )
              {
                AutoDScalar nTn = x[fvm_nodes[i]->local_offset() + node_Tn_offset];
                nTn.setADValue ( n_node_var_0+node_Tn_offset, 1.0 );
                Tn = nTn/n;
              }

              // hole temperature if required
              if ( regions[i]->get_advanced_model()->enable_Tp() )
              {
                AutoDScalar pTp = x[fvm_nodes[i]->local_offset() + node_Tp_offset];
                pTp.setADValue ( n_node_var_0+node_Tp_offset, 1.0 );
                Tp = pTp/p;
              }

              // mapping this node to material library
              Material::MaterialSemiconductor *mt = semi_region->material();
              mt->set_ad_num(adtl::AutoDScalar::numdir);
              mt->mapping ( fvm_nodes[i]->root_node(), n_data, SolverSpecify::clock );
              AutoDScalar Ec =  - ( e*V + n_data->affinity() + mt->band->EgNarrowToEc ( p, n, T ) + kb*T*log ( n_data->Nc() ) );
              AutoDScalar Ev =  - ( e*V + n_data->affinity() - mt->band->EgNarrowToEv ( p, n, T ) - kb*T*log ( n_data->Nv() ) + mt->band->Eg ( T ) );
              if ( semi_region->get_advanced_model()->Fermi )
              {
                Ec = Ec - kb*T*log ( gamma_f ( fabs ( n ) /n_data->Nc() ) );
                Ev = Ev + kb*T*log ( gamma_f ( fabs ( p ) /n_data->Nv() ) );
              }

              // the solution value of this node is equal to corresponding node value in the first semiconductor region
              AutoDScalar ff1 = V - V0;
              MatSetValues ( *jac, 1, &rows[n_node_var_0+node_psi_offset], cols.size(), &cols[0], ff1.getADValue(), ADD_VALUES );

              if ( regions[i]->get_advanced_model()->enable_Tl() )
              {
                AutoDScalar ff4 = T - T0;
                MatSetValues ( *jac, 1, &rows[n_node_var_0+node_Tl_offset], cols.size(), &cols[0], ff4.getADValue(), ADD_VALUES );
              }


              // area of out surface of control volume related with neighbor node
              PetscScalar cv_boundary = fvm_nodes[i]->outside_boundary_surface_area();

              // thermal emit current
              if ( Ec0 > Ec )
              {
                // electrons are leaving region 0
                AutoDScalar pm = mt0->band->EffecElecMass ( T ) /mt->band->EffecElecMass ( T );
                AutoDScalar Jn = 2*( mt0->band->ThermalVn ( T ) *n0 - pm*mt->band->ThermalVn ( T ) *n*exp ( - ( Ec0-Ec ) / ( kb*T ) ) ) *cv_boundary;
                MatSetValues ( *jac, 1, &rows[n_node_var_0+node_n_offset], cols.size(), &cols[0], Jn.getADValue(), ADD_VALUES );
                MatSetValues ( *jac, 1, &rows[node_n_offset], cols.size(), &cols[0], (-Jn).getADValue(), ADD_VALUES );

                AutoDScalar Sn = -2*e* ( mt0->band->ThermalVn ( T ) *n0*2.5*kb*Tn0 - pm*mt->band->ThermalVn ( T ) *n*exp ( - ( Ec0-Ec ) / ( kb*T ) ) *2.5*kb*Tn ) *cv_boundary;
                // electron energy flux
                if ( regions[i]->get_advanced_model()->enable_Tn() )
                  MatSetValues ( *jac, 1, &rows[n_node_var_0+node_Tn_offset], cols.size(), &cols[0], ( -Sn ).getADValue(), ADD_VALUES );
                if ( regions[0]->get_advanced_model()->enable_Tn() )
                  MatSetValues ( *jac, 1, &rows[node_Tn_offset], cols.size(), &cols[0], Sn.getADValue(), ADD_VALUES );
              }
              else
              {
                // electrons are leaving this region
                AutoDScalar pm = mt->band->EffecElecMass ( T ) /mt0->band->EffecElecMass ( T );
                AutoDScalar Jn = 2*( mt->band->ThermalVn ( T ) *n - pm*mt0->band->ThermalVn ( T ) *n0*exp ( - ( Ec-Ec0 ) / ( kb*T ) ) ) *cv_boundary;
                MatSetValues ( *jac, 1, &rows[n_node_var_0+node_n_offset], cols.size(), &cols[0], (-Jn).getADValue(), ADD_VALUES );
                MatSetValues ( *jac, 1, &rows[node_n_offset], cols.size(), &cols[0], Jn.getADValue(), ADD_VALUES );

                AutoDScalar Sn = -2*e* ( mt->band->ThermalVn ( T ) *n*2.5*kb*Tn - pm*mt0->band->ThermalVn ( T ) *n0*exp ( - ( Ec-Ec0 ) / ( kb*T ) ) *2.5*kb*Tn0 ) *cv_boundary;
                // electron energy flux
                if ( regions[i]->get_advanced_model()->enable_Tn() )
                  MatSetValues ( *jac, 1, &rows[n_node_var_0+node_Tn_offset], cols.size(), &cols[0], Sn.getADValue(), ADD_VALUES );
                if ( regions[0]->get_advanced_model()->enable_Tn() )
                  MatSetValues ( *jac, 1, &rows[node_Tn_offset], cols.size(), &cols[0], ( -Sn ).getADValue(), ADD_VALUES );
              }


              if ( Ev0 < Ev )
              {
                // holes are leaving region 0
                AutoDScalar pm = mt0->band->EffecHoleMass ( T ) /mt->band->EffecHoleMass ( T );
                AutoDScalar Jp = 2*( mt0->band->ThermalVp ( T ) *p0 - pm*mt->band->ThermalVp ( T ) *p*exp ( ( Ev0-Ev ) / ( kb*T ) ) )*cv_boundary;
                MatSetValues ( *jac, 1, &rows[n_node_var_0+node_p_offset], cols.size(), &cols[0], Jp.getADValue(), ADD_VALUES );
                MatSetValues ( *jac, 1, &rows[node_p_offset], cols.size(), &cols[0], ( -Jp ).getADValue(), ADD_VALUES );

                // hole energy flux
                AutoDScalar Sp = 2*e* ( mt0->band->ThermalVp ( T ) *p0*2.5*kb*Tp0 - pm*mt->band->ThermalVp ( T ) *p*exp ( ( Ev0-Ev ) / ( kb*T ) ) *2.5*kb*Tp )*cv_boundary;
                if ( regions[i]->get_advanced_model()->enable_Tp() )
                  MatSetValues ( *jac, 1, &rows[n_node_var_0+node_Tp_offset], cols.size(), &cols[0], ( -Sp ).getADValue(), ADD_VALUES );
                if ( regions[0]->get_advanced_model()->enable_Tp() )
                  MatSetValues ( *jac, 1, &rows[node_Tp_offset], cols.size(), &cols[0], ( Sp ).getADValue(), ADD_VALUES );
              }
              else
              {
                // holes are leaving this region
                AutoDScalar pm = mt->band->EffecHoleMass ( T ) /mt0->band->EffecHoleMass ( T );
                AutoDScalar Jp = 2*( mt->band->ThermalVp ( T ) *p - pm*mt0->band->ThermalVp ( T ) *p0*exp ( ( Ev-Ev0 ) / ( kb*T ) ) )*cv_boundary;
                MatSetValues ( *jac, 1, &rows[n_node_var_0+node_p_offset], cols.size(), &cols[0], (-Jp).getADValue(), ADD_VALUES );
                MatSetValues ( *jac, 1, &rows[node_p_offset], cols.size(), &cols[0], Jp.getADValue(), ADD_VALUES );

                // hole energy flux
                AutoDScalar Sp = 2*e* ( mt->band->ThermalVp ( T ) *p*2.5*kb*Tp - pm*mt0->band->ThermalVp ( T ) *p0*exp ( ( Ev-Ev0 ) / ( kb*T ) ) *2.5*kb*Tp0 )*cv_boundary;
                if ( regions[i]->get_advanced_model()->enable_Tp() )
                  MatSetValues ( *jac, 1, &rows[n_node_var_0+node_Tp_offset], cols.size(), &cols[0], Sp.getADValue(), ADD_VALUES );
                if ( regions[0]->get_advanced_model()->enable_Tp() )
                  MatSetValues ( *jac, 1, &rows[node_Tp_offset], cols.size(), &cols[0], ( -Sp ).getADValue(), ADD_VALUES );
              }

              break;
            }

            case InsulatorRegion:
            {
              unsigned int global_offset   = fvm_nodes[i]->global_offset();
              unsigned int node_psi_offset = regions[i]->ebm_variable_offset(POTENTIAL);
              unsigned int node_Tl_offset  = regions[i]->ebm_variable_offset(TEMPERATURE);

              unsigned int semiconductor_node_psi_offset = regions[0]->ebm_variable_offset(POTENTIAL);
              unsigned int semiconductor_node_Tl_offset  = regions[0]->ebm_variable_offset(TEMPERATURE);


              // psi of this node
              AutoDScalar  V = x[fvm_nodes[i]->local_offset()+node_psi_offset]; V.setADValue(0,1.0);

              // since the region is sorted, we know region[0] is semiconductor region
              // as a result, x[fvm_nodes[0]->local_offset()] is psi for corresponding semiconductor region
              //genius_assert( regions[0]->type()==SemiconductorRegion );
              AutoDScalar  V_semi = x[fvm_nodes[0]->local_offset()+semiconductor_node_psi_offset]; V_semi.setADValue(1,1.0);

              // the psi of this node is equal to corresponding psi of semiconductor node
              AutoDScalar  ff1 = V - V_semi;

              // set Jacobian of governing equation ff
              MatSetValue(*jac, global_offset+node_psi_offset, fvm_nodes[i]->global_offset()+node_psi_offset, ff1.getADValue(0), ADD_VALUES);
              MatSetValue(*jac, global_offset+node_psi_offset, fvm_nodes[0]->global_offset()+semiconductor_node_psi_offset, ff1.getADValue(1), ADD_VALUES);

              if(regions[i]->get_advanced_model()->enable_Tl())
              {
                // T of this node
                AutoDScalar  T = x[fvm_nodes[i]->local_offset()+node_Tl_offset]; T.setADValue(0,1.0);

                // T for corresponding semiconductor region
                AutoDScalar  T_semi = x[fvm_nodes[0]->local_offset()+semiconductor_node_Tl_offset]; T_semi.setADValue(1,1.0);

                // the T of this node is equal to corresponding T of semiconductor node
                // we assuming T is continuous for the interface
                AutoDScalar ff2 = T - T_semi;

                // set Jacobian of governing equation ff2
                MatSetValue(*jac, global_offset+node_Tl_offset, fvm_nodes[i]->global_offset()+node_Tl_offset, ff2.getADValue(0), ADD_VALUES);
                MatSetValue(*jac, global_offset+node_Tl_offset, fvm_nodes[0]->global_offset()+semiconductor_node_Tl_offset, ff2.getADValue(1), ADD_VALUES);
              }

              break;
            }
            default: genius_error();
        }
      }
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;
}
