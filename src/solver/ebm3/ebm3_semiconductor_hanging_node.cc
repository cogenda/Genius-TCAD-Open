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

#include "elem.h"
#include "simulation_system.h"
#include "semiconductor_region.h"
#include "jflux1.h"
#include "jflux2.h"
#include "jflux3.h"
#include "petsc_utils.h"


using PhysicalUnit::kb;
using PhysicalUnit::e;


void SemiconductorSimulationRegion::EBM3_Function_Hanging_Node(PetscScalar *x, Vec f, InsertMode &add_value_flag)
{
  if( !has_2d_hanging_node() && !has_3d_hanging_node()  ) return;

  // find the node variable offset
  unsigned int node_psi_offset = ebm_variable_offset(POTENTIAL);
  unsigned int node_n_offset   = ebm_variable_offset(ELECTRON);
  unsigned int node_p_offset   = ebm_variable_offset(HOLE);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);
  unsigned int node_Tn_offset  = ebm_variable_offset(E_TEMP);
  unsigned int node_Tp_offset  = ebm_variable_offset(H_TEMP);

  // process hanging node lies on side center
  {
    // buffer for record src and dst rows
    std::vector<PetscInt>    src_row;
    std::vector<PetscInt>    dst_row;
    std::vector<PetscScalar> alpha_buffer;

    // buffer for record position and value for insert
    std::vector<PetscInt>    insert_index;
    std::vector<PetscScalar> insert_buffer;

    hanging_node_on_elem_side_iterator  hanging_node_it = hanging_node_on_elem_side_begin();
    hanging_node_on_elem_side_iterator  hanging_node_it_end = hanging_node_on_elem_side_end();

    for(; hanging_node_it!=hanging_node_it_end; ++hanging_node_it )
    {
      const FVM_Node * fvm_node = (*hanging_node_it).first;

      // skip node not belongs to this processor
      if( fvm_node->root_node()->processor_id()!=Genius::processor_id() ) continue;

      // let the flux of hanging node flow into other "non-hanging" node on the element side averagely
      // this process ensure the global conservation of flux

      const Elem * elem = (*hanging_node_it).second.first;
      unsigned int side_index = (*hanging_node_it).second.second;

      AutoPtr<Elem>  side = elem->build_side(side_index);

      unsigned int n_side_node = side->n_nodes();
      std::vector<const FVM_Node *> side_fvm_nodes;

      for(unsigned int n=0; n < n_side_node; n++)
      {
        const Node * side_node = side->get_node(n);
        const FVM_Node * side_fvm_node = region_fvm_node(side_node);
        side_fvm_nodes.push_back(side_fvm_node);

        src_row.push_back(fvm_node->global_offset()+node_psi_offset);
        src_row.push_back(fvm_node->global_offset()+node_n_offset);
        src_row.push_back(fvm_node->global_offset()+node_p_offset);

        dst_row.push_back(side_fvm_node->global_offset()+node_psi_offset);
        dst_row.push_back(side_fvm_node->global_offset()+node_n_offset);
        dst_row.push_back(side_fvm_node->global_offset()+node_p_offset);

        alpha_buffer.push_back(1.0/n_side_node);
        alpha_buffer.push_back(1.0/n_side_node);
        alpha_buffer.push_back(1.0/n_side_node);

        // for extra Tl temperature equation
        if(get_advanced_model()->enable_Tl())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tl_offset);
          dst_row.push_back(side_fvm_node->global_offset()+node_Tl_offset);
          alpha_buffer.push_back(1.0/n_side_node);
        }

        // for extra Tn temperature equation
        if(get_advanced_model()->enable_Tn())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tn_offset);
          dst_row.push_back(side_fvm_node->global_offset()+node_Tn_offset);
          alpha_buffer.push_back(1.0/n_side_node);
        }

        // for extra Tp temperature equation
        if(get_advanced_model()->enable_Tp())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tp_offset);
          dst_row.push_back(side_fvm_node->global_offset()+node_Tp_offset);
          alpha_buffer.push_back(1.0/n_side_node);
        }
      }

      // and then, the value of hanging node is interpolated.

      const FVM_Node * interpolation_p1;
      const FVM_Node * interpolation_p2;

      // find the interpolation point
      if (side_fvm_nodes.size() == 2)
      {
        // for 2D case, the side should be an edge
        interpolation_p1 = side_fvm_nodes[0];
        interpolation_p2 = side_fvm_nodes[1];
      }
      else if (side_fvm_nodes.size() == 4)
      {
        // for 3D case, the side should be QUAD4, we use 2 point which has little psi difference as interpolation point
        PetscScalar dv1 = std::abs( x[side_fvm_nodes[0]->local_offset()] - x[side_fvm_nodes[2]->local_offset()] );
        PetscScalar dv2 = std::abs( x[side_fvm_nodes[1]->local_offset()] - x[side_fvm_nodes[3]->local_offset()] );

        if( dv1 < dv2 )
        {
          interpolation_p1 = side_fvm_nodes[0];
          interpolation_p2 = side_fvm_nodes[2];
        }
        else
        {
          interpolation_p1 = side_fvm_nodes[1];
          interpolation_p2 = side_fvm_nodes[3];
        }
      }
      else
      {
        // we should never reach here
        genius_error();
      }


      PetscScalar V = x[fvm_node->local_offset()+node_psi_offset];
      PetscScalar n = x[fvm_node->local_offset()+node_n_offset];
      PetscScalar p = x[fvm_node->local_offset()+node_p_offset];

      PetscScalar V1 = x[interpolation_p1->local_offset()+node_psi_offset];
      PetscScalar n1 = x[interpolation_p1->local_offset()+node_n_offset];
      PetscScalar p1 = x[interpolation_p1->local_offset()+node_p_offset];

      PetscScalar V2 = x[interpolation_p2->local_offset()+node_psi_offset];
      PetscScalar n2 = x[interpolation_p2->local_offset()+node_n_offset];
      PetscScalar p2 = x[interpolation_p2->local_offset()+node_p_offset];

      PetscScalar Tl = fvm_node->node_data()->T();

      // the psi is linear interpolated
      insert_index.push_back(fvm_node->global_offset()+node_psi_offset);
      insert_buffer.push_back( V - 0.5*(V1+V2) );

      if( get_advanced_model()->Jn_level() == 1 )
      {
        // the electron density is interpolated by basic S-G scheme
        insert_index.push_back(fvm_node->global_offset()+node_n_offset);
        //insert_buffer.push_back( n - 0.5*(n1+n2) );
        insert_buffer.push_back( n - nmid_dd(kb*Tl/e, V1, V2, n1, n2) );
      }

      if( get_advanced_model()->Jp_level() == 1 )
      {
        // the hole density is interpolated by basic S-G scheme
        insert_index.push_back(fvm_node->global_offset()+node_p_offset);
        //insert_buffer.push_back( p - 0.5*(p1+p2) );
        insert_buffer.push_back( p - pmid_dd(kb*Tl/e, V1, V2, p1, p2) );
      }


      if(get_advanced_model()->enable_Tl())
      {
        // the lattice temperature is linear interpolated
        PetscScalar T  = x[fvm_node->local_offset()+node_Tl_offset];
        PetscScalar T1 = x[interpolation_p1->local_offset()+node_Tl_offset];
        PetscScalar T2 = x[interpolation_p2->local_offset()+node_Tl_offset];
        insert_index.push_back(fvm_node->global_offset()+node_Tl_offset);
        insert_buffer.push_back( T - 0.5*(T1+T2) );

        if( get_advanced_model()->Jn_level() == 2 )
        {
          // the electron density is interpolated by lattice temperature corrected S-G scheme
          insert_index.push_back(fvm_node->global_offset()+node_n_offset);
          //insert_buffer.push_back( n - 0.5*(n1+n2) );
          insert_buffer.push_back( n - nmid_lt(kb, e, V2-V1, n1, n2, 0.5*(T1+T2), T2-T1) );
        }

        if( get_advanced_model()->Jp_level() == 2 )
        {
          // the hole density is interpolated by lattice temperature corrected S-G scheme
          insert_index.push_back(fvm_node->global_offset()+node_p_offset);
          //insert_buffer.push_back( p - 0.5*(p1+p2) );
          insert_buffer.push_back( p - pmid_lt(kb, e, V2-V1, p1, p2, 0.5*(T1+T2), T2-T1) );
        }
      }

      // for extra Tn temperature equation
      if(get_advanced_model()->enable_Tn())
      {
        // the electron temperature is linear interpolated
        PetscScalar Tn = x[fvm_node->local_offset()+node_Tn_offset]/n;
        PetscScalar Tn1 = x[interpolation_p1->local_offset()+node_Tn_offset]/n1;
        PetscScalar Tn2 = x[interpolation_p2->local_offset()+node_Tn_offset]/n2;
        insert_index.push_back(fvm_node->global_offset()+node_Tn_offset);
        insert_buffer.push_back( Tn - 0.5*(Tn1+Tn2) );

        if( get_advanced_model()->Jn_level() == 3 )
        {
          // the electron density is interpolated by EB Level S-G scheme
          insert_index.push_back(fvm_node->global_offset()+node_n_offset);
          //insert_buffer.push_back( n - 0.5*(n1+n2) );
          insert_buffer.push_back( n - nmid_eb(kb, e, V1, V2, n1, n2, Tn1, Tn2) );
        }
      }

      // for extra Tp temperature equation
      if(get_advanced_model()->enable_Tp())
      {
        // the hole temperature is linear interpolated
        PetscScalar Tp = x[fvm_node->local_offset()+node_Tp_offset]/p;
        PetscScalar Tp1 = x[interpolation_p1->local_offset()+node_Tp_offset]/p1;
        PetscScalar Tp2 = x[interpolation_p2->local_offset()+node_Tp_offset]/p2;
        insert_index.push_back(fvm_node->global_offset()+node_Tp_offset);
        insert_buffer.push_back( Tp - 0.5*(Tp1+Tp2) );

        if( get_advanced_model()->Jp_level() == 3 )
        {
          // the hole density is interpolated by EB Level S-G scheme
          insert_index.push_back(fvm_node->global_offset()+node_p_offset);
          //insert_buffer.push_back( p - 0.5*(p1+p2) );
          insert_buffer.push_back( p - pmid_eb(kb, e, V1, V2, p1, p2, Tp1, Tp2) );
        }
      }

    }

    PetscUtils::VecAddRowToRow(f, src_row, dst_row, alpha_buffer);

    // do INSERT_VALUES to Vec
    if(insert_index.size())
      VecSetValues(f, insert_index.size(), &insert_index[0], &insert_buffer[0], INSERT_VALUES);

  }


  // process hanging node lies on edge center
  {
    // buffer for record src and dst rows
    std::vector<PetscInt>    src_row;
    std::vector<PetscInt>    dst_row;
    std::vector<PetscScalar> alpha_buffer;

    // buffer for record position and value for insert
    std::vector<PetscInt>    insert_index;
    std::vector<PetscScalar> insert_buffer;

    hanging_node_on_elem_edge_iterator  hanging_node_it = hanging_node_on_elem_edge_begin();
    hanging_node_on_elem_edge_iterator  hanging_node_it_end = hanging_node_on_elem_edge_end();

    for(; hanging_node_it!=hanging_node_it_end; ++hanging_node_it )
    {
      const FVM_Node * fvm_node = (*hanging_node_it).first;

      // skip node not belongs to this processor
      if( fvm_node->root_node()->processor_id()!=Genius::processor_id() ) continue;

      // let the flux of hanging node flow into other "non-hanging" node on the element side averagely
      // this process ensure the global conservation of flux

      const Elem * elem = (*hanging_node_it).second.first;
      unsigned int edge_index = (*hanging_node_it).second.second;

      AutoPtr<Elem>  edge = elem->build_edge(edge_index);

      unsigned int n_edge_node = 2;
      std::vector<const FVM_Node *> edge_fvm_nodes;

      for(unsigned int n=0; n < n_edge_node; n++)
      {
        const Node * edge_node = edge->get_node(n);
        const FVM_Node * edge_fvm_node = region_fvm_node(edge_node);
        genius_assert(edge_fvm_node!=NULL);

        edge_fvm_nodes.push_back(edge_fvm_node);

        src_row.push_back(fvm_node->global_offset()+node_psi_offset);
        src_row.push_back(fvm_node->global_offset()+node_n_offset);
        src_row.push_back(fvm_node->global_offset()+node_p_offset);

        dst_row.push_back(edge_fvm_node->global_offset()+node_psi_offset);
        dst_row.push_back(edge_fvm_node->global_offset()+node_n_offset);
        dst_row.push_back(edge_fvm_node->global_offset()+node_p_offset);

        alpha_buffer.push_back(1.0/n_edge_node);
        alpha_buffer.push_back(1.0/n_edge_node);
        alpha_buffer.push_back(1.0/n_edge_node);

        // for extra Tl temperature equation
        if(get_advanced_model()->enable_Tl())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tl_offset);
          dst_row.push_back(edge_fvm_node->global_offset()+node_Tl_offset);
          alpha_buffer.push_back(1.0/n_edge_node);
        }

        // for extra Tn temperature equation
        if(get_advanced_model()->enable_Tn())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tn_offset);
          dst_row.push_back(edge_fvm_node->global_offset()+node_Tn_offset);
          alpha_buffer.push_back(1.0/n_edge_node);
        }

        // for extra Tp temperature equation
        if(get_advanced_model()->enable_Tp())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tp_offset);
          dst_row.push_back(edge_fvm_node->global_offset()+node_Tp_offset);
          alpha_buffer.push_back(1.0/n_edge_node);
        }
      }

      // and then, the value of hanging node is interpolated.
      PetscScalar V = x[fvm_node->local_offset()+node_psi_offset];
      PetscScalar n = x[fvm_node->local_offset()+node_n_offset];
      PetscScalar p = x[fvm_node->local_offset()+node_p_offset];

      PetscScalar V1 = x[edge_fvm_nodes[0]->local_offset()+node_psi_offset];
      PetscScalar n1 = x[edge_fvm_nodes[0]->local_offset()+node_n_offset];
      PetscScalar p1 = x[edge_fvm_nodes[0]->local_offset()+node_p_offset];

      PetscScalar V2 = x[edge_fvm_nodes[1]->local_offset()+node_psi_offset];
      PetscScalar n2 = x[edge_fvm_nodes[1]->local_offset()+node_n_offset];
      PetscScalar p2 = x[edge_fvm_nodes[1]->local_offset()+node_p_offset];

      PetscScalar Tl = fvm_node->node_data()->T();

      // the psi is linear interpolated
      insert_index.push_back(fvm_node->global_offset()+node_psi_offset);
      insert_buffer.push_back( V - 0.5*(V1+V2) );

      if( get_advanced_model()->Jn_level() == 1 )
      {
        // the electron density is interpolated by basic S-G scheme
        insert_index.push_back(fvm_node->global_offset()+node_n_offset);
        //insert_buffer.push_back( n - 0.5*(n1+n2) );
        insert_buffer.push_back( n - nmid_dd(kb*Tl/e, V1, V2, n1, n2) );
      }

      if( get_advanced_model()->Jp_level() == 1 )
      {
        // the hole density is interpolated by basic S-G scheme
        insert_index.push_back(fvm_node->global_offset()+node_p_offset);
        //insert_buffer.push_back( p - 0.5*(p1+p2) );
        insert_buffer.push_back( p - pmid_dd(kb*Tl/e, V1, V2, p1, p2) );
      }

      if(get_advanced_model()->enable_Tl())
      {
        // the lattice temperature is linear interpolated
        PetscScalar T  = x[fvm_node->local_offset()+node_Tl_offset];
        PetscScalar T1 = x[edge_fvm_nodes[0]->local_offset()+node_Tl_offset];
        PetscScalar T2 = x[edge_fvm_nodes[1]->local_offset()+node_Tl_offset];
        insert_index.push_back(fvm_node->global_offset()+node_Tl_offset);
        insert_buffer.push_back( T - 0.5*(T1+T2) );

        if( get_advanced_model()->Jn_level() == 2 )
        {
          // the electron density is interpolated by lattice temperature corrected S-G scheme
          insert_index.push_back(fvm_node->global_offset()+node_n_offset);
          //insert_buffer.push_back( n - 0.5*(n1+n2) );
          insert_buffer.push_back( n - nmid_lt(kb, e, V2-V1, n1, n2, 0.5*(T1+T2), T2-T1) );
        }

        if( get_advanced_model()->Jp_level() == 2 )
        {
          // the hole density is interpolated by lattice temperature corrected S-G scheme
          insert_index.push_back(fvm_node->global_offset()+node_p_offset);
          //insert_buffer.push_back( p - 0.5*(p1+p2) );
          insert_buffer.push_back( p - pmid_lt(kb, e, V2-V1, p1, p2, 0.5*(T1+T2), T2-T1) );
        }
      }

      // for extra Tn temperature equation
      if(get_advanced_model()->enable_Tn())
      {
        // the electron temperature is linear interpolated
        PetscScalar Tn = x[fvm_node->local_offset()+node_Tn_offset]/n;
        PetscScalar Tn1 = x[edge_fvm_nodes[0]->local_offset()+node_Tn_offset]/n1;
        PetscScalar Tn2 = x[edge_fvm_nodes[1]->local_offset()+node_Tn_offset]/n2;
        insert_index.push_back(fvm_node->global_offset()+node_Tn_offset);
        insert_buffer.push_back( Tn - 0.5*(Tn1+Tn2) );

        if( get_advanced_model()->Jn_level() == 3 )
        {
          // the electron density is interpolated by EB level S-G scheme
          insert_index.push_back(fvm_node->global_offset()+node_n_offset);
          //insert_buffer.push_back( n - 0.5*(n1+n2) );
          insert_buffer.push_back( n - nmid_eb(kb, e, V1, V2, n1, n2, Tn1, Tn2) );
        }
      }

      // for extra Tp temperature equation
      if(get_advanced_model()->enable_Tp())
      {
        // the hole temperature is linear interpolated
        PetscScalar Tp = x[fvm_node->local_offset()+node_Tp_offset]/p;
        PetscScalar Tp1 = x[edge_fvm_nodes[0]->local_offset()+node_Tp_offset]/p1;
        PetscScalar Tp2 = x[edge_fvm_nodes[1]->local_offset()+node_Tp_offset]/p2;
        insert_index.push_back(fvm_node->global_offset()+node_Tp_offset);
        insert_buffer.push_back( Tp - 0.5*(Tp1+Tp2) );

        if( get_advanced_model()->Jp_level() == 3 )
        {
          // the hole density is interpolated by EB level S-G scheme
          insert_index.push_back(fvm_node->global_offset()+node_p_offset);
          //insert_buffer.push_back( p - 0.5*(p1+p2) );
          insert_buffer.push_back( p - pmid_eb(kb, e, V1, V2, p1, p2, Tp1, Tp2) );
        }
      }

    }

    PetscUtils::VecAddRowToRow(f, src_row, dst_row, alpha_buffer);

    // do INSERT_VALUES to Vec
    if(insert_index.size())
      VecSetValues(f, insert_index.size(), &insert_index[0], &insert_buffer[0], INSERT_VALUES);
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  add_value_flag = INSERT_VALUES;



}





void SemiconductorSimulationRegion::EBM3_Jacobian_Hanging_Node(PetscScalar *x, Mat *jac, InsertMode &add_value_flag)
{

  if( !has_2d_hanging_node() && !has_3d_hanging_node()  ) return;

  // find the node variable offset
  unsigned int n_node_var      = ebm_n_variables();
  unsigned int node_psi_offset = ebm_variable_offset(POTENTIAL);
  unsigned int node_n_offset   = ebm_variable_offset(ELECTRON);
  unsigned int node_p_offset   = ebm_variable_offset(HOLE);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);
  unsigned int node_Tn_offset  = ebm_variable_offset(E_TEMP);
  unsigned int node_Tp_offset  = ebm_variable_offset(H_TEMP);

  // process hanging node lies on side center
  {

    // buffer for add matrix row to row
    std::vector<PetscInt>    src_row;
    std::vector<PetscInt>    dst_row;
    std::vector<PetscScalar> alpha_buffer;

    // buffer for matrix entrance
    std::vector< PetscInt >                 row_index;
    std::vector< std::vector<PetscInt> >    cols_index;
    std::vector< AutoDScalar >              ad_values;

    hanging_node_on_elem_side_iterator  hanging_node_it = hanging_node_on_elem_side_begin();
    hanging_node_on_elem_side_iterator  hanging_node_it_end = hanging_node_on_elem_side_end();

    for(; hanging_node_it!=hanging_node_it_end; ++hanging_node_it )
    {
      const FVM_Node * fvm_node = (*hanging_node_it).first;

      // skip node not belongs to this processor
      if( fvm_node->root_node()->processor_id()!=Genius::processor_id() ) continue;

      // let the flux of hanging node flow into other "non-hanging" node on the element side averagely
      // this process ensure the global conservation of flux

      const Elem * elem = (*hanging_node_it).second.first;
      unsigned int side_index = (*hanging_node_it).second.second;

      AutoPtr<Elem>  side = elem->build_side(side_index);

      unsigned int n_side_node = side->n_nodes();
      std::vector<const FVM_Node *> side_fvm_nodes;

      for(unsigned int n=0; n < n_side_node; n++)
      {
        const Node * side_node = side->get_node(n);
        const FVM_Node * side_fvm_node = region_fvm_node(side_node);
        side_fvm_nodes.push_back(side_fvm_node);

        src_row.push_back(fvm_node->global_offset()+node_psi_offset);
        src_row.push_back(fvm_node->global_offset()+node_n_offset);
        src_row.push_back(fvm_node->global_offset()+node_p_offset);

        dst_row.push_back(side_fvm_node->global_offset()+node_psi_offset);
        dst_row.push_back(side_fvm_node->global_offset()+node_n_offset);
        dst_row.push_back(side_fvm_node->global_offset()+node_p_offset);

        alpha_buffer.push_back(1.0/n_side_node);
        alpha_buffer.push_back(1.0/n_side_node);
        alpha_buffer.push_back(1.0/n_side_node);

        // for extra Tl temperature equation
        if(get_advanced_model()->enable_Tl())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tl_offset);
          dst_row.push_back(side_fvm_node->global_offset()+node_Tl_offset);
          alpha_buffer.push_back(1.0/n_side_node);
        }

        // for extra Tn temperature equation
        if(get_advanced_model()->enable_Tn())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tn_offset);
          dst_row.push_back(side_fvm_node->global_offset()+node_Tn_offset);
          alpha_buffer.push_back(1.0/n_side_node);
        }

        // for extra Tp temperature equation
        if(get_advanced_model()->enable_Tp())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tp_offset);
          dst_row.push_back(side_fvm_node->global_offset()+node_Tp_offset);
          alpha_buffer.push_back(1.0/n_side_node);
        }

      }

      // and then, the value of hanging node is interpolated.
      //the indepedent variable number, we need 3*n_node_var here
      adtl::AutoDScalar::numdir = 3*n_node_var;

      const FVM_Node * interpolation_p1;
      const FVM_Node * interpolation_p2;

      // find the interpolation point
      if (side_fvm_nodes.size() == 2)
      {
        // for 2D case, the side should be an edge
        interpolation_p1 = side_fvm_nodes[0];
        interpolation_p2 = side_fvm_nodes[1];
      }
      else if (side_fvm_nodes.size() == 4)
      {
        // for 3D case, the side should be QUAD4, we use 2 point which has little psi difference as interpolation point
        PetscScalar dv1 = std::abs( x[side_fvm_nodes[0]->local_offset()] - x[side_fvm_nodes[2]->local_offset()] );
        PetscScalar dv2 = std::abs( x[side_fvm_nodes[1]->local_offset()] - x[side_fvm_nodes[3]->local_offset()] );

        if( dv1 < dv2 )
        {
          interpolation_p1 = side_fvm_nodes[0];
          interpolation_p2 = side_fvm_nodes[2];
        }
        else
        {
          interpolation_p1 = side_fvm_nodes[1];
          interpolation_p2 = side_fvm_nodes[3];
        }
      }
      else
      {
        // we should never reach here
        genius_error();
      }

      AutoDScalar V = x[fvm_node->local_offset()+node_psi_offset]; V.setADValue(0*n_node_var + node_psi_offset, 1.0);
      AutoDScalar n = x[fvm_node->local_offset()+node_n_offset];   n.setADValue(0*n_node_var + node_n_offset, 1.0);
      AutoDScalar p = x[fvm_node->local_offset()+node_p_offset];   p.setADValue(0*n_node_var + node_p_offset, 1.0);

      AutoDScalar V1 = x[interpolation_p1->local_offset()+node_psi_offset]; V1.setADValue(1*n_node_var + node_psi_offset, 1.0);
      AutoDScalar n1 = x[interpolation_p1->local_offset()+node_n_offset];   n1.setADValue(1*n_node_var + node_n_offset, 1.0);
      AutoDScalar p1 = x[interpolation_p1->local_offset()+node_p_offset];   p1.setADValue(1*n_node_var + node_p_offset, 1.0);

      AutoDScalar V2 = x[interpolation_p2->local_offset()+node_psi_offset]; V2.setADValue(2*n_node_var + node_psi_offset, 1.0);
      AutoDScalar n2 = x[interpolation_p2->local_offset()+node_n_offset];   n2.setADValue(2*n_node_var + node_n_offset, 1.0);
      AutoDScalar p2 = x[interpolation_p2->local_offset()+node_p_offset];   p2.setADValue(2*n_node_var + node_p_offset, 1.0);

      PetscScalar Tl = fvm_node->node_data()->T();

      std::vector<PetscInt> cols(3*n_node_var);
      for(unsigned int n=0; n<n_node_var; ++n)
      {
        cols[0*n_node_var + n] = fvm_node->global_offset()+n;
        cols[1*n_node_var + n] = interpolation_p1->global_offset()+n;
        cols[2*n_node_var + n] = interpolation_p2->global_offset()+n;
      }

      // the psi is linear interpolated
      AutoDScalar ff1 = V - 0.5*(V1+V2);
      row_index.push_back(fvm_node->global_offset()+node_psi_offset);
      cols_index.push_back(cols);
      ad_values.push_back(ff1);

      if( get_advanced_model()->Jn_level() == 1 )
      {
        // the electron density is interpolated by basic S-G scheme
        AutoDScalar ff2 = n - nmid_dd(kb*Tl/e, V1, V2, n1, n2);
        row_index.push_back(fvm_node->global_offset()+node_n_offset);
        cols_index.push_back(cols);
        ad_values.push_back(ff2);
      }

      if( get_advanced_model()->Jp_level() == 1 )
      {
        // the hole density is interpolated by basic S-G scheme
        AutoDScalar ff3 =  p - pmid_dd(kb*Tl/e, V1, V2, p1, p2);
        row_index.push_back(fvm_node->global_offset()+node_p_offset);
        cols_index.push_back(cols);
        ad_values.push_back(ff3);
      }

      if(get_advanced_model()->enable_Tl())
      {
        // the lattice temperature is linear interpolated
        AutoDScalar T = x[fvm_node->local_offset()+node_Tl_offset];           T.setADValue (0*n_node_var + node_Tl_offset, 1.0);
        AutoDScalar T1 = x[interpolation_p1->local_offset()+node_Tl_offset];  T1.setADValue(1*n_node_var + node_Tl_offset, 1.0);
        AutoDScalar T2 = x[interpolation_p2->local_offset()+node_Tl_offset];  T2.setADValue(2*n_node_var + node_Tl_offset, 1.0);

        AutoDScalar ff4 = T - 0.5*(T1+T2);
        row_index.push_back(fvm_node->global_offset()+node_Tl_offset);
        cols_index.push_back(cols);
        ad_values.push_back(ff4);

        if( get_advanced_model()->Jn_level() == 2 )
        {
          // the electron density is interpolated by lattice temperature corrected S-G scheme
          AutoDScalar ff2 = n - nmid_lt(kb, e, V2-V1, n1, n2, 0.5*(T1+T2), T2-T1);
          row_index.push_back(fvm_node->global_offset()+node_n_offset);
          cols_index.push_back(cols);
          ad_values.push_back(ff2);
        }

        if( get_advanced_model()->Jp_level() == 2 )
        {
          // the hole density is interpolated by lattice temperature corrected S-G scheme
          AutoDScalar ff3 =  p - pmid_lt(kb, e, V2-V1, p1, p2, 0.5*(T1+T2), T2-T1);
          row_index.push_back(fvm_node->global_offset()+node_p_offset);
          cols_index.push_back(cols);
          ad_values.push_back(ff3);
        }
      }

      // for extra Tn temperature equation
      if(get_advanced_model()->enable_Tn())
      {
        // the electron temperature is linear interpolated
        AutoDScalar nTn = x[fvm_node->local_offset()+node_Tn_offset];           nTn.setADValue (0*n_node_var + node_Tn_offset, 1.0);
        AutoDScalar nTn1 = x[interpolation_p1->local_offset()+node_Tn_offset];  nTn1.setADValue(1*n_node_var + node_Tn_offset, 1.0);
        AutoDScalar nTn2 = x[interpolation_p2->local_offset()+node_Tn_offset];  nTn2.setADValue(2*n_node_var + node_Tn_offset, 1.0);

        AutoDScalar Tn  = nTn/n;
	AutoDScalar Tn1 = nTn1/n1;
	AutoDScalar Tn2 = nTn2/n2;

	AutoDScalar ff5 = Tn - 0.5*(Tn1+Tn2);
        row_index.push_back(fvm_node->global_offset()+node_Tn_offset);
        cols_index.push_back(cols);
        ad_values.push_back(ff5);

        if( get_advanced_model()->Jn_level() == 3 )
        {
          // the electron density is interpolated by EB Level S-G scheme
          AutoDScalar ff2 = n - nmid_eb(kb, e, V1, V2, n1, n2, Tn1, Tn2);
          row_index.push_back(fvm_node->global_offset()+node_n_offset);
          cols_index.push_back(cols);
          ad_values.push_back(ff2);
        }
      }

      // for extra Tp temperature equation
      if(get_advanced_model()->enable_Tp())
      {
        // the hole temperature is linear interpolated
        AutoDScalar pTp = x[fvm_node->local_offset()+node_Tp_offset];           pTp.setADValue (0*n_node_var + node_Tp_offset, 1.0);
        AutoDScalar pTp1 = x[interpolation_p1->local_offset()+node_Tp_offset];  pTp1.setADValue (1*n_node_var + node_Tp_offset, 1.0);
        AutoDScalar pTp2 = x[interpolation_p2->local_offset()+node_Tp_offset];  pTp2.setADValue (2*n_node_var + node_Tp_offset, 1.0);

        AutoDScalar Tp  = pTp/p;
	AutoDScalar Tp1 = pTp1/p1;
	AutoDScalar Tp2 = pTp2/p2;

	AutoDScalar ff6 = Tp - 0.5*(Tp1+Tp2);
        row_index.push_back(fvm_node->global_offset()+node_Tp_offset);
        cols_index.push_back(cols);
        ad_values.push_back(ff6);

        if( get_advanced_model()->Jp_level() == 3 )
        {
          // the hole density is interpolated by EB Level S-G scheme
          AutoDScalar ff3 = p - pmid_eb(kb, e, V1, V2, p1, p2, Tp1, Tp2);
          row_index.push_back(fvm_node->global_offset()+node_p_offset);
          cols_index.push_back(cols);
          ad_values.push_back(ff3);
        }
      }

    }


    //ok, we add source rows to destination rows
    PetscUtils::MatAddRowToRow(*jac, src_row, dst_row, alpha_buffer);

    // clear row_index
    PetscUtils::MatZeroRows(*jac, row_index.size(), row_index.empty() ? NULL : &row_index[0], 0.0);


    // insert buffered AD values to Mat
    for(unsigned int n=0; n< row_index.size(); ++n )
      MatSetValues(*jac, 1, &row_index[n], (cols_index[n]).size(), &((cols_index[n])[0]), ad_values[n].getADValue(), INSERT_VALUES);

  }


  // process hanging node lies on edge center
  {

    // buffer for add matrix row to row
    std::vector<PetscInt>    src_row;
    std::vector<PetscInt>    dst_row;
    std::vector<PetscScalar> alpha_buffer;

    // buffer for matrix entrance
    std::vector< PetscInt >                 row_index;
    std::vector< std::vector<PetscInt> >    cols_index;
    std::vector< AutoDScalar >              ad_values;

    hanging_node_on_elem_edge_iterator  hanging_node_it = hanging_node_on_elem_edge_begin();
    hanging_node_on_elem_edge_iterator  hanging_node_it_end = hanging_node_on_elem_edge_end();

    for(; hanging_node_it!=hanging_node_it_end; ++hanging_node_it )
    {
      const FVM_Node * fvm_node = (*hanging_node_it).first;

      // skip node not belongs to this processor
      if( fvm_node->root_node()->processor_id()!=Genius::processor_id() ) continue;

      // let the flux of hanging node flow into other "non-hanging" node on the element side averagely
      // this process ensure the global conservation of flux

      const Elem * elem = (*hanging_node_it).second.first;
      unsigned int edge_index = (*hanging_node_it).second.second;

      AutoPtr<Elem>  edge = elem->build_edge(edge_index);

      unsigned int n_edge_node = 2;
      std::vector<const FVM_Node *> edge_fvm_nodes;

      for(unsigned int n=0; n < n_edge_node; n++)
      {
        const Node * edge_node = edge->get_node(n);
        const FVM_Node * edge_fvm_node = region_fvm_node(edge_node);
        edge_fvm_nodes.push_back(edge_fvm_node);

        src_row.push_back(fvm_node->global_offset()+node_psi_offset);
        src_row.push_back(fvm_node->global_offset()+node_n_offset);
        src_row.push_back(fvm_node->global_offset()+node_p_offset);

        dst_row.push_back(edge_fvm_node->global_offset()+node_psi_offset);
        dst_row.push_back(edge_fvm_node->global_offset()+node_n_offset);
        dst_row.push_back(edge_fvm_node->global_offset()+node_p_offset);

        alpha_buffer.push_back(1.0/n_edge_node);
        alpha_buffer.push_back(1.0/n_edge_node);
        alpha_buffer.push_back(1.0/n_edge_node);

        // for extra Tl temperature equation
        if(get_advanced_model()->enable_Tl())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tl_offset);
          dst_row.push_back(edge_fvm_node->global_offset()+node_Tl_offset);
          alpha_buffer.push_back(1.0/n_edge_node);
        }

        // for extra Tn temperature equation
        if(get_advanced_model()->enable_Tn())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tn_offset);
          dst_row.push_back(edge_fvm_node->global_offset()+node_Tn_offset);
          alpha_buffer.push_back(1.0/n_edge_node);
        }

        // for extra Tp temperature equation
        if(get_advanced_model()->enable_Tp())
        {
          src_row.push_back(fvm_node->global_offset()+node_Tp_offset);
          dst_row.push_back(edge_fvm_node->global_offset()+node_Tp_offset);
          alpha_buffer.push_back(1.0/n_edge_node);
        }
      }

      // and then, the value of hanging node is interpolated.

      // the indepedent variable number, we need 3*n_node_var here
      adtl::AutoDScalar::numdir = 3*n_node_var;

      AutoDScalar V = x[fvm_node->local_offset()+node_psi_offset]; V.setADValue(0*n_node_var + node_psi_offset, 1.0);
      AutoDScalar n = x[fvm_node->local_offset()+node_n_offset];   n.setADValue(0*n_node_var + node_n_offset, 1.0);
      AutoDScalar p = x[fvm_node->local_offset()+node_p_offset];   p.setADValue(0*n_node_var + node_p_offset, 1.0);

      AutoDScalar V1 = x[edge_fvm_nodes[0]->local_offset()+node_psi_offset]; V1.setADValue(1*n_node_var + node_psi_offset, 1.0);
      AutoDScalar n1 = x[edge_fvm_nodes[0]->local_offset()+node_n_offset];   n1.setADValue(1*n_node_var + node_n_offset, 1.0);
      AutoDScalar p1 = x[edge_fvm_nodes[0]->local_offset()+node_p_offset];   p1.setADValue(1*n_node_var + node_p_offset, 1.0);

      AutoDScalar V2 = x[edge_fvm_nodes[1]->local_offset()+node_psi_offset]; V2.setADValue(2*n_node_var + node_psi_offset, 1.0);
      AutoDScalar n2 = x[edge_fvm_nodes[1]->local_offset()+node_n_offset];   n2.setADValue(2*n_node_var + node_n_offset, 1.0);
      AutoDScalar p2 = x[edge_fvm_nodes[1]->local_offset()+node_p_offset];   p2.setADValue(2*n_node_var + node_p_offset, 1.0);

      PetscScalar Tl = fvm_node->node_data()->T();

      std::vector<PetscInt> cols(3*n_node_var);
      for(unsigned int n=0; n<n_node_var; ++n)
      {
        cols[0*n_node_var + n] = fvm_node->global_offset()+n;
        cols[1*n_node_var + n] = edge_fvm_nodes[0]->global_offset()+n;
        cols[2*n_node_var + n] = edge_fvm_nodes[1]->global_offset()+n;
      }

      // the psi is linear interpolated
      AutoDScalar ff1 = V - 0.5*(V1+V2);
      row_index.push_back(fvm_node->global_offset()+node_psi_offset);
      cols_index.push_back(cols);
      ad_values.push_back(ff1);

      if( get_advanced_model()->Jn_level() == 1 )
      {
        // the electron density is interpolated by basic S-G scheme
        AutoDScalar ff2 = n - nmid_dd(kb*Tl/e, V1, V2, n1, n2);
        row_index.push_back(fvm_node->global_offset()+node_n_offset);
        cols_index.push_back(cols);
        ad_values.push_back(ff2);
      }

      if( get_advanced_model()->Jp_level() == 1 )
      {
        // the hole density is interpolated by basic S-G scheme
        AutoDScalar ff3 =  p - pmid_dd(kb*Tl/e, V1, V2, p1, p2);
        row_index.push_back(fvm_node->global_offset()+node_p_offset);
        cols_index.push_back(cols);
        ad_values.push_back(ff3);
      }

      if(get_advanced_model()->enable_Tl())
      {
        // the lattice temperature is linear interpolated
        AutoDScalar T = x[fvm_node->local_offset()+node_Tl_offset];            T.setADValue (0*n_node_var + node_Tl_offset, 1.0);
        AutoDScalar T1 = x[edge_fvm_nodes[0]->local_offset()+node_Tl_offset];  T1.setADValue(1*n_node_var + node_Tl_offset, 1.0);
        AutoDScalar T2 = x[edge_fvm_nodes[1]->local_offset()+node_Tl_offset];  T2.setADValue(2*n_node_var + node_Tl_offset, 1.0);

        AutoDScalar ff4 = T - 0.5*(T1+T2);
        row_index.push_back(fvm_node->global_offset()+node_Tl_offset);
        cols_index.push_back(cols);
        ad_values.push_back(ff4);

        if( get_advanced_model()->Jn_level() == 2 )
        {
          // the electron density is interpolated by lattice temperature corrected S-G scheme
          AutoDScalar ff2 = n - nmid_lt(kb, e, V2-V1, n1, n2, 0.5*(T1+T2), T2-T1);
          row_index.push_back(fvm_node->global_offset()+node_n_offset);
          cols_index.push_back(cols);
          ad_values.push_back(ff2);
        }

        if( get_advanced_model()->Jp_level() == 2 )
        {
          // the hole density is interpolated by lattice temperature corrected S-G scheme
          AutoDScalar ff3 =  p - pmid_lt(kb, e, V2-V1, p1, p2, 0.5*(T1+T2), T2-T1);
          row_index.push_back(fvm_node->global_offset()+node_p_offset);
          cols_index.push_back(cols);
          ad_values.push_back(ff3);
        }
      }

      // for extra Tn temperature equation
      if(get_advanced_model()->enable_Tn())
      {
        // the electron temperature is linear interpolated
        AutoDScalar nTn  = x[fvm_node->local_offset()+node_Tn_offset];           nTn.setADValue (0*n_node_var + node_Tn_offset, 1.0);
        AutoDScalar nTn1 = x[edge_fvm_nodes[0]->local_offset()+node_Tn_offset];  nTn1.setADValue(1*n_node_var + node_Tn_offset, 1.0);
        AutoDScalar nTn2 = x[edge_fvm_nodes[1]->local_offset()+node_Tn_offset];  nTn2.setADValue(2*n_node_var + node_Tn_offset, 1.0);

	AutoDScalar Tn  = nTn/n;
	AutoDScalar Tn1 = nTn1/n1;
	AutoDScalar Tn2 = nTn2/n2;

	AutoDScalar ff5 = Tn - 0.5*(Tn1+Tn2);
        row_index.push_back(fvm_node->global_offset()+node_Tn_offset);
        cols_index.push_back(cols);
        ad_values.push_back(ff5);

        if( get_advanced_model()->Jn_level() == 3 )
        {
          // the electron density is interpolated by EB Level S-G scheme
          AutoDScalar ff2 = n - nmid_eb(kb, e, V1, V2, n1, n2, Tn1, Tn2);
          row_index.push_back(fvm_node->global_offset()+node_n_offset);
          cols_index.push_back(cols);
          ad_values.push_back(ff2);
        }
      }

      // for extra Tp temperature equation
      if(get_advanced_model()->enable_Tp())
      {
        // the hole temperature is linear interpolated
        AutoDScalar pTp  = x[fvm_node->local_offset()+node_Tp_offset];           pTp.setADValue (0*n_node_var + node_Tp_offset, 1.0);
        AutoDScalar pTp1 = x[edge_fvm_nodes[0]->local_offset()+node_Tp_offset];  pTp1.setADValue (1*n_node_var + node_Tp_offset, 1.0);
        AutoDScalar pTp2 = x[edge_fvm_nodes[1]->local_offset()+node_Tp_offset];  pTp2.setADValue (2*n_node_var + node_Tp_offset, 1.0);

	AutoDScalar Tp  = pTp/p;
	AutoDScalar Tp1 = pTp1/p1;
	AutoDScalar Tp2 = pTp2/p2;

	AutoDScalar ff6 = Tp - 0.5*(Tp1+Tp2);
        row_index.push_back(fvm_node->global_offset()+node_Tp_offset);
        cols_index.push_back(cols);
        ad_values.push_back(ff6);

        if( get_advanced_model()->Jp_level() == 3 )
        {
          // the hole density is interpolated by EB Level S-G scheme
          AutoDScalar ff3 = p - pmid_eb(kb, e, V1, V2, p1, p2, Tp1, Tp2);
          row_index.push_back(fvm_node->global_offset()+node_p_offset);
          cols_index.push_back(cols);
          ad_values.push_back(ff3);
        }
      }

    }


    //ok, we add source rows to destination rows
    PetscUtils::MatAddRowToRow(*jac, src_row, dst_row, alpha_buffer);

    // clear row_index
    PetscUtils::MatZeroRows(*jac, row_index.size(), row_index.empty() ? NULL : &row_index[0], 0.0);


    // insert buffered AD values to Mat
    for(unsigned int n=0; n< row_index.size(); ++n )
      MatSetValues(*jac, 1, &row_index[n], (cols_index[n]).size(), &((cols_index[n])[0]), ad_values[n].getADValue(), INSERT_VALUES);

  }

  add_value_flag = INSERT_VALUES;

}





