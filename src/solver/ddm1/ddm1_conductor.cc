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

//  $Id: ddm1_conductor.cc,v 1.3 2008/07/09 05:58:16 gdiso Exp $

#include "elem.h"
#include "simulation_system.h"
#include "conductor_region.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;

///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

// for DDML1 solver, only poisson's equation should be considered in conductor region
// this file is the same as poisson_conductor.cc
// However, we can skip solving poisson's equation for a conductor...


/*---------------------------------------------------------------------
 * fill solution vector with value.
 */
void ConductorSimulationRegion::DDM1_Fill_Value(Vec x, Vec L)
{
  PetscInt n_local_dofs;
  VecGetLocalSize(x, &n_local_dofs);

  std::vector<int> ix;
  std::vector<PetscScalar> y;
  std::vector<PetscScalar> s;

  ix.reserve(n_local_dofs);
  y.reserve(n_local_dofs);
  s.reserve(n_local_dofs);

  const_node_iterator it = nodes_begin();
  const_node_iterator it_end = nodes_end();
  for(; it!=it_end; ++it)
  {
    const FVM_Node * node = (*it).second;
    //if this node NOT belongs to this processor, continue
    if( node->root_node()->processor_id() != Genius::processor_id() ) continue;

    const FVM_NodeData * node_data = node->node_data();

    ix.push_back(node->global_offset());
    y.push_back(node_data->psi());
    s.push_back(1.0/node->volume());
  }

  if( ix.size() )
  {
    VecSetValues(x, ix.size(), &ix[0], &y[0], INSERT_VALUES) ;
    VecSetValues(L, ix.size(), &ix[0], &s[0], INSERT_VALUES) ;
  }
}



/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void ConductorSimulationRegion::DDM1_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    VecAssemblyBegin(f);
    VecAssemblyEnd(f);
  }

  // set local buf here
  std::vector<int>          iy;
  std::vector<PetscScalar>  y;


  // search all the element in this region.
  // note, they are all local element, thus must be processed

  element_iterator it = elements_begin();
  element_iterator it_end = elements_end();
  for(; it!=it_end; ++it)
  {

    //search for all the Edge this cell own
    for(unsigned int ne=0; ne<(*it)->n_edges(); ++ne )
    {
      AutoPtr<Elem> edge = (*it)->build_edge (ne);
      genius_assert(edge->type()==EDGE2);

      std::vector<unsigned int > edge_nodes;
      (*it)->nodes_on_edge(ne, edge_nodes);

      // the length of this edge
      double length = edge->volume();

      // partial area associated with this edge
      double partial_area = (*it)->partial_area_with_edge(ne);

      // fvm_node of node1
      const FVM_Node * fvm_n1 = (*it)->get_fvm_node(edge_nodes[0]);   genius_assert(fvm_n1);
      // fvm_node of node2
      const FVM_Node * fvm_n2 = (*it)->get_fvm_node(edge_nodes[1]);   genius_assert(fvm_n2);

      // fvm_node_data of node1
      const FVM_NodeData * n1_data = fvm_n1->node_data();  genius_assert(n1_data);
      // fvm_node_data of node2
      const FVM_NodeData * n2_data = fvm_n2->node_data();  genius_assert(n2_data);

      unsigned int n1_local_offset = fvm_n1->local_offset();
      unsigned int n2_local_offset = fvm_n2->local_offset();

      // build poisson's equation here
      {

        PetscScalar V1   =  x[n1_local_offset];  // electrostatic potential of node1
        PetscScalar rho1 =  0;                   // free charge density of node1
        PetscScalar eps1 =  n1_data->eps();      // permittivity of node1


        PetscScalar V2   =  x[n2_local_offset];  // electrostatic potential of node2
        PetscScalar rho2 =  0;                   // free charge density of node2
        PetscScalar eps2 =  n2_data->eps();      // permittivity of node2

        PetscScalar eps = 0.5*(eps1+eps2);       // eps at mid point of the edge

        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {
          iy.push_back( fvm_n1->global_offset() );
          y.push_back ( eps*partial_area*(V2 - V1)/length + ( rho1 )*(*it)->partial_volume(edge_nodes[0]) );
        }

        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {
          iy.push_back( fvm_n2->global_offset() );
          y.push_back ( eps*partial_area*(V1 - V2)/length + ( rho2 )*(*it)->partial_volume(edge_nodes[1]) );
        }
      }

    }

  }

  if(iy.size()) VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);

  // after the first scan, every nodes are updated.
  // however, boundary condition should be processed later.

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}






/*---------------------------------------------------------------------
 * build function and its jacobian for DDML1 solver
 */
void ConductorSimulationRegion::DDM1_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  //the indepedent variable number, since we only process edges, 2 is enough
  adtl::AutoDScalar::numdir=2;

  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  // note, we will use ADD_VALUES to set values of matrix J
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  // search all the element in this region.
  // note, they are all local element, thus must be processed

  element_iterator it = elements_begin();
  element_iterator it_end = elements_end();
  for(; it!=it_end; ++it)
  {

    //search for all the Edge this cell own
    for(unsigned int ne=0; ne<(*it)->n_edges(); ++ne )
    {
      AutoPtr<Elem> edge = (*it)->build_edge (ne);
      genius_assert(edge->type()==EDGE2);

      std::vector<unsigned int > edge_nodes;
      (*it)->nodes_on_edge(ne, edge_nodes);

      // the length of this edge
      double length = edge->volume();

      // partial area associated with this edge
      double partial_area = (*it)->partial_area_with_edge(ne);

      // fvm_node of node1
      const FVM_Node * fvm_n1 = (*it)->get_fvm_node(edge_nodes[0]);   genius_assert(fvm_n1);
      // fvm_node of node2
      const FVM_Node * fvm_n2 = (*it)->get_fvm_node(edge_nodes[1]);   genius_assert(fvm_n2);

      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data() ;   genius_assert(n1_data);
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data() ;   genius_assert(n2_data);

      unsigned int n1_local_offset = fvm_n1->local_offset();
      unsigned int n2_local_offset = fvm_n2->local_offset();

      // the row/colume position of variables in the matrix
      std::vector<PetscInt> row,col;
      row.push_back(fvm_n1->global_offset());
      row.push_back(fvm_n2->global_offset());
      col = row;

      // here we use AD, however it is great overkill for such a simple problem.
      {

        // electrostatic potential, as independent variable
        AutoDScalar V1   =  x[n1_local_offset];   V1.setADValue(0,1.0);
        PetscScalar rho1 =  0;
        PetscScalar eps1 =  n1_data->eps();


        AutoDScalar V2   =  x[n2_local_offset];   V2.setADValue(1,1.0);
        PetscScalar rho2 =  0;
        PetscScalar eps2 =  n2_data->eps();

        PetscScalar eps = 0.5*(eps1+eps2); // eps at mid point of the edge


        // ignore thoese ghost nodes
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {
          AutoDScalar f =  eps*partial_area*(V2 - V1)/length + ( rho1 )*(*it)->partial_volume(edge_nodes[0]);
          MatSetValues(*jac, 1, &row[0], col.size(), &col[0], f.getADValue(), ADD_VALUES);
        }

        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {
          AutoDScalar f =  eps*partial_area*(V1 - V2)/length + ( rho2 )*(*it)->partial_volume(edge_nodes[1]);
          MatSetValues(*jac, 1, &row[1], col.size(), &col[0], f.getADValue(), ADD_VALUES);
        }

      }
    }


  }

  // boundary condition should be processed later!

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}





void ConductorSimulationRegion::DDM1_Update_Solution(PetscScalar *lxx)
{

  node_iterator it = nodes_begin();
  node_iterator it_end = nodes_end();
  for(; it!=it_end; ++it)
  {
    FVM_Node * fvm_node = (*it).second;

    // NOTE: here, solution for all the local node should be updated!
    if( !fvm_node->root_node()->on_local() ) continue;

    FVM_NodeData * node_data = fvm_node->node_data();  genius_assert(node_data!=NULL);

    //update psi
    node_data->psi_last() =  node_data->psi();
    node_data->psi() = lxx[fvm_node->local_offset()];
  }

  // addtional work: compute electrical field for all the node.
  // however, the electrical field is always zero. We needn't do anything here.

}




