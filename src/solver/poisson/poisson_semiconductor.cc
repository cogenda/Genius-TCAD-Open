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

//  $Id: poisson_semiconductor.cc,v 1.22 2008/07/09 05:58:16 gdiso Exp $

#include "elem.h"
#include "simulation_system.h"
#include "semiconductor_region.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;

///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


void SemiconductorSimulationRegion::Poissin_Fill_Value(Vec x, Vec L)
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
 * build function and its jacobian for poisson solver
 */
void SemiconductorSimulationRegion::Poissin_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of vec f
  // if the previous operator is not ADD_VALUES, we should assembly the vec first!
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

  //common used variable
  PetscScalar T   = T_external();
  PetscScalar nie = mt->band->ni(T);

  // process \nabla operator for all cells

  const_element_iterator it = elements_begin();
  const_element_iterator it_end = elements_end();
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

      // build nonlinear poisson's equation
      {

        PetscScalar V1   =  x[n1_local_offset]; // electrostatic potential of node1
        PetscScalar eps1 =  n1_data->eps();     // permittivity of node1

        PetscScalar V2   =  x[n2_local_offset]; // electrostatic potential of node2
        PetscScalar eps2 =  n2_data->eps();     // permittivity of node2

        PetscScalar eps = 0.5*(eps1+eps2);      // eps at mid point of the edge

        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {
          iy.push_back(fvm_n1->global_offset());
          y.push_back( eps*partial_area*(V2 - V1)/length );
        }

        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {
          iy.push_back( fvm_n2->global_offset() );
          y.push_back( eps*partial_area*(V1 - V2)/length );
        }
      }

    }

  }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // process \rho of poisson's equation
  const_node_iterator node_it = nodes_begin();
  const_node_iterator node_it_end = nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {

    const FVM_Node * fvm_node = (*node_it).second;

    // skip node not on this processor
    if( fvm_node->root_node()->processor_id() != Genius::processor_id() ) continue;

    const FVM_NodeData * fvm_node_data = fvm_node->node_data();  genius_assert(fvm_node_data);

    PetscScalar V   =  x[fvm_node->local_offset()];
    mt->mapping(fvm_node->root_node(), fvm_node_data, 0.0);

    // intrinsic Fermi potential.
    PetscScalar V_i =  V + fvm_node_data->affinity() + fvm_node_data->Eg()/(2*e) + kb*T*log(fvm_node_data->Nc()/fvm_node_data->Nv())/(2*e);
    // electron density
    PetscScalar n   =  nie*exp( e/(kb*T)*V_i);
    // hole density
    PetscScalar p   =  nie*exp(-e/(kb*T)*V_i);
    // total chrage density
    PetscScalar rho =  e*(fvm_node_data->Net_doping() + p - n )*fvm_node->volume();

    iy.push_back(fvm_node->global_offset());                                  // save index in the buffer
    y.push_back( rho );                                                       // save value in the buffer

  }

  if(iy.size()) VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);

  // after the first scan, every nodes are updated.
  // however, boundary condition should be processed later.

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}






/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void SemiconductorSimulationRegion::Poissin_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of matrix J
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  //common used variable
  PetscScalar T   = T_external();
  PetscScalar nie = mt->band->ni(T);

  // search all the element in this region.
  // note, they are all local element, thus must be processed

  //the indepedent variable is 2
  adtl::AutoDScalar::numdir = 2;
  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  const_element_iterator it = elements_begin();
  const_element_iterator it_end = elements_end();
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
        AutoDScalar V1   =  x[n1_local_offset];   V1.setADValue(0, 1.0);
        // eps
        PetscScalar eps1 =  n1_data->eps();


        // electrostatic potential, as independent variable
        AutoDScalar V2   =  x[n2_local_offset];   V2.setADValue(1, 1.0);
        // eps
        PetscScalar eps2 =  n2_data->eps();

        PetscScalar eps = 0.5*(eps1+eps2);

        // ignore thoese ghost nodes
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {
          AutoDScalar f =  eps*partial_area*(V2 - V1)/length ;
          MatSetValues(*jac, 1, &row[0], col.size(), &col[0], f.getADValue(), ADD_VALUES);
        }

        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {
          AutoDScalar f =  eps*partial_area*(V1 - V2)/length ;
          MatSetValues(*jac, 1, &row[1], col.size(), &col[0], f.getADValue(), ADD_VALUES);
        }

      }
    }


  }

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // process \rho of poisson's equation

  //the indepedent variable is 1
  adtl::AutoDScalar::numdir = 1;
  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  // process \rho of poisson's equation, search all the vertex of this cell
  const_node_iterator node_it = nodes_begin();
  const_node_iterator node_it_end = nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {

    const FVM_Node * fvm_node = (*node_it).second;

    // skip node not on this processor
    if( fvm_node->root_node()->processor_id() != Genius::processor_id() ) continue;

    const FVM_NodeData * fvm_node_data = fvm_node->node_data();  genius_assert(fvm_node_data);

    AutoDScalar V   =  x[fvm_node->local_offset()];  V.setADValue(0, 1.0);
    mt->mapping(fvm_node->root_node(), fvm_node_data, 0.0);

    // intrinsic Fermi potential.
    AutoDScalar V_i =  V + fvm_node_data->affinity() + fvm_node_data->Eg()/(2*e) + kb*T*log(fvm_node_data->Nc()/fvm_node_data->Nv())/(2*e);
    // electron density
    AutoDScalar n   =  nie*exp( e/(kb*T)*V_i);
    // hole density
    AutoDScalar p   =  nie*exp(-e/(kb*T)*V_i);

    // total chrage density
    AutoDScalar rho =  e*(fvm_node_data->Net_doping() + p - n)*fvm_node->volume();

    MatSetValue(*jac, fvm_node->global_offset(), fvm_node->global_offset(), rho.getADValue(0), ADD_VALUES);  // save value in the buffer

  }

  // boundary condition should be processed later!

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}



void SemiconductorSimulationRegion::Poissin_Update_Solution(PetscScalar *lxx)
{


  node_iterator it = nodes_begin();
  node_iterator it_end = nodes_end();
  for(; it!=it_end; ++it)
  {
    FVM_Node * fvm_node = (*it).second;

    // NOTE: here, potential for all the local node should be updated!
    if( !fvm_node->root_node()->on_local() ) continue;

    PetscScalar V   = lxx[fvm_node->local_offset()];

    FVM_NodeData * node_data = fvm_node->node_data();

    //update psi
    node_data->psi() = V;

    //for semiconductor region, also update n and p
    {
      PetscScalar T   = T_external();
      mt->mapping(fvm_node->root_node(), node_data, 0.0);
      PetscScalar ni = mt->band->ni(T);

      // intrinsic Fermi potential.
      PetscScalar V_i = V
                        + node_data->affinity()
                        + node_data->Eg()/(2*e)
                        + kb*T*log(node_data->Nc()/node_data->Nv())/(2*e);
      // electron density
      PetscScalar n    =  ni*exp(e/(kb*T)*V_i);
      node_data->n()   =  n;
      // hole density
      PetscScalar p    =  ni*exp(-e/(kb*T)*V_i);
      node_data->p()   =  p;

      node_data->Ec() = -(e*V + node_data->affinity());
      node_data->Ev() = -(e*V + node_data->affinity() + mt->band->Eg(T));
    }
  }

  // addtional work: compute electrical field for all the node.
  // Since this value is only used for reference.
  // It can be done simply by weighted average of cell's electrical field
  for(unsigned int n=0; n<n_cell(); ++n)
  {
    const Elem * elem = this->get_region_elem(n);
    FVM_CellData * elem_data = this->get_region_elem_Data(n);

    std::vector<PetscScalar> psi_vertex;

    for(unsigned int nd=0; nd<elem->n_nodes(); ++nd)
    {
      const FVM_Node * fvm_node = elem->get_fvm_node(nd);
      const FVM_NodeData * fvm_node_data = fvm_node->node_data();
      psi_vertex.push_back  ( fvm_node_data->psi() );
    }
    // compute the gradient in the cell
    elem_data->E()  = - elem->gradient(psi_vertex);  // E = - grad(psi)
  }

}




