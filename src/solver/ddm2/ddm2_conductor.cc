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
#include "conductor_region.h"
#include "solver_specify.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;

///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////

// for DDML2 solver, poisson's equation and lattice temperature equation
// should be both considered in insulator region


/*---------------------------------------------------------------------
 * fill solution vector with value.
 */
void ConductorSimulationRegion::DDM2_Fill_Value(Vec x, Vec L)
{
  PetscInt n_local_dofs;
  VecGetLocalSize(x, &n_local_dofs);

  std::vector<int> ix;
  std::vector<PetscScalar> y;
  std::vector<PetscScalar> s;

  ix.reserve(n_local_dofs);
  y.reserve(n_local_dofs);
  s.reserve(n_local_dofs);

  node_iterator it = nodes_begin();
  node_iterator it_end = nodes_end();
  for(; it!=it_end; ++it)
  {
    const FVM_Node * node = (*it).second;
    //if this node NOT belongs to this processor, continue
    if( node->root_node()->processor_id() != Genius::processor_id() ) continue;

    const FVM_NodeData * node_data = node->node_data();

    // psi
    ix.push_back(node->global_offset()+0);
    y.push_back(node_data->psi());
    s.push_back(1.0/node->volume());

    // lattice temperature
    ix.push_back(node->global_offset()+1);
    y.push_back(node_data->T());
    s.push_back(1.0/node->volume());
  }

  if( ix.size() )
  {
    VecSetValues(x, ix.size(), &ix[0], &y[0], INSERT_VALUES) ;
    VecSetValues(L, ix.size(), &ix[0], &s[0], INSERT_VALUES) ;
  }
}



/*---------------------------------------------------------------------
 * build function and its jacobian for DDML2 solver
 */
void ConductorSimulationRegion::DDM2_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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

      // build governing equations for conductor region
      {

        //for node 1 of the edge
        mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
        PetscScalar V1   =  x[n1_local_offset+0];             // electrostatic potential
        PetscScalar T1   =  x[n1_local_offset+1];             // lattice temperature
        PetscScalar rho1 =  0;                                // charge density
        PetscScalar eps1 =  n1_data->eps();                   // permittivity
        PetscScalar kap1 =  mt->thermal->HeatConduction(T1);

        //for node 2 of the edge
        mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
        PetscScalar V2   =  x[n2_local_offset+0];
        PetscScalar T2   =  x[n2_local_offset+1];
        PetscScalar rho2 =  0;
        PetscScalar eps2 =  n2_data->eps();
        PetscScalar kap2 =  mt->thermal->HeatConduction(T2);

        PetscScalar eps = 0.5*(eps1+eps2);       // eps at mid point of the edge
        PetscScalar kap = 0.5*(kap1+kap2);       // kapa at mid point of the edge

        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {
          // poisson's equation
          iy.push_back( fvm_n1->global_offset()+0 );
          y.push_back ( eps*partial_area*(V2 - V1)/length + ( rho1 )*(*it)->partial_volume(edge_nodes[0]) );

          // heat transport equation
          iy.push_back( fvm_n1->global_offset()+1 );
          y.push_back ( kap*partial_area*(T2 - T1)/length );
        }

        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {
          // poisson's equation
          iy.push_back( fvm_n2->global_offset()+0 );
          y.push_back ( eps*partial_area*(V1 - V2)/length + ( rho2 )*(*it)->partial_volume(edge_nodes[1]) );

          // heat transport equation
          iy.push_back( fvm_n2->global_offset()+1 );
          y.push_back ( kap*partial_area*(T1 - T2)/length );
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
 * build function and its jacobian for DDML2 solver
 */
void ConductorSimulationRegion::DDM2_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of matrix J
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }

  //the indepedent variable number, since we only process edges, 2 is enough
  adtl::AutoDScalar::numdir=2;

  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  // search all the element in this region.
  // note, they are all local element, thus must be processed

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

      // here we use AD, however it is great overkill for such a simple problem.
      {

        //for node 1 of the edge
        mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
        AutoDScalar V1   =  x[n1_local_offset+0];  V1.setADValue(0,1.0);           // electrostatic potential
        AutoDScalar T1   =  x[n1_local_offset+1];  T1.setADValue(0,1.0);           // lattice temperature
        PetscScalar rho1 =  0;                                // charge density
        PetscScalar eps1 =  n1_data->eps();                   // permittivity
        PetscScalar kap1 =  mt->thermal->HeatConduction(T1.getValue());

        //for node 2 of the edge
        mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
        AutoDScalar V2   =  x[n2_local_offset+0];  V2.setADValue(1,1.0);
        AutoDScalar T2   =  x[n2_local_offset+1];  T2.setADValue(1,1.0);
        PetscScalar rho2 =  0;
        PetscScalar eps2 =  n2_data->eps();
        PetscScalar kap2 =  mt->thermal->HeatConduction(T2.getValue());

        PetscScalar eps = 0.5*(eps1+eps2);       // eps at mid point of the edge
        PetscScalar kap = 0.5*(kap1+kap2);       // kapa at mid point of the edge


        // ignore thoese ghost nodes
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {
          AutoDScalar f1 =  eps*partial_area*(V2 - V1)/length + ( rho1 )*(*it)->partial_volume(edge_nodes[0]) ;
          MatSetValue(*jac, fvm_n1->global_offset(), fvm_n1->global_offset(), f1.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_n1->global_offset(), fvm_n2->global_offset(), f1.getADValue(1), ADD_VALUES);

          AutoDScalar f2 =  kap*partial_area*(T2 - T1)/length;
          MatSetValue(*jac, fvm_n1->global_offset()+1, fvm_n1->global_offset()+1, f2.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_n1->global_offset()+1, fvm_n2->global_offset()+1, f2.getADValue(1), ADD_VALUES);
        }

        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {

          AutoDScalar f1 =  eps*partial_area*(V1 - V2)/length + ( rho2 )*(*it)->partial_volume(edge_nodes[1]) ;
          MatSetValue(*jac, fvm_n2->global_offset(), fvm_n1->global_offset(), f1.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_n2->global_offset(), fvm_n2->global_offset(), f1.getADValue(1), ADD_VALUES);

          AutoDScalar f2 =  kap*partial_area*(T1 - T2)/length;
          MatSetValue(*jac, fvm_n2->global_offset()+1, fvm_n1->global_offset()+1, f2.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_n2->global_offset()+1, fvm_n2->global_offset()+1, f2.getADValue(1), ADD_VALUES);
        }
      }
    }


  }

  // boundary condition should be processed later!

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




void ConductorSimulationRegion::DDM2_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
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

  // process node related terms
  const_node_iterator node_it = nodes_begin();
  const_node_iterator node_it_end = nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {

    const FVM_Node * fvm_node = (*node_it).second;
    //if this node NOT belongs to this processor, continue
    if( fvm_node->root_node()->processor_id() != Genius::processor_id() ) continue;

    const FVM_NodeData * node_data = fvm_node->node_data();
    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

    PetscScalar T   =  x[fvm_node->local_offset()+1];                         // lattice temperature
    PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(T);

    // process \partial t

    iy.push_back(fvm_node->global_offset()+1);                                // save index in the buffer

    //second order
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false)
    {
      PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
      PetscScalar TT = -((2-r)/(1-r)*T - 1.0/(r*(1-r))*node_data->T() + (1-r)/r*node_data->T_last())*node_data->density()*HeatCapacity
                       / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
      y.push_back( TT );
    }
    else //first order
    {
      PetscScalar TT = -(T - node_data->T())*node_data->density()*HeatCapacity/SolverSpecify::dt*fvm_node->volume();
      y.push_back( TT );
    }
  }


  // add into petsc vector, we should prevent zero length vector add here.
  if(iy.size())  VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;


#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}




void ConductorSimulationRegion::DDM2_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // note, we will use ADD_VALUES to set values of matrix J
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( (add_value_flag != ADD_VALUES) && (add_value_flag != NOT_SET_VALUES) )
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


  //the indepedent variable number, 1 for each node
  adtl::AutoDScalar::numdir = 1;

  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  const_node_iterator node_it = nodes_begin();
  const_node_iterator node_it_end = nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {

    const FVM_Node * fvm_node = (*node_it).second;
    //if this node NOT belongs to this processor, continue
    if( fvm_node->root_node()->processor_id() != Genius::processor_id() ) continue;

    const FVM_NodeData * node_data = fvm_node->node_data();
    mt->mapping(fvm_node->root_node(), node_data, SolverSpecify::clock);

    AutoDScalar T   =  x[fvm_node->local_offset()+1];   T.setADValue(0, 1.0);              // lattice temperature
    PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(T.getValue());

    // process \partial t

    //second order
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false)
    {
      PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
      AutoDScalar TT = -((2-r)/(1-r)*T - 1.0/(r*(1-r))*node_data->T() + (1-r)/r*node_data->T_last())*node_data->density()*HeatCapacity
                       / (SolverSpecify::dt_last+SolverSpecify::dt)*fvm_node->volume();
      // ADD to Jacobian matrix,
      MatSetValue(*jac, fvm_node->global_offset()+1, fvm_node->global_offset()+1, TT.getADValue(0), ADD_VALUES);
    }
    else //first order
    {
      AutoDScalar TT = -(T - node_data->T())*node_data->density()*HeatCapacity/SolverSpecify::dt*fvm_node->volume();
      // ADD to Jacobian matrix,
      MatSetValue(*jac, fvm_node->global_offset()+1, fvm_node->global_offset()+1, TT.getADValue(0), ADD_VALUES);
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}



void ConductorSimulationRegion::DDM2_Update_Solution(PetscScalar *lxx)
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
    node_data->psi() = lxx[fvm_node->local_offset()+0];

    //update T
    node_data->T_last()   = node_data->T();
    node_data->T()   = lxx[fvm_node->local_offset()+1];
  }

  // addtional work: compute electrical field for all the node.
  // however, the electrical field is always zero. We needn't do anything here.

}




