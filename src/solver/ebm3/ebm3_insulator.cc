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
#include "insulator_region.h"
#include "solver_specify.h"

using PhysicalUnit::kb;
using PhysicalUnit::e;


/*-------------------------------------------------------------------------------------
 * get nodal variable number, which depedent on EBM Level
 */
unsigned int InsulatorSimulationRegion::ebm_n_variables() const
{
  switch(get_advanced_model()->EB_Level)
  {
  case ModelSpecify::NONE :
  case ModelSpecify::Tn   :
  case ModelSpecify::Tp   :
  case ModelSpecify::TnTp : return 1; // 1 dofs for psi

  case ModelSpecify::Tl   :
  case ModelSpecify::TnTl :
  case ModelSpecify::TpTl :
  case ModelSpecify::ALL  : return 2; // 2 dofs for psi and Tl
  }
  return 0; // prevent compiler warning.
}


/*-------------------------------------------------------------------------------------
 * get offset of nodal variable, which depedent on EBM Level
 */
unsigned int InsulatorSimulationRegion::ebm_variable_offset(SolutionVariable var) const
{
  switch(var)
  {
  case POTENTIAL     : return 0; // psi is always at offset 0
  case TEMPERATURE   : return 1; // if lattice temperature is required?
  case ELECTRON      :           // followed by electron density
  case HOLE          :           // and hole density
  case E_TEMP        :           // if electron temperature is required?
  case H_TEMP        :           // if hole temperature is required?
  default : return invalid_uint;
  }
}


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////



void InsulatorSimulationRegion::EBM3_Fill_Value(Vec x, Vec L)
{

  // find the node variable offset
  unsigned int node_psi_offset = ebm_variable_offset(POTENTIAL);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);

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

    // psi
    genius_assert(node_psi_offset!=invalid_uint);
    ix.push_back(node->global_offset()+node_psi_offset);
    y.push_back(node_data->psi());
    s.push_back(1.0/node->volume());

    // for extra Tl temperature equations
    if(get_advanced_model()->enable_Tl())
    {
      genius_assert(node_Tl_offset!=invalid_uint);
      ix.push_back(node->global_offset()+node_Tl_offset);
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



/*---------------------------------------------------------------------
 * build function and its jacobian for EBM3 solver
 */
void InsulatorSimulationRegion::EBM3_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // find the node variable offset
  unsigned int node_psi_offset = ebm_variable_offset(POTENTIAL);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);

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

      // build governing equations for insulator region
      {

        /*
         * process psi equation
         */

        //for node 1 of the edge
        PetscScalar V1   =  x[n1_local_offset + node_psi_offset];   // electrostatic potential
        PetscScalar rho1 =  0;                                      // charge density
        PetscScalar eps1 =  n1_data->eps();                         // permittivity


        //for node 2 of the edge
        PetscScalar V2   =  x[n2_local_offset + node_psi_offset];
        PetscScalar rho2 =  0;
        PetscScalar eps2 =  n2_data->eps();


        PetscScalar eps = 0.5*(eps1+eps2);       // eps at mid point of the edge

        // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {
          // poisson's equation
          iy.push_back( fvm_n1->global_offset()+node_psi_offset );
          y.push_back ( eps*partial_area*(V2 - V1)/length + ( rho1 )*(*it)->partial_volume(edge_nodes[0]) );
        }

        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {
          // poisson's equation
          iy.push_back( fvm_n2->global_offset()+node_psi_offset );
          y.push_back ( eps*partial_area*(V1 - V2)/length + ( rho2 )*(*it)->partial_volume(edge_nodes[1]) );
        }

        /*
         * process heating equation if required
         */
        if(get_advanced_model()->enable_Tl())
        {
          mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
          PetscScalar T1   =  x[n1_local_offset + node_Tl_offset];             // lattice temperature
          PetscScalar kap1 =  mt->thermal->HeatConduction(T1);

          mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
          PetscScalar T2   =  x[n2_local_offset + node_Tl_offset];
          PetscScalar kap2 =  mt->thermal->HeatConduction(T2);

          PetscScalar kap = 0.5*(kap1+kap2);       // kapa at mid point of the edge

          // ignore thoese ghost nodes (ghost nodes is local but with different processor_id())
          if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
          {
            // heat transport equation
            iy.push_back( fvm_n1->global_offset()+node_Tl_offset );
            y.push_back ( kap*partial_area*(T2 - T1)/length );

          }

          if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
          {
            // heat transport equation
            iy.push_back( fvm_n2->global_offset()+node_Tl_offset );
            y.push_back ( kap*partial_area*(T1 - T2)/length );
          }
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
 * build function and its jacobian for EBM3 solver
 */
void InsulatorSimulationRegion::EBM3_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // find the node variable offset
  unsigned int node_psi_offset = ebm_variable_offset(POTENTIAL);
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);

  //the indepedent variable number, since we only process edges, 2 (2 nodes x 2 variables_per_node) is enough
  adtl::AutoDScalar::numdir=2;

  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  // note, we will use ADD_VALUES to set values of matrix J
  // if the previous operator is not ADD_VALUES, we should flush the matrix
  if( add_value_flag != ADD_VALUES && add_value_flag != NOT_SET_VALUES)
  {
    MatAssemblyBegin(*jac, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(*jac, MAT_FLUSH_ASSEMBLY);
  }


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

        /*
         * process psi equation
         */

        //for node 1 of the edge
        AutoDScalar V1   =  x[n1_local_offset + node_psi_offset];  V1.setADValue(0, 1.0);  // electrostatic potential
        PetscScalar rho1 =  0;                                                             // charge density
        PetscScalar eps1 =  n1_data->eps();                                                // permittivity

        //for node 2 of the edge
        AutoDScalar V2   =  x[n2_local_offset + node_psi_offset];  V2.setADValue(1, 1.0);
        PetscScalar rho2 =  0;
        PetscScalar eps2 =  n2_data->eps();

        PetscScalar eps = 0.5*(eps1+eps2);       // eps at mid point of the edge

        // ignore thoese ghost nodes
        if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
        {
          AutoDScalar poisson =  eps*partial_area*(V2 - V1)/length + ( rho1 )*(*it)->partial_volume(edge_nodes[0]) ;
          MatSetValue(*jac, fvm_n1->global_offset()+node_psi_offset, fvm_n1->global_offset()+node_psi_offset, poisson.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_n1->global_offset()+node_psi_offset, fvm_n2->global_offset()+node_psi_offset, poisson.getADValue(1), ADD_VALUES);
        }

        if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
        {
          AutoDScalar poisson =  eps*partial_area*(V1 - V2)/length + ( rho2 )*(*it)->partial_volume(edge_nodes[1]) ;
          MatSetValue(*jac, fvm_n2->global_offset()+node_psi_offset, fvm_n1->global_offset()+node_psi_offset, poisson.getADValue(0), ADD_VALUES);
          MatSetValue(*jac, fvm_n2->global_offset()+node_psi_offset, fvm_n2->global_offset()+node_psi_offset, poisson.getADValue(1), ADD_VALUES);
        }

        /*
         * process heating equation if required
         */
        if(get_advanced_model()->enable_Tl())
        {
          //for node 1 of the edge
          mt->mapping(fvm_n1->root_node(), n1_data, SolverSpecify::clock);
          AutoDScalar T1   =  x[n1_local_offset+1];  T1.setADValue(0,1.0);           // lattice temperature
          PetscScalar kap1 =  mt->thermal->HeatConduction(T1.getValue());

          //for node 2 of the edge
          mt->mapping(fvm_n2->root_node(), n2_data, SolverSpecify::clock);
          AutoDScalar T2   =  x[n2_local_offset+1];  T2.setADValue(1,1.0);
          PetscScalar kap2 =  mt->thermal->HeatConduction(T2.getValue());

          PetscScalar kap = 0.5*(kap1+kap2);       // kapa at mid point of the edge


          // ignore thoese ghost nodes
          if( fvm_n1->root_node()->processor_id()==Genius::processor_id() )
          {
            AutoDScalar heat =  kap*partial_area*(T2 - T1)/length;
            MatSetValue(*jac, fvm_n1->global_offset()+node_Tl_offset, fvm_n1->global_offset()+node_Tl_offset, heat.getADValue(0), ADD_VALUES);
            MatSetValue(*jac, fvm_n1->global_offset()+node_Tl_offset, fvm_n2->global_offset()+node_Tl_offset, heat.getADValue(1), ADD_VALUES);
          }

          if( fvm_n2->root_node()->processor_id()==Genius::processor_id() )
          {
            AutoDScalar heat =  kap*partial_area*(T1 - T2)/length;
            MatSetValue(*jac, fvm_n2->global_offset()+node_Tl_offset, fvm_n1->global_offset()+node_Tl_offset, heat.getADValue(0), ADD_VALUES);
            MatSetValue(*jac, fvm_n2->global_offset()+node_Tl_offset, fvm_n2->global_offset()+node_Tl_offset, heat.getADValue(1), ADD_VALUES);
          }

        }
      }
    }


  }

  // boundary condition should be processed later!

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

}




void InsulatorSimulationRegion::EBM3_Time_Dependent_Function(PetscScalar * x, Vec f, InsertMode &add_value_flag)
{

  // find the node variable offset
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);

  // do time depedent calculation only heating equation required
  if(!get_advanced_model()->enable_Tl()) return;


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

    PetscScalar T   =  x[fvm_node->local_offset()+node_Tl_offset];                         // lattice temperature
    PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(T);

    // process \partial t

    iy.push_back(fvm_node->global_offset()+1);                                // save index in the buffer

    //second order
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false)
    {
      PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
      PetscScalar dTldt = -((2-r)/(1-r)*T - 1.0/(r*(1-r))*node_data->T() + (1-r)/r*node_data->T_last())*node_data->density()*HeatCapacity
                          / (SolverSpecify::dt_last+SolverSpecify::dt) * fvm_node->volume();
      y.push_back( dTldt );
    }
    else //first order
    {
      PetscScalar dTldt = -(T - node_data->T())*node_data->density()*HeatCapacity/SolverSpecify::dt*fvm_node->volume();
      y.push_back( dTldt );
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



void InsulatorSimulationRegion::EBM3_Time_Dependent_Jacobian(PetscScalar * x, Mat *jac, InsertMode &add_value_flag)
{

  // find the node variable offset
  unsigned int node_Tl_offset  = ebm_variable_offset(TEMPERATURE);

  // do time depedent calculation only heating equation required
  if(!get_advanced_model()->enable_Tl()) return;


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

    AutoDScalar T   =  x[fvm_node->local_offset()+node_Tl_offset];   T.setADValue(0, 1.0);              // lattice temperature
    PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(T.getValue());

    // process \partial t

    //second order
    if(SolverSpecify::TS_type==SolverSpecify::BDF2 && SolverSpecify::BDF2_restart==false)
    {
      PetscScalar r = SolverSpecify::dt_last/(SolverSpecify::dt_last + SolverSpecify::dt);
      AutoDScalar dTldt = -((2-r)/(1-r)*T - 1.0/(r*(1-r))*node_data->T() + (1-r)/r*node_data->T_last())*node_data->density()*HeatCapacity
                       / (SolverSpecify::dt_last+SolverSpecify::dt)*fvm_node->volume();
      // ADD to Jacobian matrix,
      MatSetValue(*jac, fvm_node->global_offset()+node_Tl_offset, fvm_node->global_offset()+node_Tl_offset, dTldt.getADValue(0), ADD_VALUES);
    }
    else //first order
    {
      AutoDScalar dTldt = -(T - node_data->T())*node_data->density()*HeatCapacity/SolverSpecify::dt*fvm_node->volume();
      // ADD to Jacobian matrix,
      MatSetValue(*jac, fvm_node->global_offset()+node_Tl_offset, fvm_node->global_offset()+node_Tl_offset, dTldt.getADValue(0), ADD_VALUES);
    }
  }

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#ifdef HAVE_FENV_H
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

}



void InsulatorSimulationRegion::EBM3_Update_Solution(PetscScalar *lxx)
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
    if(get_advanced_model()->enable_Tl())
    {
      node_data->T_last()   = node_data->T();
      node_data->T()   = lxx[fvm_node->local_offset()+1];
    }
  }

  // addtional work: compute electrical field for all the cell.
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


