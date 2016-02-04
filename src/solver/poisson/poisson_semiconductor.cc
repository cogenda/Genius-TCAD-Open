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

using PhysicalUnit::cm;
using PhysicalUnit::kb;
using PhysicalUnit::e;
using PhysicalUnit::eps0;


///////////////////////////////////////////////////////////////////////
//----------------Function and Jacobian evaluate---------------------//
///////////////////////////////////////////////////////////////////////


void SemiconductorSimulationRegion::Poissin_Fill_Value(Vec x, Vec L)
{
  std::vector<int> ix;
  std::vector<PetscScalar> y;
  std::vector<PetscScalar> s;

  ix.reserve(this->n_node());
  y.reserve(this->n_node());
  s.reserve(this->n_node());



  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * node_data = fvm_node->node_data();

    mt->mapping(fvm_node->root_node(), node_data, 0.0);

    ix.push_back(fvm_node->global_offset());
    y.push_back(node_data->psi());
    s.push_back(1.0/(node_data->eps()*fvm_node->volume()));
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

  // reserve memory
  iy.reserve(2*this->n_edge()+this->n_node());
  y.reserve(2*this->n_edge()+this->n_node());


  // search all the element in this region.
  // note, they are all local element, thus must be processed

  //common used variable
  PetscScalar T   = T_external();


  // process \nabla operator for all cells

  // search all the edges of this region, do integral over control volume...
  const_edge_iterator it = edges_begin();
  const_edge_iterator it_end = edges_end();
  for(; it!=it_end; ++it)
  {
    // fvm_node of node1
    const FVM_Node * fvm_n1 = (*it).first;
    // fvm_node of node2
    const FVM_Node * fvm_n2 = (*it).second;

    // fvm_node_data of node1
    const FVM_NodeData * n1_data =  fvm_n1->node_data();
    // fvm_node_data of node2
    const FVM_NodeData * n2_data =  fvm_n2->node_data();

    const unsigned int n1_local_offset = fvm_n1->local_offset();
    const unsigned int n2_local_offset = fvm_n2->local_offset();

    {
      // electrostatic potential, as independent variable
      PetscScalar V1   =  x[n1_local_offset];
      PetscScalar eps1 =  n1_data->eps();


      PetscScalar V2   =  x[n2_local_offset];
      PetscScalar eps2 =  n2_data->eps();

      PetscScalar eps = 0.5*(eps1+eps2);

      // "flux" from node 2 to node 1
      PetscScalar f =  eps*fvm_n1->cv_surface_area(fvm_n2)*(V2 - V1)/fvm_n1->distance(fvm_n2) ;

      // ignore thoese ghost nodes
      if( fvm_n1->on_processor() )
      {
        iy.push_back(fvm_n1->global_offset());
        y.push_back(f);
      }

      if( fvm_n2->on_processor() )
      {
        iy.push_back(fvm_n2->global_offset());
        y.push_back(-f);
      }
    }
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // process \rho of poisson's equation
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * fvm_node_data = fvm_node->node_data();

    PetscScalar V   =  x[fvm_node->local_offset()];
    mt->mapping(fvm_node->root_node(), fvm_node_data, 0.0);

    // intrinsic Fermi potential.
    PetscScalar V_i =  V + fvm_node_data->affinity() + fvm_node_data->Eg()/(2*e) + kb*T*log(fvm_node_data->Nc()/fvm_node_data->Nv())/(2*e);

    PetscScalar ni = mt->band->ni(T);
    // approx electron density
    PetscScalar n   =  ni*exp( e/(kb*T)*V_i);
    // approx hole density
    PetscScalar p   =  ni*exp(-e/(kb*T)*V_i);
#if 0
    // bandgap narrowing
    PetscScalar dEg =  mt->band->EgNarrow(p, n, T);
    n = n*exp(dEg);
    p = p*exp(dEg);
#endif
    // charge density
    PetscScalar doping = fvm_node_data->Net_doping();
    if(get_advanced_model()->IncompleteIonization)
      doping = mt->band->Nd_II(n, T, get_advanced_model()->Fermi) - mt->band->Na_II(p, T, get_advanced_model()->Fermi);
    PetscScalar rho = e*( doping + p - n)*fvm_node->volume(); // the charge density

    iy.push_back(fvm_node->global_offset());                                  // save index in the buffer
    y.push_back( rho );                                                       // save value in the buffer

  }

  if(iy.size()) VecSetValues(f, iy.size(), &iy[0], &y[0], ADD_VALUES);

  // after the first scan, every nodes are updated.
  // however, boundary condition should be processed later.

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}






/*---------------------------------------------------------------------
 * build function and its jacobian for poisson solver
 */
void SemiconductorSimulationRegion::Poissin_Jacobian(PetscScalar * x, SparseMatrix<PetscScalar> *jac, InsertMode &add_value_flag)
{

  //common used variable
  PetscScalar T   = T_external();


  // search all the element in this region.
  // note, they are all local element, thus must be processed

  //the indepedent variable is 2
  adtl::AutoDScalar::numdir = 2;
  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

 // search all the edges of this region, do integral over control volume...
  const_edge_iterator it = edges_begin();
  const_edge_iterator it_end = edges_end();
  for(; it!=it_end; ++it)
  {
    // fvm_node of node1
    const FVM_Node * fvm_n1 = (*it).first;
    // fvm_node of node2
    const FVM_Node * fvm_n2 = (*it).second;

    // fvm_node_data of node1
    const FVM_NodeData * n1_data =  fvm_n1->node_data();
    // fvm_node_data of node2
    const FVM_NodeData * n2_data =  fvm_n2->node_data();

    const unsigned int n1_local_offset = fvm_n1->local_offset();
    const unsigned int n2_local_offset = fvm_n2->local_offset();

    // the row/colume position of variables in the matrix
    unsigned int row[2],col[2];
    row[0] = col[0] = fvm_n1->global_offset();
    row[1] = col[1] = fvm_n2->global_offset();

    // here we use AD, however it is great overkill for such a simple problem.
    {
      // electrostatic potential, as independent variable
      AutoDScalar V1   =  x[n1_local_offset];   V1.setADValue(0,1.0);
      PetscScalar eps1 =  n1_data->eps();


      AutoDScalar V2   =  x[n2_local_offset];   V2.setADValue(1,1.0);
      PetscScalar eps2 =  n2_data->eps();

      PetscScalar eps = 0.5*(eps1+eps2);

      AutoDScalar f =  eps*fvm_n1->cv_surface_area(fvm_n2)*(V2 - V1)/fvm_n1->distance(fvm_n2) ;

      // ignore thoese ghost nodes
      if( fvm_n1->on_processor() )
      {
        jac->add_row(row[0], 2, &col[0], f.getADValue());
      }

      if( fvm_n2->on_processor() )
      {
        jac->add_row(row[1], 2, &col[0], (-f).getADValue());
      }
    }
  }

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif

  // process \rho of poisson's equation

  //the indepedent variable is 1
  adtl::AutoDScalar::numdir = 1;
  //synchronize with material database
  mt->set_ad_num(adtl::AutoDScalar::numdir);

  // process \rho of poisson's equation, search all the vertex of this cell
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;
    const FVM_NodeData * fvm_node_data = fvm_node->node_data();

    AutoDScalar V   =  x[fvm_node->local_offset()];  V.setADValue(0, 1.0);
    mt->mapping(fvm_node->root_node(), fvm_node_data, 0.0);

    PetscScalar ni = mt->band->ni(T);

    // intrinsic Fermi potential.
    AutoDScalar V_i =  V + fvm_node_data->affinity() + fvm_node_data->Eg()/(2*e) + kb*T*log(fvm_node_data->Nc()/fvm_node_data->Nv())/(2*e);
    // electron density
    AutoDScalar n   =  ni*exp( e/(kb*T)*V_i);
    // hole density
    AutoDScalar p   =  ni*exp(-e/(kb*T)*V_i);

#if 0
    // bandgap narrowing
    AutoDScalar dEg =  mt->band->EgNarrow(p, n, T);
    n = n*exp(dEg);
    p = p*exp(dEg);
#endif

    // charge density
    AutoDScalar doping = fvm_node_data->Net_doping();
    if(get_advanced_model()->IncompleteIonization)
      doping = mt->band->Nd_II(n, T, get_advanced_model()->Fermi) - mt->band->Na_II(p, T, get_advanced_model()->Fermi);
    AutoDScalar rho = e*( doping + p - n)*fvm_node->volume(); // the charge density

    jac->add(fvm_node->global_offset(), fvm_node->global_offset(), rho.getADValue(0));  // save value in the buffer

  }

  // boundary condition should be processed later!

  // the last operator is ADD_VALUES
  add_value_flag = ADD_VALUES;

#if defined(HAVE_FENV_H) && defined(DEBUG)
  genius_assert( !fetestexcept(FE_INVALID) );
#endif
}



void SemiconductorSimulationRegion::Poissin_Update_Solution(PetscScalar *lxx)
{
  //common used variable
  const PetscScalar T   = T_external();

  local_node_iterator node_it = on_local_nodes_begin();
  local_node_iterator node_it_end = on_local_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    FVM_Node * fvm_node = *node_it;
    FVM_NodeData * node_data = fvm_node->node_data();

    PetscScalar V   = lxx[fvm_node->local_offset()];


    //for semiconductor region, also update n and p
    {
      PetscScalar T   = T_external();
      mt->mapping(fvm_node->root_node(), node_data, 0.0);

      PetscScalar ni = mt->band->ni(T);

      // intrinsic Fermi potential.
      PetscScalar V_i = V
                        + node_data->affinity()/e
                        + node_data->Eg()/(2*e)
                        + kb*T*log(node_data->Nc()/node_data->Nv())/(2*e);
      // electron density
      PetscScalar n    =  ni*exp(e/(kb*T)*V_i);
      // hole density
      PetscScalar p    =  ni*exp(-e/(kb*T)*V_i);
#if 0
      // consider bandgap narrow
      PetscScalar dEg = mt->band->EgNarrow(p, n, T);
      node_data->n()   =  n*exp(dEg/(2*kb*T));
      node_data->p()   =  p*exp(dEg/(2*kb*T));
#endif
      //update psi
      PetscScalar Eg = mt->band->Eg(T);
      node_data->psi() = V;
      node_data->Eg() = Eg;
      node_data->Ec() = -(e*V + node_data->affinity());
      node_data->Ev() = -(e*V + node_data->affinity()+ Eg);

      node_data->qFn() = node_data->Ec() + log(fabs(n/node_data->Nc()))*kb*T;
      node_data->qFp() = node_data->Ev() - log(fabs(p/node_data->Nv()))*kb*T;
    }
  }

  // addtional work: compute electrical field for all the node.
  // Since this value is only used for reference.
  // It can be done simply by weighted average of cell's electrical field
  for(unsigned int n=0; n<n_cell(); ++n)
  {
    const Elem * elem = this->get_region_elem(n);
    FVM_CellData * elem_data = this->get_region_elem_data(n);

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




