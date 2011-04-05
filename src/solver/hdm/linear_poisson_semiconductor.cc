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
#include "solver_specify.h"

using PhysicalUnit::cm;
using PhysicalUnit::kb;
using PhysicalUnit::e;


void SemiconductorSimulationRegion::LinearPoissin_RHS(Vec b, InsertMode &add_value_flag)
{
  PetscInt n_local_dofs;
  VecGetLocalSize(b, &n_local_dofs);

  // set local buf here
  std::vector<int>          iy;
  std::vector<PetscScalar>  y;
  iy.reserve(n_local_dofs);
  y.reserve(n_local_dofs);

  // process \rho of poisson's equation
  const_processor_node_iterator node_it = on_processor_nodes_begin();
  const_processor_node_iterator node_it_end = on_processor_nodes_end();
  for(; node_it!=node_it_end; ++node_it)
  {
    const FVM_Node * fvm_node = *node_it;

    const FVM_NodeData * fvm_node_data = fvm_node->node_data();

    iy.push_back(fvm_node->global_offset());        // save index in the buffer
    y.push_back( -fvm_node_data->rho()*fvm_node->volume() );            // save value in the buffer
  }

  if(iy.size())
    VecSetValues(b, iy.size(), &iy[0], &y[0], INSERT_VALUES);

  add_value_flag = INSERT_VALUES;
}


void SemiconductorSimulationRegion::LinearPoissin_Matrix(Mat A, InsertMode &add_value_flag)
{
  //the indepedent variable number, since we only process edges, 2 is enough
  adtl::AutoDScalar::numdir=2;

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
      const FVM_Node * fvm_n1 = (*it)->get_fvm_node(edge_nodes[0]);
      // fvm_node of node2
      const FVM_Node * fvm_n2 = (*it)->get_fvm_node(edge_nodes[1]);


      // fvm_node_data of node1
      const FVM_NodeData * n1_data =  fvm_n1->node_data();
      // fvm_node_data of node2
      const FVM_NodeData * n2_data =  fvm_n2->node_data();

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
        AutoDScalar V1   =  0.0;   V1.setADValue(0, 1.0);
        // eps
        PetscScalar eps1 =  n1_data->eps();


        // electrostatic potential, as independent variable
        AutoDScalar V2   =  0.0;   V2.setADValue(1, 1.0);
        // eps
        PetscScalar eps2 =  n2_data->eps();

        PetscScalar eps = 0.5*(eps1+eps2);

        // ignore thoese ghost nodes
        if( fvm_n1->on_processor() )
        {
          AutoDScalar f =  eps*partial_area*(V2 - V1)/length ;
          MatSetValues(A, 1, &row[0], col.size(), &col[0], f.getADValue(), ADD_VALUES);
        }

        if( fvm_n2->on_processor() )
        {
          AutoDScalar f =  eps*partial_area*(V1 - V2)/length ;
          MatSetValues(A, 1, &row[1], col.size(), &col[0], f.getADValue(), ADD_VALUES);
        }

      }
    }

  }

  add_value_flag = ADD_VALUES;
}

#include <TNT/jama_lu.h>
void SemiconductorSimulationRegion::LinearPoissin_Update_Solution(const PetscScalar * x)
{
  PetscScalar damping_density = 1e17*std::pow(cm, -3);
  PetscScalar dV_max = 0.0; // the max changes of psi

  {
    const_processor_node_iterator node_it = on_processor_nodes_begin();
    const_processor_node_iterator node_it_end = on_processor_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      const FVM_Node * fvm_node = *node_it;
      const FVM_NodeData * node_data = fvm_node->node_data();
      dV_max = std::max(dV_max, fabs(x[fvm_node->local_offset()] - node_data->psi()));
    }
  }

  // for parallel situation, we should find the dv_max in global.
  //Parallel::max( dV_max );
  PetscScalar Vt = kb*T_external()/e;
  PetscScalar f = log(1+dV_max/Vt)/(dV_max/Vt);

  {
    local_node_iterator node_it = on_local_nodes_begin();
    local_node_iterator node_it_end = on_local_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      FVM_Node * fvm_node = *node_it;
      FVM_NodeData * node_data = fvm_node->node_data();

      //update psi
      node_data->psi_last() =  node_data->psi();
      node_data->psi()      =  node_data->psi() + 1.0*(x[fvm_node->local_offset()] - node_data->psi()) ;
    }
  }

  // calculate E with least squares
  {
    TNT::Array2D<Real> A (3, 3, 0.0);
    TNT::Array1D<Real> r (3, 0.0);

    processor_node_iterator node_it = on_processor_nodes_begin();
    processor_node_iterator node_it_end = on_processor_nodes_end();
    for(; node_it!=node_it_end; ++node_it)
    {
      FVM_Node * fvm_node = *node_it;
      FVM_NodeData * node_data = fvm_node->node_data();

      PetscScalar a=0,b=0,c=0,d=0,e=0,f=0;
      PetscScalar r1=0,r2=0,r3=0;
      FVM_Node::fvm_neighbor_node_iterator neighbor_it = fvm_node->neighbor_node_begin();
      FVM_Node::fvm_neighbor_node_iterator neighbor_it_end =  fvm_node->neighbor_node_end();
      for(;neighbor_it!=neighbor_it_end; ++neighbor_it)
      {
        const FVM_Node * neighbor_node = neighbor_it->second;
        const FVM_NodeData * neighbor_node_data = neighbor_node->node_data();

        PetscScalar w = 1.0/ fvm_node->distance(neighbor_node);
        PetscScalar dx = neighbor_node->root_node()->x() - fvm_node->root_node()->x();
        PetscScalar dy = neighbor_node->root_node()->y() - fvm_node->root_node()->y();
        PetscScalar dz = neighbor_node->root_node()->z() - fvm_node->root_node()->z();
        PetscScalar dphi = x[neighbor_node->local_offset()] - x[fvm_node->local_offset()];
        a += w*dx*dx;
        b += w*dx*dy;
        c += w*dx*dz;
        d += w*dy*dy;
        e += w*dy*dz;
        f += w*dz*dz;

        r1 += w*dx*dphi;
        r2 += w*dy*dphi;
        r3 += w*dz*dphi;
      }


      // for 2D case, f is 0, which breaks LU solver
      // so we set it to 1.0
      if(fabs(f)<1e-6) f=1.0;

      A[0][0] = a; A[0][1] = b; A[0][2] = c;
      A[1][0] = b; A[1][1] = d; A[1][2] = e;
      A[2][0] = c; A[2][1] = e; A[2][2] = f;

      r[0] = r1;
      r[1] = r2;
      r[2] = r3;

      JAMA::LU<Real> solver(A);
      TNT::Array1D<Real> dphi = solver.solve(r);

      //PetscScalar carrier = node_data->n() + node_data->p();
      //PetscScalar damping = carrier > damping_density ? damping_density/carrier : 1.0;

      node_data->E() = -VectorValue<PetscScalar>(dphi[0], dphi[1], dphi[2]);

    }
  }
}


