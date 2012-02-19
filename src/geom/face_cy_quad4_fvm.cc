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


#include "face_cy_quad4_fvm.h"
#include "edge_edge2_fvm.h"
#include <TNT/jama_lu.h>

//#define DEBUG

AutoPtr<Elem> Quad4_CY_FVM::build_fvm_side (const unsigned int i, bool proxy) const
{
  assert (i < this->n_sides());

  if (proxy)
  {
    return Quad4::build_side(i, proxy);
  }
  else
  {
    Edge2_FVM* edge = new Edge2_FVM;
    switch (i)
    {
        case 0:
        {
          edge->set_node(0) = this->get_node(0);
          edge->set_node(1) = this->get_node(1);
          edge->subdomain_id() = this->subdomain_id();
          edge->hold_fvm_node(0, this->get_fvm_node(0));
          edge->hold_fvm_node(1, this->get_fvm_node(1));
          edge->prepare_for_fvm();
          AutoPtr<Elem> ap(edge);  return ap;
        }
        case 1:
        {
          edge->set_node(0) = this->get_node(1);
          edge->set_node(1) = this->get_node(2);
          edge->subdomain_id() = this->subdomain_id();
          edge->hold_fvm_node(0, this->get_fvm_node(1));
          edge->hold_fvm_node(1, this->get_fvm_node(2));
          edge->prepare_for_fvm();
          AutoPtr<Elem> ap(edge);  return ap;
        }
        case 2:
        {
          edge->set_node(0) = this->get_node(2);
          edge->set_node(1) = this->get_node(3);
          edge->subdomain_id() = this->subdomain_id();
          edge->hold_fvm_node(0, this->get_fvm_node(2));
          edge->hold_fvm_node(1, this->get_fvm_node(3));
          edge->prepare_for_fvm();
          AutoPtr<Elem> ap(edge);  return ap;
        }
        case 3:
        {
          edge->set_node(0) = this->get_node(3);
          edge->set_node(1) = this->get_node(0);
          edge->subdomain_id() = this->subdomain_id();
          edge->hold_fvm_node(0, this->get_fvm_node(3));
          edge->hold_fvm_node(1, this->get_fvm_node(0));
          edge->prepare_for_fvm();
          AutoPtr<Elem> ap(edge);  return ap;
        }
        default:
        {
          genius_error();
        }
    }
  }

  // We will never get here...
  AutoPtr<Elem> ap(NULL);  return ap;
}


// build the Geom information here
void Quad4_CY_FVM::prepare_for_fvm()
{
  const Real pi = 3.14159265358979323846;
  // The A,B,C,D naming scheme here corresponds exactly to the
  // libmesh counter-clockwise numbering scheme.

  //        3           2        D           C
  // QUAD4: o-----------o    o-----------o
  //        |           |    |           |
  //        |           |    |           |
  //        |           |    |           |
  //        |           |    |           |
  //        |           |    |           |
  //        o-----------o    o-----------o
  //        0           1    A           B

  Point circumcircle_center;
  {
    Point v12 = this->point(0) - this->point(1);
    Point v23 = this->point(1) - this->point(2);
    Point v13 = this->point(0) - this->point(2);

    Real ccdet = v12.cross(v23).size_sq();
    Real alpha = v23.size_sq() * v12.dot(v13) /2.0/ccdet;
    Real beta  = - v13.size_sq() * v12.dot(v23) /2.0/ccdet;
    Real gamma = v12.size_sq() * v13.dot(v23) /2.0/ccdet;

    circumcircle_center = alpha * this->point(0) +
                          beta  * this->point(1) +
                          gamma * this->point(2);
  }

  Point weight_center = this->centroid();

  //store this vlaue for efficiency reason
  vol = Quad4::volume()*2*pi*weight_center[0]; //Pappus's centroid theorem

  // clear partitial volume
  v[0] = v[1] = v[2] = v[3] = 0;

  // the side (edge) center
  Point side_centers[4];
  for( unsigned int i=0; i<4; i++ )
  {
    unsigned int node1 = side_nodes_map[i][0];
    unsigned int node2 = side_nodes_map[i][1];
    Point p1 =  this->point(node1) ;
    Point p2 =  this->point(node2) ;
    Point p3 =  this->point((2+i)%4) ;
    // the side (edge) center
    side_centers[i] = 0.5*(p1+p2);
    // the side (edge) length
    l[i] = (p1-p2).size();
    // the distance from circumcircle center to edge
    if( (p1-p3).cos_angle(p2-p3) < 0 )
      d[i] = -(side_centers[i] - circumcircle_center).size();
    else
      d[i] =  (side_centers[i] - circumcircle_center).size();

    Real L  = 2*pi*(0.5*(side_centers[i] + circumcircle_center)(0));
    Real K1 = 2*pi*(1.0/3.0*(side_centers[i] + circumcircle_center + p1)(0));
    Real K2 = 2*pi*(1.0/3.0*(side_centers[i] + circumcircle_center + p2)(0));

    v[node1] +=  0.5 * 0.5 * l[i] * d[i] * K1;
    v[node2] +=  0.5 * 0.5 * l[i] * d[i] * K2;
    d[i] *= L;
  }

  // self test
#ifdef DEBUG
  Real V = 0;
  for( unsigned int i=0; i<4; i++ )
  {
    V += v[i];
  }

  genius_assert( (vol > 1e-10 ? (std::abs(V-vol) < 1e-3*vol) : (std::abs(V-vol) < std::max(1e-13, 1e-2*vol))) );
#endif

  // set data for least square fitting
  prepare_for_least_squares();

  // set least square fitting data for vector reconstruction
  prepare_for_vector_reconstruct();
}



void Quad4_CY_FVM::prepare_for_least_squares()
{
  TNT::Array2D<Real> A (n_nodes(), 3, 0.0);
  TNT::Array2D<Real> AT(3, n_nodes(), 0.0);

  for( unsigned int i=0; i<n_nodes(); i++ )
  {
    A[i][0] = 1.0;
    A[i][1] = (this->point(i))(0);
    A[i][2] = (this->point(i))(1);

    AT[0][i] = 1.0;
    AT[1][i] = (this->point(i))(0);
    AT[2][i] = (this->point(i))(1);
  }

  TNT::Array2D<Real> ATA = TNT::matmult(AT, A);

  JAMA::LU<Real> solver(ATA);

  TNT::Array2D<Real> inv_ATA = solver.inv();

  least_squares_gradient_matrix = TNT::matmult(inv_ATA, AT);

}



VectorValue<PetscScalar> Quad4_CY_FVM::gradient( const std::vector<PetscScalar> & var) const
{
  // FIXME we assume the quad on xy plane here. not a general case
  genius_assert( (this->point(0))(2) == (this->point(1))(2) );
  genius_assert( (this->point(0))(2) == (this->point(2))(2) );
  genius_assert( (this->point(0))(2) == (this->point(3))(2) );
  genius_assert( var.size()==n_nodes() );

  PetscScalar dx=0, dy=0, dz=0;
  for( unsigned int i=0; i<n_nodes(); i++ )
  {
    dx += least_squares_gradient_matrix[1][i] * var[i];
    dy += least_squares_gradient_matrix[2][i] * var[i];
  }

  return VectorValue<PetscScalar>(dx, dy, dz);

}


VectorValue<Complex> Quad4_CY_FVM::gradient( const std::vector<Complex> & var) const
{
  // FIXME we assume the quad on xy plane here. not a general case
  genius_assert( (this->point(0))(2) == (this->point(1))(2) );
  genius_assert( (this->point(0))(2) == (this->point(2))(2) );
  genius_assert( (this->point(0))(2) == (this->point(3))(2) );
  genius_assert( var.size()==n_nodes() );

  Complex dx=0, dy=0, dz=0;
  for( unsigned int i=0; i<n_nodes(); i++ )
  {
    dx += least_squares_gradient_matrix[1][i] * var[i];
    dy += least_squares_gradient_matrix[2][i] * var[i];
  }

  return VectorValue<Complex>(dx, dy, dz);

}


VectorValue<AutoDScalar> Quad4_CY_FVM::gradient( const std::vector<AutoDScalar> & var) const
{
  // FIXME we assume the quad on xy plane here. not a general case
  genius_assert( (this->point(0))(2) == (this->point(1))(2) );
  genius_assert( (this->point(0))(2) == (this->point(2))(2) );
  genius_assert( (this->point(0))(2) == (this->point(3))(2) );
  genius_assert( var.size()==n_nodes() );

  AutoDScalar dx=0, dy=0, dz=0;

  for( unsigned int i=0; i<n_nodes(); i++ )
  {
    dx += least_squares_gradient_matrix[1][i] * var[i];
    dy += least_squares_gradient_matrix[2][i] * var[i];
  }

  return VectorValue<AutoDScalar>(dx, dy, dz);
}


void Quad4_CY_FVM::prepare_for_vector_reconstruct()
{
  TNT::Array2D<Real> A (n_edges(), 2, 0.0);
  TNT::Array2D<Real> AT(2, n_edges(), 0.0);

  for( unsigned int e=0; e<n_edges(); e++ )
  {
    AutoPtr<Elem> edge = this->build_edge (e);
    VectorValue<double> dir = (edge->point(1) - edge->point(0)).unit(); // unit direction of the edge

    A[e][0] = dir(0);
    A[e][1] = dir(1);

    AT[0][e] = dir(0);
    AT[1][e] = dir(1);
  }

  TNT::Array2D<Real> ATA = TNT::matmult(AT, A);

  JAMA::LU<Real> solver(ATA);

  TNT::Array2D<Real> inv_ATA = solver.inv();

  least_squares_vector_reconstruct_matrix = TNT::matmult(inv_ATA, AT);
}


VectorValue<PetscScalar> Quad4_CY_FVM::reconstruct_vector( const std::vector<PetscScalar> & projects) const
{
  assert(projects.size() == n_edges());
  PetscScalar Vx=0, Vy=0, Vz=0;
  for( unsigned int e=0; e<n_edges(); e++ )
  {
    Vx += least_squares_vector_reconstruct_matrix[0][e] * projects[e];
    Vy += least_squares_vector_reconstruct_matrix[1][e] * projects[e];
  }

  return VectorValue<PetscScalar>(Vx, Vy, Vz);
}


VectorValue<AutoDScalar> Quad4_CY_FVM::reconstruct_vector( const std::vector<AutoDScalar> & projects) const
{
  assert(projects.size() == n_edges());
  AutoDScalar Vx=0, Vy=0, Vz=0;
  for( unsigned int e=0; e<n_edges(); e++ )
  {
    Vx += least_squares_vector_reconstruct_matrix[0][e] * projects[e];
    Vy += least_squares_vector_reconstruct_matrix[1][e] * projects[e];
  }

  return VectorValue<AutoDScalar>(Vx, Vy, Vz);
}


