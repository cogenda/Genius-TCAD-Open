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

//  $Id: face_quad4_fvm.cc,v 1.6 2008/07/10 09:39:38 gdiso Exp $

#include "face_quad4_fvm.h"

#include <TNT/jama_lu.h>

// build the Geom information here
void Quad4_FVM::prepare_for_fvm()
{
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

  // Vector pointing from A to C
  Point AC ( this->point(2) - this->point(0) );

  // Vector pointing from A to B
  Point AB ( this->point(1) - this->point(0) );

  // Vector pointing from B to C
  Point BC ( this->point(2) - this->point(1) );

  // Vector normal to plan ABC
  Point n = AB.cross(AC);

  // Vector in plan ABC and normal to AB
  Point a = n.cross(AB);

  // the angle between AB and AC
  Real angle_a = AB.angle(AC);

  // the angle between CA and CB
  Real cos_angle_c = (-AC).cos_angle(-BC);

  // the radius of circumcircle
  Real R = BC.size() / (2*sin(angle_a));

  // the circle centre of ABC
  Point circumcircle_center;
  if(cos_angle_c > 0) // for sharp angle, the circle centre is inside of triangle ABC
    circumcircle_center = this->point(0) + 0.5*AB + a.unit() * sqrt( std::abs(R*R - 0.25*AB.size()*AB.size()) );
  else                 // oh, it is an obtuse angle, the circle centre is out of triangle ABC
    circumcircle_center = this->point(0) + 0.5*AB - a.unit() * sqrt( std::abs(R*R - 0.25*AB.size()*AB.size()) );

  // the side (edge) center
  Point side_centers[4];
  for( unsigned int i=0; i<4; i++ )
  {
    Point p1 =  this->point(side_nodes_map[i][0]) ;
    Point p2 =  this->point(side_nodes_map[i][1]) ;
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
  }

  //store this vlaue for efficiency reason
  vol = Quad4::volume();

  // clear partitial volume
  v[0] = v[1] = v[2] = v[3] = 0;
  // add volume
  for( unsigned int i=0; i<4; i++ ) // has 4 side (edge)
   {
     unsigned int node1 = side_nodes_map[i][0];
     unsigned int node2 = side_nodes_map[i][1];
     v[node1] +=  0.5 * 0.5 * l[i] * d[i];
     v[node2] +=  0.5 * 0.5 * l[i] * d[i];
   }

  // self test
  Real V = 0;
  for( unsigned int i=0; i<4; i++ )
  {
     V += v[i];
  }

  genius_assert( std::abs(V-vol) < 1e-4*vol );

  // set data for least square fitting
  prepare_for_least_squares();

  // set least square fitting data for vector reconstruction
  prepare_for_vector_reconstruct();
}


void Quad4_FVM::prepare_for_least_squares()
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



VectorValue<PetscScalar> Quad4_FVM::gradient( const std::vector<PetscScalar> & var) const
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


VectorValue<Complex> Quad4_FVM::gradient( const std::vector<Complex> & var) const
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


VectorValue<AutoDScalar> Quad4_FVM::gradient( const std::vector<AutoDScalar> & var) const
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


void Quad4_FVM::prepare_for_vector_reconstruct()
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


VectorValue<PetscScalar> Quad4_FVM::reconstruct_vector( const std::vector<PetscScalar> & projects) const
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


VectorValue<AutoDScalar> Quad4_FVM::reconstruct_vector( const std::vector<AutoDScalar> & projects) const
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

