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

//  $Id: face_tri3_fvm.cc,v 1.6 2008/07/10 09:39:38 gdiso Exp $

#include "face_tri3_fvm.h"

#include <TNT/jama_lu.h>

/*
 *   TRI3:  2               TRI3:  C
 *          o                      o
 *         / \                    / \
 *        /   \                  /   \
 *       /     \                /     \
 *      /       \              /       \
 *     /         \            /         \
 *    o-----------o          o-----------o
 *    0           1          A           B
 */


// build the Geom information here
void Tri3_FVM::prepare_for_fvm()
{
  // the circle centre of ABC
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

  Point side_centers[3];
  unsigned int obtuse = invalid_uint;
  for( unsigned int i=0; i<3; i++ )
  {
    Point p1 =  this->point(side_nodes_map[i][0]) ;
    Point p2 =  this->point(side_nodes_map[i][1]) ;
    Point p3 =  this->point((2+i)%3) ;
    // the side (edge) center
    side_centers[i] = 0.5*(p1+p2);
    // the side (edge) length
    l[i] = (p1-p2).size();
    // the distance from circumcircle center to edge
    if( (p1-p3).cos_angle(p2-p3) < 0 )
      obtuse = i; //we need special process to obtuse angle
    else
      d[i] =  (side_centers[i] - circumcircle_center).size();
  }

  // special process to obtuse angle: Truncated negtive length
  if(obtuse!=invalid_uint)
  {
    Point p1 =  this->point(side_nodes_map[obtuse][0]) ;
    Point p2 =  this->point(side_nodes_map[obtuse][1]) ;
    Point p3 =  this->point((2+obtuse)%3) ;
    double a1 = (p1-p2).angle(p1-p3);
    double a2 = (p2-p1).angle(p2-p3);

    unsigned int pre = (obtuse == 0 ? 2 :  obtuse-1);
    unsigned int pos = (obtuse == 2 ? 0 :  obtuse+1);

    d[obtuse] = 0;
    d[pre] = 0.5*(p3-p1).size()*tan(a1);
    d[pos] = 0.5*(p3-p2).size()*tan(a2);
  }

  //store this vlaue for efficiency reason
  vol = Tri3::volume();

  // clear partitial volume
  v[0] = v[1] = v[2] = 0;
  // add volume
  for( unsigned int i=0; i<3; i++ ) // has 3 side (edge)
  {
    unsigned int node1 = side_nodes_map[i][0];
    unsigned int node2 = side_nodes_map[i][1];
    v[node1] +=  0.5 * 0.5 * l[i] * d[i];
    v[node2] +=  0.5 * 0.5 * l[i] * d[i];
  }

  if(obtuse!=invalid_uint)
  {
    v[obtuse] += vol - v[0] - v[1] - v[2];
  }

  // self test
  Real V = 0;
  for( unsigned int i=0; i<3; i++ )
  {
    V += v[i];
  }

  genius_assert( vol > 1e-10 ? (std::abs(V-vol) < 1e-3*vol) : (std::abs(V-vol) < 1e-2*vol) );

  prepare_for_vector_reconstruct();
}


#if 0
// use weight center instead of circumcircle center?
void Tri3_FVM::prepare_for_fvm()
{
  //weight center
  Point c = (this->point(0) + this->point(1) + this->point(2))/3;

  Point side_centers[3];
  for( unsigned int i=0; i<3; i++ )
  {
    Point p1 =  this->point(side_nodes_map[i][0]) ;
    Point p2 =  this->point(side_nodes_map[i][1]) ;
    Point p3 =  this->point((2+i)%3) ;
    // the side (edge) center
    side_centers[i] = 0.5*(p1+p2);
    double angle = (p1-side_centers[i]).angle(c-side_centers[i]);
    // the side (edge) length
    l[i] = (p1-p2).size();
    d[i] =  (side_centers[i] - c).size()*sin(angle);
  }

  //store this vlaue for efficiency reason
  vol = Tri3::volume();

  // clear partitial volume
  v[0] = v[1] = v[2] = 0;
  // add volume
  for( unsigned int i=0; i<3; i++ ) // has 3 side (edge)
  {
    unsigned int node1 = side_nodes_map[i][0];
    unsigned int node2 = side_nodes_map[i][1];
    v[node1] +=  0.5 * 0.5 * l[i] * d[i];
    v[node2] +=  0.5 * 0.5 * l[i] * d[i];
  }

  // self test
  Real V = 0;
  for( unsigned int i=0; i<3; i++ )
  {
    V += v[i];
  }

  genius_assert( vol > 1e-10 ? (std::abs(V-vol) < 1e-3*vol) : (std::abs(V-vol) < 1e-2*vol) );
}
#endif



VectorValue<PetscScalar> Tri3_FVM::gradient( const std::vector<PetscScalar> & var) const
{
  // FIXME we assume the triangle on xy plane here. not a general case
  genius_assert( (this->point(0))(2) == (this->point(1))(2) );
  genius_assert( (this->point(0))(2) == (this->point(2))(2) );

  Real xa = (this->point(0))(0);
  Real xb = (this->point(1))(0);
  Real xc = (this->point(2))(0);
  Real ya = (this->point(0))(1);
  Real yb = (this->point(1))(1);
  Real yc = (this->point(2))(1);

  PetscScalar dx = ((yb-yc)*var[0] + (yc-ya)*var[1] +(ya-yb)*var[2])/(2*vol);
  PetscScalar dy = ((xc-xb)*var[0] + (xa-xc)*var[1] +(xb-xa)*var[2])/(2*vol);
  return VectorValue<PetscScalar>(dx, dy, 0.0);
}


VectorValue<Complex> Tri3_FVM::gradient( const std::vector<Complex> & var) const
{
  // FIXME we assume the triangle on xy plane here. not a general case
  genius_assert( (this->point(0))(2) == (this->point(1))(2) );
  genius_assert( (this->point(0))(2) == (this->point(2))(2) );

  Real xa = (this->point(0))(0);
  Real xb = (this->point(1))(0);
  Real xc = (this->point(2))(0);
  Real ya = (this->point(0))(1);
  Real yb = (this->point(1))(1);
  Real yc = (this->point(2))(1);

  Complex dx = ((yb-yc)*var[0] + (yc-ya)*var[1] +(ya-yb)*var[2])/(2*vol);
  Complex dy = ((xc-xb)*var[0] + (xa-xc)*var[1] +(xb-xa)*var[2])/(2*vol);
  return VectorValue<Complex>(dx, dy, 0.0);
}


VectorValue<AutoDScalar> Tri3_FVM::gradient( const std::vector<AutoDScalar> & var) const
{
  // FIXME we assume the triangle on xy plane here. not a general case
  genius_assert( (this->point(0))(2) == (this->point(1))(2) );
  genius_assert( (this->point(0))(2) == (this->point(2))(2) );

  Real xa = (this->point(0))(0);
  Real xb = (this->point(1))(0);
  Real xc = (this->point(2))(0);
  Real ya = (this->point(0))(1);
  Real yb = (this->point(1))(1);
  Real yc = (this->point(2))(1);

  AutoDScalar dx = ((yb-yc)*var[0] + (yc-ya)*var[1] +(ya-yb)*var[2])/(2*vol);
  AutoDScalar dy = ((xc-xb)*var[0] + (xa-xc)*var[1] +(xb-xa)*var[2])/(2*vol);
  return VectorValue<AutoDScalar>(dx, dy, 0.0);
}


void Tri3_FVM::prepare_for_vector_reconstruct()
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


VectorValue<PetscScalar> Tri3_FVM::reconstruct_vector( const std::vector<PetscScalar> & projects) const
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

VectorValue<AutoDScalar> Tri3_FVM::reconstruct_vector( const std::vector<AutoDScalar> & projects) const
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


