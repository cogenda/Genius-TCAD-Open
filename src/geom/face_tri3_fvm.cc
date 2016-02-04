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
#include "edge_edge2_fvm.h"

// TNT matrix-vector library
#include <TNT/tnt.h>
#include <TNT/jama_lu.h>

//#define DEBUG

/*
 *   TRI3:  2               TRI3:  C
 *          o                      o
 *         / \                    / \
 *        /   \                  /   \
 *    e2 /     \ e1             /     \
 *      /       \              /       \
 *     /    e0   \            /         \
 *    o-----------o          o-----------o
 *    0           1          A           B
 */

AutoPtr<Elem> Tri3_FVM::build_fvm_side (const unsigned int i, bool proxy) const
{
  assert (i < this->n_sides());

  if (proxy)
  {
    return Tri3::build_side(i, proxy);
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
        edge->set_node(1) = this->get_node(0);
        edge->subdomain_id() = this->subdomain_id();
        edge->hold_fvm_node(0, this->get_fvm_node(2));
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

  // We will never get here...  Look at the code above.
  genius_error();
  AutoPtr<Elem> ap(NULL);  return ap;
}


// build the Geom information here
void Tri3_FVM::prepare_for_fvm()
{
  //store this vlaue for efficiency reason
  vol = Tri3::volume();

  // clear partitial volume
  v[0] = v[1] = v[2] = 0;
  vt[0] = vt[1] = vt[2] = 0;

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
  unsigned int obtuse_edge = invalid_uint;
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
    {
      d[i] = -(side_centers[i] - circumcircle_center).size();
      obtuse_edge = i; //we need special process to obtuse angle
    }
    else
      d[i] =  (side_centers[i] - circumcircle_center).size();

    dt[i] = d[i];
  }

  // special process to obtuse angle: Truncated negative length
  if(obtuse_edge!=invalid_uint)
  {
   /*
    *     TRI3:              p3
    *                        o
    *                     *     *
    *                  *          *
    *               *               *
    *            *   \            /   *
    *         * a1    \          /    a2*
    *      o-----------o--------o---------o
    *      p1          m1      m2         p2
    */
    unsigned int obtuse_node = (2+obtuse_edge)%3;
    Point p1 =  this->point(side_nodes_map[obtuse_edge][0]) ;
    Point p2 =  this->point(side_nodes_map[obtuse_edge][1]) ;
    Point p3 =  this->point(obtuse_node) ;
    Real  a1 = (p1-p2).angle(p1-p3);
    Real  a2 = (p2-p1).angle(p2-p3);

    Point pre_edge_center = 0.5*(p1 + p3);
    Point pos_edge_center = 0.5*(p2 + p3);
    Point m1 = p1 + (p2-p1).unit()*(pre_edge_center-p1).size()/cos(a1);
    Point m2 = p2 + (p1-p2).unit()*(pos_edge_center-p2).size()/cos(a2);

    unsigned int pre_edge = (obtuse_edge + 3 - 1)%3;
    unsigned int pos_edge = (obtuse_edge + 3 + 1)%3;

    dt[obtuse_edge] = 0;
    dt[pre_edge] = (pre_edge_center-m1).size();
    dt[pos_edge] = (pos_edge_center-m2).size();
  }


  // compute partial volume
  for( unsigned int i=0; i<3; i++ ) // has 3 side (edge)
  {
    unsigned int node1 = side_nodes_map[i][0];
    unsigned int node2 = side_nodes_map[i][1];
    v[node1] +=  0.5 * 0.5 * l[i] * d[i];
    v[node2] +=  0.5 * 0.5 * l[i] * d[i];
    vt[node1] +=  0.5 * 0.5 * l[i] * dt[i];
    vt[node2] +=  0.5 * 0.5 * l[i] * dt[i];
  }
  
  // truncated partial volum consistant with element volume
  if(obtuse_edge!=invalid_uint)
  {
    unsigned int obtuse_node = (2+obtuse_edge)%3;
    unsigned int n1 = side_nodes_map[obtuse_edge][0];
    unsigned int n2 = side_nodes_map[obtuse_edge][1];
    vt[obtuse_node] = vol - vt[n1] - vt[n2];
  }
  
  


  // self test
#ifdef DEBUG
  Real V = 0;
  for( unsigned int i=0; i<3; i++ )
    V += v[i];
  genius_assert( (vol > 1e-10 ? (std::abs(V-vol) < 1e-3*vol) : (std::abs(V-vol) < std::max(1e-13, 1e-2*vol))) );
#endif

  prepare_for_vector_reconstruct();
}



#if 0
// build the Geom information here
void Tri3_FVM::prepare_for_fvm()
{
  //store this vlaue for efficiency reason
  vol = Tri3::volume();

  // clear partitial volume
  v[0] = v[1] = v[2] = 0;

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
  unsigned int obtuse_edge = invalid_uint;
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
      d[i] = -(side_centers[i] - circumcircle_center).size(); //we need special process to obtuse angle
    else
      d[i] =  (side_centers[i] - circumcircle_center).size();
  }


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
    V += v[i];
  genius_assert( vol > 1e-10 ? (std::abs(V-vol) < 1e-3*vol) : (std::abs(V-vol) < 1e-2*vol) );


  prepare_for_vector_reconstruct();
}
#endif


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
   //FIXME NOT work for 3D triangle
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

   if( solver.isNonsingular() )
   {
     TNT::Array2D<Real> inv_ATA = solver.inv();
     TNT::Array2D<Real>  M = TNT::matmult(inv_ATA, AT);

     for(unsigned int m=0; m<M.dim1(); m++)
       for( unsigned int e=0; e<M.dim2(); e++ )
         least_squares_vector_reconstruct_matrix[m][e] = M[m][e];
   }
   else
   {
     for(unsigned int m=0; m<2; m++)
       for( unsigned int e=0; e<3; e++ )
         least_squares_vector_reconstruct_matrix[m][e] = 0.0;
   }
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

