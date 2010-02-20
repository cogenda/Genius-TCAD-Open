// $Id: sphere.cc,v 1.2 2008/04/18 14:10:43 gdiso Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// C++ includes
#include <cmath> // for std::abs


// Local includes
#include "sphere.h"
#include "tensor_value.h"


// ------------------------------------------------------------
// Sphere class member functions
Sphere::Sphere () :
  _rad(-1.)
{
}



Sphere::Sphere (const Point& c,
        const Real   r)
{
  assert (r > 0.);

  this->create_from_center_radius (c, r);
}


Sphere::Sphere (const Point& pa, const Point& pb, const Point& pc, const Point& pd)
{
  Point pad = pa - pd;
  Point pbd = pb - pd;
  Point pcd = pc - pd;

  TensorValue<Real> T(pad,pbd,pcd);

  // the det
  Real D = T.det();
  assert (std::abs(D) > 0);

  Real e = 0.5*(pa.size_sq() - pd.size_sq());
  Real f = 0.5*(pb.size_sq() - pd.size_sq());
  Real g = 0.5*(pc.size_sq() - pd.size_sq());

  TensorValue<Real> T1(e,pad(1),pad(2),
                       f,pbd(1),pbd(2),
                       g,pcd(1),pcd(2));
  Real  sx = T1.det()/D;

  TensorValue<Real> T2(pad(0),e,pad(2),
                       pbd(0),f,pbd(2),
                       pcd(0),g,pcd(2));
  Real  sy = T2.det()/D;

  TensorValue<Real> T3(pad(0),pad(1),e,
                       pbd(0),pbd(1),f,
                       pcd(0),pcd(1),g);
  Real  sz = T3.det()/D;

  Point c(sx,sy,sz);
  Real  r = (c-pa).size();

  this->create_from_center_radius (c, r);

/*
    double arrayD[3][3]={p1.x-p4.x,p1.y-p4.y,p1.z-p4.z,
                 p2.x-p4.x,p2.y-p4.y,p2.z-p4.z,
                 p3.x-p4.x,p3.y-p4.y,p3.z-p4.z};
    D=det3(arrayD);
    E=(pow(p1.x,2)-pow(p4.x,2)+pow(p1.y,2)-pow(p4.y,2)+pow(p1.z,2)-pow(p4.z,2))/2;
    F=(pow(p2.x,2)-pow(p4.x,2)+pow(p2.y,2)-pow(p4.y,2)+pow(p2.z,2)-pow(p4.z,2))/2;
    G=(pow(p3.x,2)-pow(p4.x,2)+pow(p3.y,2)-pow(p4.y,2)+pow(p3.z,2)-pow(p4.z,2))/2;

    double arrayPX[3][3]={E,p1.y-p4.y,p1.z-p4.z,
                  F,p2.y-p4.y,p2.z-p4.z,
                  G,p3.y-p4.y,p3.z-p4.z};
    sc.x=det3(arrayPX)/D;

    double arrayPY[3][3]={p1.x-p4.x,E,p1.z-p4.z,
                  p2.x-p4.x,F,p2.z-p4.z,
                  p3.x-p4.x,G,p3.z-p4.z};
    sc.y=det3(arrayPY)/D;

    double arrayPZ[3][3]={p1.x-p4.x,p1.y-p4.y,E,
                  p2.x-p4.x,p2.y-p4.y,F,
                  p3.x-p4.x,p3.y-p4.y,G};

    sc.z=det3(arrayPZ)/D;
    sr=sqrt(pow(p1.x-sc.x,2)+pow(p1.y-sc.y,2)+pow(p1.z-sc.z,2));
*/

}

Sphere::Sphere (const Sphere& other_sphere) :
  Surface()
{
  this->create_from_center_radius (other_sphere.center(),
                   other_sphere.radius());
}



Sphere::~Sphere ()
{
}



void Sphere::create_from_center_radius (const Point& c, const Real r)
{
  this->center() = c;
  this->radius() = r;

  assert (this->radius() > 0.);
}



bool Sphere::intersects (const Sphere& other_sphere) const
{
  assert ( this->radius() > 0. );
  assert ( other_sphere.radius() > 0. );

  const Real distance = (this->center() - other_sphere.center()).size();

  if (distance < (this->radius() + other_sphere.radius()) )
    return true;

  return false;
}



bool Sphere::above_surface (const Point& p) const
{
  assert (this->radius() > 0.);

  // create a vector from the center to the point.
  const Point w = p - this->center();

  if (w.size() > this->radius())
    return true;

  return false;
}



bool Sphere::below_surface (const Point& p) const
{
  assert (this->radius() > 0.);

  return ( !this->above_surface (p) );
}



bool Sphere::on_surface (const Point& p) const
{
  assert (this->radius() > 0.);

  // Create a vector from the center to the point.
  const Point w = p - this->center();

  // if the size of that vector is the same as the radius() then
  // the point is on the surface.
  if (std::abs(w.size() - this->radius()) < 1.e-10)
    return true;

  return false;
}



Point Sphere::closest_point (const Point& p) const
{
  assert (this->radius() > 0.);

  // get the normal from the surface in the direction
  // of p
  Point normal = this->unit_normal (p);

  // The closest point on the sphere is in the direction
  // of the normal a distance r from the center.
  const Point cp = this->center() + normal*this->radius();

  return cp;
}



Point Sphere::unit_normal (const Point& p) const
{
  assert (this->radius() > 0.);

  assert ( !(p == this->center()) );

  // Create a vector from the center to the point
  Point n = p - this->center();

  return n.unit();
}
