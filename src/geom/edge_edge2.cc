// $Id: edge_edge2.cc,v 1.5 2008/05/25 02:42:09 gdiso Exp $

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



// Local includes
#include "edge_edge2.h"
#include "tensor_value.h"

// ------------------------------------------------------------
// Edge2 class static member initializations
const unsigned int Edge2::node_connect_graph[2][2] =
  {
    // the 2 node is linked to each other
    {1, 1},
    {1, 1}
  };



#ifdef ENABLE_AMR

const float Edge2::_embedding_matrix[2][2][2] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2
      {1.0, 0.0}, // 0
      {0.5, 0.5}  // 1
    },

    // embedding matrix for child 1
    {
      // 0    1    2
      {0.5, 0.5}, // 0
      {0.0, 1.0}  // 1
    }
  };

#endif

bool Edge2::is_vertex(const unsigned int) const
{
  return true;
}

bool Edge2::is_edge(const unsigned int) const
{
  return false;
}

bool Edge2::is_face(const unsigned int) const
{
  return false;
}

bool Edge2::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  genius_assert(s < 2);
  return (s == n);
}

bool Edge2::is_node_on_edge(const unsigned int,
                            const unsigned int e) const
{
  genius_assert(e == 0);
  return true;
}

bool Edge2::is_edge_on_side(const unsigned int e,
                            const unsigned int s) const
{
  return true;
}

void Edge2::nodes_on_edge (const unsigned int ,
                           std::vector<unsigned int> & nodes ) const
{
  nodes.clear();
  nodes.push_back(0);
  nodes.push_back(1);
}



bool Edge2::contains_point (const Point& p) const
{
  Point p0 = this->point(0);
  Point p1 = this->point(1);

  Point v10 ( p - p0 );
  Point v20 ( p - p1 );

  // three point should co-line
  if ( 0.5*(v10.cross(v20)).size() > TOLERANCE ) return false;

  // p should lies on p0-p1
  if( std::abs(v10.size() + v20.size() - (p0-p1).size()) > TOLERANCE ) return false;

  return true;

}



void Edge2::connectivity(const unsigned int sc,
                         const IOPackage iop,
                         std::vector<unsigned int>& conn) const
{
  assert (sc == 0);
  assert (sc < this->n_sub_elem());
  assert (iop != INVALID_IO_PACKAGE);

  // Create storage
  conn.resize(2);

  switch (iop)
  {
  case TECPLOT:
    {
      conn[0] = this->node(0)+1;
      conn[1] = this->node(1)+1;
      return;
    }

  case VTK:
    {
      conn[0] = this->node(0);
      conn[1] = this->node(1);
      return;
    }

  default:
    {
      genius_error();
    }
  }

  genius_error();
}


void Edge2::ray_hit(const Point &p , const Point &d , IntersectionResult &result, unsigned int) const
{
  Point p0 = this->point(0);
  Point d0 = (this->point(1)-this->point(0)).unit();

  double h = (d.cross(d0)).size();

  // two lines are parallel
  if(h<1e-10)
  {
    Point p0_proj = p + (p0-p).dot(d)*d;
    double dist = (p0 - p0_proj).size();
    // two lines are overlap
    if(dist<1e-10)
    {
      result.state = Overlap_Edge;

      result.hit_points.resize(2);

      result.hit_points[0].p = this->point(0);
      result.hit_points[0].t = d.dot((this->point(0) - p));
      result.hit_points[0].point_location = on_vertex;
      result.hit_points[0].mark = 0;

      result.hit_points[1].p = this->point(1);
      result.hit_points[1].t = d.dot((this->point(1) - p));
      result.hit_points[1].point_location = on_vertex;
      result.hit_points[1].mark = 1;

      if(result.hit_points[0].t > result.hit_points[1].t)
        std::swap(result.hit_points[0], result.hit_points[1]);

      return;
    }
    else // no intersection
    {
      result.state = Missed;
      return;
    }
  }

  TensorValue<double> T0(p-p0, d,  d0.cross(d));
  TensorValue<double> T (p-p0, d0, d0.cross(d));

  double t0 = T0.det()/(h*h);
  double t  = T.det() /(h*h);

  Point near_point = p+t*d;
  double length = (this->point(1)-this->point(0)).size();

  // has one intersection
  if( ((p0+t0*d0) - (near_point)).size()<1e-10 && t0>-1e-10 && t0<length+1e-10 )
  {
    result.state = Intersect_Body;

    result.hit_points.resize(1);

    result.hit_points[0].p = near_point;
    result.hit_points[0].t = t;
    if( std::abs(t0) <1e-10)
    {
      result.hit_points[0].point_location = on_vertex;
      result.hit_points[0].p = this->point(0);
      result.hit_points[0].mark = 0;
    }
    else if( std::abs(t0-length)<1e-10)
    {
      result.hit_points[0].point_location = on_vertex;
      result.hit_points[0].p = this->point(1);
      result.hit_points[0].mark = 1;
    }
    else
    {
      result.hit_points[0].point_location = on_edge;
      result.hit_points[0].mark = 0;
    }
    return;
  }

  result.state = Missed;
  return;

}


Real Edge2::volume () const
{
  // OK, so this is probably overkill, since it is equivalent to
  // Elem::hmax() for the Edge2, but here it is nonetheless...
  return (this->point(1) - this->point(0)).size();
}
