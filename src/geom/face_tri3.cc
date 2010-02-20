// $Id: face_tri3.cc,v 1.4 2008/05/24 12:32:04 gdiso Exp $

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
#include <map>

// Local includes
#include "side.h"
#include "edge_edge2.h"
#include "face_tri3.h"



// ------------------------------------------------------------
// Tri3 class static member initializations
const unsigned int Tri3::side_nodes_map[3][2] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {2, 0}  // Side 2
  };

const unsigned int Tri3::node_connect_graph[3][3] =
  {
    // the 3 node is linked to each other
    {1, 1, 1},
    {1, 1, 1},
    {1, 1, 1}
  };

#ifdef ENABLE_AMR

const float Tri3::_embedding_matrix[4][3][3] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2
      {1.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0}, // 1
      {0.5, 0.0, 0.5}  // 2
    },

    // embedding matrix for child 1
    {
      // 0    1    2
      {0.5, 0.5, 0.0}, // 0
      {0.0, 1.0, 0.0}, // 1
      {0.0, 0.5, 0.5}  // 2
    },

    // embedding matrix for child 2
    {
      // 0    1    2
      {0.5, 0.0, 0.5}, // 0
      {0.0, 0.5, 0.5}, // 1
      {0.0, 0.0, 1.0}  // 2
    },

    // embedding matrix for child 3
    {
      // 0    1    2
      {0.5, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.5}, // 1
      {0.5, 0.0, 0.5}  // 2
    }
  };

#endif



// ------------------------------------------------------------
// Tri3 class member functions

bool Tri3::is_vertex(const unsigned int) const
{
  return true;
}

bool Tri3::is_edge(const unsigned int) const
{
  return false;
}

bool Tri3::is_face(const unsigned int) const
{
  return false;
}

bool Tri3::is_node_on_side(const unsigned int n,
                           const unsigned int s) const
{
  genius_assert(s < n_sides());
  for (unsigned int i = 0; i != 2; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}


void Tri3::nodes_on_edge (const unsigned int e,
                          std::vector<unsigned int> & nodes ) const
{
  genius_assert(e<3);

  nodes.clear();
  nodes.push_back(side_nodes_map[e][0]);
  nodes.push_back(side_nodes_map[e][1]);
}


AutoPtr<Elem> Tri3::build_side (const unsigned int i,
                                bool proxy) const
{
  assert (i < this->n_sides());

  if (proxy)
  {
    AutoPtr<Elem> ap(new Side<Edge2,Tri3>(this,i));
    return ap;
  }

  else
  {
    Edge2* edge = new Edge2;

    switch (i)
    {
    case 0:
      {
        edge->set_node(0) = this->get_node(0);
        edge->set_node(1) = this->get_node(1);

        AutoPtr<Elem> ap(edge);  return ap;
      }
    case 1:
      {
        edge->set_node(0) = this->get_node(1);
        edge->set_node(1) = this->get_node(2);

        AutoPtr<Elem> ap(edge);  return ap;
      }
    case 2:
      {
        edge->set_node(0) = this->get_node(2);
        edge->set_node(1) = this->get_node(0);

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



const Node * Tri3::is_hanging_node_on_side(unsigned int s) const
{

  // Only constrain active and ancestor elements
  if (this->subactive())
    return NULL;

  // a boundary edge ? no hanging node
  const Elem * neighbor_elem = this->neighbor(s);

  if( neighbor_elem == NULL )    return NULL;

  // the active neighbor elem is equal or higher than me, no hanging node exist
  if( neighbor_elem->active() && neighbor_elem->level() <= this->level() ) return NULL;

  // neighbor_elem->neighbor(this_side) == this
  unsigned int this_side = neighbor_elem->which_neighbor_am_i(this);


  std::set<const Node *> hanging_nodes;

  // for 2D elem, the side is identical to edge
  AutoPtr<Elem> edge = this->build_side(s);


  std::vector<const Elem*> family;
  neighbor_elem->active_family_tree(family);


  for(unsigned int el=0; el<family.size(); ++el)
  {
    genius_assert( family[el]->level() > this->level() );

    AutoPtr<Elem> child_edge = family[el]->build_side(this_side);

    if(child_edge->node(0) == edge->node(0) ||  child_edge->node(0) == edge->node(1))
      hanging_nodes.insert(child_edge->get_node(1));

    if(child_edge->node(1) == edge->node(0) ||  child_edge->node(1) == edge->node(1))
      hanging_nodes.insert(child_edge->get_node(0));
  }

  // should have only one hanging node
  genius_assert(hanging_nodes.size()==1);

  return *(hanging_nodes.begin());

}



void Tri3::connectivity(const unsigned int sf,
                        const IOPackage iop,
                        std::vector<unsigned int>& conn) const
{
  assert (sf <this->n_sub_elem());
  assert (iop != INVALID_IO_PACKAGE);

  switch (iop)
  {
  case TECPLOT:
    {
      conn.resize(4);
      conn[0] = this->node(0)+1;
      conn[1] = this->node(1)+1;
      conn[2] = this->node(2)+1;
      conn[3] = this->node(2)+1;
      return;
    }

  case VTK:
    {
      conn.resize(3);
      conn[0] = this->node(0);
      conn[1] = this->node(1);
      conn[2] = this->node(2);
      return;
    }

  default:
    genius_error();
  }

  genius_error();
}


Real Tri3::volume () const
{
  // 3-node triangles have the following formula for computing the area
  Point v10 ( *(this->get_node(1)) - *(this->get_node(0)) );

  Point v20 ( *(this->get_node(2)) - *(this->get_node(0)) );

  return 0.5 * (v10.cross(v20)).size() ;
}



bool Tri3::contains_point (const Point & p) const
{
  //volume method to judging point p is in this 3-node triangles or not
  //a Tri3 element are splitted into 3 triangulars because of point p
  static const unsigned char sub_tri[3][2] =
    {
      {0, 1},
      {1, 2},
      {2, 0}
    };

  Real tri_volume = this->volume();
  Real pt_volume = 0.0;

  //temporary storage for Nodes which form the base of the subelements
  std::vector<const Node*> base(2);

  //sum the sub triangulars volumes
  for (unsigned int n=0; n<3; ++n)
  {
    // Set the nodes of the triangulars base
    for (unsigned int i=0; i<2; ++i)
      base[i] = this->get_node(sub_tri[n][i]);

    //computing diff vectors
    Point v0p ( p - *base[0] );
    Point v21 ( *base[1] - *base[0] );

    //sub volume of this small triangular
    Real sub_volume = 0.5 * (v0p.cross(v21)).size();

    assert (sub_volume >= 0.0);

    pt_volume += sub_volume;
  }

  if ( std::abs(tri_volume-pt_volume) < 1.0e-6*tri_volume )
    return true;
  else
    return false;
}


void Tri3::ray_hit(const Point &p , const Point &d , IntersectionResult &result, unsigned int dim) const
{
  Point e1 = this->point(1) - this->point(0);
  Point e2 = this->point(2) - this->point(0);

  Point h = d.cross(e2);

  double a = e1.dot(h);

  // the ray and the triangle is co-plane
  if (std::abs(a)<1e-10)
  {
    IntersectionResult edge_intersections;

    for(unsigned int s=0; s<this->n_sides(); ++s)
    {
      AutoPtr<Elem> edge = this->build_side(s);
      edge->ray_hit(p, d, edge_intersections);

      if(edge_intersections.state==Missed)
        continue;

      if(edge_intersections.state==Intersect_Body)
      {
        if(dim==2)
          result.state=Intersect_Body;
        else
          result.state=On_Face;

        if(edge_intersections.hit_points[0].point_location == on_edge)
        {
          //for 2D mesh, on edge equals to on side
          if(dim==2) edge_intersections.hit_points[0].point_location = on_side;
          edge_intersections.hit_points[0].mark = s;
          result.hit_points.push_back(edge_intersections.hit_points[0]);
        }
        else if(edge_intersections.hit_points[0].point_location == on_vertex)
        {
          edge_intersections.hit_points[0].mark = this->side_nodes_map[s][edge_intersections.hit_points[0].mark];
          if(!result.is_mark_exist(on_vertex, edge_intersections.hit_points[0].mark))
            result.hit_points.push_back(edge_intersections.hit_points[0]);
        }
      }

      if(edge_intersections.state==Overlap_Edge)
      {
        // collect result
        result = edge_intersections;
        result.mark = s;
        result.hit_points[0].mark = this->side_nodes_map[s][result.hit_points[0].mark];
        result.hit_points[1].mark = this->side_nodes_map[s][result.hit_points[1].mark];
        return;
      }
    }

    if(result.hit_points.size()==0)
    {
      result.state = Missed;
      return;
    }

    if(result.hit_points.size()==1)
    {
      return;
    }

    genius_assert(result.hit_points.size()==2);
    if(result.hit_points[0].t > result.hit_points[1].t)
      std::swap(result.hit_points[0], result.hit_points[1]);

    return;
  }

  double f = 1/a;
  Point s = p - this->point(0);

  double u = f * (s.dot(h));
  if (u < -1e-10 || u > 1.0+1e-10)
  {
    result.state = Missed;
    return;
  }

  Point q = s.cross(e1);

  double v = f * d.dot(q);
  if (v < -1e-10 || u + v > 1.0+1e-10)
  {
    result.state = Missed;
    return;
  }

  // at this stage we can compute t to find out where
  // the intersection point is on the line
  double t = f * e2.dot(q);

  // this means that there is a line intersection but not a ray intersection
  // however, the ray start point on elem is ok
  if ( t < -1e-10 )
  {
    result.state = Missed;
    return;
  }

  result.state=Intersect_Body;
  result.hit_points.resize(1);
  result.hit_points[0].p = p+t*d;
  result.hit_points[0].t = t;

  // the intersection point location
  if(std::abs(u+v)<1e-10) //
  {
    result.hit_points[0].point_location = on_vertex;
    result.hit_points[0].p = this->point(0);
    result.hit_points[0].mark = 0;
    return;
  }

  if(std::abs(1-u)<1e-10) //
  {
    result.hit_points[0].point_location = on_vertex;
    result.hit_points[0].p = this->point(1);
    result.hit_points[0].mark = 1;
    return;
  }

  if(std::abs(1-v)<1e-10) //
  {
    result.hit_points[0].point_location = on_vertex;
    result.hit_points[0].p = this->point(2);
    result.hit_points[0].mark = 2;
    return;
  }

  if(std::abs(1-u-v)<1e-10)
  {
    result.hit_points[0].point_location = on_edge;
    result.hit_points[0].mark = 1;
    return;
  }

  if(std::abs(u)<1e-10)
  {
    result.hit_points[0].point_location = on_edge;
    result.hit_points[0].mark = 2;
    return;
  }

  if(std::abs(v)<1e-10)
  {
    result.hit_points[0].point_location = on_edge;
    result.hit_points[0].mark = 0;
    return;
  }

  result.hit_points[0].point_location = on_face;
  result.hit_points[0].mark = 0;
  return;

}



Point Tri3::outside_unit_normal(unsigned short int side_id) const
{
  genius_assert(side_id < this->n_sides());

  Point *p1 = this->get_node(Tri3::side_nodes_map[side_id][0]);
  Point *p2 = this->get_node(Tri3::side_nodes_map[side_id][1]);
  Point p3 = this->centroid();

  Point cross_3 = (*p2 - *p1).cross(p3 - *p2);
  return   ((*p2 - *p1).cross(cross_3)).unit();
}


std::pair<Real, Real> Tri3::min_and_max_angle() const
{
  Point v10 ( this->point(1) - this->point(0) );
  Point v20 ( this->point(2) - this->point(0) );
  Point v21 ( this->point(2) - this->point(1) );

  const Real
  len_10=v10.size(),
         len_20=v20.size(),
                len_21=v21.size()
                       ;

  const Real
  theta0=std::acos(( v10*v20)/len_10/len_20),
         theta1=std::acos((-v10*v21)/len_10/len_21),
                theta2=M_PI - theta0 - theta1
                       ;

  genius_assert(theta0 > 0.);
  genius_assert(theta1 > 0.);
  genius_assert(theta2 > 0.);

  return std::make_pair(std::min(theta0, std::min(theta1,theta2)),
                        std::max(theta0, std::max(theta1,theta2)));
}
