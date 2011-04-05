// $Id: face_quad4.cc,v 1.8 2008/05/24 12:32:04 gdiso Exp $

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
#include "face_quad4.h"




// ------------------------------------------------------------
// Quad class static member initialization
const unsigned int Quad4::side_nodes_map[4][2] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {2, 3}, // Side 2
    {3, 0}  // Side 3
  };

const unsigned int Quad4::node_connect_graph[4][4] =
  {
    //0  1  2  3
    {1, 1, 0, 1}, //0
    {1, 1, 1, 0}, //1
    {0, 1, 1, 1}, //2
    {3, 0, 1, 1}  //3
  };


#ifdef ENABLE_AMR

const float Quad4::_embedding_matrix[4][4][4] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2    3
      {1.0, 0.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0, 0.0}, // 1
      {.25, .25, .25, .25}, // 2
      {0.5, 0.0, 0.0, 0.5}  // 3
    },

    // embedding matrix for child 1
    {
      // 0    1    2    3
      {0.5, 0.5, 0.0, 0.0}, // 0
      {0.0, 1.0, 0.0, 0.0}, // 1
      {0.0, 0.5, 0.5, 0.0}, // 2
      {.25, .25, .25, .25}  // 3
    },

    // embedding matrix for child 2
    {
      // 0    1    2    3
      {.25, .25, .25, .25}, // 0
      {0.0, 0.5, 0.5, 0.0}, // 1
      {0.0, 0.0, 1.0, 0.0}, // 2
      {0.0, 0.0, 0.5, 0.5}  // 3
    },

    // embedding matrix for child 3
    {
      // 0    1    2    3
      {0.5, 0.0, 0.0, 0.5}, // 0
      {.25, .25, .25, .25}, // 1
      {0.0, 0.0, 0.5, 0.5}, // 2
      {0.0, 0.0, 0.0, 1.0}  // 3
    }

  };

#endif





// ------------------------------------------------------------
// Quad4 class member functions

bool Quad4::is_vertex(const unsigned int) const
{
  return true;
}

bool Quad4::is_edge(const unsigned int) const
{
  return false;
}

bool Quad4::is_face(const unsigned int) const
{
  return false;
}

bool Quad4::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  genius_assert(s < n_sides());
  for (unsigned int i = 0; i != 2; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}


void Quad4::nodes_on_edge (const unsigned int e,
                           std::vector<unsigned int> & nodes ) const
{
  genius_assert(e<4);

  nodes.clear();
  nodes.push_back(side_nodes_map[e][0]);
  nodes.push_back(side_nodes_map[e][1]);
}


void Quad4::nodes_on_edge (const unsigned int e,
                          std::pair<unsigned int, unsigned int> & nodes ) const
{
  nodes.first = side_nodes_map[e][0];
  nodes.second = side_nodes_map[e][1];
}


Real Quad4::edge_length(const unsigned int e) const
{
  return (point(side_nodes_map[e][0]) - point(side_nodes_map[e][1])).size();
}


bool Quad4::has_affine_map() const
{
  Point v = this->point(3) - this->point(0);
  return (v.relative_fuzzy_equals(this->point(2) - this->point(1)));
}



AutoPtr<Elem> Quad4::build_side (const unsigned int i,
                                 bool proxy) const
{
  assert (i < this->n_sides());

  if (proxy)
  {
    AutoPtr<Elem> ap(new Side<Edge2,Quad4>(this,i));
    return ap;
  }

  else
  {

    switch (i)
    {
    case 0:
      {
        Edge2* edge = new Edge2;

        edge->set_node(0) = this->get_node(0);
        edge->set_node(1) = this->get_node(1);
        edge->subdomain_id() = this->subdomain_id();
        AutoPtr<Elem> ap(edge);  return ap;
      }
    case 1:
      {
        Edge2* edge = new Edge2;

        edge->set_node(0) = this->get_node(1);
        edge->set_node(1) = this->get_node(2);
        edge->subdomain_id() = this->subdomain_id();
        AutoPtr<Elem> ap(edge);  return ap;
      }
    case 2:
      {
        Edge2* edge = new Edge2;

        edge->set_node(0) = this->get_node(2);
        edge->set_node(1) = this->get_node(3);
        edge->subdomain_id() = this->subdomain_id();
        AutoPtr<Elem> ap(edge);  return ap;
      }
    case 3:
      {
        Edge2* edge = new Edge2;

        edge->set_node(0) = this->get_node(3);
        edge->set_node(1) = this->get_node(0);
        edge->subdomain_id() = this->subdomain_id();
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



const Node * Quad4::is_hanging_node_on_side(unsigned int s) const
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



bool Quad4::in_plane() const
{
  // The A,B,C,D naming scheme here corresponds exactly to the
  // libmesh counter-clockwise numbering scheme.

  //        3           2        D           C
  // QUAD4: o-----------o        o-----------o
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        o-----------o        o-----------o
  //        0           1        A           B

  // Vector pointing from A to C
  Point AC ( this->point(2) - this->point(0) );

  // Vector pointing from A to B
  Point AB ( this->point(1) - this->point(0) );

  // Vector pointing from A to D
  Point AD ( this->point(3) - this->point(0) );

  // Vector normal to plan ABC
  Point n1 = AB.cross(AC);

  // Vector normal to plan ACD
  Point n2 = AC.cross(AD);

  // n1 cross n2
  Point n = n1.cross(n2);

  return (n == Point(0.,0.,0.));
}


bool Quad4::in_circle() const
{
  // make sure the 4 points co-plane first.
  if (in_plane() == false) return false;

  // The A,B,C,D naming scheme here corresponds exactly to the
  // libmesh counter-clockwise numbering scheme.

  //        3           2        D           C
  // QUAD4: o-----------o        o-----------o
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        o-----------o        o-----------o
  //        0           1        A           B

  Point circumcircle_center;
  Real  radii;
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
    radii = (circumcircle_center - this->point(0)).size();
  }

  return ( std::abs((circumcircle_center - this->point(3)).size() - radii) < TOLERANCE ) ;

}


void Quad4::connectivity(const unsigned int sf,
                         const IOPackage iop,
                         std::vector<unsigned int>& conn) const
{
  assert (sf < this->n_sub_elem());
  assert (iop != INVALID_IO_PACKAGE);

  // Create storage.
  conn.resize(4);

  switch (iop)
  {
  case TECPLOT:
    {
      conn[0] = this->node(0)+1;
      conn[1] = this->node(1)+1;
      conn[2] = this->node(2)+1;
      conn[3] = this->node(3)+1;
      return;
    }

  case VTK:
    {
      conn[0] = this->node(0);
      conn[1] = this->node(1);
      conn[2] = this->node(2);
      conn[3] = this->node(3);
      return;
    }

  default:
    genius_error();
  }

  genius_error();
}


void Quad4::side_order( const IOPackage iop, std::vector<unsigned int>& order) const
{
  order.resize(4);
  switch (iop)
  {
    case ISE:
    {
      order[0] = 0;
      order[1] = 1;
      order[2] = 2;
      order[3] = 3;
      return;
    }
    default:
      genius_error();
  }
}



Real Quad4::volume () const
{
  // The A,B,C,D naming scheme here corresponds exactly to the
  // libmesh counter-clockwise numbering scheme.

  //        3           2        D           C
  // QUAD4: o-----------o        o-----------o
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        o-----------o        o-----------o
  //        0           1        A           B

  // Vector pointing from A to C
  Point AC ( this->point(2) - this->point(0) );

  // Vector pointing from A to B
  Point AB ( this->point(1) - this->point(0) );

  // Vector pointing from A to D
  Point AD ( this->point(3) - this->point(0) );

  // The diagonal vector minus the side vectors
  Point AC_AB_AD (AC - AB - AD);

  // Check for quick return for planar QUAD4.  This will
  // be the most common case, occuring for all purely 2D meshes.
  if (AC_AB_AD == Point(0.,0.,0.))
    return AB.cross(AD).size();

  else
  {
    // Use 2x2 quadrature to approximate the surface area.  (The
    // true integral is too difficult to compute analytically.)  The
    // accuracy here is exactly the same as would be obtained via a
    // call to Elem::volume(), however it is a bit more optimized to
    // do it this way.  The technique used is to integrate the magnitude
    // of the normal vector over the whole area.  See for example,
    //
    // Y. Zhang, C. Bajaj, G. Xu. Surface Smoothing and Quality
    // Improvement of Quadrilateral/Hexahedral Meshes with Geometric
    // Flow. The special issue of the Journal Communications in
    // Numerical Methods in Engineering (CNME), submitted as an
    // invited paper, 2006.
    // http://www.ices.utexas.edu/~jessica/paper/quadhexgf/quadhex_geomflow_CNM.pdf

    // 4-point rule
    const Real q[2] =
      {
        0.5 - std::sqrt(3.) / 6.,
        0.5 + std::sqrt(3.) / 6.
      };

    Real vol=0.;
    for (unsigned int i=0; i<2; ++i)
      for (unsigned int j=0; j<2; ++j)
        vol += (AB + q[i]*AC_AB_AD).cross(AD + q[j]*AC_AB_AD).size();

    return 0.25*vol;
  }
}


bool Quad4::contains_point (const Point& p) const
{
  //volume method for judging point p is in this Quad4 element or not
  //A Quad4 element is splitted into 4 triangulars because of point p
  static const unsigned char sub_tri[4][2] =
    {
      {0, 1},
      {1, 2},
      {2, 3},
      {3, 0}
    };

  Real Quad4_volume = this->volume();
  Real pt_volume = 0.0;

  //temporary storage for Nodes which form the base of the subelements
  Node* base[2];

  //Add upp the sub triangulars volumes
  for (unsigned int n=0; n<4; ++n)
  {
    // Set the nodes of the triangulars base
    for (unsigned int i=0; i<2; ++i)
      base[i] = this->get_node(sub_tri[n][i]);

    //computing diff vectors
    Point v0p ( p - *base[0] );
    Point v21 ( *base[1] - *base[0] );

    //sub volume of this small triangular
    Real sub_volume = 0.5 * (v0p.cross(v21)).size();

    assert (sub_volume >= 0.);

    pt_volume += sub_volume;
  }

  if (fabs(Quad4_volume-pt_volume)<=1.0e-6*Quad4_volume)
    return true;
  else
    return false;
}



void Quad4::ray_hit(const Point &p , const Point &d , IntersectionResult &result, unsigned int dim) const
{
  // the ray and the quad is co-plane
  if (std::abs((this->point(1) - this->point(0)).dot(d.cross(this->point(2) - this->point(0))))<1e-10)
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

        Hit_Point & hit_point = edge_intersections.hit_points[0];
        if(hit_point.point_location == on_edge)
        {
          //for 2D mesh, on edge equals to on side
          if(dim==2) hit_point.point_location = on_side;
          hit_point.mark = s;
          if(!result.is_point_overlap_exist(hit_point))
            result.hit_points.push_back(hit_point);
        }
        else if(hit_point.point_location == on_vertex)
        {
          hit_point.mark = this->side_nodes_map[s][hit_point.mark];
          if(!result.is_mark_exist(on_vertex, hit_point.mark) && !result.is_point_overlap_exist(hit_point))
            result.hit_points.push_back(hit_point);
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

  // on triangle 0,1,2
  {

    Point e1 = this->point(1) - this->point(0);
    Point e2 = this->point(2) - this->point(0);

    Point h = d.cross(e2);

    double a = e1.dot(h);

    double f = 1/a;
    Point s = p - this->point(0);
    double u = f * (s.dot(h));
    Point q = s.cross(e1);
    double v = f * d.dot(q);

    if (u > -1e-10 && u < 1.0+1e-10 && v > -1e-10 && u + v < 1.0+1e-10)
    {
      double t = f * e2.dot(q);

      // this means that there is a line intersection
      // but not a ray intersection
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

  }

  // on triangle 0 2 3
  {

    Point e1 = this->point(2) - this->point(0);
    Point e2 = this->point(3) - this->point(0);

    Point h = d.cross(e2);

    double a = e1.dot(h);

    double f = 1/a;
    Point s = p - this->point(0);
    double u = f * (s.dot(h));
    Point q = s.cross(e1);
    double v = f * d.dot(q);

    if (u > -1e-10 && u < 1.0+1e-10 && v > -1e-10 && u + v < 1.0+1e-10)
    {
      double t = f * e2.dot(q);

      // this means that there is a line intersection
      // but not a ray intersection
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
        result.hit_points[0].p = this->point(2);
        result.hit_points[0].mark = 2;
        return;
      }

      if(std::abs(1-v)<1e-10) //
      {
        result.hit_points[0].point_location = on_vertex;
        result.hit_points[0].p = this->point(3);
        result.hit_points[0].mark = 3;
        return;
      }

      if(std::abs(1-u-v)<1e-10)
      {
        result.hit_points[0].point_location = on_edge;
        result.hit_points[0].mark = 2;
        return;
      }

      if(std::abs(u)<1e-10)
      {
        result.hit_points[0].point_location = on_edge;
        result.hit_points[0].mark = 3;
        return;
      }

      result.hit_points[0].point_location = on_face;
      result.hit_points[0].mark = 0;
      return;
    }

  }

  result.state = Missed;
  return;
}


Point Quad4::nearest_point(const Point &p, Real * dist) const
{
  // Just use p0 for the point.
  const Point plane_point = this->point(0);
  const Point e0 = this->point(1) - this->point(0);
  const Point e1 = this->point(2) - this->point(0);
  const Point plane_normal = e0.cross(e1).unit();

  // Create a vector from the surface to point p;
  const Point w = p - plane_point;

  // The closest point in the plane to point p
  // is in the negative normal direction
  // a distance w (dot) p.
  const Point cp = p - plane_normal*(w*plane_normal);

  if(this->contains_point(cp))
  {
    if(dist) *dist = (cp-p).size();
    return cp;
  }

  Point np;
  Real distance = 1e30;
  for(unsigned int s=0; s<n_sides(); ++s)
  {
    Real d;
    AutoPtr<Elem> side = Quad4::build_side(s, false);
    Point n = side->nearest_point(p, &d);
    if( d < distance )
    {
      np = n;
      distance = d;
    }
  }

  if(dist) *dist = distance;
  return np;
}


Point Quad4::outside_unit_normal(unsigned short int side_id) const
{
  genius_assert(side_id < this->n_sides());

  Point *p1 = this->get_node(Quad4::side_nodes_map[side_id][0]);
  Point *p2 = this->get_node(Quad4::side_nodes_map[side_id][1]);
  Point p3 = this->centroid();

  Point cross_3 = (*p2 - *p1).cross(p3 - *p2);
  return   ((*p2 - *p1).cross(cross_3)).unit();
}
