// $Id: elem.cc,v 1.6 2008/05/25 02:42:09 gdiso Exp $

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
#include <numeric>
#include <algorithm> // for std::sort
#include <iterator>  // for std::ostream_iterator

// Local includes
#include "elem.h"
#include "edge_edge2.h"
#include "edge_edge2_fvm.h"
#include "face_tri3.h"
#include "face_tri3_fvm.h"
#include "face_quad4.h"
#include "face_quad4_fvm.h"

#ifdef COGENDA_COMMERCIAL_PRODUCT
  #include "cell_tet4.h"
  #include "cell_pyramid5.h"
  #include "cell_prism6.h"
  #include "cell_hex8.h"
  #include "cell_tet4_fvm.h"
  #include "cell_pyramid5_fvm.h"
  #include "cell_prism6_fvm.h"
  #include "cell_hex8_fvm.h"
#endif

#include "elem_clone.h"


// Initialize static member variables
const unsigned int Elem::_bp1 = 65449;
const unsigned int Elem::_bp2 = 48661;

// ------------------------------------------------------------
// Elem class member funcions
AutoPtr<Elem> Elem::build(const ElemType type,
                          Elem* p)
{
  Elem* elem = NULL;

  switch (type)
  {
      // 1D elements
      case EDGE2:
      {
        elem = new Edge2(p);
        break;
      }

      case EDGE2_FVM:
      {
        elem = new Edge2_FVM(p);
        break;
      }

      // 2D elements
      case TRI3:
      {
        elem = new Tri3(p);
        break;
      }

      case TRI3_FVM:
      {
        elem = new Tri3_FVM(p);
        break;
      }

      case QUAD4:
      {
        elem = new Quad4(p);
        break;
      }

      case QUAD4_FVM:
      {
        elem = new Quad4_FVM(p);
        break;
      }
#ifdef COGENDA_COMMERCIAL_PRODUCT
      // 3D elements
      case TET4:
      {
        elem = new Tet4(p);
        break;
      }
      case TET4_FVM:
      {
        elem = new Tet4_FVM(p);
        break;
      }

      case HEX8:
      {
        elem = new Hex8(p);
        break;
      }
      case HEX8_FVM:
      {
        elem = new Hex8_FVM(p);
        break;
      }

      case PRISM6:
      {
        elem = new Prism6(p);
        break;
      }
      case PRISM6_FVM:
      {
        elem = new Prism6_FVM(p);
        break;
      }

      case PYRAMID5:
      {
        elem = new Pyramid5(p);
        break;
      }
      case PYRAMID5_FVM:
      {
        elem = new Pyramid5_FVM(p);
        break;
      }

#else
      case TET4:
      case TET4_FVM:
      case HEX8:
      case HEX8_FVM:
      case PRISM6:
      case PRISM6_FVM:
      case PYRAMID5:
      case PYRAMID5_FVM:
      {
        std::cerr << "ERROR: 3D elements are not supported by Open Source Version." << std::endl;
        genius_error();
      }
#endif

      default:
      {
        std::cerr << "ERROR: Unsupported element type!." << std::endl;
        genius_error();
      }
  }

  AutoPtr<Elem> ap(elem);

  if (p) ap->subdomain_id() = p->subdomain_id();

  return ap;
}



AutoPtr<Elem> Elem::build_clone(const ElemType type, Elem* p)
{
  Elem* elem = NULL;

  switch (type)
  {
      // 1D elements
      case EDGE2:
      {
        elem = new ElemClone<Edge2>(p);
        break;
      }

      case EDGE2_FVM:
      {
        elem = new ElemClone<Edge2_FVM>(p);
        break;
      }

      // 2D elements
      case TRI3:
      {
        elem = new ElemClone<Tri3>(p);
        break;
      }

      case TRI3_FVM:
      {
        elem = new ElemClone<Tri3_FVM>(p);
        break;
      }

      case QUAD4:
      {
        elem = new ElemClone<Quad4>(p);
        break;
      }

      case QUAD4_FVM:
      {
        elem = new ElemClone<Quad4_FVM>(p);
        break;
      }
#ifdef COGENDA_COMMERCIAL_PRODUCT
      // 3D elements
      case TET4:
      {
        elem = new ElemClone<Tet4>(p);
        break;
      }
      case TET4_FVM:
      {
        elem = new ElemClone<Tet4_FVM>(p);
        break;
      }

      case HEX8:
      {
        elem = new ElemClone<Hex8>(p);
        break;
      }
      case HEX8_FVM:
      {
        elem = new ElemClone<Hex8_FVM>(p);
        break;
      }

      case PRISM6:
      {
        elem = new ElemClone<Prism6>(p);
        break;
      }
      case PRISM6_FVM:
      {
        elem = new ElemClone<Prism6_FVM>(p);
        break;
      }

      case PYRAMID5:
      {
        elem = new ElemClone<Pyramid5>(p);
        break;
      }
      case PYRAMID5_FVM:
      {
        elem = new ElemClone<Pyramid5_FVM>(p);
        break;
      }

#else
      case TET4:
      case TET4_FVM:
      case HEX8:
      case HEX8_FVM:
      case PRISM6:
      case PRISM6_FVM:
      case PYRAMID5:
      case PYRAMID5_FVM:
      {
        std::cerr << "ERROR: 3D elements are not supported by Open Source Version." << std::endl;
        genius_error();
      }
#endif

      default:
      {
        std::cerr << "ERROR: Unsupported element type!." << std::endl;
        genius_error();
      }
  }

  AutoPtr<Elem> ap(elem);

  if (p) ap->subdomain_id() = p->subdomain_id();

  return ap;
}



unsigned int Elem::key() const
{
  switch (this->type())
  {
      case EDGE2     :
      case EDGE2_FVM :
      return this->compute_key(this->node(0), this->node(1));
      case TRI3      :
      case TRI3_FVM  :
      return this->compute_key(this->node(0), this->node(1), this->node(2));
      case QUAD4     :
      case QUAD4_FVM :
      return this->compute_key(this->node(0), this->node(1), this->node(2), this->node(3));
      case TET4         :
      case TET4_FVM     :
      case PYRAMID5     :
      case PYRAMID5_FVM :
      case PRISM6       :
      case PRISM6_FVM   :
      case HEX8         :
      case HEX8_FVM     :
      {
        unsigned int _key=0;
        for(unsigned int s=0; s<this->n_sides(); ++s)
          _key += this->key(s);
        return _key;
      }
      default: genius_error();
  }

  genius_error();

  return 0;
}



Point Elem::centroid() const
{
  Point cp;

  for (unsigned int n=0; n<this->n_vertices(); n++)
    cp.add (this->point(n));

  return (cp /= static_cast<Real>(this->n_vertices()));
}



Real Elem::hmin() const
{
  Real h_min=1.e30;

  for (unsigned int n_outer=0; n_outer<this->n_vertices(); n_outer++)
    for (unsigned int n_inner=n_outer+1; n_inner<this->n_vertices(); n_inner++)
    {
      const Point diff = (this->point(n_outer) - this->point(n_inner));

      h_min = std::min(h_min,diff.size());
    }

  return h_min;
}



Real Elem::hmax() const
{
  Real h_max=0;

  for (unsigned int n_outer=0; n_outer<this->n_vertices(); n_outer++)
    for (unsigned int n_inner=n_outer+1; n_inner<this->n_vertices(); n_inner++)
    {
      const Point diff = (this->point(n_outer) - this->point(n_inner));

      h_max = std::max(h_max,diff.size());
    }

  return h_max;
}



Real Elem::length(const unsigned int n1,
                  const unsigned int n2) const
{
  assert ( n1 < this->n_vertices() );
  assert ( n2 < this->n_vertices() );

  return (this->point(n1) - this->point(n2)).size();
}



bool Elem::operator == (const DofObject& rhs) const
{

  // Cast rhs to an Elem*
  const Elem* rhs_elem = dynamic_cast<const Elem*>(&rhs);

  // If we cannot cast to an Elem*, rhs must be a Node
  if(rhs_elem == static_cast<const Elem*>(NULL))
    return false;

  //   assert (n_nodes());
  //   assert (rhs.n_nodes());

  //   // Elements can only be equal if they
  //   // contain the same number of nodes.
  //   if (this->n_nodes() == rhs.n_nodes())
  //     {
  //       // Create a set that contains our global
  //       // node numbers and those of our neighbor.
  //       // If the set is the same size as the number
  //       // of nodes in both elements then they must
  //       // be connected to the same nodes.
  //       std::set<unsigned int> nodes_set;

  //       for (unsigned int n=0; n<this->n_nodes(); n++)
  //         {
  //           nodes_set.insert(this->node(n));
  //           nodes_set.insert(rhs.node(n));
  //         }

  //       // If this passes the elements are connected
  //       // to the same global nodes
  //       if (nodes_set.size() == this->n_nodes())
  //         return true;
  //     }

  //   // If we get here it is because the elements either
  //   // do not have the same number of nodes or they are
  //   // connected to different nodes.  Either way they
  //   // are not the same element
  //   return false;

  // Useful typedefs
  typedef std::vector<unsigned int>::iterator iterator;


  // Elements can only be equal if they
  // contain the same number of nodes.
  // However, we will only test the vertices,
  // which is sufficient & cheaper
  if (this->n_nodes() == rhs_elem->n_nodes())
  {
    // The number of nodes in the element
    const unsigned int nn = this->n_nodes();

    // Create a vector that contains our global
    // node numbers and those of our neighbor.
    // If the sorted, unique vector is the same size
    // as the number of nodes in both elements then
    // they must be connected to the same nodes.
    //
    // The vector will be no larger than 2*n_nodes(),
    // so we might as well reserve the space.
    std::vector<unsigned int> common_nodes;
    common_nodes.reserve (2*nn);

    // Add the global indices of the nodes
    for (unsigned int n=0; n<nn; n++)
    {
      common_nodes.push_back (this->node(n));
      common_nodes.push_back (rhs_elem->node(n));
    }

    // Sort the vector and find out how long
    // the sorted vector is.
    std::sort (common_nodes.begin(), common_nodes.end());

    iterator new_end = std::unique (common_nodes.begin(),
                                    common_nodes.end());

    const int new_size = std::distance (common_nodes.begin(),
                                        new_end);

    // If this passes the elements are connected
    // to the same global vertex nodes
    if (new_size == static_cast<int>(nn))
      return true;
  }

  // If we get here it is because the elements either
  // do not have the same number of nodes or they are
  // connected to different nodes.  Either way they
  // are not the same element
  return false;
}



bool Elem::contains_vertex_of(const Elem *e) const
{
  // Our vertices are the first numbered nodes
  for (unsigned int n = 0; n != e->n_vertices(); ++n)
    for (unsigned int i = 0; i != this->n_vertices(); ++i)
      if( e->node(n) == this->node(i) )
        return true;
  return false;
}



bool Elem::contains_all_vertex_of(const Elem *e) const
{
  genius_assert(e);

  unsigned int match_count=0;

  for (unsigned int n = 0; n < e->n_vertices(); ++n)
    for (unsigned int i = 0; i < this->n_vertices(); ++i)
      if( e->node(n) == this->node(i) )
        ++match_count;

  return (match_count == e->n_vertices() );
}



void Elem::find_point_neighbors(std::set<const Elem *> &neighbor_set) const
  {
    neighbor_set.clear();
    neighbor_set.insert(this);

    unsigned int old_size;
    do
    {
      old_size = neighbor_set.size();

      // Loop over all the elements in the patch
      std::set<const Elem*>::const_iterator       it  = neighbor_set.begin();
      const std::set<const Elem*>::const_iterator end = neighbor_set.end();

      for (; it != end; ++it)
      {
        const Elem* elem = *it;

        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor(s) != NULL)           // we have a neighbor on this side
          {
            const Elem* neighbor = elem->neighbor(s);

            if (neighbor->active())                // ... if it is active
            {
              if (this->contains_vertex_of(neighbor) // ... and touches us
                  || neighbor->contains_vertex_of(this))
                neighbor_set.insert (neighbor);  // ... then add it
            }
#ifdef ENABLE_AMR
            else                                 // ... the neighbor is *not* active,
            {                                  // ... so add *all* neighboring
              // active children
              std::vector<const Elem*> active_neighbor_children;

              neighbor->active_family_tree_by_neighbor
              (active_neighbor_children, elem);

              std::vector<const Elem*>::const_iterator
              child_it = active_neighbor_children.begin();
              const std::vector<const Elem*>::const_iterator
              child_end = active_neighbor_children.end();
              for (; child_it != child_end; ++child_it)
                if (this->contains_vertex_of(*child_it) ||
                    (*child_it)->contains_vertex_of(this))
                  neighbor_set.insert (*child_it);
            }
#endif // #ifdef ENABLE_AMR
          }
      }
    }
    while (old_size != neighbor_set.size());
  }



void Elem::write_connectivity (std::ostream& out,
                               const IOPackage iop) const
{
  assert (out.good());
  assert (_nodes != NULL);
  assert (iop != INVALID_IO_PACKAGE);

  switch (iop)
  {
      case TECPLOT:
      {
        // This connectivity vector will be used repeatedly instead
        // of being reconstructed inside the loop.
        std::vector<unsigned int> conn;
        for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
        {
          this->connectivity(sc, TECPLOT, conn);

          std::copy(conn.begin(),
                    conn.end(),
                    std::ostream_iterator<unsigned int>(out, " "));

          out << '\n';
        }
        return;
      }

      case UCD:
      {
        for (unsigned int i=0; i<this->n_nodes(); i++)
          out << this->node(i)+1 << "\t";

        out << '\n';
        return;
      }

      default:
      genius_error();
  }

  genius_error();
}


// void Elem::write_tecplot_connectivity(std::ostream& out) const
// {
//   assert (!out.bad());
//   assert (_nodes != NULL);

//   // This connectivity vector will be used repeatedly instead
//   // of being reconstructed inside the loop.
//   std::vector<unsigned int> conn;
//   for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
//     {
//       this->connectivity(sc, TECPLOT, conn);

//       std::copy(conn.begin(),
//              conn.end(),
//              std::ostream_iterator<unsigned int>(out, " "));

//       out << std::endl;
//     }
// }



// void Elem::write_ucd_connectivity(std::ostream &out) const
// {
//   assert (out);
//   assert (_nodes != NULL);

//   for (unsigned int i=0; i<this->n_nodes(); i++)
//     out << this->node(i)+1 << "\t";

//   out << std::endl;
// }



Real Elem::quality (const ElemQuality q) const
{
  switch (q)
  {
      /**
       * I don't know what to do for this metric.
       */
      default:
      {
        genius_here();

        std::cerr << "ERROR:  unknown quality metric: "
        << q
        << std::endl
        << "Cowardly returning 1."
        << std::endl;

        return 1.;
      }
  }


  // Will never get here...
  genius_error();
  return 0.;
}



/**
 * return the interpolated value at given point
 */
PetscScalar Elem::interpolation( const std::vector<PetscScalar> & value, const Point &p) const
{
  genius_assert(value.size() == this->n_nodes());

  PetscScalar v = 0.0;
  std::vector<Real> w(this->n_nodes());
  for(unsigned int n=0; n<this->n_nodes(); ++n)
  {
    w[n] = 1.0/((p - this->point(n)).size_sq()+1e-6);
    v += w[n]*value[n];
  }

  return v/std::accumulate(w.begin(), w.end(), 0.0);
}


#ifdef ENABLE_AMR



bool Elem::ancestor() const
{
#ifdef ENABLE_AMR

  if (this->active())
    return false;

  if (!this->has_children())
    return false;
  if (this->child(0)->active())
    return true;

  return this->child(0)->ancestor();
#else
  return false;
#endif
}



void Elem::add_child (Elem* elem)
{
  if(_children == NULL)
  {
    _children = new Elem*[this->n_children()];

    for (unsigned int c=0; c<this->n_children(); c++)
      _children[c] = NULL;
  }

  for (unsigned int c=0; c<this->n_children(); c++)
  {
    if(_children[c] == NULL)
    {
      _children[c] = elem;
      return;
    }
  }

  std::cerr << "Error: Tried to add a child to an element with full children array"
  << std::endl;
  genius_error();
}



void Elem::add_child (Elem* elem, unsigned int pos)
{
  if(_children == NULL)
  {
    _children = new Elem*[this->n_children()];

    for (unsigned int c=0; c<this->n_children(); c++)
      _children[c] = NULL;
  }

  if(_children[pos] == NULL)
  {
    _children[pos] = elem;
    return;
  }

  std::cerr << "Error: Tried to add a child to an element with full children array"
  << std::endl;
  genius_error();
}


void Elem::delete_child (Elem* elem)
{
  for (unsigned int c=0; c<this->n_children(); c++)
  {
    if(_children[c] == elem)
    {
      _children[c] = NULL;
      return;
    }
  }

  std::cerr << "Error: Tried to delete a child to an element with NULL children array"
  << std::endl;
  genius_error();
}


bool Elem::is_child_on_edge(const unsigned int c,
                            const unsigned int e) const
{
  assert (c < this->n_children());
  assert (e < this->n_edges());

  AutoPtr<Elem> my_edge = this->build_edge(e);
  AutoPtr<Elem> child_edge = this->build_edge(e);

  // We're assuming that an overlapping child edge has the same
  // number and orientation as its parent
  return (child_edge->node(0) == my_edge->node(0) ||
          child_edge->node(1) == my_edge->node(1));
}


bool Elem::is_child_on_edge(const Elem * child,
                            const unsigned int e) const
{
  assert (child->parent() == this);
  assert (e < this->n_edges());

  AutoPtr<Elem> my_edge = this->build_edge(e);
  AutoPtr<Elem> child_edge = child->build_edge(e);

  // We're assuming that an overlapping child edge has the same
  // number and orientation as its parent
  return (child_edge->node(0) == my_edge->node(0) ||
          child_edge->node(1) == my_edge->node(1));
}



bool Elem::is_child_on_side(const Elem * child,
                            const unsigned int s) const
{
  unsigned int c = this->which_child_am_i(child);

  return this->is_child_on_side(c, s);
}




void Elem::family_tree (std::vector<const Elem*>& family,
                        const bool reset) const
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      this->child(c)->family_tree (family, false);
}



void Elem::active_family_tree (std::vector<const Elem*>& active_family,
                               const bool reset) const
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    active_family.clear();

  // Add this element to the family tree if it is active
  if (this->active())
    active_family.push_back(this);

  // Otherwise recurse into the element's children.
  // Do not clear the vector any more.
  else
    for (unsigned int c=0; c<this->n_children(); c++)
      this->child(c)->active_family_tree (active_family, false);

}



void Elem::family_tree_by_neighbor (std::vector<const Elem*>& family,
                                    const Elem* neighbor,
                                    const bool reset) const
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
  assert (this->is_neighbor(neighbor));

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      if (this->child(c)->is_neighbor(neighbor))
        this->child(c)->family_tree_by_neighbor (family, neighbor, false);
}



void Elem::active_family_tree_by_neighbor (std::vector<const Elem*>& family,
    const Elem* neighbor,
    const bool reset) const
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
  if (this->level() >= neighbor->level())
    assert (this->is_neighbor(neighbor));

  // Add an active element to the family tree.
  if (this->active())
    family.push_back(this);

  // Or recurse into an ancestor element's children.
  // Do not clear the vector any more.
  else
    for (unsigned int c=0; c<this->n_children(); c++)
      if (this->child(c)->is_neighbor(neighbor))
        this->child(c)->active_family_tree_by_neighbor (family, neighbor, false);
}




void Elem::family_tree_by_side (std::vector<const Elem*>& family,
                                const unsigned int s,
                                const bool reset)  const
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  genius_assert( s<this->n_sides() );

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      if (this->child(c)->is_child_on_side(c, s))
        this->child(c)->family_tree_by_side (family, s, false);

}



void Elem::active_family_tree_by_side (std::vector<const Elem*>& family,
                                       const unsigned int s,
                                       const bool reset) const
{

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  genius_assert( s<this->n_sides() );

  // Add an active element to the family tree.
  if (this->active())
    family.push_back(this);

  // Or recurse into an ancestor element's children.
  // Do not clear the vector any more.
  else
    for (unsigned int c=0; c<this->n_children(); c++)
      if (this->child(c)->is_child_on_side(c, s))
        this->child(c)->active_family_tree_by_side (family, s, false);
}



unsigned int Elem::min_p_level_by_neighbor(const Elem* neighbor,
    unsigned int current_min) const
{
  genius_assert(!this->subactive());
  genius_assert(neighbor->active());

  // If we're an active element this is simple
  if (this->active())
    return std::min(current_min, this->p_level());

  genius_assert(is_neighbor(neighbor));

  // The p_level() of an ancestor element is already the minimum
  // p_level() of its children - so if that's high enough, we don't
  // need to examine any children.
  if (current_min <= this->p_level())
    return current_min;

  unsigned int min_p_level = current_min;

  for (unsigned int c=0; c<this->n_children(); c++)
  {
    const Elem* const child = this->child(c);
    if (child->is_neighbor(neighbor))
      min_p_level =
        child->min_p_level_by_neighbor(neighbor,
                                       min_p_level);
  }

  return min_p_level;
}


unsigned int Elem::min_new_p_level_by_neighbor(const Elem* neighbor,
    unsigned int current_min) const
{
  genius_assert(!this->subactive());
  genius_assert(neighbor->active());

  // If we're an active element this is simple
  if (this->active())
  {
    unsigned int new_p_level = this->p_level();
    if (this->p_refinement_flag() == Elem::REFINE)
      new_p_level += 1;
    if (this->p_refinement_flag() == Elem::COARSEN)
    {
      assert (new_p_level > 0);
      new_p_level -= 1;
    }
    return std::min(current_min, new_p_level);
  }

  genius_assert(is_neighbor(neighbor));

  unsigned int min_p_level = current_min;

  for (unsigned int c=0; c<this->n_children(); c++)
  {
    const Elem* const child = this->child(c);
    if (child->is_neighbor(neighbor))
      min_p_level =
        child->min_new_p_level_by_neighbor(neighbor,
                                           min_p_level);
  }

  return min_p_level;
}

#endif // #ifdef ENABLE_AMR



bool Elem::contains_point (const Point & p) const
{
  genius_assert(&p);
  // we should never reach here
  // derived class should override it.
  genius_error();
  return true;
}



void Elem::nullify_neighbors ()
{
  // Tell any of my neighbors about my death...
  // Looks strange, huh?
  for (unsigned int n=0; n<this->n_neighbors(); n++)
    if (this->neighbor(n) != NULL)
    {
      Elem* neighbor = this->neighbor(n);

      // Note:  it is possible that I see the neighbor
      // (which is coarser than me)
      // but they don't see me, so avoid that case.
      if (neighbor->level() == this->level())
      {
        const unsigned int w_n_a_i = neighbor->which_neighbor_am_i(this);
        neighbor->set_neighbor(w_n_a_i, NULL);
        this->set_neighbor(n, NULL);
      }
    }
}



unsigned int Elem::n_second_order_adjacent_vertices (const unsigned int) const
{
  // for linear elements, always return 0
  return 0;
}



unsigned short int Elem::second_order_adjacent_vertex (const unsigned int,
    const unsigned int) const
{
  // for linear elements, always return 0
  return 0;
}



std::pair<unsigned short int, unsigned short int>
Elem::second_order_child_vertex (const unsigned int) const
{
  // for linear elements, always return 0
  return std::pair<unsigned short int, unsigned short int>(0,0);
}



ElemType Elem::first_order_equivalent_type (const ElemType et)
{
  switch (et)
  {
      case EDGE2:
      case EDGE3:
      case EDGE4:
      return EDGE2;
      case EDGE2_FVM:
      return EDGE2_FVM;
      case TRI3:
      case TRI6:
      return TRI3;
      case TRI3_FVM:
      return TRI3_FVM;
      case QUAD4:
      case QUAD8:
      case QUAD9:
      return QUAD4;
      case QUAD4_FVM:
      return QUAD4_FVM;
      case TET4:
      case TET10:
      return TET4;
      case TET4_FVM:
      return TET4_FVM;
      case HEX8:
      case HEX27:
      case HEX20:
      return HEX8;
      case HEX8_FVM:
      return HEX8_FVM;
      case PRISM6:
      case PRISM15:
      case PRISM18:
      return PRISM6;
      case PRISM6_FVM:
      return PRISM6_FVM;
      case PYRAMID5:
      return PYRAMID5;
      case PYRAMID5_FVM:
      return PYRAMID5_FVM;

      default:
      // unknown element
      return INVALID_ELEM;
  }
}



ElemType Elem::second_order_equivalent_type (const ElemType et,
    const bool full_ordered)
{
  /* for second-order elements, always return \p INVALID_ELEM
   * since second-order elements should not be converted
   * into something else.  Only linear elements should
   * return something sensible here
   */
  switch (et)
  {
      case EDGE2:
      case EDGE2_FVM:
      {
        // full_ordered not relevant
        return EDGE3;
      }

      case TRI3:
      case TRI3_FVM:
      {
        // full_ordered not relevant
        return TRI6;
      }

      case QUAD4:
      case QUAD4_FVM:
      {
        if (full_ordered)
          return QUAD9;
        else
          return QUAD8;
      }

      case TET4:
      case TET4_FVM:
      {
        // full_ordered not relevant
        return TET10;
      }

      case HEX8:
      case HEX8_FVM:
      {
        // see below how this correlates with INFHEX8
        if (full_ordered)
          return HEX27;
        else
          return HEX20;
      }

      case PRISM6:
      case PRISM6_FVM:
      {
        if (full_ordered)
          return PRISM18;
        else
          return PRISM15;
      }

      case PYRAMID5:
      case PYRAMID5_FVM:
      {
        genius_error();
        return INVALID_ELEM;
      }

      default:
      {
        // second-order element
        return INVALID_ELEM;
      }
  }
}


ElemType Elem::fvm_compatible_type (const ElemType et)
{
  switch (et)
  {
      case EDGE2:
      case EDGE2_FVM:
      case EDGE3:
      case EDGE4:
      return EDGE2_FVM;
      case TRI3:
      case TRI3_FVM:
      case TRI6:
      return TRI3_FVM;
      case QUAD4:
      case QUAD4_FVM:
      case QUAD8:
      case QUAD9:
      return QUAD4_FVM;
      case TET4:
      case TET4_FVM:
      case TET10:
      return TET4_FVM;
      case HEX8:
      case HEX8_FVM:
      case HEX27:
      case HEX20:
      return HEX8_FVM;
      case PRISM6:
      case PRISM6_FVM:
      case PRISM15:
      case PRISM18:
      return PRISM6_FVM;
      case PYRAMID5:
      case PYRAMID5_FVM:
      return PYRAMID5_FVM;

      default:
      // unknown element
      return INVALID_ELEM;
  }
}


unsigned int Elem::dim (const ElemType et)
{
  switch (et)
  {
      case EDGE2:
      case EDGE2_FVM:
      case EDGE3:
      case EDGE4:
      return 1;
      case TRI3:
      case TRI3_FVM:
      case TRI6:
      return 2;
      case QUAD4:
      case QUAD4_FVM:
      case QUAD8:
      case QUAD9:
      return 2;
      case TET4:
      case TET4_FVM:
      case TET10:
      return 3;
      case HEX8:
      case HEX8_FVM:
      case HEX27:
      case HEX20:
      return 3;
      case PRISM6:
      case PRISM6_FVM:
      case PRISM15:
      case PRISM18:
      return 3;
      case PYRAMID5:
      case PYRAMID5_FVM:
      return 3;

      default:
      // unknown element
      return 0;
  }
}

Elem::side_iterator Elem::boundary_sides_begin()
{
  Predicates::BoundarySide<SideIter> bsp;
  return side_iterator(this->_first_side(), this->_last_side(), bsp);
}




Elem::side_iterator Elem::boundary_sides_end()
{
  Predicates::BoundarySide<SideIter> bsp;
  return side_iterator(this->_last_side(), this->_last_side(), bsp);
}


Point Elem::outside_unit_normal(unsigned short int side_id) const
{
  genius_assert(side_id < this->n_sides());

  genius_error();

  return Point(0.0, 0.0, 0.0);
}



std::vector<unsigned short int> Elem::get_terminate_side(const Point & terminate_point) const
{
  std::vector<unsigned short int> side_contain_point;

  for (unsigned short int side=0; side<this->n_sides(); ++side)
  {
    AutoPtr<Elem> side_elem = this->build_side(side);
    if (side_elem->contains_point(terminate_point))
      side_contain_point.push_back(side);
  }
  return side_contain_point;
}


// Pack all this information into one communication to avoid two latency hits
// For each element it is of the form
// [ level p_level r_flag p_flag etype subdomain_id
//   self_ID parent_ID which_child node_0 node_1 ... node_n]
// We cannot use unsigned int because parent_ID can be negative
void Elem::pack_element (std::vector<int> &conn) const
{
  conn.push_back (static_cast<int>(this->type()));
  conn.push_back (static_cast<int>(this->processor_id()));
  conn.push_back (static_cast<int>(this->subdomain_id()));
  conn.push_back (this->id());

#ifdef ENABLE_AMR
  conn.push_back (static_cast<int>(this->level()));
  conn.push_back (static_cast<int>(this->p_level()));
  conn.push_back (static_cast<int>(this->refinement_flag()));
  conn.push_back (static_cast<int>(this->p_refinement_flag()));

  // use parent_ID of -1 to indicate a level 0 element
  if (this->level() == 0)
  {
    conn.push_back(-1);
    conn.push_back(-1);
  }
  else
  {
    conn.push_back(this->parent()->id());
    conn.push_back(this->parent()->which_child_am_i(this));
  }
#endif

  for (unsigned int n=0; n<this->n_nodes(); n++)
    conn.push_back (this->node(n));

  for (unsigned int n=0; n<this->n_sides(); n++)
  {
    const Elem * neighbor = this->neighbor(n);
    int neighbor_id = neighbor ? static_cast<int>(neighbor->id()) : -1;
    conn.push_back (neighbor_id);
  }

}


unsigned int Elem::pack_size( ElemType t )
{
#ifdef ENABLE_AMR
  unsigned int header = 10;
#else
  unsigned int header = 4;
#endif

  switch (t)
  {
      case EDGE2        :
      case EDGE2_FVM    :      return header + 2 + 2;
      case TRI3         :
      case TRI3_FVM     :      return header + 3 + 3;
      case QUAD4        :
      case QUAD4_FVM    :      return header + 4 + 4;
      case TET4         :
      case TET4_FVM     :      return header + 4 + 4;
      case PYRAMID5     :
      case PYRAMID5_FVM :      return header + 5 + 5;
      case PRISM6       :
      case PRISM6_FVM   :      return header + 6 + 5;
      case HEX8         :
      case HEX8_FVM     :      return header + 8 + 6;
      default: genius_error();
  }
  genius_error();
  return 0;
}

