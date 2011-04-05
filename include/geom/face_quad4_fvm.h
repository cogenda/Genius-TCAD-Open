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

//  $Id: face_quad4_fvm.h,v 1.6 2008/07/10 09:39:38 gdiso Exp $



#ifndef __quad4_fvm_h__
#define __quad4_fvm_h__


// C++ includes


// Local includes
#include "genius_common.h"
#include "face_quad4.h"


// TNT matrix-vector library
#include <TNT/tnt.h>

// Forward declarations



/**
 * The \p QUAD4_FVM is QUAD4 elem with more Geom information for
 * FVM usage. It should satisfy in circle test.
 */

// ------------------------------------------------------------
// Quad4_FVM class definition
class Quad4_FVM : public Quad4
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Quad4_FVM (Elem* p=NULL) :
  Quad4(p) {}

  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  Quad4_FVM (const unsigned int nn,
             const unsigned int ns,
             Elem* p) :
  Quad4(nn, ns, p) {}

  /**
   * @returns \p QUAD4_FVM
   */
  virtual ElemType type () const { return QUAD4_FVM; }

  /**
   * @returns the pointer to local \p FVM_Node \p i.
   */
  virtual FVM_Node * get_fvm_node(const unsigned int i) const
    { return _fvm_node[i]; }

  /**
   * @returns the pointer to local \p FVM_Node \p i on side \p s.
   */
  virtual FVM_Node * get_side_fvm_node(const unsigned int s, const unsigned int i) const
  { return _fvm_node[side_nodes_map[s][i]]; }

  /**
   * set the \p ith FVM_Node pointer.
   */
  virtual void hold_fvm_node(const unsigned int i, FVM_Node *pn)
  { _fvm_node[i] = pn; }

  /**
   * @returns a proxy element coincident with side \p i.
   */
  virtual AutoPtr<Elem> build_fvm_side (const unsigned int i, bool proxy=true) const;

  /**
   * @return the gradient of input variable in the cell
   */
  virtual VectorValue<PetscScalar> gradient( const std::vector<PetscScalar> & var) const;

  /**
   * @return the gradient of input \p complex variable in the cell,
   */
  virtual VectorValue<Complex> gradient( const std::vector<Complex> & var) const;

  /**
   * @return the gradient of input \p AD variable in the cell
   */
  virtual VectorValue<AutoDScalar> gradient( const std::vector<AutoDScalar> & var) const;

  /**
   * when we know the projection of vector V to each edge of the cell, use least-squares method to
   * reconstruct vector V
   */
  virtual VectorValue<PetscScalar> reconstruct_vector( const std::vector<PetscScalar> & ) const;

  /**
   * when we know the projection of vector V to each edge of the cell, use least-squares method to
   * reconstruct vector V
   */
  virtual VectorValue<AutoDScalar> reconstruct_vector( const std::vector<AutoDScalar> & ) const;

  /**
   * @return precomputed volume
   */
  virtual Real volume () const
    { return vol; }

  /**
   * @return the partial (length/area/volume) of the geometric element by local index.
   */
  virtual Real partial_volume (unsigned int i) const
  { genius_assert( i<n_nodes() ); return v[i]; }

  /**
   * @return the edge associated partial (length/area) of the geometric element with local edge index.
   */
  virtual Real partial_area_with_edge(unsigned int e) const
    { genius_assert( e<n_edges() ); return d[e]; }

  /**
   * @return the edge associated partial (area/volume) of the geometric element with local edge index.
   */
  virtual Real partial_volume_with_edge(unsigned int e) const
  { return 0.5 * partial_area_with_edge(e) * l[e]; }

  /**
   * calculate geom information for fvm usage
   */
  virtual void prepare_for_fvm();

  // For FVM usage, we need more Geom information of a QUAD4
private:

  /**
   * hold fvm_node pointer
   */
  FVM_Node * _fvm_node[4];

  /**
   * the circumcircle center of QUAD4_FVM
   */
  //Point circumcircle_center;

  /**
   *  the projection point of circumcircle center to 4 edges
   */
  //Point side_centers[4];

  /**
   * the distance of circumcircle center to each side (between center and side center)
   */
  Real d[4]; // 4- side index

  /**
   * the length of the edge
   */
  Real l[4]; // 4- side index

  /**
   * partial volumn of each region seperated by segment made up of
   * circumcircle center and edge center
   */
  Real v[4]; // 4- node index

  /**
   * the volume of QUAD4, store this vlaue for efficiency reason
   */
  Real vol;

  /**
   * precomputed matrix inv[A^T.A].A^T
   */
  TNT::Array2D<Real> least_squares_gradient_matrix;

  /**
   * precomput matrix inv[A^T.A].A^T
   */
  void prepare_for_least_squares();


  /**
   * precomputed matrix inv[A^T.A].A^T for fast vector reconstruct computation
   */
  TNT::Array2D<Real> least_squares_vector_reconstruct_matrix;


  void prepare_for_vector_reconstruct();

};


#endif
