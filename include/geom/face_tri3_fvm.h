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

//  $Id: face_tri3_fvm.h,v 1.6 2008/07/10 09:39:38 gdiso Exp $



#ifndef __tri3_fvm_h__
#define __tri3_fvm_h__


// C++ includes


// Local includes
#include "genius_common.h"
#include "face_tri3.h"



// Forward declarations

/**
 * The \p TRI3_FVM is TRI3 elem with more Geom information for
 * FVM usage.
 */

// ------------------------------------------------------------
// Tri3 class definition
class Tri3_FVM : public Tri3
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Tri3_FVM (Elem* p=NULL) :
  Tri3(p) {}

  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  Tri3_FVM (const unsigned int nn,
            const unsigned int ns,
            Elem* p) :
  Tri3(nn, ns, p) {}

  /**
   * @returns \p TRI3_FVM
   */
  virtual ElemType type () const { return TRI3_FVM; }

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
    { return v[i]; }

  /**
   * @return the truncated partial (length/area/volume) of the geometric element by local index.
   */
  virtual Real partial_volume_truncated (unsigned int i) const
    { return vt[i]; }


  /**
   * @return the edge associated partial (length/area) of the geometric element with local edge index.
   */
  virtual Real partial_area_with_edge(unsigned int e) const
    { return d[e]; }

  /**
   * @return the edge associated volume (area/volume) of the geometric element with local edge index.
   */
  virtual Real partial_volume_with_edge(unsigned int e) const
    { return 0.5*d[e]*l[e]; }


  /**
   * @return the edge associated truncated partial (length/area) of the geometric element with local edge index.
   */
  virtual Real partial_area_with_edge_truncated(unsigned int e) const
  { return dt[e]; }

  /**
   * @return the edge associated truncated partial (area/volume) of the geometric element with local edge index.
   */
  virtual Real partial_volume_with_edge_truncated(unsigned int e) const
  { return 0.5*dt[e]*l[e]; }


  /**
   * calculate geom information for fvm usage
   */
  virtual void prepare_for_fvm();

  // For FVM usage, we need more Geom information of a TRI3
  // memory storage is not a problem for 2D elem, here we buffer as many data as possible
private:

  /**
   * hold fvm_node pointer
   */
  FVM_Node * _fvm_node[3];

  /**
   * the circumcircle center of TRI3_FVM
   */
  //Point circumcircle_center;

  /**
   *  the projection point of circumcircle center to 3 edges
   */
  //Point side_centers[3];

  /**
   * the distance of circumcircle center to each side (between center and side center)
   */
  Real d[3]; // 3- side index

  /**
   * the truncated distance of circumcircle center to each side (between center and side center)
   */
  Real dt[3]; // 3- side index

  /**
   * the edge length
   */
  Real l[3]; // 3- side index

  /**
   * partial volumn of each region seperated by segment made up of
   * circumcircle center and edge center
   */
  Real v[3]; // 3- node index

  /**
   * truncated partial volumn of each region seperated by segment made up of
   * circumcircle center and edge center
   */
  Real vt[3]; // 3- node index

  /**
   * the volume of TRI3, store this vlaue for efficiency reason
   */
  Real vol;

  /**
   * precomputed matrix inv[A^T.A].A^T for fast vector reconstruct computation
   */
  Real least_squares_vector_reconstruct_matrix[2][3];

  /**
   * precomput matrix inv[A^T.A].A^T for fast vector reconstruct computation
   */
  void prepare_for_vector_reconstruct();


};


#endif
