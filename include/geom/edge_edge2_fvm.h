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

//  $Id: edge_edge2_fvm.h,v 1.2 2008/07/10 07:16:23 gdiso Exp $



#ifndef __edge2_fvm_h__
#define __edge2_fvm_h__

// C++ includes


// Local includes
#include "genius_common.h"
#include "edge_edge2.h"


/**
 * The \p Edge2 is an element in 1D composed of 2 nodes. It is numbered
 * like this:
 *
   \verbatim
    EDGE2: o--------o
           0        1
   \endverbatim
 */

// ------------------------------------------------------------
// Edge class definition
class Edge2_FVM : public Edge2
{
public:

  /**
   * Constructor.  By default this element has no parent.
   */
  Edge2_FVM (Elem* p=NULL) :
  Edge2(p) {}

  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  Edge2_FVM (const unsigned int nn,
             const unsigned int ns,
             Elem* p) :
  Edge2(nn, ns, p) {}

  /**
   * @returns \p Edge2_FVM
   */
  virtual ElemType type () const { return EDGE2_FVM; }

  /**
   * @returns the pointer to local \p FVM_Node \p i.
   */
  virtual FVM_Node * get_fvm_node(const unsigned int i) const
    { return _fvm_node[i]; }

  /**
   * set the \p ith FVM_Node pointer.
   */
  virtual void hold_fvm_node(const unsigned int i, FVM_Node *pn)
  { _fvm_node[i] = pn; }


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
   * @return precomputed volume
   */
  virtual Real volume () const
    { return vol; }

  /**
   * @return the partial (length/area/volume) of the geometric element by local index.
   */
  virtual Real partial_volume (unsigned int) const
    { return 0.5*vol; }

  /**
   * @return the edge associated partial (length/area) of the geometric element with local edge index.
   * for an 1D element, return 0 here.
   */
  virtual Real partial_area_with_edge(unsigned int) const
    { return 0; }

  /**
   * @return the edge associated partial (area/volume) of the geometric element with local edge index.
   */
  virtual Real partial_volume_with_edge(unsigned int ) const
    { return 0;}

  // For FVM usage, we need more Geometry information of an Edge
private:

  /**
   * hold fvm_node pointer
   */
  FVM_Node * _fvm_node[2];

  /**
   * the volume of Edge2, store this vlaue for efficiency reason
   */
  Real vol;

  /**
   * calculate geom information for fvm usage
   */
  virtual void prepare_for_fvm();
};


#endif
