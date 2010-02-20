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

//  $Id: edge_edge2_fvm.cc,v 1.1 2008/07/10 04:03:45 gdiso Exp $


// Local includes
#include "edge_edge2_fvm.h"

// build the Geom information here
void Edge2_FVM::prepare_for_fvm()
{
  //store this vlaue for efficiency reason
  vol = Edge2::volume();
}

/**
 * return the gradient of input variable in the cell
 */
VectorValue<PetscScalar> Edge2_FVM::gradient( const std::vector<PetscScalar> & var) const
{
  VectorValue<Real> p = (point(1)-point(0));
  PetscScalar dv = var[1]-var[0];
  return VectorValue<PetscScalar>(dv*p(0), dv*p(1), dv*p(2));
}

/**
 * return the gradient of input variable in the cell
 */
VectorValue<Complex> Edge2_FVM::gradient( const std::vector<Complex> & var) const
{
  VectorValue<Real> p = (point(1)-point(0));
  Complex dv = var[1]-var[0];
  return VectorValue<Complex>(dv*p(0), dv*p(1), dv*p(2));
}


/**
 * return the gradient of input \p AD variable in the cell
 */
VectorValue<AutoDScalar> Edge2_FVM::gradient( const std::vector<AutoDScalar> & var) const
{
  VectorValue<Real> p = (point(1)-point(0));
  AutoDScalar dv = var[1]-var[0];
  return VectorValue<AutoDScalar>(dv*p(0), dv*p(1), dv*p(2));
}
