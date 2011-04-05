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

//  $Id: point_solver.h,v 1.3 2008/07/09 05:58:16 gdiso Exp $


#ifndef __point_solver_h__
#define __point_solver_h__

#include "solver_base.h"



/**
 * The point solver is a most simple solver which node solution only related with
 * the node itself.
 */
class PointSolver :  public  SolverBase
{

public:
  /**
   * constructor
   */
  PointSolver(SimulationSystem & system) : SolverBase(system) {}

  /**
   * destructor
   */
  ~PointSolver() {}

};

#endif //#ifndef __point_solver_h__

