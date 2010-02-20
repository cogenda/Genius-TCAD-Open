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


#ifndef __partical_solver_h__
#define __partical_solver_h__

#include <vector>

#include "enum_petsc_type.h"
#include "solver_base.h"
//#include "petscis.h"
//#include "petscvec.h"


 
/**
 * The partical solver contex. This solver is used for partical simulation, i.e.
 * ray-tracing for optical problem, M.C. simulation and so on. 
 */
class ParticalSolver : public SolverBase
{
public:  
  
  /**
   * construct the partical solver contex:
   */
  ParticalSolver(SimulationSystem & system);
  
  /**
   * free all the contex
   */
  virtual ~ParticalSolver();


  // the main problem for parallel partical simulation is the exchange of particals in different processors
  //  
  
}; 


#endif // #define __partical_solver_h__
 
