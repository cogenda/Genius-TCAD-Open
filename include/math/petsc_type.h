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

//  $Id: petsc_type.h,v 1.3 2008/07/09 05:58:16 gdiso Exp $

#ifndef __petsc_type_h__
#define __petsc_type_h__

// convert string to enum

namespace SolverSpecify
{

  extern NonLinearSolverType nonlinear_solver_type(const std::string & ns);

  extern LinearSolverType linear_solver_type(const std::string & ls);
 
  extern PreconditionerType preconditioner_type(const std::string & pc);

}


#endif


