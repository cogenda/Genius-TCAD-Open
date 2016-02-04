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

#include <numeric>

#include "genius_common.h"
#include "fvm_flex_pde_solver.h"

#ifdef COGENDA_COMMERCIAL_PRODUCT
  #include "fvm_flex_parallel_dof_map.h"
#else
  #include "fvm_flex_serial_dof_map.h"
#endif



void FVM_FlexPDESolver::build_dof_map()
{
#ifdef COGENDA_COMMERCIAL_PRODUCT
    // for commercial version
    set_parallel_dof_map();
#else
    // for open source version
    set_serial_dof_map();
#endif
}









