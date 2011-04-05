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

// TIF file parser

#include <cstdio>
#include <cstring>

#include "tif_data.h"

// avoid isatty() problem of Bison 2.3
#include "config.h"
#define YY_NEVER_INTERACTIVE 1
#ifdef CYGWIN
  #define YY_NO_UNISTD_H 1
extern "C"
{
extern int isatty (int );
}
#endif

namespace TIF
{
  #include "tif_lex.yy.c"
  #include "tif_parser.tab.c"


  void clear()
  {
    node_array.clear();
    edge_array.clear();
    tri_array.clear();
    quad_array.clear();
    region_array.clear();
    interface_array.clear();
    component_array.clear();
    parameter_array.clear();
    sol_head.clear();
    sol_data.clear();
    tran_sol.clear();
  }



}

