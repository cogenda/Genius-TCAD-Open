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

//  $Id: adolc_init.cc,v 1.2 2008/07/09 05:58:16 gdiso Exp $

#include "adolc.h"

unsigned int adtl::AutoDScalar::numdir = 12;


extern "C"
{
  DLL_EXPORT_DECLARE  void  set_ad_number(const unsigned int p)
  {
    adtl::AutoDScalar::numdir = p;
  }
}
