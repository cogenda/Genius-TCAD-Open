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

//  $Id: enum_region.h,v 1.2 2008/07/09 05:58:15 gdiso Exp $



#ifndef __enum_region_h__
#define __enum_region_h__


/**
 * the region type genius support
 */
enum  SimulationRegionType 
{ 
  SemiconductorRegion=0, 
  InsulatorRegion, 
  ConductorRegion, 
  VacuumRegion, 
  PMLRegion, 
  InvalidRegion
};

#endif // #define __enum_region_h__




