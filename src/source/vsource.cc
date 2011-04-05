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

#include <cassert>

#include "vsource.h"

#ifdef CYGWIN
  #include <Windows.h>
  #undef max
  #undef min
#else
  #include <dlfcn.h>
#endif


#ifdef CYGWIN
VSHELL::VSHELL(const std::string & s, HINSTANCE dp, void * fp, double s_t,double s_V): VSource(s)
#else
VSHELL::VSHELL(const std::string & s, void * dp, void * fp, double s_t,double s_V): VSource(s)
#endif
{
  dll = dp;
  Vapp_Shell = (double (*)(double)) fp;
  assert(Vapp_Shell);
  scale_t = s_t;
  scale_V = s_V;
}


VSHELL::~VSHELL()
{
     //delete Vapp_Shell;
#ifdef CYGWIN
    FreeLibrary(dll);
#else
    if ( dll )
      dlclose( dll );
#endif
}