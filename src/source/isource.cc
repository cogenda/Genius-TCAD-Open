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

#include "isource.h"

#ifdef WINDOWS
  #include <Windows.h>
  #undef max
  #undef min
#else
  #include <dlfcn.h>
#endif


#ifdef WINDOWS
ISHELL::ISHELL(const std::string & s, HINSTANCE dp, void * fp, double s_t,double s_A): ISource(s)
#else
ISHELL::ISHELL(const std::string & s, void * dp, void * fp, double s_t,double s_A): ISource(s)
#endif
{
  dll = dp;
  Iapp_Shell = (double (*)(double)) fp;
  scale_t = s_t;
  scale_A = s_A;
}


ISHELL::~ISHELL()
{
  //delete Iapp_Shell;
#ifdef WINDOWS
  FreeLibrary(dll);
#else
  if ( dll )
    dlclose( dll );
#endif
}
