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

//  $Id: dlhook.cc,v 1.5 2008/07/09 05:58:16 gdiso Exp $

#include "dlhook.h"
#include "solver_base.h"

#ifndef CYGWIN

#include <dlfcn.h>

DllHook::DllHook(SolverBase & solver, const std::string & name, void * fun_data)
:Hook(solver, name), dll_handle(0), hook(0)
{
  char * error;

  std::string genius_dir(Genius::genius_dir());
  std::string filename =  genius_dir + "/lib/" + name + ".so";

  // Get 'dlopen()' handle to the extension function
  dll_handle = dlopen(filename.c_str(),RTLD_LAZY);
  error = dlerror();
  if(error)
  {
    std::cout<< "Load hook faild: "<< error << std::endl;
    dll_handle = NULL;
    return;
  }

  // get the address of function get_hook in the dll file
  get_hook = (GET_HOOK *) dlsym(dll_handle, "get_hook");
  genius_assert( !dlerror() );

  // call it to get the real hook pointer
  hook = (*get_hook)(_solver, _name, fun_data);
  genius_assert(hook);

}



DllHook::~DllHook()
{
  if(hook) delete hook;
  if ( dll_handle ) dlclose( dll_handle );
}

#endif
