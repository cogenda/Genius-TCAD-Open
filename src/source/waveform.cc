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

#include "waveform.h"

#ifdef WINDOWS
  #include <Windows.h>
  #undef max
  #undef min
#else
  #include <dlfcn.h>
#endif


WaveformShell::WaveformShell(const std::string & s, const std::string & dll_file, const std::string & function, double s_t):Waveform(s)
{

#ifdef WINDOWS

  dll = LoadLibrary(dll_file.c_str());
  void *fp = GetProcAddress(dll, function.c_str());

#else

#ifdef RTLD_DEEPBIND
  dll = dlopen(dll_file.c_str(), RTLD_LAZY|RTLD_DEEPBIND);  assert(dll);
#else
  dll = dlopen(dll_file.c_str(), RTLD_LAZY);  assert(dll);
#endif

  void *fp = dlsym(dll, function.c_str());  assert(fp);

#endif

  assert(fp);
  Waveform_Shell = (double (*)(double)) fp;

  scale_t = s_t;
}

  /**
 * destructor, free the pointer
   */
WaveformShell::~WaveformShell()
{
  //delete Vapp_Shell;
#ifdef WINDOWS
  FreeLibrary(dll);
#else
  dlclose( dll );
#endif
}
