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

#include <fstream>
#include <cassert>

#include "waveform.h"

#ifdef WINDOWS
  #include <Windows.h>
  #undef max
  #undef min
#else
  #include <dlfcn.h>
#endif

#include "parallel.h"
#include "physical_unit.h"


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



//-------------------------------------------------------------------------------
#include "monot_cubic_interpolator.h"

WaveformFile::WaveformFile(const std::string & s, const std::string & filename):Waveform(s)
{
  if(Genius::processor_id()==0)
  {
    std::ifstream in;
    in.open(filename.c_str());
    genius_assert(in.good());

    while(!in.eof())
    {
      double t, v;
      in >> t >> v;
      _time.push_back(t*PhysicalUnit::s);
      _wave.push_back(v);
    }

    in.close();
  }

  Parallel::broadcast(_time);
  Parallel::broadcast(_wave);

  _int = new MonotCubicInterpolator(_time, _wave);
}


WaveformFile::~WaveformFile()
{
  delete _int;
}



double WaveformFile::waveform(double t)
{
  if(_time.empty()) return 0.0;
  if( t <_int->getMinimumX().first ) return 0.0;
  if( t >_int->getMaximumX().first ) return 0.0;
  return _int->evaluate(t);
}




