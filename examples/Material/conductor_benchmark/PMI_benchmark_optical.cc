#include <iostream>
#include <cmath>

#include "PMI_benchmark_optical.h"

#include <dlfcn.h>
#define LDFUN dlsym


PMI_Benchmark_Optical::PMI_Benchmark_Optical(const std::string &path, const std::string &material, const std::string &model)
  :dll_file(0)
{
  cm = 1e6;
  s  = 1e12;
  V  = 1.0;
  C  = 1.0/1.602176462e-19;
  K  = 1.0/300;
  um = 1e-4*cm;

  // load material lib
  std::string filename =  path + "/lib" + material + ".so";
  dll_file = dlopen(filename.c_str(), RTLD_LAZY|RTLD_DEEPBIND);
  if(!dll_file)
  {
    std::cerr<<"Open material file "<< filename <<" error." << '\n';
    std::cerr<<"Error code: " << dlerror() << '\n';
    throw(1);
  }

  // open mobility model
  std::string model_fun_name = "PMIC_" + material + "_Optical_" + model;
  PMIC_Optical*      (*woptical)      (const PMI_Environment& env);
  woptical  =  (PMIC_Optical* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());

  PMI_Environment env(100*cm, s, V, C, K);
  optical = woptical(env);
}


PMI_Benchmark_Optical::~PMI_Benchmark_Optical()
{
  delete optical;
  dlclose( dll_file );
}



double PMI_Benchmark_Optical::n(double lambda, double T)
{
  return optical->RefractionIndex(lambda*um, T*K).real();
}


double PMI_Benchmark_Optical::k(double lambda, double T)
{
  return optical->RefractionIndex(lambda*um, T*K).imag();
}


double PMI_Benchmark_Optical::alpha(double lambda, double T)
{
  double k= optical->RefractionIndex(lambda*um, T*K).imag();
  double alpha = 12.5663706144*k/(lambda*um)*cm;
  return alpha;
}
