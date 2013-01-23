#include <iostream>
#include <cmath>

#include "PMI_benchmark_optical.h"

#include <dlfcn.h>
#define LDFUN dlsym


PMI_Benckmark_Optical::PMI_Benckmark_Optical(const std::string &path, const std::string &material, const std::string &model)
  :dll_file(0)
{
  cm = 1e6;
  s  = 1e12;
  V  = 1.0;
  C  = 1.0/1.602176462e-19;
  K  = 1.0/300;

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
  std::string model_fun_name = "PMIS_" + material + "_Optical_" + model;
  PMIS_Optical*      (*woptical)      (const PMI_Environment& env);
  woptical  =  (PMIS_Optical* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());

  PMI_Environment env(100*cm, s, V, C, K);
  optical = woptical(env);
}


PMI_Benckmark_Optical::~PMI_Benckmark_Optical()
{
  delete optical;
  dlclose( dll_file );
}



void PMI_Benckmark_Optical::set_doping(double Na, double Nd)
{
  optical->SetFakeDopingEnvironment(Na*pow(cm,-3), Nd*pow(cm,-3));
}


void PMI_Benckmark_Optical::set_mole(double mole_x, double mole_y)
{
  optical->SetFakeMoleEnvironment(mole_x, mole_y);
}


double PMI_Benckmark_Optical::n(double lambda)
{
  return optical->RefractionIndex(lambda*1e-4*cm, 300*K).real();
}


double PMI_Benckmark_Optical::k(double lambda)
{
  return optical->RefractionIndex(lambda*1e-4*cm, 300*K).imag();
}



