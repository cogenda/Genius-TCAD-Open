#include <iostream>
#include <cmath>

#include "PMI_benchmark_avalanche.h"

#include <dlfcn.h>
#define LDFUN dlsym


PMI_Benchmark_ImpactIonization::PMI_Benchmark_ImpactIonization(const std::string &path, const std::string &material, const std::string &model)
    :dll_file(0)
{
  cm = 1e6;
  s  = 1e12;
  V  = 1.0;
  C  = 1.0/1.602176462e-19;
  K  = 1.0/300;
  eV = 1.0;
  m  = 1e2*cm;
  J  = C*V;
  kg = J/(m*m)*s*s;
  g  = 1e-3*kg;

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
  std::string model_fun_name = "PMIS_" + material + "_Avalanche_" + model;
  PMIS_Avalanche*     (*wgen)      (const PMI_Environment& env);
  wgen = (PMIS_Avalanche* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());

  PMI_Environment env(100*cm, s, V, C, K);
  agen = wgen(env);
}


PMI_Benchmark_ImpactIonization::~PMI_Benchmark_ImpactIonization()
{
  delete agen;
  dlclose( dll_file );
}


bool PMI_Benchmark_ImpactIonization::calibrate_real_parameter(const std::string & var_name, double var_value)
{
  if(agen->calibrate_real_parameter(var_name, var_value)) return true;
  return false;
}

bool PMI_Benchmark_ImpactIonization::calibrate_string_parameter(const std::string & var_name, const std::string &var_value)
{
  if(agen->calibrate_string_parameter(var_name, var_value)) return true;
  return false;
}

void PMI_Benchmark_ImpactIonization::set_doping(double Na, double Nd)
{
  agen->SetFakeDopingEnvironment(Na*pow(cm,-3), Nd*pow(cm,-3));
}


void PMI_Benchmark_ImpactIonization::set_mole(double mole_x, double mole_y)
{
  agen->SetFakeMoleEnvironment(mole_x, mole_y);
}


/**
 * @return the electron generation rate for DDM simulation
 */
double  PMI_Benchmark_ImpactIonization::ElecGenRate (const double Tl,const double Ep,const double Eg) const
{
  return agen->ElecGenRate(Tl*K, Ep*V/cm, Eg*eV)*cm;
}

/**
 * @return the hole generation rate for DDM simulation
 */
double PMI_Benchmark_ImpactIonization::HoleGenRate (const double Tl,const double Ep,const double Eg) const
{
  return agen->HoleGenRate(Tl*K, Ep*V/cm, Eg*eV)*cm;
}


