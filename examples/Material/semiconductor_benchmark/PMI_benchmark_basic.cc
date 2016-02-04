#include <iostream>
#include <cmath>

#include "PMI_benchmark_basic.h"

#include <dlfcn.h>
#define LDFUN dlsym


PMI_Benchmark_Basic::PMI_Benchmark_Basic(const std::string &path, const std::string &material, const std::string &model)
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
  std::string model_fun_name = "PMIS_" + material + "_BasicParameter_" + model;
  PMIS_BasicParameter*(*wbasic)    (const PMI_Environment& env);
  wbasic = (PMIS_BasicParameter* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());

  PMI_Environment env(100*cm, s, V, C, K);
  basic = wbasic(env);
}


PMI_Benchmark_Basic::~PMI_Benchmark_Basic()
{
  delete basic;
  dlclose( dll_file );
}


bool PMI_Benchmark_Basic::calibrate_real_parameter(const std::string & var_name, double var_value)
{
  if(basic->calibrate_real_parameter(var_name, var_value)) return true;
  return false;
}

bool PMI_Benchmark_Basic::calibrate_string_parameter(const std::string & var_name, const std::string &var_value)
{
  if(basic->calibrate_string_parameter(var_name, var_value)) return true;
  return false;
}

void PMI_Benchmark_Basic::set_doping(double Na, double Nd)
{
  basic->SetFakeDopingEnvironment(Na*pow(cm,-3), Nd*pow(cm,-3));
}


void PMI_Benchmark_Basic::set_mole(double mole_x, double mole_y)
{
  basic->SetFakeMoleEnvironment(mole_x, mole_y);
}



/**
 * @return the mass density [g cm^-3] of material
 */
double PMI_Benchmark_Basic::Density       (const double Tl) const
{ return basic->Density(Tl*K)/(g/(cm*cm*cm)); }


/**
 * @return the \p relative \p permittivity of material
 */
double PMI_Benchmark_Basic::Permittivity  () const
{ return basic->Permittivity(); }

/**
 * @return the \p relative \p permeability of material
 */
double PMI_Benchmark_Basic::Permeability  () const
{ return basic->Permeability(); }

/**
 * @return the affinity energy [eV] of material
 */
double PMI_Benchmark_Basic::Affinity      (const double Tl) const
{ return basic->Affinity(Tl*K)/eV; }



