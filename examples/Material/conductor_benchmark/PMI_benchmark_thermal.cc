#include <iostream>
#include <cmath>

#include "PMI_benchmark_thermal.h"

#include <dlfcn.h>
#define LDFUN dlsym


PMI_Benchmark_Thermal::PMI_Benchmark_Thermal(const std::string &path, const std::string &material, const std::string &model)
  :dll_file(0)
{
  cm = 1e6;
  s  = 1e12;
  V  = 1.0;
  C  = 1.0/1.602176462e-19;
  K  = 1.0/300;


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

  // open thermal model
  std::string model_fun_name = "PMIC_" + material + "_Thermal_Default";
  PMIC_Thermal*       (*wthermal)  (const PMI_Environment& env);
  wthermal  = (PMIC_Thermal* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());

  PMI_Environment env(100*cm, s, V, C, K);
  thermal = wthermal(env);
}


PMI_Benchmark_Thermal::~PMI_Benchmark_Thermal()
{
  delete thermal;
  dlclose( dll_file );
}


bool PMI_Benchmark_Thermal::calibrate_real_parameter(const std::string & var_name, double var_value)
{
  if(thermal->calibrate_real_parameter(var_name, var_value)) return true;
  return false;
}

bool PMI_Benchmark_Thermal::calibrate_string_parameter(const std::string & var_name, const std::string &var_value)
{
  if(thermal->calibrate_string_parameter(var_name, var_value)) return true;
  return false;
}


/**
 * @return the heat capacity [J/g/K] of the material
 */
double PMI_Benchmark_Thermal::HeatCapacity  (const double Tl) const
{
  return thermal->HeatCapacity(Tl*K)/(J/g/K);
}

/**
 * @return the heat conduction [W/cm/K] of the material
 */
double PMI_Benchmark_Thermal::HeatConduction(const double Tl) const
{
  return thermal->HeatConduction(Tl*K)/(J/s/cm/K);
}




