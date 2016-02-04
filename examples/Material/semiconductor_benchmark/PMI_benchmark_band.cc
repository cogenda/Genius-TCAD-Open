#include <iostream>
#include <cmath>

#include "PMI_benchmark_band.h"

#include <dlfcn.h>
#define LDFUN dlsym


PMI_Benchmark_Band::PMI_Benchmark_Band(const std::string &path, const std::string &material, const std::string &model)
    :dll_file(0)
{
  cm = 1e6;
  s  = 1e12;
  V  = 1.0;
  C  = 1.0/1.602176462e-19;
  K  = 1.0/300;
  eV = 1.0;

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
  std::string model_fun_name = "PMIS_" + material + "_BandStructure_" + model;
  PMIS_BandStructure*      (*wband)      (const PMI_Environment& env);
  wband  =  (PMIS_BandStructure* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());

  PMI_Environment env(100*cm, s, V, C, K);
  band = wband(env);
}


PMI_Benchmark_Band::~PMI_Benchmark_Band()
{
  delete band;
  dlclose( dll_file );
}


bool PMI_Benchmark_Band::calibrate_real_parameter(const std::string & var_name, double var_value)
{
  if(band->calibrate_real_parameter(var_name, var_value)) return true;
  return false;
}

bool PMI_Benchmark_Band::calibrate_string_parameter(const std::string & var_name, const std::string &var_value)
{
  if(band->calibrate_string_parameter(var_name, var_value)) return true;
  return false;
}


void PMI_Benchmark_Band::set_doping(double Na, double Nd)
{
  band->SetFakeDopingEnvironment(Na*pow(cm,-3), Nd*pow(cm,-3));
}


void PMI_Benchmark_Band::set_mole(double mole_x, double mole_y)
{
  band->SetFakeMoleEnvironment(mole_x, mole_y);
}


/**
 * @return band gap of semiconductor
 */
double   PMI_Benchmark_Band::Eg(const double Tl)
{ return band->Eg(Tl*K)/eV; }

/**
 * @return band gap narrowing due to heavy doping
 */
double PMI_Benchmark_Band::EgNarrow(const double p, const double n, const double Tl)
{ return band->EgNarrow(p*pow(cm,-3), n*pow(cm,-3), Tl*K)/eV; }

/**
 * @return effective density of states in the conduction band
 */
double PMI_Benchmark_Band::Nc(const double Tl)
{ return band->Nc(Tl*K)/pow(cm,-3); }

/**
 * @return effective density of states in the valence band
 */
double PMI_Benchmark_Band::Nv(const double Tl)
{ return band->Nv(Tl*K)/pow(cm,-3); }


/**
 * @return intrinsic carrier concentration
 */
double PMI_Benchmark_Band::ni(const double Tl)
{ return band->ni(Tl*K)/pow(cm,-3); }

/**
 * @return effective intrinsic carrier concentration
 */
double PMI_Benchmark_Band::nie(const double p, const double n, const double Tl)
{ return band->nie(p*pow(cm,-3), n*pow(cm,-3), Tl*K)/pow(cm,-3); }


/**
 * @return direct Recombination rate
 */
double PMI_Benchmark_Band::CDIR           (const double Tl)
{
  return band->CDIR(Tl*K)/(pow(cm,3)/s);
}


/**
 * @return electron lift time in SHR Recombination
 */
double PMI_Benchmark_Band::TAUN(const double Tl)
{ return band->TAUN(Tl*K)/s; }

/**
 * @return hole lift time in SHR Recombination
 */
double PMI_Benchmark_Band::TAUP(const double Tl)
{ return band->TAUP(Tl*K)/s; }


/**
 * @return electron Auger Recombination rate
 */
double PMI_Benchmark_Band::AUGERN(const double p, const double n, const double Tl)
{
  return band->AUGERN(p*pow(cm,-3), n*pow(cm,-3), Tl*K)/(pow(cm,6)/s);
}

/**
 * @return hole Auger Recombination rate
 */
double PMI_Benchmark_Band::AUGERP(const double p, const double n, const double Tl)
{
  return band->AUGERP(p*pow(cm,-3), n*pow(cm,-3), Tl*K)/(pow(cm,6)/s);
}

