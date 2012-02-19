#include <iostream>
#include <cmath>

#include "PMI_benchmark_mob.h"

#include <dlfcn.h>
#define LDFUN dlsym


PMI_Benckmark_Mob::PMI_Benckmark_Mob(const std::string &material, const std::string &model)
  :dll_file(0)
{
  cm = 1e6;
  s  = 1e12;
  V  = 1.0;
  C  = 1.0/1.602176462e-19;
  K  = 1.0/300;

  // load material lib
  std::string filename =  std::string(getenv("GENIUS_DIR")) + "/lib/lib" + material + ".so";
  dll_file = dlopen(filename.c_str(), RTLD_LAZY|RTLD_DEEPBIND);
  if(!dll_file)
  {
    std::cerr<<"Open material file lib"<< material <<".so error." << '\n';
    std::cerr<<"Error code: " << dlerror() << '\n';
    throw(1);
  }

  // open mobility model
  std::string model_fun_name = "PMIS_" + material + "_Mob_" + model;
  PMIS_Mobility*      (*wmob)      (const PMI_Environment& env);
  wmob  =  (PMIS_Mobility* (*) (const PMI_Environment& env))LDFUN(dll_file, model_fun_name.c_str());

  PMI_Environment env(100*cm, s, V, C, K);
  mob = wmob(env);
}


PMI_Benckmark_Mob::~PMI_Benckmark_Mob()
{
  delete mob;
  dlclose( dll_file );
}



void PMI_Benckmark_Mob::set_doping(double Na, double Nd)
{
  mob->SetFakeDopingEnvironment(Na*pow(cm,-3), Nd*pow(cm,-3));
}


void PMI_Benckmark_Mob::set_mole(double mole_x, double mole_y)
{
  mob->SetFakeMoleEnvironment(mole_x, mole_y);
}


double PMI_Benckmark_Mob::mob_electron(double p, double n, double Ep, double Et, double T)
{
  return mob->ElecMob(p*pow(cm,-3), n*pow(cm,-3), T*K, Ep*(V/cm), Et*(V/cm), T*K)/(cm*cm/V/s);
}


double PMI_Benckmark_Mob::mob_hole(double p, double n, double Ep, double Et, double T)
{
  return mob->HoleMob(p*pow(cm,-3), n*pow(cm,-3), T*K, Ep*(V/cm), Et*(V/cm), T*K)/(cm*cm/V/s);
}



