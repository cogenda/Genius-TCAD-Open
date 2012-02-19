#include <iostream>
#include "PMI_benchmark_mob.h"


int main(int argc, char ** args)
{
  PMI_Benckmark_Mob * mob = new PMI_Benckmark_Mob("Si", "Lucent");
  mob->set_doping(1e18, 0.0);
  std::cout<<mob->mob_electron(1e18, 0, 0, 0, 300)<<std::endl;

  delete mob;
}

