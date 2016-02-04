import os
from PMI_Benchmark import *


mob = PMI_Benchmark_Mob("/home/gdiso/develop/genius/lib/", "Si", "Lucent")
mob.calibrate_real_parameter("NRFP.UM", 3e17)

mob.set_doping(1e18, 0.0)
print mob.mob_electron(1e18, 0, 0, 0, 300)

