import os
from PMI_Benchmark import *


optical = PMI_Benchmark_Optical("/home/gdiso/develop/genius/lib/", "Si", "Green")
#optical.set_doping(1e18, 0.0)
#optical.set_mole(0.71, 0.0)
T=249

for i in range(0, 400):
  wave = 0.2 + i/300.0
  print wave, optical.n(wave, T), optical.alpha(wave, T)
  #optical.n(wave, T)

