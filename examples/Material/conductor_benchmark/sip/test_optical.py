import os
from PMI_Benchmark import *


optical = PMI_Benchmark_Optical("/home/gdiso/develop/genius/lib/", "Cu", "Default")
T=249

for i in range(0, 400):
  wave = 0.2 + i/300.0
  print wave, optical.n(wave, T), optical.alpha(wave, T)
  #optical.n(wave, T)

