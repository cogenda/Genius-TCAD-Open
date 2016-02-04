import os
from PMI_Benchmark import *


band = PMI_Benchmark_Band("/home/gdiso/develop/genius/lib/", "AlGaAs", "Default")

band.set_doping(0.0, 1e19)

for mole_x in range(0, 100):
  mole_x = mole_x*0.01
  band.set_mole(mole_x, 0.0)
  print mole_x, band.Eg(300)

