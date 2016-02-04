import os
import numpy as np
from PMI_Benchmark import *


basic = PMI_Benchmark_Basic("/home/gdiso/develop/genius/lib/", "ORPS", "Default")
#print basic.CurrentDensity(0.5/6e-4,250)


for V in np.linspace(0, 150, 150):
  print V, basic.CurrentDensity(V/6e-4,250), basic.CurrentDensity(V/6e-4,300), basic.CurrentDensity(V/6e-4,323), basic.CurrentDensity(V/6e-4,423), basic.CurrentDensity(V/6e-4,473)
  