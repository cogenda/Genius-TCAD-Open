import os
from PMI_Benchmark import *


basic = PMI_Benchmark_Basic("/home/gdiso/develop/genius/lib/", "AlGaAs", "Default")
print basic.Density(300)
print basic.Affinity(300)
