import os
from PMI_Benchmark import *


thermal = PMI_Benchmark_Thermal("/home/gdiso/develop/genius/lib/", "Cu", "Default")
print thermal.HeatCapacity(300)
print thermal.HeatConduction(300)
