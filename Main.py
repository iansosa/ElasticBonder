import structures
from mdhandler import Handler as MDH
import numpy as np
import matplotlib.pyplot as plt
from bondcalc import Bonds
import time
import sys
import utils
import numpy as np
from scipy.optimize import curve_fit
from functools import partial

Nat=10
R0=2.4

# Chain1 = structures.Chain(Nat,R0)
# Chain1.SaveGeometry()
# Chain1.RunOptimize()
# Chain1.LoadGeometry()

# Chain2 = structures.Chain(Nat,R0)
# Chain2.SaveGeometry()
# Chain2.RunOptimize()
# Chain2.LoadGeometry()

# Chain2.MoveAll([0,0,30])

# Chain1.add(Chain2)

# Chain1.SaveGeometry()
# Chain1.RunStatic()

# F = Chain1.GetForces()
# F_lower = 0
# for k in range(10):
#     F_lower = F_lower + F[k][2]
# F_upper = 0

# for k in range(10,20):
#     print(k)
#     F_upper = F_upper + F[k][2]
# print(F_upper,F_lower)


for i in range(10):
	utils.StaticOverEvolve(i*30+50)



# Chain = structures.Chain(Nat,R0)

# md = MDH(Chain,False)
# md.RunMD(100,300,[0,9])
# md.LoadEvolution()
