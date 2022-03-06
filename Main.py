import structures
from mdhandler import Handler as MDH
import numpy as np
import matplotlib.pyplot as plt
from bondcalc import Bonds
import time
import sys
import imp
import numpy as np
from scipy.optimize import curve_fit
from functools import partial

Nat=10
R0=2.4

Sphere = structures.Chain(Nat,R0)
# Sphere.LoadGeometry("C76-Td.cc1")
Sphere.RunOptimize()
Sphere.SaveGeometry()
Sphere.UpdateR0s()
print(Sphere.R0)
Sphere.ShowStruct()
Sphere.RunHessian()

# for i in range(4):
# 	#imp.StaticOverEvolve_Chains(i*100+50,40,40)
# 	imp.CorrelationOverEvolve_Chains(i*100+50,20,20)


# Chain = structures.Chain(Nat,R0)

# md = MDH(Chain,False)
# md.RunMD(100,300,[0,9])
# md.LoadEvolution()
