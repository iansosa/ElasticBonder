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


Sphere = structures.Ring(Nat,R0)
Sphere.LoadGeometry("Ring_C10.gen")

# print(Chain.R0)
# Chain.ShowStruct()
# Chain.SaveGeometry()
# Chain.RunOptimize()
# Chain.LoadGeometry()
# print(Chain.R0)
# Chain.ShowStruct()

# Bonder = Bonds(Chain,False)
# Bonder.FitEnergy(1000)


Chain = structures.Chain(Nat,R0)

Chain.MoveAll([11.403,0,0])
Sphere.add(Chain)

md = MDH(Sphere,False)
md.RunMD(50000,10)
