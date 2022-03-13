import structures
from mdhandler import Handler as MDH
from bondcalc import Bonds
import imp
import numpy as np
import vdw
import numpy.linalg as LA

imp.bucklingtest("PW")

# Nat = 5
# R0 = 2.4
# chain = structures.Sphere(Nat,R0)
# chain.LoadGeometry("geom_MBD_-0.3.gen")
# chain.ShowStruct()