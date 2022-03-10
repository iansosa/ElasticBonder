import structures
from mdhandler import Handler as MDH
from bondcalc import Bonds
import imp
import numpy as np
import vdw

Nat=50
R0=2.4
	


Chain = structures.Sphere(Nat,R0)
Chain.LoadGeometry("nanotube.txt")
Chain.SaveGeometry()
Chain.RunOptimize()
Chain.LoadGeometry()
Chain.SaveGeometry("_nanotube")

bonder = Bonds(Chain,False,False)
bonder.FitEnergy(1000,"three")


