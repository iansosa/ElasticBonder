import structures
from mdhandler import Handler as MDH
from bondcalc import Bonds
import imp
import numpy as np
import vdw

from pymbd import mbd_energy as MBDcalc_Py, from_volumes

Nat=50
R0=2.4
	


Chain1 = structures.Sphere(Nat,R0)
Chain1.LoadGeometry("C320.cc1")
# Chain1 = structures.Ring(Nat,R0)
bonder = Bonds(Chain1,False,False)
bonder.CalcSaveHessianCompOur("C320")
