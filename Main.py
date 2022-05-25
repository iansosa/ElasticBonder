import structures
from mdhandler import Handler as MDH
from bondcalc import Bonds
import imp
import numpy as np
import vdw
import numpy.linalg as LA

imp.buckling_stats()




# Nat = 100
# R0 = 2.4
# chain = structures.Chain(Nat,R0)
# chain.SaveGeometry()
# chain.RunOptimize("PW")
# chain.LoadGeometry()

# chain2 = structures.Chain(Nat,R0)
# chain2.SaveGeometry()
# chain2.RunOptimize("PW")
# chain2.LoadGeometry()

# chain2.MoveAll([0,0,37.7945197720])

# chain.add(chain2)
# chain.SaveGeometry()
# chain.ShowStruct()
# chain.RunOptimize("PW",[0,99,100,199])
# chain.LoadGeometry()
# chain.SaveGeometry("chain_PW")

# imp.UHMWPE_comp_stats()

#########################################################################################

# Nat = 100
# R0 = 2.4
# chain = structures.Sphere(Nat,R0)

# chain.LoadGeometry("UHMWPE/UHMWPE_MBD.gen")
# chain.SaveGeometry()
# chain.RunOptimize("MBD")
# chain.LoadGeometry("geo_end.gen")
# chain.SaveGeometry("_UHMWPE_MBD")

# chain.LoadGeometry("UHMWPE/UHMWPE.gen")
# chain.SaveGeometry()
# chain.RunOptimize()
# chain.LoadGeometry("geo_end.gen")
# chain.SaveGeometry("_UHMWPE_novdw")

# chain.LoadGeometry("UHMWPE/UHMWPE_PW.gen")
# chain.SaveGeometry()
# chain.RunOptimize("PW")
# chain.LoadGeometry("geo_end.gen")
# chain.SaveGeometry("_UHMWPE_PW")
