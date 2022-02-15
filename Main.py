import structures
import numpy as np
import matplotlib.pyplot as plt
from bondcalc import Bonds
import time
import sys
import utils
import numpy as np
from scipy.optimize import curve_fit
from functools import partial

Nat=4
R0=2.4

Chain = structures.Sphere(Nat,R0)
Chain.LoadGeometry("Graphene-C92.sdf")
Bonder = Bonds(Chain,True)
Bonder.FitEnergy(1000)
