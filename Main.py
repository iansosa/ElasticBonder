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


# Sphere = structures.Ring(Nat,R0)
# Sphere.LoadGeometry("Ring_C10.gen")


Chain = structures.Chain(Nat,R0)

md = MDH(Chain,False)
# md.RunMD(100,300)
md.LoadEvolution()

total_F = 0
for i in range(Nat):
	total_F = total_F + md.forcespmass[50][i][1]
print(total_F)

# stats = []
# for i in range(1):
# 	md = MDH(Chain,False)
# 	md.RunMD(10,i*10+30)
# 	md.LoadForces()
# 	averages, dispersions, averager = md.ComputeTempDispersions()
# 	r = 0
# 	d = 0
# 	for j in range(len(dispersions)):
# 		r = r + averager[j]
# 		d = d + dispersions[j]
# 	r = r/len(averager)
# 	d = d/len(dispersions)
# 	stats.append([r,d])
# 	print("r",r)
# 	print("i",i)
# 	print("d",d)

# with open('out/TempStats_both.txt', 'w') as f:
#     for i in range(len(stats)):
#         f.write(str(i*10+30)+' '+str(stats[i][0])+' '+str(stats[i][1])+'\n')