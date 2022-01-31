import structures
import numpy as np
import matplotlib.pyplot as plt
from bondcalc import Bonds


Nat=50
R0=2.4

Chain = structures.Chain(Nat,R0)
calculator = Bonds(Chain)
calculator.CalcSaveBondMatrix()

# bonds = []
# for i in range(Nat):
# 	column = []
# 	for k in range(Nat):
# 		column.append(calculator.GetK(k,i))
# 	bonds.append(column)

# with open('out/Forces.txt', 'w') as f:
# 	for i in range(len(bonds)):
# 		f.write(str(i))
# 		for k in range(len(bonds)):
# 			f.write(' '+str(bonds[i][k]))
# 		f.write('\n')