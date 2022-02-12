import structures
import numpy as np
import matplotlib.pyplot as plt
from bondcalc import Bonds
import time
import sys

# Nat=5
# R0=2.455

# Chain = structures.Chain(Nat,R0)
# BondCalc = Bonds(Chain,True)
# BondCalc.CalcSaveBondMatrix()

# def getForces(dy):

# 	K1 = 0.717808095025
# 	K2 = 0.0873187293
# 	K3 = -0.04184426625
# 	K4 = 0.010964463525

# 	p0 = [0,dy,0]
# 	p0 = np.array(p0)
# 	p1 = [2.42258207,0,0]
# 	p1 = np.array(p1)
# 	p2 = [4.89002449,0,0]
# 	p2 = np.array(p2)
# 	p3 = [7.357468,0,0]
# 	p3 = np.array(p3)
# 	p4 = [9.7800508,0,0]
# 	p4 = np.array(p4)


# 	versor1 = p0 - p1
# 	norm1 = np.sqrt((versor1[0])**2+(versor1[1])**2+(versor1[2])**2)
# 	versor1 = versor1/norm1
# 	d1 = p0 - p1
# 	d1 = np.sqrt((d1[0])**2+(d1[1])**2+(d1[2])**2)

# 	versor2 = p0 - p2
# 	norm2 = np.sqrt((versor2[0])**2+(versor2[1])**2+(versor2[2])**2)
# 	versor2 = versor2/norm2
# 	d2 = p0 - p2
# 	d2 = np.sqrt((d2[0])**2+(d2[1])**2+(d2[2])**2)

# 	versor3 = p0 - p3
# 	norm3 = np.sqrt((versor3[0])**2+(versor3[1])**2+(versor3[2])**2)
# 	versor3 = versor3/norm3
# 	d3 = p0 - p3
# 	d3 = np.sqrt((d3[0])**2+(d3[1])**2+(d3[2])**2)

# 	versor4 = p0 - p4
# 	norm4 = np.sqrt((versor4[0])**2+(versor4[1])**2+(versor4[2])**2)
# 	versor4 = versor4/norm4
# 	d4 = p0 - p4
# 	d4 = np.sqrt((d4[0])**2+(d4[1])**2+(d4[2])**2)


# 	F = -K1*(d1-2.42258207)*versor1 - K2*(d2-4.89002449)*versor2 - K3*(d3-7.357468)*versor3 - K4*(d4-9.7800508426)*versor4
# 	return F

def getForces(dy):

	K1 = 0.71780690794
	K2 = 0.5943861564
	K3 = 0.120717892
	K4 = -0.041834921

	p0 = [0,0,0]
	p0 = np.array(p0)
	p1 = [2.42258207,dy,0]
	p1 = np.array(p1)
	p2 = [4.89002449,0,0]
	p2 = np.array(p2)
	p3 = [7.357468,0,0]
	p3 = np.array(p3)
	p4 = [9.7800508,0,0]
	p4 = np.array(p4)


	versor1 = p1 - p0
	norm1 = np.sqrt((versor1[0])**2+(versor1[1])**2+(versor1[2])**2)
	versor1 = versor1/norm1
	d1 = p1 - p0
	d1 = np.sqrt((d1[0])**2+(d1[1])**2+(d1[2])**2)

	versor2 = p1 - p2
	norm2 = np.sqrt((versor2[0])**2+(versor2[1])**2+(versor2[2])**2)
	versor2 = versor2/norm2
	d2 = p1 - p2
	d2 = np.sqrt((d2[0])**2+(d2[1])**2+(d2[2])**2)

	versor3 = p1 - p3
	norm3 = np.sqrt((versor3[0])**2+(versor3[1])**2+(versor3[2])**2)
	versor3 = versor3/norm3
	d3 = p1 - p3
	d3 = np.sqrt((d3[0])**2+(d3[1])**2+(d3[2])**2)

	versor4 = p1 - p4
	norm4 = np.sqrt((versor4[0])**2+(versor4[1])**2+(versor4[2])**2)
	versor4 = versor4/norm4
	d4 = p1 - p4
	d4 = np.sqrt((d4[0])**2+(d4[1])**2+(d4[2])**2)


	F = -K1*(d1-2.42258207)*versor1 - K2*(d2-2.46744241644)*versor2 - K3*(d3-4.934885)*versor3 - K4*(d4-7.35746876)*versor4

	return F


Nat=10
R0=2.4

Chain = structures.Sphere(Nat,R0)
Chain.SaveGeometry()
Chain.RunOptimize()
Chain.LoadGeometry()
Chain.CalcBondAngles()
Chain.ShowStruct()

sys.exit("bye")

F0 = Chain.GetForces()[1]
# F1 = Chain.GetForces()[1]

# Chain.MoveBond(0,"y",dv=0.04)
# Chain.SaveGeometry()
# Chain.RunStatic()

# F0y = Chain.GetForces()[0]
# F1y = Chain.GetForces()[1]

# Chain.MoveBond(0,"y",dv=-0.04)
# Chain.SaveGeometry()

# Chain.MoveBond(0,"x",dv=0.02)
# Chain.SaveGeometry()
# Chain.RunStatic()

# F0x = Chain.GetForces()[0]
# F1x = Chain.GetForces()[1]

# F0y = F0y-F0

with open('out/ForcesY.txt', 'w') as f:
	for i in range(1,100):
		Chain.MoveBond(1,"y",dv=0.02)
		Chain.SaveGeometry()
		Chain.RunStatic()
		F0y = Chain.GetForces()[1]
		F0y = F0y-F0
		F0yT = getForces(i*0.02)
		f.write(str(i*0.02)+" "+str(F0y[0])+" "+str(F0y[1])+" "+str(F0yT[0])+" "+str(F0yT[1])+"\n")


# print(F0x-F0, F0y-F0)
# print(F1x-F1, F1y-F1)

# print(F)


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