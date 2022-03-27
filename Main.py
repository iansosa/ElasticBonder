import structures
from mdhandler import Handler as MDH
from bondcalc import Bonds
import imp
import numpy as np
import vdw
import numpy.linalg as LA


Nat = 400
R0 = 2.4


du = 0.05
u=0
displ = []
displ.append('0.00')
for i in range(1,100):
	u=i*du
	ru = str(round(u,2))
	if len(ru) == 3:
		ru = ru + '0'
	displ.append(ru)

Fz_PW = np.zeros(len(displ))
for i in range(len(displ)):
	name = "/buckling/PW/geom_PW_-"+displ[i]+".gen"
	print(name)
	chain = structures.Sphere(Nat,R0)
	chain.LoadGeometry(name)
	chain.SaveGeometry()

	chain.RunStatic("PW")
	E = chain.GetEnergy()
	Fz_PW[i] = E
	print(Fz_PW)

with open('out/Buckling_PW_E.txt', 'w') as f:
    for i in range(len(Fz_PW)):
        f.write(displ[i]+' '+str(Fz_PW[i])+'\n')



Fz_MBD = np.zeros(len(displ))
for i in range(len(displ)):
	name = "/buckling/MBD/geom_MBD_-"+displ[i]+".gen"
	print(name)
	chain = structures.Sphere(Nat,R0)
	chain.LoadGeometry(name)
	chain.SaveGeometry()

	chain.RunStatic("MBD")
	E = chain.GetEnergy()
	Fz_MBD[i] = E
	print(Fz_MBD)

with open('out/Buckling_MBD_E.txt', 'w') as f:
    for i in range(len(Fz_MBD)):
        f.write(displ[i]+' '+str(Fz_MBD[i])+'\n')



Fz = np.zeros(len(displ))
for i in range(len(displ)):
	name = "/buckling/novdw/geom_-"+displ[i]+".gen"
	print(name)
	chain = structures.Sphere(Nat,R0)
	chain.LoadGeometry(name)
	chain.SaveGeometry()

	chain.RunStatic()
	E = chain.GetEnergy()
	Fz[i] = E
	print(Fz)

with open('out/Buckling_E.txt', 'w') as f:
    for i in range(len(Fz)):
        f.write(displ[i]+' '+str(Fz[i])+'\n')