import structures
import numpy as np
import matplotlib.pyplot as plt



Nat=22
R0=2.73

Sphere = structures.Sphere(Nat,R0)
#Sphere.SaveGeometry()
#Sphere.Optimize()
#Sphere.LoadGeometry()
#Sphere.ShowStruct()
#Sphere.ShowR0s()
#Sphere.SaveDistances()

#print(Sphere.R0)
R0vec = []

for i in range(69,101):
	Sphere.SetPos(i,R0)
	Sphere.SaveGeometry()
	Sphere.Optimize()
	Sphere.LoadGeometry()
	R0vec.append(Sphere.R0)

with open('R0.txt', 'w') as f:
    for i in range(len(R0vec)):
        f.write(str(i+69)+' '+str(R0vec[i])+'\n')
