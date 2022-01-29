import structures
import numpy as np
import matplotlib.pyplot as plt



Nat=17
R0=2.7

Sphere = structures.Ring(Nat,R0)
Sphere.SetPos(Nat,R0)
print(Sphere.R0)
Sphere.ShowStruct()
Sphere.ShowR0s()
Sphere.SaveGeometry()
Sphere.Optimize()
Sphere.LoadGeometry()
Sphere.ShowStruct()
Sphere.ShowR0s()
