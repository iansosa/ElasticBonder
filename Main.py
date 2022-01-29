import structures
import numpy as np
import matplotlib.pyplot as plt



Nat=60
R0=2.72

Sphere = structures.Sphere(Nat,R0)

Sphere.LoadGeometry("Fullerene-C60.sdf")
Sphere.ShowStruct()
Sphere.ShowR0s()
Sphere.SaveGeometry()
Sphere.Optimize()
Sphere.LoadGeometry()
Sphere.ShowStruct()
Sphere.ShowR0s()

