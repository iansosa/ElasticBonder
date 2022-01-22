import structures
import numpy as np
import matplotlib.pyplot as plt


Nat=50
R0=1
Sphere = structures.Sphere(Nat,1)
Sphere.SaveGeometry()
Sphere.ShowR0s()
Sphere.LoadGeometry()
Sphere.ShowR0s()
Sphere.SaveGeometry()
Sphere.ShowStruct()