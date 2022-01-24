import structures
import numpy as np
import matplotlib.pyplot as plt


Nat=64
R0=2.8
Sphere = structures.Sphere(Nat,R0)
Sphere.ShowStruct()
Sphere.ShowR0s()
Sphere.SaveGeometry()
Sphere.LoadGeometry()
Sphere.ShowStruct()
Sphere.ShowR0s()

