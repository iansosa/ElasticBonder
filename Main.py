import structures
import numpy as np
import matplotlib.pyplot as plt


Nat=64
R0=2.74
Sphere = structures.Sphere(Nat,R0)
Sphere.LoadGeometry()
Sphere.ShowStruct()
Sphere.SaveGeometry()
Sphere.SaveDistances()


