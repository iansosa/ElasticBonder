import structures
import numpy as np
import matplotlib.pyplot as plt



Nat=64
R0=2.7

Chain = structures.Sphere(Nat,R0)
Chain.LoadGeometry()
print(Chain.R0)
Chain.ShowStruct()


# for i in range(3,101):
# 	Chain.SetPos(i,R0)
# 	Chain.SaveGeometry()
# 	Chain.Optimize()
# 	Chain.LoadGeometry()
# 	print("sphere",i,Chain.R0)
# 	with open('R0_sphere.txt', 'a') as f:
# 		f.write(str(i)+' '+str(Chain.R0)+'\n')


# Chain = structures.Ring(Nat,R0)

# for i in range(3,101):
# 	Chain.SetPos(i,R0)
# 	Chain.SaveGeometry()
# 	Chain.Optimize()
# 	Chain.LoadGeometry()
# 	print("ring",i,Chain.R0)
# 	with open('R0_ring.txt', 'a') as f:
# 		f.write(str(i)+' '+str(Chain.R0)+'\n')


# Chain = structures.Chain(Nat,R0)

# for i in range(3,101):
# 	Chain.SetPos(i,R0)
# 	Chain.SaveGeometry()
# 	Chain.Optimize()
# 	Chain.LoadGeometry()
# 	print("chain",i,Chain.R0)
# 	with open('R0_chain.txt', 'a') as f:
# 		f.write(str(i)+' '+str(Chain.R0)+'\n')



