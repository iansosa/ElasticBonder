from numpy import pi, cos, sin, arccos, arange, sqrt
import numpy as np
import matplotlib.pyplot as plot

class Sphere():

    def __init__(self,Nat,R0):
        self.Nat = Nat
        self.R0s = None
        self.R0 = R0
        self.widths = None
        self.SetPos(self.Nat)

    def SetPos(self,Nat):
        self.Nat = Nat

        indices = arange(0, Nat, dtype=float) + 0.5
        phi = arccos(1 - 2*indices/Nat)
        theta = pi * (1 + 5**0.5) * indices
        self.x, self.y, self.z = cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi);
        self.R0s, self.widths = self.GetR0s(3)
        norm=np.median(self.R0s)/self.R0

        self.x = self.x/norm
        self.y = self.y/norm
        self.z = self.z/norm
        self.R0s = self.R0s/norm
        self.widths = self.widths/norm

    def Pos(self):
        return self.x, self.y, self.z

    def ShowStruct(self):
        plot.figure().add_subplot(111, projection='3d').scatter(self.x,self.y,self.z);
        plot.show()

    def ShowR0s(self):
        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, self.R0s,color="black")
        plot.show()

    def ShowWidths(self):
        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, self.widths,color="black")
        plot.show()

    def Distances(self,idx):
        dist = []
        for i in range(self.Nat):
            dist.append(sqrt((self.x[idx]-self.x[i])**2+(self.y[idx]-self.y[i])**2+(self.z[idx]-self.z[i])**2))
        return dist

    def GetR0s(self,Nneighbours):
        R0= []
        width= []
        for i in range(self.Nat):
            dist=Sphere.Distances(self,i)
            dist.sort()
            median = 0
            for k in range(Nneighbours):
                median= median + dist[k+1]
            median=median/Nneighbours
            R0.append(median)
            width.append(dist[Nneighbours]-dist[1])
        return R0, width

