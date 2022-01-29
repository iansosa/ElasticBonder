import numpy as np
from numpy import random
from geohandler import Handler

class Sphere(Handler):

    def __init__(self, Nat,R0):
        super().__init__(Nat,R0)

    def SetPos(self,Nat,R0): #creates a sphere with N equidistant atoms and interatomic distance R0
        self.Nat = Nat
        self.R0 = R0

        indices = np.arange(0, Nat, dtype=float) + 0.5
        phi = np.arccos(1 - 2*indices/Nat)
        theta = np.pi * (1 + 5**0.5) * indices
        self.x, self.y, self.z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi);
        self.R0s, self.widths = self.GetR0s(self.R0neighbours)
        norm=np.mean(self.R0s)/self.R0

        self.x = self.x/norm
        self.y = self.y/norm
        self.z = self.z/norm
        self.R0s = self.R0s/norm
        self.widths = self.widths/norm


class Ring(Handler):

    def __init__(self, Nat,R0):
        super().__init__(Nat,R0)

    def SetPos(self,Nat,R0): #creates a ring with N equidistant atoms and interatomic distance R0
        self.Nat = Nat
        self.R0 = R0

        indices = np.arange(0, Nat, dtype=float)
        for i in range(len(indices)):
            indices[i]=indices[i]+random.rand()/(10*self.R0)

        Dtheta = 2*np.pi /self.Nat

        distance = np.sqrt(2-2*np.cos(Dtheta))
        Radius = self.R0/distance
        self.x, self.y, self.z = Radius*np.cos(indices*Dtheta) ,Radius*np.sin(indices*Dtheta), indices*0;
        self.R0s, self.widths = self.GetR0s(self.R0neighbours)

    def GetR0s(self,Nneighbours): #returns a list of R0 estimations and errors from every atom considering Nneighbours closest neighbours
        print("Calculating R0s..")
        R0= []
        width= []
        for i in range(self.Nat-1):
            R0.append(np.sqrt((self.x[i+1]-self.x[i])**2+(self.y[i+1]-self.y[i])**2+(self.z[i+1]-self.z[i])**2))
        R0.append(np.sqrt((self.x[self.Nat-1]-self.x[0])**2+(self.y[self.Nat-1]-self.y[0])**2+(self.z[self.Nat-1]-self.z[0])**2))
        return R0, width


class Chain(Handler):

    def __init__(self, Nat,R0):
        super().__init__(Nat,R0)

    def SetPos(self,Nat,R0): #creates a chain with N equidistant atoms and interatomic distance R0
        self.Nat = Nat
        self.R0 = R0

        indices = np.arange(0, self.Nat, dtype=float)
        dx = []
        dx = R0*(0.98+0.02*((indices-float(self.Nat)/2.0)/(float(self.Nat)/2.0))*((indices-float(self.Nat)/2.0)/(float(self.Nat)/2.0)))

        self.x = []

        for i in range(len(dx)):
            cummulative = 0
            for k in range(i):
                cummulative = cummulative + dx[k]
            self.x.append(cummulative)


        self.y, self.z =0*indices, 0*indices;
        self.R0s, self.widths = self.GetR0s(self.R0neighbours)

    def GetR0s(self,Nneighbours): #returns a list of R0 estimations and errors from every atom considering Nneighbours closest neighbours
        print("Calculating R0s..")
        R0= []
        width= []
        for i in range(self.Nat-1):
            R0.append(self.x[i+1]-self.x[i])
        R0.append(self.x[self.Nat-1]-self.x[self.Nat-2])
        return R0, width