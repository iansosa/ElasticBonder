import numpy as np
import matplotlib.pyplot as plot
import sys

class Sphere():

    def __init__(self,Nat,R0): #initializes equally distributed points on the surface of a sphere with interatomic distance R0 (Bohr)
        self.Nat = Nat
        self.R0s = None
        self.R0 = R0 #interatomic distance in bohr
        self.widths = None
        self.R0neighbours = 3 #number of closest neighbours to calculate the R0 approximation
        if self.Nat <= self.R0neighbours:
            self.R0neighbours=self.Nat-1
        self.SetPos(self.Nat)

    def SetPos(self,Nat): #creates the sphere with N equidistant atoms and interatomic distance R0
        self.Nat = Nat

        indices = np.arange(0, Nat, dtype=float) + 0.5
        phi = np.arccos(1 - 2*indices/Nat)
        theta = np.pi * (1 + 5**0.5) * indices
        self.x, self.y, self.z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi);
        self.R0s, self.widths = self.GetR0s(self.R0neighbours)
        norm=np.median(self.R0s)/self.R0

        self.x = self.x/norm
        self.y = self.y/norm
        self.z = self.z/norm
        self.R0s = self.R0s/norm
        self.widths = self.widths/norm

    def Pos(self): #returns the positions of every atom
        return self.x, self.y, self.z

    def ShowStruct(self): #displays the 3D structure of the sphere
        plot.figure().add_subplot(111, projection='3d').scatter(self.x,self.y,self.z);
        plot.show()

    def ShowR0s(self): #displays a list of the R0 estimations for every atom
        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, self.R0s,color="black")
        plot.show()

    def ShowWidths(self): #displays a list of the errors in the R0 estimations for every atom
        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, self.widths,color="black")
        plot.show()

    def Distances(self,idx): #returns a list of all the interatomic distances from atom idx
        dist = []
        for i in range(self.Nat):
            dist.append(np.sqrt((self.x[idx]-self.x[i])**2+(self.y[idx]-self.y[i])**2+(self.z[idx]-self.z[i])**2))
        return dist

    def GetR0s(self,Nneighbours): #returns a list of R0 estimations and errors from every atom considering Nneighbours closest neighbours
        print("Calculating R0s..")
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

    def SaveGeometry(self): #saves the geometry to a gen file in angstroms
        print("Saving geometry..")
        angstrom = 0.529177249
        with open('sphere.gen', 'w') as f:
            f.write(str(self.Nat)+' C\n')
            f.write('  C\n')
            f.write('\n')
            for i in range(self.Nat):
                f.write('  '+str(i+1)+' 1  '+str(angstrom*self.x[i])+' '+str(angstrom*self.y[i])+' '+str(angstrom*self.z[i])+'\n')

    def LoadGeometry(self,xyz=True): #loads the geometry from a gen or xyz file in angstroms and transforms it into Bohr
        print("Loading geometry..")
        angstrom = 0.529177249

        if xyz==False:
            try:
                file = open("sphere.gen", "r+")
            except OSError:
                print ("Could not find file: sphere.gen")
                sys.exit()

            lines = file.readlines()

            aux = lines[0].split(' ')
            aux = list(filter(lambda x: x != '', aux))
            self.Nat = int(aux[0])
            if self.Nat <= self.R0neighbours:
                self.R0neighbours=self.Nat-1
            lines = lines[3:]

            geometry = []
            for i in range(len(lines)):
                a = lines[i].split(' ')
                a = list(filter(lambda x: x != '', a))
                a = list(map(float, a[2:]))
                geometry.append(a)

            arr_t = np.array(geometry).T/angstrom
            geometry = arr_t.tolist()

            self.x = geometry[0]
            self.y = geometry[1]
            self.z = geometry[2]
            self.R0s, self.widths = self.GetR0s(self.R0neighbours)
        else:
            try:
                file = open("geom.out.xyz", "r+")
            except OSError:
                print ("Could not find file: geom.out.xyz")
                sys.exit()

            lines = file.readlines()

            aux = lines[0].split(' ')
            aux = list(filter(lambda x: x != '', aux))
            self.Nat = int(aux[0])
            if self.Nat <= self.R0neighbours:
                self.R0neighbours=self.Nat-1
            lines = lines[2:]

            geometry = []
            for i in range(len(lines)):
                a = lines[i].split(' ')
                a = list(filter(lambda x: x != '', a))
                a = list(map(float, a[1:-1]))
                geometry.append(a)

            self.x = geometry[0]
            self.y = geometry[1]
            self.z = geometry[2]
            self.R0s, self.widths = self.GetR0s(self.R0neighbours)


