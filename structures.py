import numpy as np
import matplotlib.pyplot as plot
import sys
import subprocess

from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d import proj3d


class Sphere():

    def __init__(self,Nat,R0): #initializes equally distributed points on the surface of a sphere with interatomic distance R0 (Bohr)
        self.Nat = Nat
        self.R0s = None
        self.R0 = R0 #interatomic distance in bohr
        self.widths = None
        self.R0neighbours = 2 #number of closest neighbours to calculate the R0 approximation
        if self.Nat <= self.R0neighbours:
            self.R0neighbours=self.Nat-1
        self.SetPos(self.Nat,self.R0)

    def SetPos(self,Nat,R0): #creates the sphere with N equidistant atoms and interatomic distance R0
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

    def Pos(self): #returns the positions of every atom
        return self.x, self.y, self.z

    def ShowStruct(self): #displays the 3D structure of the sphere

        fig = plot.figure(figsize=(5,5))
        ax = fig.gca(projection='3d')


        ###scaling
        x_scale=1
        y_scale=1
        z_scale=1.3

        scale=np.diag([x_scale, y_scale, z_scale, 1.0])
        scale=scale*(1.0/scale.max())
        scale[3,3]=1.0

        def short_proj():
          return np.dot(Axes3D.get_proj(ax), scale)

        ax.get_proj=short_proj
        ###end scaling

        ax.scatter(self.x,self.y,self.z);
        plot.show()

    def ShowR0s(self): #displays a list of the R0 estimations for every atom
        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, self.R0s,color="black")
        plot.show()

    def ShowWidths(self): #displays a list of the errors in the R0 estimations for every atom
        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, self.widths,color="black")
        plot.show()

    def ShowDistances(self,idx): #displays a list of the errors in the R0 estimations for every atom
        dist=self.Distances(idx)
        dist.sort()

        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, dist,color="black")
        plot.show()

    def SaveR0s(self): #Saves R0s to file
        print("Saving R0s..")
        with open('out/R0s.txt', 'w') as f:
            for i in range(len(self.R0s)):
                f.write(str(i)+' '+str(self.R0s[i])+'\n')

    def Distances(self,idx): #returns a list of all the interatomic distances from atom idx
        dist = []
        for i in range(self.Nat):
            dist.append(np.sqrt((self.x[idx]-self.x[i])**2+(self.y[idx]-self.y[i])**2+(self.z[idx]-self.z[i])**2))
        return dist

    def SaveDistances(self): #Saves all the distances (Bohr) for each atom to a file
        dist = []
        for i in range(self.Nat):
            distances = self.Distances(i)
            distances.sort()
            dist.append(distances)

        print("Saving distances..")
        with open('Distances.txt', 'w') as f:
            for i in range(self.Nat):
                f.write(str(i)+' ')
                for k in range(self.Nat):
                    f.write(str(dist[k][i])+' ')
                f.write('\n')


    def GetR0s(self,Nneighbours): #returns a list of R0 estimations and errors from every atom considering Nneighbours closest neighbours
        print("Calculating R0s..")
        R0= []
        width= []
        for i in range(self.Nat):
            dist=self.Distances(i)
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
        with open('DFTB+/sphere.gen', 'w') as f:
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
                file = open("DFTB+/sphere.gen", "r+")
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
            self.R0 = np.median(self.R0s)
        else:
            try:
                file = open("DFTB+/geom.out.xyz", "r+")
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

            arr_t = np.array(geometry).T/angstrom
            geometry = arr_t.tolist()

            self.x = geometry[0]
            self.y = geometry[1]
            self.z = geometry[2]
            self.R0s, self.widths = self.GetR0s(self.R0neighbours)
            self.R0 = np.mean(self.R0s)

    def Optimize(self):
        subprocess.run("./dftbOpt.sh", shell=True)



class Ring():

    def __init__(self,Nat,R0): #initializes equally distributed points on the surface of a sphere with interatomic distance R0 (Bohr)
        self.Nat = Nat
        self.R0s = None
        self.R0 = R0 #interatomic distance in bohr
        self.widths = None
        self.R0neighbours = 2 #number of closest neighbours to calculate the R0 approximation
        if self.Nat <= self.R0neighbours:
            self.R0neighbours=self.Nat-1
        self.SetPos(self.Nat,self.R0)

    def SetPos(self,Nat,R0): #creates the sphere with N equidistant atoms and interatomic distance R0
        self.Nat = Nat
        self.R0 = R0

        indices = np.arange(0, Nat, dtype=float)

        Dtheta = 2*np.pi /self.Nat

        distance = np.sqrt(2-2*np.cos(Dtheta))
        Radius = self.R0/distance
        self.x, self.y, self.z = Radius*np.cos(indices*Dtheta) ,Radius*np.sin(indices*Dtheta), indices*0;
        self.R0s, self.widths = self.GetR0s(self.R0neighbours)

    def Pos(self): #returns the positions of every atom
        return self.x, self.y, self.z

    def ShowStruct(self): #displays the 3D structure of the sphere

        fig = plot.figure(figsize=(5,5))
        ax = fig.gca(projection='3d')


        ###scaling
        x_scale=1
        y_scale=1
        z_scale=1.3

        scale=np.diag([x_scale, y_scale, z_scale, 1.0])
        scale=scale*(1.0/scale.max())
        scale[3,3]=1.0

        def short_proj():
          return np.dot(Axes3D.get_proj(ax), scale)

        ax.get_proj=short_proj
        ###end scaling

        ax.scatter(self.x,self.y,self.z);
        plot.show()

    def ShowR0s(self): #displays a list of the R0 estimations for every atom
        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, self.R0s,color="black")
        plot.show()

    def ShowWidths(self): #displays a list of the errors in the R0 estimations for every atom
        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, self.widths,color="black")
        plot.show()

    def ShowDistances(self,idx): #displays a list of the errors in the R0 estimations for every atom
        dist=self.Distances(idx)
        dist.sort()

        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, dist,color="black")
        plot.show()

    def SaveR0s(self): #Saves R0s to file
        print("Saving R0s..")
        with open('out/R0s.txt', 'w') as f:
            for i in range(len(self.R0s)):
                f.write(str(i)+' '+str(self.R0s[i])+'\n')

    def Distances(self,idx): #returns a list of all the interatomic distances from atom idx
        dist = []
        for i in range(self.Nat):
            dist.append(np.sqrt((self.x[idx]-self.x[i])**2+(self.y[idx]-self.y[i])**2+(self.z[idx]-self.z[i])**2))
        return dist

    def SaveDistances(self): #Saves all the distances (Bohr) for each atom to a file
        dist = []
        for i in range(self.Nat):
            distances = self.Distances(i)
            distances.sort()
            dist.append(distances)

        print("Saving distances..")
        with open('Distances.txt', 'w') as f:
            for i in range(self.Nat):
                f.write(str(i)+' ')
                for k in range(self.Nat):
                    f.write(str(dist[k][i])+' ')
                f.write('\n')


    def GetR0s(self,Nneighbours): #returns a list of R0 estimations and errors from every atom considering Nneighbours closest neighbours
        print("Calculating R0s..")
        R0= []
        width= []
        for i in range(self.Nat):
            dist=self.Distances(i)
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
        with open('DFTB+/sphere.gen', 'w') as f:
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
                file = open("DFTB+/sphere.gen", "r+")
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
            self.R0 = np.median(self.R0s)
        else:
            try:
                file = open("DFTB+/geom.out.xyz", "r+")
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

            arr_t = np.array(geometry).T/angstrom
            geometry = arr_t.tolist()

            self.x = geometry[0]
            self.y = geometry[1]
            self.z = geometry[2]
            self.R0s, self.widths = self.GetR0s(self.R0neighbours)
            self.R0 = np.mean(self.R0s)

    def Optimize(self):
        subprocess.run("./dftbOpt.sh", shell=True)


class Chain():

    def __init__(self,Nat,R0): #initializes equally distributed points on the surface of a sphere with interatomic distance R0 (Bohr)
        self.Nat = Nat
        self.R0s = None
        self.R0 = R0 #interatomic distance in bohr
        self.widths = None
        self.R0neighbours = 1 #number of closest neighbours to calculate the R0 approximation
        if self.Nat <= self.R0neighbours:
            self.R0neighbours=self.Nat-1
        self.SetPos(self.Nat,self.R0)

    def SetPos(self,Nat,R0): #creates the sphere with N equidistant atoms and interatomic distance R0
        self.Nat = Nat
        self.R0 = R0

        indices = np.arange(0, self.Nat, dtype=float)
        dx = []
        dx = R0*(0.98+0.02*((indices-float(self.Nat)/2.0)/(float(self.Nat)/2.0))*((indices-float(self.Nat)/2.0)/(float(self.Nat)/2.0)))

        #for i in range(self.Nat):
        #    if (i%2 == 1):
        #        dx.append(self.R0)
        #    else:
        #        dx.append(self.R0*0.98)

        self.x = []

        for i in range(len(dx)):
            cummulative = 0
            for k in range(i):
                cummulative = cummulative + dx[k]
            self.x.append(cummulative)


        self.y, self.z =0*indices, 0*indices;
        self.R0s, self.widths = self.GetR0s(self.R0neighbours)

    def Pos(self): #returns the positions of every atom
        return self.x, self.y, self.z

    def ShowStruct(self): #displays the 3D structure of the sphere

        fig = plot.figure(figsize=(5,5))
        ax = fig.gca(projection='3d')


        ###scaling
        x_scale=1
        y_scale=1
        z_scale=1.3

        scale=np.diag([x_scale, y_scale, z_scale, 1.0])
        scale=scale*(1.0/scale.max())
        scale[3,3]=1.0

        def short_proj():
          return np.dot(Axes3D.get_proj(ax), scale)

        ax.get_proj=short_proj
        ###end scaling

        ax.scatter(self.x,self.y,self.z);
        plot.show()

    def ShowR0s(self): #displays a list of the R0 estimations for every atom
        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, self.R0s,color="black")
        plot.show()

    def ShowWidths(self): #displays a list of the errors in the R0 estimations for every atom
        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, self.widths,color="black")
        plot.show()

    def ShowDistances(self,idx): #displays a list of the errors in the R0 estimations for every atom
        dist=self.Distances(idx)
        dist.sort()

        x = np.linspace(1,self.Nat,self.Nat)
        plot.scatter(x, dist,color="black")
        plot.show()

    def SaveR0s(self): #Saves R0s to file
        print("Saving R0s..")
        with open('out/R0s.txt', 'w') as f:
            for i in range(len(self.R0s)):
                f.write(str(i)+' '+str(self.R0s[i])+'\n')

    def Distances(self,idx): #returns a list of all the interatomic distances from atom idx
        dist = []
        for i in range(self.Nat):
            dist.append(np.sqrt((self.x[idx]-self.x[i])**2+(self.y[idx]-self.y[i])**2+(self.z[idx]-self.z[i])**2))
        return dist

    def SaveDistances(self): #Saves all the distances (Bohr) for each atom to a file
        dist = []
        for i in range(self.Nat):
            distances = self.Distances(i)
            distances.sort()
            dist.append(distances)

        print("Saving distances..")
        with open('Distances.txt', 'w') as f:
            for i in range(self.Nat):
                f.write(str(i)+' ')
                for k in range(self.Nat):
                    f.write(str(dist[k][i])+' ')
                f.write('\n')


    def GetR0s(self,Nneighbours): #returns a list of R0 estimations and errors from every atom considering Nneighbours closest neighbours
        print("Calculating R0s..")
        R0= []
        width= []
        for i in range(self.Nat-1):
            R0.append(self.x[i+1]-self.x[i])
        R0.append(self.x[self.Nat-1]-self.x[self.Nat-2])
        return R0, width

    def SaveGeometry(self): #saves the geometry to a gen file in angstroms
        print("Saving geometry..")
        angstrom = 0.529177249
        with open('DFTB+/sphere.gen', 'w') as f:
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
                file = open("DFTB+/sphere.gen", "r+")
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
            self.R0 = np.median(self.R0s)
        else:
            try:
                file = open("DFTB+/geom.out.xyz", "r+")
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

            arr_t = np.array(geometry).T/angstrom
            geometry = arr_t.tolist()

            self.x = geometry[0]
            self.y = geometry[1]
            self.z = geometry[2]
            self.R0s, self.widths = self.GetR0s(self.R0neighbours)
            self.R0 = np.mean(self.R0s)

    def Optimize(self):
        subprocess.run("./dftbOpt.sh", shell=True)