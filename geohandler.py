import numpy as np
import matplotlib.pyplot as plot
import sys
import subprocess
import os.path
import filetypes
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

class Handler():

    def __init__(self,Nat,R0): #initialize the desired geometry with interatomic distance R0 (Bohr)
        self.Nat = Nat
        self.R0s = None
        self.R0 = R0 #interatomic distance in bohr
        self.widths = None
        self.R0neighbours = 2 #number of closest neighbours to calculate the R0 approximation
        if self.Nat <= self.R0neighbours:
            self.R0neighbours=self.Nat-1
        self.SetPos(self.Nat,self.R0)

    def SetPos(self,Nat,R0): #set the position of the geometry
        print ("SetPos Unimplemented")
        sys.exit()

    def Pos(self): #returns the positions of every atom
        return self.x, self.y, self.z

    def ShowStruct(self): #displays the 3D structure

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
        with open('out/Distances.txt', 'w') as f:
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
        with open('DFTB+/geom.gen', 'w') as f:
            f.write(str(self.Nat)+' C\n')
            f.write('  C\n')
            for i in range(self.Nat):
                f.write('  '+str(i+1)+' 1  '+str(angstrom*self.x[i])+' '+str(angstrom*self.y[i])+' '+str(angstrom*self.z[i])+'\n')

    def LoadGeometry(self,path="geom.out.xyz"): #loads the geometry from a gen, xyz or sdf file in angstroms and converts it into Bohr
        print("Loading geometry..")
        angstrom = 0.529177249
        extension = path[-3:]
        recognized = False

        if extension == "sdf":
            recognized = True
            self.Nat, geometry = filetypes.Loadsdf("SavedStructures/"+path,angstrom)
                
        if extension == "gen":
            recognized = True
            self.Nat, geometry = filetypes.Loadgen("DFTB+/"+path,angstrom)

        if extension == "xyz":
            recognized = True
            self.Nat, geometry = filetypes.Loadxyz("DFTB+/"+path,angstrom)

        if recognized == False:
            print ("Extension not recognized")
            sys.exit()

        if self.Nat <= self.R0neighbours:
            self.R0neighbours=self.Nat-1

        self.x = geometry[0]
        self.y = geometry[1]
        self.z = geometry[2]
        self.R0s, self.widths = self.GetR0s(self.R0neighbours)
        self.R0 = np.mean(self.R0s)


    def Optimize(self):
        subprocess.run("./dftbOpt.sh", shell=True)