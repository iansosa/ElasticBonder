import numpy as np
import matplotlib.pyplot as plot
import math
import sys
import subprocess
import shutil
import os.path
import copy
import filetypes
import numpy.linalg as LA
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from itertools import permutations 
from scipy.spatial.transform import Rotation as R

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
        self.Cutoffneighbours = self.R0*1.1 #cutoffdistance for bonded neighbours
        self.bonds = None
        self.CalcBondedNeighbours(self.Cutoffneighbours)
        self.offplane = None
        self.angles = None #every angle formed by 3 atoms. angles[i] holds all of the angles related to atom i, None if there are none.  angles[i][k][0],angles[i][k][1] indices of two atoms forming angle k. angles[i][k][2] the angle in radians. angles[i][k][3] a perpendicular vector to angles[i][k][0],angles[i][k][1]
        self.distances = None
        self.CalcBondDistances()
        self.CalcBondAngles()
        self.CalcBondOffPlane()

    def add(self,structure):
        self.Nat = self.Nat+structure.Nat
        # self.x = self.x.tolist()
        # self.y = self.y.tolist()
        # self.z = self.z.tolist()
        for i in range(structure.Nat):
            self.x.append(structure.x[i])
            self.y.append(structure.y[i])
            self.z.append(structure.z[i])
        self.CalcBondedNeighbours(self.Cutoffneighbours)
        self.CalcBondDistances()
        self.CalcBondAngles()
        self.CalcBondOffPlane()

    def SetPos(self,Nat,R0): #set the position of the geometry
        print ("SetPos Unimplemented")
        sys.exit()

    def Pos(self): #returns the positions of every atom
        return self.x, self.y, self.z
    def PosAsList(self): #returns the positions of every atom as a list
        pos = []
        for i in range(self.Nat):
            pos.append([self.x[i],self.y[i],self.z[i]])
        return pos

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
        print(dist)

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

    def Distance(self,i,j): #returns the distance between atom i and j
        return np.sqrt((self.x[i]-self.x[j])**2+(self.y[i]-self.y[j])**2+(self.z[i]-self.z[j])**2)

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

    def UpdateR0s(self): #returns a list of R0 estimations and errors from every atom considering Nneighbours closest neighbours
        self.R0s, self.widths = self.GetR0s(self.R0neighbours)
        self.R0 = np.mean(self.R0s)
        self.CalcBondedNeighbours(self.R0*1.1)
        self.CalcBondDistances()
        self.CalcBondAngles()
        self.CalcBondOffPlane()

    def SaveGeometry(self,decour="",path=None): #saves the geometry to a gen file in angstroms
        print("Saving geometry..")
        angstrom = 0.529177249
        name='DFTB+/geom'+decour+'.gen'
        if path != None:
            name = 'DFTB+/'+path+'/geom'+decour+'.gen'
        with open(name, 'w') as f:
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
            if path=="Graphene-C92.sdf":
                angstrom = angstrom * 0.71121438
            self.Nat, geometry = filetypes.Loadsdf("SavedStructures/"+path,angstrom)
            print(str(self.Nat)+" atoms loaded")
                
        if extension == "gen":
            recognized = True
            if path != "geom.out.gen":
                self.Nat, geometry = filetypes.Loadgen("SavedStructures/"+path,angstrom)
                print(str(self.Nat)+" atoms loaded")
            else:
                self.Nat, geometry = filetypes.Loadgen("DFTB+/"+path,angstrom)

        if extension == "xyz":
            recognized = True
            self.Nat, geometry = filetypes.Loadxyz_single("DFTB+/"+path,angstrom)

        if extension == "cc1":
            recognized = True
            self.Nat, geometry = filetypes.Loadcc1("SavedStructures/"+path,angstrom)

        if extension == "txt":
            recognized = True
            self.Nat, geometry = filetypes.Loadtxt("SavedStructures/"+path,angstrom)
            print(str(self.Nat)+" atoms loaded")

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
        self.CalcBondedNeighbours(self.R0*1.1)
        self.CalcBondDistances()
        self.CalcBondAngles()
        self.CalcBondOffPlane()

    def RunOptimize(self,vdw=None,static=None,read_charges=False):
        if vdw == None:
            shutil.copyfile('DFTB+/optimize.hsd', 'DFTB+/dftb_in.hsd')
        elif vdw == "MBD":
            shutil.copyfile('DFTB+/optimize_mbd.hsd', 'DFTB+/dftb_in.hsd')
        elif vdw == "PW":
            shutil.copyfile('DFTB+/optimize_pw.hsd', 'DFTB+/dftb_in.hsd')
        else:
            print ("Dispersion type not recognized")
            sys.exit()
        try:
            file = open("DFTB+/dftb_in.hsd", "r+")
        except OSError:
            print ("Could not open detailed.out file")
            sys.exit()

        lines = file.readlines()
        file.close()

        if read_charges == True:
            idx = -1
            for i in range(len(lines)):
                if lines[i].find("ReadInitialCharges") != -1:
                    idx = i
            targetline = "  ReadInitialCharges = Yes\n"
            lines[idx] = targetline

        if static != None:
            idx = -1
            for i in range(len(lines)):
                if lines[i].find("MovedAtoms") != -1:
                    idx = i
            targetline = "  MovedAtoms = !("
            for j in range(len(static)):
                targetline = targetline + str(static[j]+1) + " "
            targetline = targetline + ")"+"\n"
            lines[idx] = targetline

        with open('DFTB+/dftb_in.hsd', 'w') as f:
            for i in range(len(lines)):
                f.write(lines[i])

        subprocess.run("./dftbOpt.sh", shell=True)

    def RunStatic(self,vdw=None):
        if vdw == None:
            shutil.copyfile('DFTB+/static_calc.hsd', 'DFTB+/dftb_in.hsd')
        elif vdw == "MBD":
            shutil.copyfile('DFTB+/static_calc_mbd.hsd', 'DFTB+/dftb_in.hsd')
        elif vdw == "PW":
            shutil.copyfile('DFTB+/static_calc_pw.hsd', 'DFTB+/dftb_in.hsd')
        else:
            print ("Dispersion type not recognized")
            sys.exit()
        subprocess.run("./dftbOpt.sh", shell=True)

    def RunHessian(self,vdw=None):
        if vdw == None:
            shutil.copyfile('DFTB+/Hessian_calc.hsd', 'DFTB+/dftb_in.hsd')
        elif vdw == "MBD":
            shutil.copyfile('DFTB+/Hessian_calc_mbd.hsd', 'DFTB+/dftb_in.hsd')
        elif vdw == "PW":
            shutil.copyfile('DFTB+/Hessian_calc_pw.hsd', 'DFTB+/dftb_in.hsd')
        else:
            print ("Dispersion type not recognized")
            sys.exit()
        subprocess.run("./dftbOpt.sh", shell=True)

    def Displace(self,i,dv): #displaces atom i a dv vector distance (Bohr)
        self.x[i]=self.x[i]+dv[0]
        self.y[i]=self.y[i]+dv[1]
        self.z[i]=self.z[i]+dv[2]

    def GetVersor(self,i,j): #returns versor that points from i to j
        versor = [self.x[j]-self.x[i],self.y[j]-self.y[i],self.z[j]-self.z[i]]
        norm = np.sqrt((self.x[i]-self.x[j])**2+(self.y[i]-self.y[j])**2+(self.z[i]-self.z[j])**2)
        return np.array(versor)/norm

    def GetVector(self,i,j): #returns vector that points from i to j
        versor = [self.x[j]-self.x[i],self.y[j]-self.y[i],self.z[j]-self.z[i]]
        return np.array(versor)

    def PullBond(self,i,j,dv=0.001): #pulls atom j away from i a dv distance (Bohr)
        versor = self.GetVersor(i,j)
        versor = versor*dv

        self.x[j]=self.x[j]+versor[0]
        self.y[j]=self.y[j]+versor[1]
        self.z[j]=self.z[j]+versor[2]

    def RotateBond(self,i,j,k,rotvec,da=0.01): #rotates atom k towards j an angle da pivoting on i

        inner = np.inner(self.GetVector(i,k),self.GetVector(i,j)) / (LA.norm(self.GetVector(i,k))* LA.norm(self.GetVector(i,j)))
        rad_in = np.arccos(np.clip(inner, -1.0, 1.0))

        rotation_vector = da * rotvec
        rotation = R.from_rotvec(rotation_vector)
        v1 = rotation.apply(self.GetVector(i,k))

        inner = np.inner(v1,self.GetVector(i,j)) / (LA.norm(v1)* LA.norm(self.GetVector(i,j)))
        rad_out = np.arccos(np.clip(inner, -1.0, 1.0))

        if rad_out > rad_in:
            rotation_vector = -da * rotvec
            rotation = R.from_rotvec(rotation_vector)
            v1 = rotation.apply(self.GetVector(i,k))

        v1 = v1 + np.array([self.x[i],self.y[i],self.z[i]])
        self.x[k] = v1[0]
        self.y[k] = v1[1]
        self.z[k] = v1[2]

    def MoveBond(self,i,xyz,dv=0.001): #moves atom i in a given direction

        if xyz=="x":
            self.x[i]=self.x[i]+dv
        if xyz=="y":
            self.y[i]=self.y[i]+dv
        if xyz=="z":
            self.z[i]=self.z[i]+dv

    def MoveAll(self,dv): #moves all atoms in a given direction

        for i in range(self.Nat):
            self.x[i] = self.x[i] + dv[0]
            self.y[i] = self.y[i] + dv[1]
            self.z[i] = self.z[i] + dv[2]

    def SetDisplacements(self,dm): #moves all atoms using a displacement matrix

        for i in range(self.Nat):
            self.x[i]=self.x[i]+dm[i][0]
            self.y[i]=self.y[i]+dm[i][1]
            self.z[i]=self.z[i]+dm[i][2]

    def GetForces(self): #load all of the total forces from DFTB
        try:
            file = open("DFTB+/detailed.out", "r+")
        except OSError:
            print ("Could not open detailed.out file")
            sys.exit()

        lines = file.readlines()
        forceindex = -1
        for i in range(len(lines)):
            if lines[i].find("Total Forces") != -1:
                forceindex = i
        lines = lines[forceindex+1:forceindex+self.Nat+1]

        Forces = []
        for i in range(len(lines)):
            a = lines[i].split(' ')
            a = list(filter(lambda x: x != '', a))
            a = list(map(float, a[1:]))
            Forces.append(a)
        return np.array(Forces)

    def GetHessian(self,condensed=False): #load the hessian matrix from dftb
        try:
            file = open("DFTB+/hessian.out", "r+")
        except OSError:
            print ("Could not open detailed.out file")
            sys.exit()

        ceil = math.ceil(float(3*self.Nat)/float(4))
        lines = file.readlines()
        Hessian = []
        for i in range(int(len(lines)/ceil)):
            row = []
            for j in range(ceil):
                a = lines[i*ceil+j].split(' ')
                a = list(filter(lambda x: x != '', a))
                a = list(map(float, a))
                row = row + a
            Hessian.append(np.array(row))

        if condensed==True:
            HessianCond=[]
            for i in range(self.Nat*3):
                HessianCondColumn = []
                for j in range(self.Nat):
                    cond = np.sqrt(Hessian[i][j*3]**2+Hessian[i][j*3+1]**2+Hessian[i][j*3+2]**2)
                    HessianCondColumn.append(cond)
                HessianCond.append(np.array(HessianCondColumn))
            return np.array(HessianCond)
        return np.array(Hessian)

    def _condenseHessian(self,Hessian): #calculates the condensed Hessian for a given Hessian
        HessianCond=[]
        for i in range(self.Nat*3):
            HessianCondColumn = []
            for j in range(self.Nat):
                cond = np.sqrt(Hessian[i][j*3]**2+Hessian[i][j*3+1]**2+Hessian[i][j*3+2]**2)
                HessianCondColumn.append(cond)
            HessianCond.append(np.array(HessianCondColumn))
        return np.array(HessianCond)

    def GetEnergy(self): #load the total energy from DFTB
        try:
            file = open("DFTB+/detailed.out", "r+")
        except OSError:
            print ("Could not open detailed.out file")
            sys.exit()

        lines = file.readlines()
        forceindex = -1
        for i in range(len(lines)):
            if lines[i].find("Total energy:") != -1:
                forceindex = i
        lines = lines[forceindex]

        a = lines.split(' ')
        a = list(filter(lambda x: x != '', a))
        a = a[2]
        return float(a)
        

    def CalcBondedNeighbours(self,Cutoffneighbours): #calculates the list of bonds in the whole structure, bonds are defined using Cutoffneighbours

        self.bonds = []
        for i in range(self.Nat):
            bondidx = []
            distances = self.Distances(i)
            for k in range(self.Nat):
                if distances[k] <= Cutoffneighbours and k != i:
                    bondidx.append(k)
            self.bonds.append(bondidx)

    def CalcBondAngles(self): #calculates the list of bond angles in the whole structure, bonds are defined using Cutoffneighbours
        self.angles = []
        for i in range(self.Nat):
            p = permutations(self.bonds[i]) 
            p = list(p)
            if len(p) == 1:
                p = None
            else:  
                for k in range(len(p)):
                    p[k]=list(p[k][:2])
                for k in range(len(p)):
                    if k>=len(p):
                        break
                    for j in range(len(p)- 1, -1, -1):
                        if p[k][0] == p[j][0] and p[k][1] == p[j][1] and k!= j:
                            p.pop(j)
                    for j in range(len(p)- 1, -1, -1):
                        if p[k][0] == p[j][1] and p[k][1] == p[j][0] and k!= j:
                            p.pop(j)
            self.angles.append(p)

        for i in range(len(self.angles)): 
            if self.angles[i] != None:
                for k in range(len(self.angles[i])):
                    v1 = self.GetVersor(i,self.angles[i][k][0])
                    v2 = self.GetVersor(i,self.angles[i][k][1])
                    inner = np.inner(v1,v2)
                    rad = np.arccos(np.clip(inner, -1.0, 1.0))
                    self.angles[i][k].append(rad)
                    cross = np.cross(v1,v2)
                    if rad <= 3.14159265359 and rad > 3.14159265358:
                        inner = np.inner(v1,v2+np.array([0,0,0.001]))
                        rad = np.arccos(np.clip(inner, -1.0, 1.0))
                        cross = np.cross(v1,v2+np.array([0,0,0.001]))
                        if rad <= rad <= 3.14159265359 and rad > 3.14159265358:
                            inner = np.inner(v1,v2+np.array([0,0.001,0]))
                            rad = np.arccos(np.clip(inner, -1.0, 1.0))
                            cross = np.cross(v1,v2+np.array([0,0.001,0]))

                    cross = cross/LA.norm(cross)
                    self.angles[i][k].append(cross)


    def CalcBondOffPlane(self): #calculates the list of offplane angles in the whole structure, bonds are defined using Cutoffneighbours
        self.offplane = []
        for i in range(self.Nat):

            if len(self.bonds[i])>2:
                self.offplane.append(self.bonds[i][:3])
            else:
                self.offplane.append(None)

        for i in range(len(self.offplane)): 
            if self.offplane[i] != None:
                v1 = self.GetVersor(i,self.offplane[i][0])
                v2 = self.GetVersor(i,self.offplane[i][1])
                v3 = self.GetVersor(i,self.offplane[i][2])
                v = np.cross(v1,v2) + np.cross(v2,v3) + np.cross(v3,v1)
                h = np.inner(v/LA.norm(v), v1)
                self.offplane[i].append(np.arcsin(h))
                # v1 = self.GetVersor(self.offplane[i][0],self.offplane[i][1])
                # v2 = self.GetVersor(self.offplane[i][0],self.offplane[i][2])
                # cross = np.cross(v1,v2)
                # cross = cross/LA.norm(cross)
                # v = self.GetVector(i,self.offplane[i][1])
                # distance = LA.norm(v)
                # inner = np.inner(v,cross)
                # self.offplane[i].append(np.arcsin(inner/distance))


    def CalcBondDistances(self): #calculates the list of bond angles in the whole structure, bonds are defined using Cutoffneighbours
        self.distances = []
        for i in range(self.Nat):
            row = []
            for j in range(self.Nat):
                row.append(self.Distance(i,j))
            self.distances.append(row)


                

