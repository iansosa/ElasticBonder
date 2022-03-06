import numpy as np
import copy
import sys
import numpy.linalg as LA
from numpy import random
from scipy.optimize import curve_fit
import utils
from functools import partial

class Bonds():

    def __init__(self,structure,optimize=True):
        self.structure_eq = copy.deepcopy(structure)

        if optimize == True:
            self.structure_eq.SaveGeometry()
            self.structure_eq.RunOptimize()
            self.structure_eq.LoadGeometry()
        else:
            self.structure_eq.SaveGeometry()
            self.structure_eq.RunStatic()
        self.BondMatrix = None
        self.EqDistanceMatrix = self.structure_eq.distances
        self.EqPos = []
        x,y,z = self.structure_eq.Pos()
        for i in range(self.structure_eq.Nat):
            pos = [x[i],y[i],z[i]]
            self.EqPos.append(np.array(pos))
        self.EqAngles = self.structure_eq.angles

    def SaveForces(self,idx,variance=0.5,db=0.02):
        eqForces = self.structure_eq.GetForces()
        structure = copy.deepcopy(self.structure_eq)
        points = int(variance/db)

        for i in range(structure.Nat):
            if i != idx:
                with open('out/Forces_'+str(idx)+'_'+str(i)+'.txt', 'w') as f:
                    versor = structure.GetVersor(idx,i)
                    forward = []
                    forward.append([0,0])
                    for k in range(1,points):

                        structure.PullBond(idx,i,db)
                        structure.SaveGeometry()
                        structure.RunStatic()
                        newForces = structure.GetForces()
                        forward.append([k*db,np.inner(newForces[idx],versor)-np.inner(eqForces[idx],versor)])
                    structure = copy.deepcopy(self.structure_eq)
                    backward = []
                    for k in range(1,points):

                        structure.PullBond(idx,i,-db)
                        structure.SaveGeometry()
                        structure.RunStatic()
                        newForces = structure.GetForces()
                        backward.append([-k*db,np.inner(newForces[idx],versor)-np.inner(eqForces[idx],versor)])
                    backward.reverse()

                    total = backward + forward

                    for k in range(len(total)):
                        f.write(str(total[k][0])+' '+str(total[k][1])+'\n')
                    structure = copy.deepcopy(self.structure_eq)


    def GetK(self,i,j,db=0.02):
        if i == j:
            return 0
        self.structure_eq.SaveGeometry()
        self.structure_eq.RunStatic()
        eqForces = self.structure_eq.GetForces()
        structure = copy.deepcopy(self.structure_eq)
        versor = structure.GetVersor(i,j)

        structure.PullBond(i,j,db)
        structure.SaveGeometry()
        structure.RunStatic()
        newForces = structure.GetForces()
        ForceForward = np.inner(newForces[i],versor)-np.inner(eqForces[i],versor)

        return ForceForward/db

    def GetForces(self,dv_matrix):
        self.structure_eq.SaveGeometry()
        self.structure_eq.RunStatic()
        eqForces = self.structure_eq.GetForces()
        structure = copy.deepcopy(self.structure_eq)
        for i in range(structure.Nat):
            structure.x[i] = structure.x[i] + dv_matrix[i][0]
            structure.y[i] = structure.y[i] + dv_matrix[i][1]
            structure.z[i] = structure.z[i] + dv_matrix[i][2]
        structure.SaveGeometry()
        structure.RunStatic()
        newForces = structure.GetForces()
        return newForces - eqForces

    def GetKangle(self,i,j,k,rotvec,da=0.02): #gets coupling strenght from atom k towards j rotanting k an angle da pivoting on i
        if i == j or i == k or k == j:
            return 0
        self.structure_eq.SaveGeometry()
        self.structure_eq.RunStatic()
        eqForces = self.structure_eq.GetForces()
        structure = copy.deepcopy(self.structure_eq)

        v1_in = structure.GetVersor(i,k)
        v2_in = structure.GetVersor(i,j)
        inner_in = np.inner(v1_in,v2_in)

        #################################################################
        structure.RotateBond(i,j,k,rotvec,da)
        v2 = structure.GetVersor(i,k)
        v1 = structure.GetVersor(i,j)
        inner = np.inner(v1,v2)
        versor = (v2 - v1*inner)
        versor = versor/LA.norm(versor)

        structure.SaveGeometry()
        structure.RunStatic()
        newForces = structure.GetForces()


        K1 = - (np.inner(newForces[j],versor)-np.inner(eqForces[j],versor))*structure.Distance(i,j)/(inner-inner_in)


        ###################################################################
        structure = copy.deepcopy(self.structure_eq)
        ################################################################3333
        structure.RotateBond(i,k,j,rotvec,da)
        v2 = structure.GetVersor(i,k)
        v1 = structure.GetVersor(i,j)
        inner = np.inner(v1,v2)
        versor = (v1 - v2*inner)
        versor = versor/LA.norm(versor)

        structure.SaveGeometry()
        structure.RunStatic()
        newForces = structure.GetForces()

        K2 = - (np.inner(newForces[k],versor)-np.inner(eqForces[k],versor))*structure.Distance(i,k)/(inner-inner_in)

        return (K1+K2)/2
        ###################################################################
        # structure = copy.deepcopy(self.structure_eq)
        ################################################################3333
        # structure.RotateBond(i,k,j,rotvec,da/2.0)
        # structure.RotateBond(i,j,k,rotvec,da/2.0)
        
        # v1 = structure.GetVersor(i,k)
        # v2 = structure.GetVersor(i,j)
        # inner = np.inner(v1,v2)

        # structure.SaveGeometry()
        # structure.RunStatic()
        # newForces = structure.GetForces()

        # vector = (v1* inner - v2)/structure.Distance(i,k) + (v2* inner - v1)/structure.Distance(i,j)
        # magnitude = LA.norm(vector)
        # versor = vector/magnitude

        # pos = np.array([structure.x[i],structure.y[i],structure.z[i]])

        # K3 = -magnitude*(np.inner(newForces[i],versor)-np.inner(eqForces[i],versor))/(inner-inner_in)
        # return -magnitude*(np.inner(newForces[i],versor)-np.inner(eqForces[i],versor))/(inner-inner_in)

    def CalcSaveEnergyPerturbation(self,H0):

        dm = []
        for i in range(self.structure_eq.Nat):
            dm.append([(random.rand()-0.5)/5,(random.rand()-0.5)/5,(random.rand()-0.5)/5])

        structure = copy.deepcopy(self.structure_eq)
        structure.SetDisplacements(dm)
        structure.SaveGeometry()
        structure.RunStatic()
        structure.CalcBondDistances()
        structure.CalcBondAngles()
        structure.CalcBondOffPlane()

        bonds_eq = []
        alreadyconsidered = []
        for i in range(len(self.structure_eq.bonds)):
            for j in range(len(self.structure_eq.bonds[i])):
                if ([i,self.structure_eq.bonds[i][j]] in alreadyconsidered) == False:
                    alreadyconsidered.append([self.structure_eq.bonds[i][j],i])
                    bonds_eq.append(self.structure_eq.Distance(i,self.structure_eq.bonds[i][j]))
        angles_eq = []
        for i in range(len(self.structure_eq.angles)):
            if self.structure_eq.angles[i] != None: 
                for j in range(len(self.structure_eq.angles[i])):
                    angles_eq.append(self.structure_eq.angles[i][j][2])

        offplane_eq = []
        for i in range(len(self.structure_eq.offplane)):
            if self.structure_eq.offplane[i] != None: 
                    offplane_eq.append(self.structure_eq.offplane[i][3])

        bonds = []
        alreadyconsidered = []
        for i in range(len(structure.bonds)):
            for j in range(len(structure.bonds[i])):
                if ([i,structure.bonds[i][j]] in alreadyconsidered) == False:
                    alreadyconsidered.append([structure.bonds[i][j],i])
                    bonds.append(structure.Distance(i,structure.bonds[i][j]))
        angles = []
        for i in range(len(structure.angles)):
            if structure.angles[i] != None: 
                for j in range(len(structure.angles[i])):
                    angles.append(structure.angles[i][j][2])

        offplane = []
        for i in range(len(structure.offplane)):
            if structure.offplane[i] != None: 
                    offplane.append(structure.offplane[i][3])

        for i in range(len(bonds)):
            bonds[i] = bonds[i] - bonds_eq[i]

        for i in range(len(angles)):
            angles[i] = np.cos(angles[i]) - np.cos(angles_eq[i])

        for i in range(len(offplane)):
            offplane[i] = offplane[i] - offplane_eq[i]

        return structure.GetEnergy() - H0, bonds , angles , offplane

    def FitEnergy(self,iters=100,type_opt="two"):

        self.structure_eq.SaveGeometry()
        self.structure_eq.RunStatic()
        H0 = self.structure_eq.GetEnergy()

        H = []
        x_in = []
        distances = []
        angles = []
        offplane = []
        for i in range(iters):
            print(str(i)+"/"+str(iters))
            results = self.CalcSaveEnergyPerturbation(H0)
            H.append(results[0])
            distances.append(results[1])
            angles.append(results[2])
            offplane.append(results[3])
            x_in.append(results[1]+results[2]+results[3]) 
        Nbonds = len(distances[0])
        Nangs = len(angles[0])
        Noffplane = len(offplane[0])
        numparams = Nbonds + Nangs + Noffplane
        x_in = np.array(x_in).T
        if type_opt == "two":
            Energy = partial(functions.EbondsTwo,Nbonds,Nangs)
            pars, cov = curve_fit(f=Energy, xdata=x_in, ydata=H,p0=[0.6,0.3])
        if type_opt == "three":
            Energy = partial(functions.EbondsThree,Nbonds,Nangs,Noffplane)
            pars, cov = curve_fit(f=Energy, xdata=x_in, ydata=H,p0=[0.6,0.3,0.3])
        elif type_opt == "all":
            Energy = partial(functions.Ebonds,Nbonds,Nangs)
            pars, cov = curve_fit(f=Energy, xdata=x_in, ydata=H,p0=[0]*numparams)
        print(pars)
        print(cov)

    def CalcSaveBondMatrix(self,type_calc="restricted"):
        bonds = []
        if type_calc == "all":
            for i in range(self.structure_eq.Nat):
                column = []
                for k in range(self.structure_eq.Nat):
                    column.append(self.GetK(k,i))
                bonds.append(column)
        elif type_calc == "restricted":
            bonded = self.structure_eq.bonds
            for i in range(self.structure_eq.Nat):
                column = []
                for k in range(self.structure_eq.Nat):
                    if (k in bonded[i]):
                        column.append(self.GetK(k,i))
                    else:
                        column.append(0)
                bonds.append(column)
        else:
            print ("type_calc is either restricted or all")
            sys.exit()
        self.BondMatrix = bonds
        self.SaveBondMatrix()

    def CalcAngleBond(self):
        for i in range(len(self.EqAngles)):
            if self.EqAngles[i] != None:
                for k in range(len(self.EqAngles[i])):
                    self.EqAngles[i][k].append(self.GetKangle(i,self.EqAngles[i][k][0],self.EqAngles[i][k][1],self.EqAngles[i][k][3],da=0.1))


    def SaveDistanceMatrix(self):
        with open('out/EqDistances_.txt', 'w') as f:
            for i in range(self.structure_eq.Nat):
                for k in range(self.structure_eq.Nat):
                    f.write(str(self.EqDistanceMatrix[i][k])+' ')
                f.write('\n')

    def SaveBondMatrix(self):
        with open('out/Bonds.txt', 'w') as f:
            for i in range(len(self.BondMatrix)):
                for k in range(len(self.BondMatrix)):
                    f.write(str(self.BondMatrix[i][k])+' ')
                f.write('\n')

    def SaveBondsOverDistance(self,idx):
        bonds = []
        for k in range(self.structure_eq.Nat):
            bonds.append(self.GetK(idx,k))
        with open('out/BondsOverDistance_'+str(idx)+'.txt', 'w') as f:
            for k in range(len(bonds)):
                f.write(str(k)+' '+str(bonds[k])+'\n')

