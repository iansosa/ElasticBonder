import numpy as np
import copy
import sys

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
        eqForces = self.structure_eq.GetForces()
        structure = copy.deepcopy(self.structure_eq)
        versor = structure.GetVersor(i,j)

        structure.PullBond(i,j,db)
        structure.SaveGeometry()
        structure.RunStatic()
        newForces = structure.GetForces()
        ForceForward = np.inner(newForces[i],versor)-np.inner(eqForces[i],versor)

        structure.PullBond(i,j,-2*db)
        structure.SaveGeometry()
        structure.RunStatic()
        newForces = structure.GetForces()
        ForceBackwards= np.inner(newForces[i],versor)-np.inner(eqForces[i],versor)

        return (ForceForward-ForceBackwards)/(2*db)

    def CalcSaveBondMatrix(self):
        bonds = []
        for i in range(self.structure_eq.Nat):
            column = []
            for k in range(self.structure_eq.Nat):
                column.append(self.GetK(k,i))
            bonds.append(column)
        self.BondMatrix = bonds
        self.SaveBondMatrix()

    def SaveBondMatrix(self):
        with open('out/Bonds.txt', 'w') as f:
            for i in range(len(self.BondMatrix)):
                f.write(str(i))
                for k in range(len(self.BondMatrix)):
                    f.write(' '+str(self.BondMatrix[i][k]))
                f.write('\n')

    def SaveBondsOverDistance(self,idx):
        bonds = []
        for k in range(self.structure_eq.Nat):
            bonds.append(self.GetK(idx,k))
        with open('out/BondsOverDistance_'+str(idx)+'.txt', 'w') as f:
            for k in range(len(bonds)):
                f.write(str(k)+' '+str(bonds[k])+'\n')

