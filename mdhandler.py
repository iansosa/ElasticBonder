import copy
import numpy as np
import shutil
import subprocess
import filetypes
import sys

class Handler():

    def __init__(self,structure,optimize=True):
        self.structure_eq = copy.deepcopy(structure)

        if optimize == True:
            self.structure_eq.SaveGeometry()
            self.structure_eq.RunOptimize()
            self.structure_eq.LoadGeometry()
        else:
            self.structure_eq.SaveGeometry()
        self.evolution = None
        self.acelerations = None #acceleration in a.u.

    def RunMD(self,steps,temp=400,vdw=None,static=None,save_steps=1):
        if vdw == None:
            shutil.copyfile('DFTB+/md.hsd', 'DFTB+/dftb_in.hsd')
        elif vdw == "MBD":
            shutil.copyfile('DFTB+/md_mbd.hsd', 'DFTB+/dftb_in.hsd')
        elif vdw == "PW":
            shutil.copyfile('DFTB+/md_pw.hsd', 'DFTB+/dftb_in.hsd')
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
        idx = -1
        for i in range(len(lines)):
            if lines[i].find("Steps") != -1:
                idx = i
        targetline = "  Steps = " + str(steps) +"\n"
        lines[idx] = targetline

        idx = -1
        for i in range(len(lines)):
            if lines[i].find("MDRestartFrequency") != -1:
                idx = i
        targetline = "  MDRestartFrequency = " + str(save_steps) +"\n"
        lines[idx] = targetline

        idx = -1
        for i in range(len(lines)):
            if lines[i].find("Temperature [Kelvin]") != -1:
                idx = i
        targetline = "    Temperature [Kelvin] = " + str(temp) +"\n"
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

    def LoadEvolution(self):
        angstrom = 0.529177249
        femtosecond = 41.341374575751
        Nat, Niter, self.evolution = filetypes.Loadxyz("DFTB+/geo_end.xyz",angstrom)
        self.acelerations = []
        for i in range(len(self.evolution)-1):
            self.forcespmassiter = []
            for j in range(len(self.evolution[0])):
                self.forcespmassiter.append([(self.evolution[i+1][j][3]-self.evolution[i][j][3])/femtosecond,(self.evolution[i+1][j][4]-self.evolution[i][j][4])/femtosecond,(self.evolution[i+1][j][5]-self.evolution[i][j][5])/femtosecond])
            self.acelerations.append(self.forcespmassiter)

    def ComputeTempDispersions(self):
        if self.evolution == None:
            print ("Evolution not loaded")
            sys.exit()
        if len(self.evolution[0]) != self.structure_eq.Nat:
            print ("Evolution and equilibrium structure are not the same")
            sys.exit()

        averages = []
        for i in range(len(self.evolution[0])):
            averagex = 0
            averagey = 0
            averagez = 0
            for j in range(len(self.evolution)):
                averagex = averagex + self.evolution[j][i][0]
                averagey = averagey + self.evolution[j][i][1]
                averagez = averagez + self.evolution[j][i][2]
            averages.append([averagex/len(self.evolution),averagey/len(self.evolution),averagez/len(self.evolution)])
        dispersions = []
        for i in range(len(self.evolution[0])):
            dispersionx = 0
            dispersiony = 0
            dispersionz = 0
            for j in range(len(self.evolution)):
                dispersionx = dispersionx + (self.evolution[j][i][0]-averages[i][0])*(self.evolution[j][i][0]-averages[i][0])
                dispersiony = dispersiony + (self.evolution[j][i][1]-averages[i][1])*(self.evolution[j][i][1]-averages[i][1])
                dispersionz = dispersionz + (self.evolution[j][i][2]-averages[i][2])*(self.evolution[j][i][2]-averages[i][2])
            dispersions.append(np.sqrt((dispersionx+dispersiony+dispersionz))/len(self.evolution))

        averager = []
        for i in range(len(self.evolution[0])):
            r = 0
            for j in range(len(self.evolution)):
                diffx = (self.evolution[j][i][0]-averages[i][0])
                diffy = (self.evolution[j][i][1]-averages[i][1])
                diffz = (self.evolution[j][i][2]-averages[i][2])
                r = r + np.sqrt(diffx*diffx+diffy*diffy+diffz*diffz)
            averager.append(r/len(self.evolution))


        return averages, dispersions, averager

    def ComputeBondDispersions(self):
        if self.evolution == None:
            print ("Evolution not loaded")
            sys.exit()
        if len(self.evolution[0]) != self.structure_eq.Nat:
            print ("Evolution and equilibrium structure are not the same")
            sys.exit()

        bonds = self.structure_eq.bonds

        averages = []
        for i in range(len(self.evolution[0])):
            averages_in = []
            for j in range(len(bonds[i])):
                bx = 0
                by = 0
                bz = 0
                br = 0
                for k in range(len(self.evolution)):
                    bx = self.evolution[k][i][0]-self.evolution[k][bonds[i][j]][0]
                    by = self.evolution[k][i][1]-self.evolution[k][bonds[i][j]][1]
                    bz = self.evolution[k][i][2]-self.evolution[k][bonds[i][j]][2]
                    br = br + np.sqrt(bx*bx+by*by+bz*bz)
                averages_in.append(br/len(self.evolution))
            averages.append(averages_in)

        dispersions = []
        for i in range(len(self.evolution[0])):
            for j in range(len(bonds[i])):
                bx = 0
                by = 0
                bz = 0
                br = 0
                dispersion = 0
                for k in range(len(self.evolution)):
                    bx = self.evolution[k][i][0]-self.evolution[k][bonds[i][j]][0]
                    by = self.evolution[k][i][1]-self.evolution[k][bonds[i][j]][1]
                    bz = self.evolution[k][i][2]-self.evolution[k][bonds[i][j]][2]
                    br = np.sqrt(bx*bx+by*by+bz*bz)
                    dispersion = dispersion + (br-averages[i][j])*(br-averages[i][j])
                dispersions.append(np.sqrt(dispersion/len(self.evolution)))

        return averages, dispersions

    def GetForcesSOE(self,vdw=None):
        structure = copy.deepcopy(self.structure_eq)
        Forces=[]
        self.LoadEvolution()

        for i in range(len(self.evolution)):
            print(str(i)+"/"+str(len(self.evolution))+"\n")
            for j in range(structure.Nat):
                structure.x[j]=self.evolution[i][j][0]
                structure.y[j]=self.evolution[i][j][1]
                structure.z[j]=self.evolution[i][j][2]
            structure.SaveGeometry()
            structure.RunStatic(vdw)
            Forces.append(structure.GetForces())
        return Forces

