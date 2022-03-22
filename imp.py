import numpy as np
import numpy.linalg as LA
from mdhandler import Handler as MDH
import structures
import sys

def StaticOverEvolve_Chains(temp,Length,distance):
    if distance <0:
        print("distance should be positive")
        sys.exit()
    Nat=Length
    R0=2.4
    Chain1 = structures.Chain(Nat,R0)
    Chain1.SaveGeometry()
    Chain1.RunOptimize()
    Chain1.LoadGeometry()

    Chain2 = structures.Chain(Nat,R0)
    Chain2.SaveGeometry()
    Chain2.RunOptimize()
    Chain2.LoadGeometry()

    Chain2.MoveAll([0,0,distance])

    Chain1.add(Chain2)

    md = MDH(Chain1,False)
    md.RunMD(5000,temp,[0,Length-1,Length,2*Length-1])
    ForcesMBD = md.GetForcesSOE("MBD")
    ForcesPW = md.GetForcesSOE("PW")
    ForcesShort = md.GetForcesSOE()


    with open('out/SOE_'+str(temp)+'.txt', 'w') as f:
        for i in range(len(ForcesShort)):
            F_lower_MBD = 0
            F_lower_PW = 0
            F_lower_Short = 0
            for k in range(Length):
                F_lower_MBD = F_lower_MBD + ForcesMBD[i][k][2]
                F_lower_PW = F_lower_PW + ForcesPW[i][k][2]
                F_lower_Short = F_lower_Short + ForcesShort[i][k][2]
            F_upper_MBD = 0
            F_upper_PW = 0
            F_upper_Short = 0
            for k in range(Length,2*Length):
                F_upper_MBD = F_upper_MBD + ForcesMBD[i][k][2]
                F_upper_PW = F_upper_PW + ForcesPW[i][k][2]
                F_upper_Short = F_upper_Short + ForcesShort[i][k][2]
            f.write(str(i)+ " " +str(F_upper_MBD-F_upper_Short)+ " "+ str(F_lower_MBD-F_lower_Short)+ " " +str(F_upper_PW-F_upper_Short)+ " "+ str(F_lower_PW-F_lower_Short)+ " " +str(F_upper_Short)+ " "+ str(F_lower_Short)+"\n")

def CorrelationOverEvolve_Chains(temp,Length,distance):
    if distance <0:
        print("distance should be positive")
        sys.exit()
    Nat=Length
    R0=2.4
    Chain1 = structures.Chain(Nat,R0)
    Chain1.SaveGeometry()
    Chain1.RunOptimize()
    Chain1.LoadGeometry()

    Chain2 = structures.Chain(Nat,R0)
    Chain2.SaveGeometry()
    Chain2.RunOptimize()
    Chain2.LoadGeometry()

    Chain2.MoveAll([0,0,distance])
    Chain1.add(Chain2)

    Navg=100
    corrxx = np.zeros(Length)
    corryy = np.zeros(Length)
    corrzz = np.zeros(Length)

    corrxy = np.zeros(Length)
    corrxz = np.zeros(Length)
    corryz = np.zeros(Length)

    corrxx_abs = np.zeros(Length)
    corryy_abs = np.zeros(Length)
    corrzz_abs = np.zeros(Length)

    corrxy_abs = np.zeros(Length)
    corrxz_abs = np.zeros(Length)
    corryz_abs = np.zeros(Length)

    for k in range(Navg):
        print(str(k)+"/"+str(Navg))
        md = MDH(Chain1,False)
        md.RunMD(5000,temp,[0,Length-1,Length,2*Length-1])
        md.LoadEvolution()
        evolution = md.evolution
        for j in range(Length):
            v_lower_x = []
            v_lower_y = []
            v_lower_z = []
            v_lower_x_abs = []
            v_lower_y_abs = []
            v_lower_z_abs = []
            for i in range(len(evolution)):
                v_lower_x.append(evolution[i][j][0])
                v_lower_y.append(evolution[i][j][1])
                v_lower_z.append(evolution[i][j][2])
                v_lower_x_abs.append(np.abs(evolution[i][j][0]))
                v_lower_y_abs.append(np.abs(evolution[i][j][1]))
                v_lower_z_abs.append(np.abs(evolution[i][j][2]))
            v_lower_x=np.array(v_lower_x)
            v_lower_y=np.array(v_lower_y)
            v_lower_z=np.array(v_lower_z)
            v_lower_x_abs=np.array(v_lower_x_abs)
            v_lower_y_abs=np.array(v_lower_y_abs)
            v_lower_z_abs=np.array(v_lower_z_abs)

            v_upper_x = []
            v_upper_y = []
            v_upper_z = []
            v_upper_x_abs = []
            v_upper_y_abs = []
            v_upper_z_abs = []
            for i in range(len(evolution)):
                v_upper_x.append(evolution[i][j+Length][0])
                v_upper_y.append(evolution[i][j+Length][1])
                v_upper_z.append(evolution[i][j+Length][2])
                v_upper_x_abs.append(np.abs(evolution[i][j+Length][0]))
                v_upper_y_abs.append(np.abs(evolution[i][j+Length][1]))
                v_upper_z_abs.append(np.abs(evolution[i][j+Length][2]))
            v_upper_x=np.array(v_upper_x)
            v_upper_y=np.array(v_upper_y)
            v_upper_z=np.array(v_upper_z)
            v_upper_x_abs=np.array(v_upper_x_abs)
            v_upper_y_abs=np.array(v_upper_y_abs)
            v_upper_z_abs=np.array(v_upper_z_abs)

            v_lower = []
            v_lower.append(v_lower_x)
            v_lower.append(v_lower_y)
            v_lower.append(v_lower_z)
            v_lower.append(v_lower_x_abs)
            v_lower.append(v_lower_y_abs)
            v_lower.append(v_lower_z_abs)
            v_upper = []
            v_upper.append(v_upper_x)
            v_upper.append(v_upper_y)
            v_upper.append(v_upper_z)
            v_upper.append(v_upper_x_abs)
            v_upper.append(v_upper_y_abs)
            v_upper.append(v_upper_z_abs)

            corrxx[j]=corrxx[j]+np.corrcoef(v_lower,v_upper)[0,6]/Navg
            corryy[j]=corryy[j]+np.corrcoef(v_lower,v_upper)[1,7]/Navg
            corrzz[j]=corrzz[j]+np.corrcoef(v_lower,v_upper)[2,8]/Navg

            corrxy[j]=corrxy[j]+np.corrcoef(v_lower,v_upper)[0,7]/Navg
            corrxz[j]=corrxz[j]+np.corrcoef(v_lower,v_upper)[0,8]/Navg
            corryz[j]=corryz[j]+np.corrcoef(v_lower,v_upper)[1,8]/Navg

            corrxx_abs[j]=corrxx_abs[j]+np.corrcoef(v_lower,v_upper)[3,9]/Navg
            corryy_abs[j]=corryy_abs[j]+np.corrcoef(v_lower,v_upper)[4,10]/Navg
            corrzz_abs[j]=corrzz_abs[j]+np.corrcoef(v_lower,v_upper)[5,11]/Navg

            corrxy_abs[j]=corrxy_abs[j]+np.corrcoef(v_lower,v_upper)[3,10]/Navg
            corrxz_abs[j]=corrxz_abs[j]+np.corrcoef(v_lower,v_upper)[3,11]/Navg
            corryz_abs[j]=corryz_abs[j]+np.corrcoef(v_lower,v_upper)[4,11]/Navg

    with open('out/corr_'+str(temp)+'.txt', 'w') as f:
        for i in range(len(corrxx)):
            f.write(str(i)+ " " +str(corrxx[i])+" " +str(corryy[i])+" " +str(corrzz[i])+" " +str(corrxy[i])+" " +str(corrxz[i])+" " +str(corryz[i])+" " +str(corrxx_abs[i])+" " +str(corryy_abs[i])+" " +str(corrzz_abs[i])+" " +str(corrxy_abs[i])+" " +str(corrxz_abs[i])+" " +str(corryz_abs[i])+"\n")

def CorrelationOverEvolve_Chain(temp,Length):
    Nat=Length
    R0=2.4
    Chain1 = structures.Chain(Nat,R0)
    Chain1.SaveGeometry()
    Chain1.RunOptimize()
    Chain1.LoadGeometry()

    Navg=100
    corrxx = np.zeros(Length)
    corryy = np.zeros(Length)
    corrzz = np.zeros(Length)

    corrxy = np.zeros(Length)
    corrxz = np.zeros(Length)
    corryz = np.zeros(Length)

    corrxx_abs = np.zeros(Length)
    corryy_abs = np.zeros(Length)
    corrzz_abs = np.zeros(Length)

    corrxy_abs = np.zeros(Length)
    corrxz_abs = np.zeros(Length)
    corryz_abs = np.zeros(Length)

    for k in range(Navg):
        print(str(k)+"/"+str(Navg))
        md = MDH(Chain1,False)
        md.RunMD(5000,temp,[0,Length-1])
        md.LoadEvolution()
        evolution = md.evolution
        for j in range(Length-1):
            v_lower_x = []
            v_lower_y = []
            v_lower_z = []
            v_lower_x_abs = []
            v_lower_y_abs = []
            v_lower_z_abs = []
            for i in range(len(evolution)):
                v_lower_x.append(evolution[i][j][0])
                v_lower_y.append(evolution[i][j][1])
                v_lower_z.append(evolution[i][j][2])
                v_lower_x_abs.append(np.abs(evolution[i][j][0]))
                v_lower_y_abs.append(np.abs(evolution[i][j][1]))
                v_lower_z_abs.append(np.abs(evolution[i][j][2]))
            v_lower_x=np.array(v_lower_x)
            v_lower_y=np.array(v_lower_y)
            v_lower_z=np.array(v_lower_z)
            v_lower_x_abs=np.array(v_lower_x_abs)
            v_lower_y_abs=np.array(v_lower_y_abs)
            v_lower_z_abs=np.array(v_lower_z_abs)

            v_upper_x = []
            v_upper_y = []
            v_upper_z = []
            v_upper_x_abs = []
            v_upper_y_abs = []
            v_upper_z_abs = []
            for i in range(len(evolution)):
                v_upper_x.append(evolution[i][j+1][0])
                v_upper_y.append(evolution[i][j+1][1])
                v_upper_z.append(evolution[i][j+1][2])
                v_upper_x_abs.append(np.abs(evolution[i][j+1][0]))
                v_upper_y_abs.append(np.abs(evolution[i][j+1][1]))
                v_upper_z_abs.append(np.abs(evolution[i][j+1][2]))
            v_upper_x=np.array(v_upper_x)
            v_upper_y=np.array(v_upper_y)
            v_upper_z=np.array(v_upper_z)
            v_upper_x_abs=np.array(v_upper_x_abs)
            v_upper_y_abs=np.array(v_upper_y_abs)
            v_upper_z_abs=np.array(v_upper_z_abs)

            v_lower = []
            v_lower.append(v_lower_x)
            v_lower.append(v_lower_y)
            v_lower.append(v_lower_z)
            v_lower.append(v_lower_x_abs)
            v_lower.append(v_lower_y_abs)
            v_lower.append(v_lower_z_abs)
            v_upper = []
            v_upper.append(v_upper_x)
            v_upper.append(v_upper_y)
            v_upper.append(v_upper_z)
            v_upper.append(v_upper_x_abs)
            v_upper.append(v_upper_y_abs)
            v_upper.append(v_upper_z_abs)

            corrxx[j]=corrxx[j]+np.corrcoef(v_lower,v_upper)[0,6]/Navg
            corryy[j]=corryy[j]+np.corrcoef(v_lower,v_upper)[1,7]/Navg
            corrzz[j]=corrzz[j]+np.corrcoef(v_lower,v_upper)[2,8]/Navg

            corrxy[j]=corrxy[j]+np.corrcoef(v_lower,v_upper)[0,7]/Navg
            corrxz[j]=corrxz[j]+np.corrcoef(v_lower,v_upper)[0,8]/Navg
            corryz[j]=corryz[j]+np.corrcoef(v_lower,v_upper)[1,8]/Navg

            corrxx_abs[j]=corrxx_abs[j]+np.corrcoef(v_lower,v_upper)[3,9]/Navg
            corryy_abs[j]=corryy_abs[j]+np.corrcoef(v_lower,v_upper)[4,10]/Navg
            corrzz_abs[j]=corrzz_abs[j]+np.corrcoef(v_lower,v_upper)[5,11]/Navg

            corrxy_abs[j]=corrxy_abs[j]+np.corrcoef(v_lower,v_upper)[3,10]/Navg
            corrxz_abs[j]=corrxz_abs[j]+np.corrcoef(v_lower,v_upper)[3,11]/Navg
            corryz_abs[j]=corryz_abs[j]+np.corrcoef(v_lower,v_upper)[4,11]/Navg

    with open('out/corr_'+str(temp)+'.txt', 'w') as f:
        for i in range(len(corrxx)):
            f.write(str(i)+ " " +str(corrxx[i])+" " +str(corryy[i])+" " +str(corrzz[i])+" " +str(corrxy[i])+" " +str(corrxz[i])+" " +str(corryz[i])+" " +str(corrxx_abs[i])+" " +str(corryy_abs[i])+" " +str(corrzz_abs[i])+" " +str(corrxy_abs[i])+" " +str(corrxz_abs[i])+" " +str(corryz_abs[i])+"\n")

def fixnanotube(name):

    Nat=50
    R0=2.4

    Chain = structures.Sphere(Nat,R0)
    Chain.LoadGeometry(name)
    pos = Chain.PosAsList()

    repeat = []
    for x in range(6*16):
        repeat.append(pos[16*10+x])
    groups = []
    for i in range(int(len(repeat)/16)):
        z = 0
        for j in range(16):
            z = z + repeat[i*16+j][2]
        groups.append(z/16)
    shift = groups[2]-groups[1]

    custom1 = structures.Custom(repeat)
    for x in range(5):
        custom2 = structures.Custom(repeat)
        groups = []
        pos = custom1.PosAsList()
        for i in range(int(len(pos)/16)):
            z = 0
            for j in range(16):
                z = z + pos[i*16+j][2]
            groups.append(z/16)
        custom2.MoveAll([0,0,groups[len(groups)-1]-groups[0]+shift])
        custom1.add(custom2)
    repeat = repeat[0:64]
    custom2 = structures.Custom(repeat)
    groups = []
    pos = custom1.PosAsList()
    for i in range(int(len(pos)/16)):
        z = 0
        for j in range(16):
            z = z + pos[i*16+j][2]
        groups.append(z/16)
    custom2.MoveAll([0,0,groups[len(groups)-1]-groups[0]+shift])
    custom1.add(custom2)

    pos = custom1.PosAsList()
    groups = []
    for i in range(int(len(pos)/16)):
        z = 0
        for j in range(16):
            z = z + pos[i*16+j][2]
        groups.append(z/16)
    for i in range(len(groups)-1):
        print(groups[i+1]-groups[i])

    custom1.SaveGeometry("nanotube_fixed")
    custom1.ShowStruct()

def bucklingtest(vdw=None):

    Nat = 5
    R0 = 2.4

    u = 0
    du = -5.0/100.0

    if vdw == None:
        struct_name = "NanotubeC640.gen"
    elif vdw == "MBD":
        struct_name = "NanotubeC640_MBD.gen"
    elif vdw == "PW":
        struct_name = "NanotubeC640_PW.gen"

    chain = structures.Sphere(Nat,R0)
    chain.LoadGeometry(struct_name)
    chain.SaveGeometry()
    chain.RunOptimize(vdw=vdw,static=None,read_charges=False)
    chain.LoadGeometry()
    if vdw == None:
        chain.SaveGeometry("_"+str(u),"buckling")
    else:
        chain.SaveGeometry("_"+vdw+"_"+str(u),"buckling")


    Nat = chain.Nat

    for i in range(1,100):
        for j in range(16):
            chain.Displace(Nat-1-j,[0,0,du])
        chain.SaveGeometry()
        chain.RunOptimize(vdw=vdw,static=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,Nat-1,Nat-2,Nat-3,Nat-4,Nat-5,Nat-6,Nat-7,Nat-8,Nat-9,Nat-10,Nat-11,Nat-12,Nat-13,Nat-14,Nat-15,Nat-16],read_charges=True)
        chain.LoadGeometry()
        u = i*du
        if vdw == None:
            chain.SaveGeometry("_"+str(u),"buckling")
        else:
            chain.SaveGeometry("_"+vdw+"_"+str(u),"buckling")