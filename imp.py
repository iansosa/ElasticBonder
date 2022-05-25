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

    displ = []
    # displ.append("0.00")
    for i in range(0,100):
        u=i*du-4.95
        ru = str(round(u,2))
        if len(ru) == 4:
            ru = ru + '0'
        elif len(ru) == 3:
            ru = ru + '00'
        displ.append(ru)
    print(displ)

    if vdw == None:
        struct_name = "/buckling/novdw/geom_-4.95.gen"
    elif vdw == "MBD":
        struct_name = "/buckling/MBD/geom_MBD_-4.95.gen"
    elif vdw == "PW":
        struct_name = "/buckling/PW/geom_PW_-4.95.gen"

    chain = structures.Sphere(Nat,R0)
    chain.LoadGeometry(struct_name)
    Nat = chain.Nat
    chain.SaveGeometry()
    chain.RunOptimize(vdw=vdw,static=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,Nat-1,Nat-2,Nat-3,Nat-4,Nat-5,Nat-6,Nat-7,Nat-8,Nat-9,Nat-10,Nat-11,Nat-12,Nat-13,Nat-14,Nat-15,Nat-16],read_charges=False)
    chain.LoadGeometry()
    if vdw == None:
        chain.SaveGeometry("_"+displ[0],"buckling")
    else:
        chain.SaveGeometry("_"+vdw+"_"+displ[0],"buckling")


    

    for i in range(1,100):
        for j in range(16):
            chain.Displace(Nat-1-j,[0,0,du])
        chain.SaveGeometry()
        chain.RunOptimize(vdw=vdw,static=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,Nat-1,Nat-2,Nat-3,Nat-4,Nat-5,Nat-6,Nat-7,Nat-8,Nat-9,Nat-10,Nat-11,Nat-12,Nat-13,Nat-14,Nat-15,Nat-16],read_charges=True)
        chain.LoadGeometry()
        u = i*du
        if vdw == None:
            chain.SaveGeometry("_"+displ[i],"buckling")
        else:
            chain.SaveGeometry("_"+vdw+"_"+displ[i],"buckling")

def compressiontest(vdw=None):

    Nat = 5
    R0 = 2.4

    u = 0
    du = -5.0/100.0

    displ = []
    # displ.append("0.00")
    for i in range(0,100):
        u=i*du
        ru = str(round(u,2))
        if len(ru) == 4:
            ru = ru + '0'
        elif len(ru) == 3:
            ru = ru + '00'
        displ.append(ru)
    print(displ)

    if vdw == None:
        struct_name = "/Nanotubes/geom_5-10.gen"
    elif vdw == "MBD":
        struct_name = "/Nanotubes/geom_5-10_MBD.gen"
    elif vdw == "PW":
        struct_name = "/Nanotubes/geom_5-10_PW.gen"

    chain = structures.Sphere(Nat,R0)
    chain.LoadGeometry(struct_name)
    Nat = chain.Nat
    chain.SaveGeometry()
    chain.RunOptimize(vdw=vdw,static=[2,3,30,31,68,69,96,97,124,125,138,139,100,101,72,73,44,45,16,17],read_charges=False)
    chain.LoadGeometry()
    if vdw == None:
        chain.SaveGeometry("_"+displ[0],"nanotube")
    else:
        chain.SaveGeometry("_"+vdw+"_"+displ[0],"nanotube")


    for i in range(1,100):
        chain.Displace(2,[0,du,0])
        chain.Displace(3,[0,du,0])
        chain.Displace(30,[0,du,0])
        chain.Displace(31,[0,du,0])
        chain.Displace(68,[0,du,0])
        chain.Displace(69,[0,du,0])
        chain.Displace(96,[0,du,0])
        chain.Displace(97,[0,du,0])
        chain.Displace(124,[0,du,0])
        chain.Displace(125,[0,du,0])
        chain.Displace(138,[0,-du,0])
        chain.Displace(139,[0,-du,0])
        chain.Displace(100,[0,-du,0])
        chain.Displace(101,[0,-du,0])
        chain.Displace(72,[0,-du,0])
        chain.Displace(73,[0,-du,0])
        chain.Displace(44,[0,-du,0])
        chain.Displace(45,[0,-du,0])
        chain.Displace(16,[0,-du,0])
        chain.Displace(17,[0,-du,0])


        chain.SaveGeometry()
        chain.RunOptimize(vdw=vdw,static=[2,3,30,31,68,69,96,97,124,125,138,139,100,101,72,73,44,45,16,17],read_charges=True)
        chain.LoadGeometry()
        u = i*du
        if vdw == None:
            chain.SaveGeometry("_"+displ[i],"nanotube")
        else:
            chain.SaveGeometry("_"+vdw+"_"+displ[i],"nanotube")

def compressiontest_big(vdw=None):

    Nat = 5
    R0 = 2.4

    u = 0
    du = -5.0/100.0

    displ = []
    # displ.append("0.00")
    for i in range(0,100):
        u=i*du
        ru = str(round(u,2))
        if len(ru) == 4:
            ru = ru + '0'
        elif len(ru) == 3:
            ru = ru + '00'
        displ.append(ru)
    print(displ)

    if vdw == None:
        struct_name = "/Nanotubes/geom_25-50.gen"
    elif vdw == "MBD":
        struct_name = "/Nanotubes/geom_25-50_MBD.gen"
    elif vdw == "PW":
        struct_name = "/Nanotubes/geom_25-50_PW.gen"

    chain = structures.Sphere(Nat,R0)
    chain.LoadGeometry(struct_name)
    Nat = chain.Nat
    chain.SaveGeometry()
    chain.RunOptimize(vdw=vdw,static=[605,604,457,456,309,308,161,160,13,12,86,87,234,235,382,383,530,531,678,679],read_charges=False)
    chain.LoadGeometry()
    if vdw == None:
        chain.SaveGeometry("_"+displ[0],"nanotube")
    else:
        chain.SaveGeometry("_"+vdw+"_"+displ[0],"nanotube")


    for i in range(1,100):
        chain.Displace(605,[0,du,0])
        chain.Displace(604,[0,du,0])
        chain.Displace(457,[0,du,0])
        chain.Displace(456,[0,du,0])
        chain.Displace(309,[0,du,0])
        chain.Displace(308,[0,du,0])
        chain.Displace(161,[0,du,0])
        chain.Displace(160,[0,du,0])
        chain.Displace(13,[0,du,0])
        chain.Displace(12,[0,du,0])
        chain.Displace(86,[0,-du,0])
        chain.Displace(87,[0,-du,0])
        chain.Displace(234,[0,-du,0])
        chain.Displace(235,[0,-du,0])
        chain.Displace(382,[0,-du,0])
        chain.Displace(383,[0,-du,0])
        chain.Displace(530,[0,-du,0])
        chain.Displace(531,[0,-du,0])
        chain.Displace(678,[0,-du,0])
        chain.Displace(679,[0,-du,0])


        chain.SaveGeometry()
        chain.RunOptimize(vdw=vdw,static=[605,604,457,456,309,308,161,160,13,12,86,87,234,235,382,383,530,531,678,679],read_charges=True)
        chain.LoadGeometry()
        u = i*du
        if vdw == None:
            chain.SaveGeometry("_"+displ[i],"nanotube")
        else:
            chain.SaveGeometry("_"+vdw+"_"+displ[i],"nanotube")

def compressiontest_UHMWPE(vdw=None):

    Nat = 5
    R0 = 2.4

    u = 0
    du = 5.0/100.0

    displ = []
    # displ.append("0.00")
    for i in range(0,100):
        u=i*du
        ru = str(round(u,2))
        if len(ru) == 4:
            ru = ru + '0'
        elif len(ru) == 3:
            ru = ru + '00'
        displ.append(ru)
    print(displ)

    if vdw == None:
        struct_name = "/UHMWPE/UHMWPE.gen"
    elif vdw == "MBD":
        struct_name = "/UHMWPE/UHMWPE_MBD.gen"
    elif vdw == "PW":
        struct_name = "/UHMWPE/UHMWPE_PW.gen"

    chain = structures.Sphere(Nat,R0)
    chain.LoadGeometry(struct_name)
    Nat = chain.Nat
    chain.SaveGeometry()
    chain.RunOptimize(vdw=vdw,static=None,read_charges=False)
    chain.LoadGeometry("geo_end.gen")
    if vdw == None:
        chain.SaveGeometry("_"+displ[0],"UHMWPE")
    else:
        chain.SaveGeometry("_"+vdw+"_"+displ[0],"UHMWPE")


    for i in range(1,100):
        chain.Displace_UC([0,du,0])
        chain.SaveGeometry()
        chain.RunOptimize(vdw=vdw,static=None,read_charges=True)
        chain.LoadGeometry("geo_end.gen")
        u = i*du
        if vdw == None:
            chain.SaveGeometry("_"+displ[i],"UHMWPE")
        else:
            chain.SaveGeometry("_"+vdw+"_"+displ[i],"UHMWPE")

def buckling_stats():

    Nat = 400
    R0 = 2.4


    du = 0.05
    u=0
    displ = []
    # displ.append('0.00')
    for i in range(0,34):
        u=i*du+4.95
        ru = str(round(u,2))
        if len(ru) == 3:
            ru = ru + '0'
        displ.append(ru)
    print(displ)

    Fz_PW = np.zeros(len(displ))
    for i in range(len(displ)):
        name = "/buckling/PW/geom_PW_-"+displ[i]+".gen"
        print(name)
        chain = structures.Sphere(Nat,R0)
        chain.LoadGeometry(name)
        chain.SaveGeometry()

        chain.RunStatic("PW")
        F = chain.GetForces()
        for j in range(chain.Nat-16,chain.Nat):
            Fz_PW[i] = Fz_PW[i] + F[j][2]
        print(Fz_PW)

    with open('out/Buckling_PW.txt', 'w') as f:
        for i in range(len(Fz_PW)):
            f.write(displ[i]+' '+str(Fz_PW[i])+'\n')



    # Fz_MBD = np.zeros(len(displ))
    # for i in range(len(displ)):
    #     name = "/buckling/MBD/geom_MBD_-"+displ[i]+".gen"
    #     print(name)
    #     chain = structures.Sphere(Nat,R0)
    #     chain.LoadGeometry(name)
    #     chain.SaveGeometry()

    #     chain.RunStatic("MBD")
    #     F = chain.GetForces()
    #     for j in range(chain.Nat-16,chain.Nat):
    #         Fz_MBD[i] = Fz_MBD[i] + F[j][2]
    #     print(Fz_MBD)

    # with open('out/Buckling_MBD.txt', 'w') as f:
    #     for i in range(len(Fz_MBD)):
    #         f.write(displ[i]+' '+str(Fz_MBD[i])+'\n')



    # Fz = np.zeros(len(displ))
    # for i in range(len(displ)):
    #     name = "/buckling/novdw/geom_-"+displ[i]+".gen"
    #     print(name)
    #     chain = structures.Sphere(Nat,R0)
    #     chain.LoadGeometry(name)
    #     chain.SaveGeometry()

    #     chain.RunStatic()
    #     F = chain.GetForces()
    #     for j in range(chain.Nat-16,chain.Nat):
    #         Fz[i] = Fz[i] + F[j][2]
    #     print(Fz)

    # with open('out/Buckling.txt', 'w') as f:
    #     for i in range(len(Fz)):
    #         f.write(displ[i]+' '+str(Fz[i])+'\n')

def UHWPE():

    Nat = 40
    R0 = 2.4
    chain = structures.Chain(Nat,R0)
    chain.LoadGeometry("UHWPE.gen")
    chain.SaveGeometry()

    E = []
    u_start = 9/10
    u_end = 11/10
    du=(u_end-u_start)/100
    for i in range(100):
        for j in range(chain.Nat):
            chain.y[j]=chain.y[j]*u_start
        chain.unit_cell[2][1]=chain.unit_cell[2][1]*u_start
        
        chain.SaveGeometry()
        chain.RunStatic()
        E.append(chain.GetEnergy())
        u_start = u_start + du
        chain.LoadGeometry("UHWPE.gen")
        print(E)

    u_start = 9/10
    u_end = 11/10
    du=(u_end-u_start)/100
    with open('out/UHWPE_E_y.txt', 'w') as f:
        for i in range(len(E)):
            f.write(str(u_start)+' '+str(E[i])+'\n')
            u_start = u_start + du

    Nat = 40
    R0 = 2.4
    chain = structures.Chain(Nat,R0)
    chain.LoadGeometry("UHWPE_MBD.gen")
    chain.SaveGeometry()

    E = []
    u_start = 9/10
    u_end = 11/10
    du=(u_end-u_start)/100
    for i in range(100):
        for j in range(chain.Nat):
            chain.y[j]=chain.y[j]*u_start
        chain.unit_cell[2][1]=chain.unit_cell[2][1]*u_start
        
        chain.SaveGeometry()
        chain.RunStatic("MBD")
        E.append(chain.GetEnergy())
        u_start = u_start + du
        chain.LoadGeometry("UHWPE_MBD.gen")
        print(E)

    u_start = 9/10
    u_end = 11/10
    du=(u_end-u_start)/100
    with open('out/UHWPE_E_MBD_y.txt', 'w') as f:
        for i in range(len(E)):
            f.write(str(u_start)+' '+str(E[i])+'\n')
            u_start = u_start + du

    Nat = 40
    R0 = 2.4
    chain = structures.Chain(Nat,R0)
    chain.LoadGeometry("UHWPE_PW.gen")
    chain.SaveGeometry()

    E = []
    u_start = 9/10
    u_end = 11/10
    du=(u_end-u_start)/100
    for i in range(100):
        for j in range(chain.Nat):
            chain.y[j]=chain.y[j]*u_start
        chain.unit_cell[2][1]=chain.unit_cell[2][1]*u_start
        
        chain.SaveGeometry()
        chain.RunStatic("PW")
        E.append(chain.GetEnergy())
        u_start = u_start + du
        chain.LoadGeometry("UHWPE_PW.gen")
        print(E)

    u_start = 9/10
    u_end = 11/10
    du=(u_end-u_start)/100
    with open('out/UHWPE_E_PW_y.txt', 'w') as f:
        for i in range(len(E)):
            f.write(str(u_start)+' '+str(E[i])+'\n')
            u_start = u_start + du

def UHMWPE_comp_stats():

    Nat = 40
    R0 = 2.4
    chain = structures.Sphere(Nat,R0)

    u = 0
    du = 5.0/100.0

    displ = []
    # displ.append("0.00")
    for i in range(0,100):
        u=i*du
        ru = str(round(u,2))
        if len(ru) == 4:
            ru = ru + '0'
        elif len(ru) == 3:
            ru = ru + '00'
        displ.append(ru)
    print(displ)

    E = []

    for i in range(100):
        
        chain.LoadGeometry("/UHMWPE/compression/UHMWPE-nowdv/geom_"+displ[i]+".gen")
        chain.SaveGeometry()
        chain.RunStatic()
        E.append(chain.GetEnergy())
        print(E)

    with open('out/UHWPE_E_y.txt', 'w') as f:
        for i in range(len(E)):
            f.write(displ[i]+' '+str(E[i])+'\n')

    E = []

    for i in range(100):
        
        chain.LoadGeometry("/UHMWPE/compression/UHMWPE-PW/geom_PW_"+displ[i]+".gen")
        chain.SaveGeometry()
        chain.RunStatic("PW")
        E.append(chain.GetEnergy())
        print(E)

    with open('out/UHWPE_E_y-PW.txt', 'w') as f:
        for i in range(len(E)):
            f.write(displ[i]+' '+str(E[i])+'\n')

    E = []

    for i in range(100):
        
        chain.LoadGeometry("/UHMWPE/compression/UHMWPE-MBD/geom_MBD_"+displ[i]+".gen")
        chain.SaveGeometry()
        chain.RunStatic("MBD")
        E.append(chain.GetEnergy())
        print(E)

    with open('out/UHWPE_E_y-MBD.txt', 'w') as f:
        for i in range(len(E)):
            f.write(displ[i]+' '+str(E[i])+'\n')

def nanotube_comp_stats():

    Nat = 40
    R0 = 2.4
    chain = structures.Sphere(Nat,R0)

    u = 0
    du = -5.0/100.0

    displ = []
    # displ.append("0.00")
    for i in range(0,100):
        u=i*du
        ru = str(round(u,2))
        if len(ru) == 4:
            ru = ru + '0'
        elif len(ru) == 3:
            ru = ru + '00'
        displ.append(ru)
    print(displ)

    indexes = [2,3,30,31,68,69,96,97,124,125]

    Fz_PW = np.zeros(len(displ))
    for i in range(len(displ)):
        chain.LoadGeometry("/Nanotubes/compression_5-10/geom_PW_"+displ[i]+".gen")
        chain.SaveGeometry()
        chain.RunStatic("PW")
        F = chain.GetForces()
        for j in indexes:
            Fz_PW[i] = Fz_PW[i] + F[j][1]
        print(Fz_PW)

    with open('out/Nanotuve_comp_PW.txt', 'w') as f:
        for i in range(len(Fz_PW)):
            f.write(displ[i]+' '+str(Fz_PW[i])+'\n')

    Fz_MBD = np.zeros(len(displ))
    for i in range(len(displ)):
        chain.LoadGeometry("/Nanotubes/compression_5-10/geom_MBD_"+displ[i]+".gen")
        chain.SaveGeometry()
        chain.RunStatic("MBD")
        F = chain.GetForces()
        for j in indexes:
            Fz_MBD[i] = Fz_MBD[i] + F[j][1]
        print(Fz_MBD)

    with open('out/Nanotuve_comp_MBD.txt', 'w') as f:
        for i in range(len(Fz_MBD)):
            f.write(displ[i]+' '+str(Fz_MBD[i])+'\n')

    Fz = np.zeros(len(displ))
    for i in range(len(displ)):
        chain.LoadGeometry("/Nanotubes/compression_5-10/geom_"+displ[i]+".gen")
        chain.SaveGeometry()
        chain.RunStatic()
        F = chain.GetForces()
        for j in indexes:
            Fz[i] = Fz[i] + F[j][1]
        print(Fz)

    with open('out/Nanotuve_comp_MBD.txt', 'w') as f:
        for i in range(len(Fz)):
            f.write(displ[i]+' '+str(Fz[i])+'\n')
