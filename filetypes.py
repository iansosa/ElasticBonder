import numpy as np
import sys

def Loadsdf(path,conversion):
    try:
        file = open(path, "r+")
    except OSError:
        print ("Could not open sdf file")
        sys.exit()
    lines = file.readlines()
    aux = lines[3].split(' ')
    aux = list(filter(lambda x: x != '', aux))
    print(int(aux[0]))
    Nat = int(aux[0])
    lines = lines[4:]

    geometry = []
    noC = 0
    for i in range(Nat):
        a = lines[i].split(' ')
        a = list(filter(lambda x: x != '', a))
        if a[3] == 'C':
            a = list(map(float, a[0:3]))
            geometry.append(a)
        else:
            noC = noC + 1

    arr_t = np.array(geometry).T/conversion
    geometry = arr_t.tolist()
    Nat = Nat - noC
    return Nat, geometry

def Loadgen(path,conversion):
    try:
        file = open(path, "r+")
    except OSError:
        print ("Could not open gen file")
        sys.exit()
    lines = file.readlines()

    aux = lines[0].split(' ')
    aux = list(filter(lambda x: x != '', aux))
    Nat = int(aux[0])
    lines = lines[2:]

    geometry = []
    for i in range(len(lines)):
        a = lines[i].split(' ')
        a = list(filter(lambda x: x != '', a))
        a = list(map(float, a[2:]))
        geometry.append(a)

    arr_t = np.array(geometry).T/conversion
    geometry = arr_t.tolist()

    return Nat, geometry

def Loadxyz_single(path,conversion):
    try:
        file = open(path, "r+")
    except OSError:
        print ("Could not open xyz file")
        sys.exit()
    lines = file.readlines()

    aux = lines[0].split(' ')
    aux = list(filter(lambda x: x != '', aux))
    Nat = int(aux[0])
    lines = lines[2:]

    geometry = []
    for i in range(len(lines)):
        a = lines[i].split(' ')
        a = list(filter(lambda x: x != '', a))
        a = list(map(float, a[1:-1]))
        geometry.append(a)

    arr_t = np.array(geometry).T/conversion
    geometry = arr_t.tolist()

    return Nat, geometry

def Loadxyz(path,conversion): #dftb+ specifies velocities in A/ps, velocities are returned in bohr/a.u(time)
    femtosecond = 41.341374575751
    try:
        file = open(path, "r+")
    except OSError:
        print ("Could not open xyz file")
        sys.exit()
    lines = file.readlines()

    aux = lines[0].split(' ')
    aux = list(filter(lambda x: x != '', aux))
    Nat = int(aux[0])
    Niter = int(len(lines)/(2+Nat))
    print(Niter)

    evolution = []
    for i in range(Niter):
        instant = []
        for j in range(Nat):
            a = lines[2+j+i*(2+Nat)].split(' ')
            a = list(filter(lambda x: x != '', a))
            a = list(map(float, a[1:8]))
            a = [a[0]/conversion,a[1]/conversion,a[2]/conversion,a[4]/(femtosecond*conversion*1000),a[5]/(femtosecond*conversion*1000),a[6]/(femtosecond*conversion*1000)]
            instant.append(a)
        evolution.append(instant)
        
    return Nat, Niter, evolution