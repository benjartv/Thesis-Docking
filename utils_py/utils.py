from objects_py.vector import *
import numpy as np
import math
import commands
import string
import sys
import getopt
import os
import shutil

def calcRotM(vector, theta):
        vector=vector.copy()
        vector.normalize()
        c=np.cos(theta)
        s=np.sin(theta)
        t=1-c
        x, y, z=vector.get_array()
        rot=np.zeros((3, 3))
        # 1st row
        rot[0, 0]=t*x*x+c
        rot[0, 1]=t*x*y-s*z
        rot[0, 2]=t*x*z+s*y
        # 2nd row
        rot[1, 0]=t*x*y+s*z
        rot[1, 1]=t*y*y+c
        rot[1, 2]=t*y*z-s*x
        # 3rd row
        rot[2, 0]=t*x*z-s*y
        rot[2, 1]=t*y*z+s*x
        rot[2, 2]=t*z*z+c

        return rot

def getRMSD(ligand1, ligand2):
    dist = 0
    for i in range(len(ligand1.x)):
        dist += math.sqrt( (float(ligand1.x[i])-float(ligand2.x[i]))**2 + (float(ligand1.y[i])-float(ligand2.y[i]))**2 + (float(ligand1.z[i])-float(ligand2.z[i]))**2 )
    rmsd = dist / len(ligand1.x)
    return rmsd

def spherePoint(radius, theta, phi):
    xAux = radius * math.cos(theta) * math.sin(phi)
    yAux = radius * math.sin(theta) * math.sin(phi)
    zAux = radius * math.cos(phi)
    return [xAux, yAux, zAux]

def crossProduct(A, B):
    linesA = len(A)
    columnsA = len(A[0])

    linesB = len(B)
    columnsB = len(B[0])

    if columnsA == linesB:
        M = [[sum(A[m][n] * B[n][p] for n in range(columnsA)) \
              for p in range(columnsB)] for m in range(linesA)]
        return M
    else:
        return -1

def calculateFreeEnergy():
    a = commands.getstatusoutput('./energy/vina --config energy/config.txt --score_only')
    energiaVina = string.split(string.split(a[1],"Affinity:")[1])[0]
    #print energiaVina
    return float(energiaVina)


def convertTime(time):
    strTime = str(time).split()
    strTime = strTime[0]+"_"+strTime[1]
    strTime = strTime.split(".")[0]
    strTime = strTime.split(":")
    strTime = strTime[0].split("-") + strTime[1:]
    newTime = ""
    for i in strTime:
        newTime += i
    newTime = newTime.strip("-")
    return newTime

def initConfig():
    ligand = ""
    protein = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'l:p:h', ['ligand=', 'protein=', 'help'])
        if not opts:
          print ' Missing parameters'
          usage()
          sys.exit(2)
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-l', '--ligand'):
            ligand = arg
        elif opt in ('-p', '--protein'):
            protein = arg
        else:
            usage()
            sys.exit(2)
    return [ligand, protein]

def setPathLig(name):
    path = "molecules/ligand/"+name.lower()+"_ligand.pdbqt"
    return path

def setPathMLig(name):
    path = "molecules/m-ligand/"+name.lower()+"-modify.pdbqt"
    return path

def impProtein(name):
    shutil.copy2("molecules/protein/"+name.upper()+".pdbqt", "temp/protein.pdbqt")

def cleanTemp():
    os.remove("temp/ligand.pdbqt")
    os.remove("temp/protein.pdbqt")

def usage():
    print "\n Memetic algorithm for Molecular Docking"
    print " -l: ligand name\n -p: complex name"
    print "\t -l    -p"
    print "\t NAD - 1ENY"
    print "\t KNI - 1HPX"
    print "\t NMB - 1AJV"
    print "\t U02 - 2UPJ"
    print "\t XV6 - 1BV9\n"

def setAnglesPath(lig):
    path = "molecules/anglesKB/"+str(lig)+"/"
    return path

def configParameters(name):
    params = open("ligands-config.txt", "r")
    content = params.readlines()
    params.close()
    flag = 0
    data = []
    for line in content:
        lig = line.strip().split("_")
        lig[1] = lig[1].split(";")
        data.append(lig)
    for ligand in data:
        if name == ligand[0]: 
            center_x = ligand[1][0]
            center_y = ligand[1][1]
            center_z = ligand[1][2]
            flag = 1
            break
    if flag == 1:
        dataW = "receptor = temp/protein.pdbqt\n"
        dataW += "ligand = temp/ligand.pdbqt\n"
        dataW += "center_x = "+str(center_x)+"\n"
        dataW += "center_y = "+str(center_y)+"\n"
        dataW += "center_z = "+str(center_z)+"\n"
        dataW += "size_x = 20.0\nsize_y = 20.0\nsize_z = 20.0\n"
        dataW += "out = result.pdbqt"
        configFile = open("energy/config.txt", "w")
        configFile.write(dataW)
        configFile.close()
    else:
        print "ERROR: no config parameters for ligand", name
        exit(2)




