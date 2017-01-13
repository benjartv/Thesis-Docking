from vector import *
import numpy as np
import math
import commands
import string

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
    a = commands.getstatusoutput('./vina --config config.txt --score_only')
    energiaVina = string.split(string.split(a[1],"Affinity:")[1])[0]
    #print energiaVina
    return float(energiaVina)