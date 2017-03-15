import numpy as np
from vector import *
from decimal import *
from utils_py.utils import *
import math
import random
import copy
import datetime
import os
import sys

class Molecule(object):

    #-------------------------------- Init ------------------------------
    def __init__(self, recordName = ""):
        self.recordName = recordName
        self.root = []
        self.idatm = []
        self.atm = []
        self.residueName = []
        self.chainNAME = []
        self.chainId = []
        self.x = []
        self.y = []
        self.z = []
        self.occupancy = []
        self.tempFactor = []
        self.charge = []
        self.atmType = []
        self.branch = []
        self.branchSegment = []
        self.torsdof = 0
        self.data = []

    def rotateAtomsBranch(self, idbranch, theta):
        idatm1 = int(self.branch[idbranch][0])
        idatm2 = int(self.branch[idbranch][1])
        vect_atom1 = [float(self.x[idatm1-1]), float(self.y[idatm1-1]), float(self.z[idatm1-1])]
        vect_atom2 = [float(self.x[idatm2-1]), float(self.y[idatm2-1]), float(self.z[idatm2-1])]

        
        delthaX = float(self.x[idatm1-1])
        delthaY = float(self.y[idatm1-1])
        delthaZ = float(self.z[idatm1-1])
        
        for j in xrange(len(self.x)):
            self.x[j] = float(self.x[j]) - delthaX
            self.y[j] = float(self.y[j]) - delthaY
            self.z[j] = float(self.z[j]) - delthaZ
        
        refVect = Vector(vect_atom2[0]-vect_atom1[0],vect_atom2[1]-vect_atom1[1],vect_atom2[2]-vect_atom1[2])
        matrixRot = calcRotM(refVect, theta)

        segment = self.branchSegment[idbranch][:]
        for atm in segment:
            pointRef = [[float(self.x[atm-1])], [float(self.y[atm-1])], [float(self.z[atm-1])]]
            mresult = crossProduct(matrixRot, pointRef)
            self.x[atm-1] = mresult[0][0]
            self.y[atm-1] = mresult[1][0]
            self.z[atm-1] = mresult[2][0]
        
        for j in xrange(len(self.x)):
            self.x[j] = float(self.x[j]) + delthaX
            self.y[j] = float(self.y[j]) + delthaY
            self.z[j] = float(self.z[j]) + delthaZ

    def translateToPoint(self, point):
        center = self.findCenter()
        delta = [point[0]-center[0], point[1]-center[1], point[2]-center[2]]
        for i in range(len(self.x)):
            self.x[i] = float(self.x[i]) + delta[0]
            self.y[i] = float(self.y[i]) + delta[1]
            self.z[i] = float(self.z[i]) + delta[2]

    def findCenter(self):
        sumX = 0
        sumY = 0
        sumZ = 0
        for i in range(len(self.x)):
            sumX += float(self.x[i])
            sumY += float(self.y[i])
            sumZ += float(self.z[i])
        centerX = sumX / len(self.x)
        centerY = sumY / len(self.y)
        centerZ = sumZ / len(self.z)
        
        return [centerX, centerY, centerZ]

    def rotateByVector(self, vector, theta):
        #vector = self.validateNormCero(vector)
        originPoint = self.findCenter()

        unitX = vector[0] + originPoint[0]
        unitY = vector[1] + originPoint[1]
        unitZ = vector[2] + originPoint[2]

        for j in xrange(len(self.x)):
            self.x[j] = float(self.x[j]) - originPoint[0]
            self.y[j] = float(self.y[j]) - originPoint[1]
            self.z[j] = float(self.z[j]) - originPoint[2]

        p1 = [originPoint[0], originPoint[1], originPoint[2]]
        p2 = [unitX, unitY, unitZ]

        vectorRef = Vector(p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
        matrixRot = calcRotM(vectorRef, theta)

        for i in xrange(len(self.x)):
            pointRef = [[float(self.x[i])],[float(self.y[i])],[float(self.z[i])]]
            mresult = crossProduct(matrixRot, pointRef)
            self.x[i] = mresult[0][0] + originPoint[0]
            self.y[i] = mresult[1][0] + originPoint[1]
            self.z[i] = mresult[2][0] + originPoint[2]

    def readPDBQT(self, name):
        if os.path.isfile(name):  
            file = open(name, "r")
            data = file.readlines()
            file.close()
            content = []
            for line in data:
                content.append(line.strip().split())
            i = 1
            while True:
                if content[i][0] == "ENDROOT":
                    break
                self.idatm.append(content[i][1])
                self.atm.append(content[i][2])
                self.residueName.append(content[i][3])
                self.chainNAME.append(content[i][4])
                self.chainId.append(content[i][5])
                self.x.append(content[i][6])
                self.y.append(content[i][7])
                self.z.append(content[i][8])
                self.occupancy.append(content[i][9])
                self.tempFactor.append(content[i][10])
                self.charge.append(content[i][11])
                self.atmType.append(content[i][12])
                self.root.append(content[i][1])
                i+=1
            while content[i][0] != "TORSDOF":
                if content[i][0] == "BRANCH":
                    self.branch.append([content[i][1],content[i][2]])
                    self.data.append(["B",content[i][1],content[i][2]])
                elif content[i][0] == "HETATM":
                    self.idatm.append(content[i][1])
                    self.atm.append(content[i][2])
                    self.residueName.append(content[i][3])
                    self.chainNAME.append(content[i][4])
                    self.chainId.append(content[i][5])
                    self.x.append(content[i][6])
                    self.y.append(content[i][7])
                    self.z.append(content[i][8])
                    self.occupancy.append(content[i][9])
                    self.tempFactor.append(content[i][10])
                    self.charge.append(content[i][11])
                    self.atmType.append(content[i][12])
                    self.data.append(["H",content[i][1]])
                    #self.root.append(content[i][1])
                elif content[i][0] == "ENDBRANCH":
                    self.data.append(["E",content[i][1], content[i][2]])
                i+=1
            self.torsdof = content[i][1]
        else:
            print "ERROR: ligand file not found."
            sys.exit(2)

    def writePDBQT(self, name):
        file = open(name, 'w')
        i = 0
        file.write("ROOT\n")
        while i < len(self.root):
            file.write("{:<6}".format("HETATM"))
            file.write("{:>5}".format(str(self.idatm[i])))
            file.write("{:>4}".format(str(self.atm[i])))
            file.write("{:>5}".format(self.residueName[i]))
            file.write("{:>2}".format(self.chainNAME[i]))
            file.write("{:>4}".format(self.chainId[i]))
            file.write("{:>12}".format(str(Decimal(self.x[i]).quantize(Decimal('1.000')))))
            file.write("{:>8}".format(str(Decimal(self.y[i]).quantize(Decimal('1.000')))))
            file.write("{:>8}".format(str(Decimal(self.z[i]).quantize(Decimal('1.000')))))
            file.write("{:>6}".format(self.occupancy[i]))
            file.write("{:>6}".format(self.tempFactor[i]))
            file.write("{:>10}".format(self.charge[i]))
            file.write("{:>4}".format(self.atmType[i]))
            file.write("\n")
            i += 1
        file.write("ENDROOT\n")
        bcount = 0
        hcount = int(self.root[-1])
        ecount = 0
        for info in self.data:
            if info[0] == "B":
                file.write("BRANCH")
                file.write("{:>4}".format(str(info[1])))
                file.write("{:>4}".format(str(info[2])))
                file.write("\n")
                bcount += 1
                ecount += 1
            elif info[0] == "H":
                file.write("{:<6}".format("HETATM"))
                file.write("{:>5}".format(str(self.idatm[hcount])))
                file.write("{:>4}".format(str(self.atm[hcount])))
                file.write("{:>5}".format(self.residueName[hcount]))
                file.write("{:>2}".format(self.chainNAME[hcount]))
                file.write("{:>4}".format(self.chainId[hcount]))
                file.write("{:>12}".format(str(Decimal(self.x[hcount]).quantize(Decimal('1.000')))))
                file.write("{:>8}".format(str(Decimal(self.y[hcount]).quantize(Decimal('1.000')))))
                file.write("{:>8}".format(str(Decimal(self.z[hcount]).quantize(Decimal('1.000')))))
                file.write("{:>6}".format(self.occupancy[hcount]))
                file.write("{:>6}".format(self.tempFactor[hcount]))
                file.write("{:>10}".format(self.charge[hcount]))
                file.write("{:>4}".format(self.atmType[hcount]))
                file.write("\n")
                hcount += 1
            elif info[0] == "E":
                file.write("ENDBRANCH")
                file.write("{:>4}".format(str(info[1])))
                file.write("{:>4}".format(str(info[2])))
                ecount -= 1
                file.write("\n")

        
        file.write("TORSDOF")
        file.write("{:>3}".format(str(self.torsdof)))
        #print "done"
        file.close()

    def generateRandomCenter(self, center):
        newpoint = center[:]
        axis = random.randint(0,2)
        if axis == 0:
            newpoint[0] = random.uniform(-2,2)+center[0]
        elif axis == 1:
            newpoint[1] = random.uniform(-2,2)+center[1]
        else:
            newpoint[2] = random.uniform(-2,2)+center[2]
        self.translateToPoint(newpoint)
        return newpoint[:]

    def generateRandomRotation(self, vectorRef):
        alpha = random.randint(0,2)
        vector = vectorRef[:]
        if alpha == 0:
            vector[0] = random.uniform(0,2)*math.pi
        elif alpha == 1:
            vector[1] = random.uniform(0,1)*math.pi
        else:
            vector[2] = random.uniform(0,2)*math.pi
        sphVect = spherePoint(1, vector[0], vector[1])
        self.rotateByVector(sphVect, vector[2])
        return vector[:]

    def generateRandomAngles(self):
        branch = len(self.branch) - 1
        r = random.randint(0, branch)
        self.rotateAtomsBranch(r, random.uniform(0,2*math.pi))

   

    #Calculate the segment of each branch, save the data on rotationSegment list
    #each position fit with the position of the branch list
    def calculateSegment(self):
        rootbranch = []
        for bond in self.branch:
            if bond[0] in self.root:
                rootbranch.append(int(bond[1]))
        rootbranch.sort()
        rootbranch.append(-1)
        i = 0
        while i < len(self.branch):
            if int(self.branch[i][0]) in rootbranch:
                pos = rootbranch.index(int(self.branch[i][0]))
                if rootbranch[pos+1] != -1:
                    segment = range(int(self.branch[i][1]), int(rootbranch[pos+1]))
                else:
                    segment = range(int(self.branch[i][1]), int(self.idatm[-1])+1)
                self.branchSegment.append(segment)
            else:
                for r in rootbranch:
                    if int(self.branch[i][1]) < r:
                        segment = range(int(self.branch[i][1]),r)
                        self.branchSegment.append(segment)
                        break
                    elif r == -1:
                        segment = range(int(self.branch[i][1]), int(self.idatm[-1])+1)
                        self.branchSegment.append(segment)
                        break
            i+=1

    def validateNormCero(self,vector):
        res = 0.0
        res += vector[0] * vector[0]
        res += vector[1] * vector[1]
        res += vector[2] * vector[2]
        res = math.sqrt(res)

        if(res<0.2):
            vector = randomSpherePoint(1)
        return vector



























