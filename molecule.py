import numpy as np
from vector import *
from decimal import *
import math
import random
import commands
import string
import copy

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

    def rotatateBranch(self, idatm1, idatm2, idbranch, theta):
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
        

    def readPDBQT(self, name):
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
            file.write("{:>12}".format(str(self.x[i])))
            file.write("{:>8}".format(str(self.y[i])))
            file.write("{:>8}".format(str(self.z[i])))
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


    def generateRandom(self):
        r = random.randint(0,10)
        if r == 0:
            self.rotatateBranch(1,4,0,random.uniform(0,2*math.pi))
        elif r == 1:
            self.rotatateBranch(4,5,1,random.uniform(0,2*math.pi))
        elif r == 2:
            self.rotatateBranch(5,6,2,random.uniform(0,2*math.pi))
        elif r == 3:
            self.rotatateBranch(10,13,3,random.uniform(0,2*math.pi))
        elif r == 4:
            self.rotatateBranch(1,23,4,random.uniform(0,2*math.pi))
        elif r == 5:
            self.rotatateBranch(23,24,5,random.uniform(0,2*math.pi))
        elif r == 6:
            self.rotatateBranch(24,27,6,random.uniform(0,2*math.pi))
        elif r == 7:
            self.rotatateBranch(27,28,7,random.uniform(0,2*math.pi))
        elif r == 8:
            self.rotatateBranch(28,29,8,random.uniform(0,2*math.pi))
        elif r == 9:
            self.rotatateBranch(31,36,9,random.uniform(0,2*math.pi))
        elif r == 10:
            self.rotatateBranch(38,42,10,random.uniform(0,2*math.pi))


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

ligand = Molecule("NAD")
ligand.readPDBQT("lig.pdbqt")
ligand.calculateSegment()

#for seg in ligand.branchSegment:
 #   print seg

print ligand.branch

ligand.writePDBQT("ligand.pdbqt")
original_score = calculateFreeEnergy()
print original_score

#ligand.rotatateBranch(1,4,0,math.pi)
ligand.rotatateBranch(1,4,0,random.uniform(0,2*math.pi))
ligand.rotatateBranch(4,5,1,random.uniform(0,2*math.pi))
ligand.rotatateBranch(5,6,2,random.uniform(0,2*math.pi))
ligand.rotatateBranch(10,13,3,random.uniform(0,2*math.pi))
ligand.rotatateBranch(1,23,4,random.uniform(0,2*math.pi))
ligand.rotatateBranch(23,24,5,random.uniform(0,2*math.pi))
ligand.rotatateBranch(24,27,6,random.uniform(0,2*math.pi))
ligand.rotatateBranch(27,28,7,random.uniform(0,2*math.pi))
ligand.rotatateBranch(28,29,8,random.uniform(0,2*math.pi))
ligand.rotatateBranch(31,36,9,random.uniform(0,2*math.pi))
ligand.rotatateBranch(38,42,10,random.uniform(0,2*math.pi))

best_score = 1000.0
best_ligand = copy.deepcopy(ligand)
gen = 0

for i in range(0,10000):
    if i%100 == 0:
        print "generation: ", i
    new_ligand = copy.deepcopy(ligand)
    new_ligand.generateRandom()
    new_ligand.writePDBQT("ligand.pdbqt")
    new_score = calculateFreeEnergy()
    if new_score < best_score:
        best_score = new_score
        best_ligand = copy.deepcopy(new_ligand)
        gen = i


best_ligand.writePDBQT("best.pdbqt")
print "Best score: ", best_score
print "Gen: ",gen

#print ligand.data
#ligand.writePDBQT("ligand.pdbqt")
#print ligand.torsdof
#print ligand.branch[7]
#print ligand.branchSegment[7]

#print ligand.idatm
#print ligand.root
'''
print ligand.idatm
print ligand.atm
'''

























