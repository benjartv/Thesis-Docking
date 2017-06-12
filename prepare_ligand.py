from algorithm_py.Memetic import *
from objects_py.molecule import *
from utils_py.params import *
from utils_py.utils import *
import random
import sys
import getopt
import shutil
import copy
import math

__namemolecule = "AH1"
__searchSpace = 10

if __name__ == '__main__':
	ormolecule = Molecule()
	path = setPathLig(__namemolecule)
	ormolecule.readPDBQT(path)
	ormolecule.calculateSegment()
	#print ormolecule.findCenter()
	for i in ormolecule.branchSegment:
		print i
	
	solution = Gene()
	center = ormolecule.findCenter()
	newcenter = [random.uniform(-10,10)+center[i] for i in range(3)]
	solution.randomCell(int(ormolecule.torsdof), 3)

	for i in range(len(ormolecule.branch)):
		ormolecule.rotateAtomsBranch(i, solution.rotateBonds[i])
	ormolecule.translateToPoint([newcenter[0]+solution.x, 
								newcenter[1]+solution.y, 
								newcenter[2]+solution.z])
	sphVect = spherePoint(1, solution.sph_theta, solution.sph_phi)
	ormolecule.rotateByVector(sphVect, solution.theta)
	
	ormolecule.writePDBQT(__namemolecule.lower()+"_mod.pdbqt")
	