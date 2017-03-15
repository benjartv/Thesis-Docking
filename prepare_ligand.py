from algorithm_py.Memetic import *
from objects_py.molecule import *
from utils_py.params import *
from utils_py.utils import *
import random
import sys
import getopt
import shutil
import copy

__namemolecule = "KNI"
__searchSpace = 2

if __name__ == '__main__':
	ormolecule = Molecule()
	path = setPathLig(__namemolecule)
	ormolecule.readPDBQT(path)
	ormolecule.calculateSegment()

	solution = Gene()
	solution.randomCell(int(ormolecule.torsdof), __searchSpace)
	aux = copy.deepcopy(ormolecule)
	center = aux.findCenter()
	aux.translateToPoint([center[0]+solution.x,
						center[1]+solution.y,
						center[2]+solution.z,])
	sphVect = spherePoint(1, solution.sph_theta, solution.sph_phi)
	aux.rotateByVector(sphVect, solution.theta)
	for i in range(len(aux.branch)):
		aux.rotateAtomsBranch(i, solution.rotateBonds[i])

	aux.writePDBQT(__namemolecule.lower()+"_mod.pdbqt")
