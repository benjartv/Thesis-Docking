from algorithm_py.Memetic import *
from objects_py.molecule import *
from utils_py.params import *
import random


if __name__ == "__main__":
	print "Preparing molecules..."
	NADOR = Molecule("NAD")
	NADOR.readPDBQT("lig.pdbqt")
	NADOR.calculateSegment()

	ligand = Molecule("NAD")
	ligand.readPDBQT("modify.pdbqt")
	ligand.calculateSegment()

	spaceCenter = NADOR.findCenter()
	__searchSpace = 2
	newSearchSpace = [random.uniform(-__searchSpace,__searchSpace)+spaceCenter[i] for i in range(3)]
	__generations = 2
	__pocketSize = 5
	__treeNodes = 13
	__mutProbability = 0.2
	__isLocalSearch = False
	__typeCO = 0
	__typeLS = 0
	__distanceCriteria = 2.0
	__nodeByTree = 3
	__tempLS = 1000.0
	__minTemp = 1.0
	__alphaTemp = 0.9
	__numberIteration = 1 #by Local Search loop


	parameters = params(__searchSpace,
						newSearchSpace,
						__generations,
						__pocketSize,
						__treeNodes,
						__mutProbability,
						__isLocalSearch,
						__typeLS,
						__typeCO,
						__distanceCriteria,
						__nodeByTree,
						__tempLS,
						__minTemp,
						__alphaTemp,
						__numberIteration)
	print "Init memetic algorithm..."
	Memetic(parameters, ligand, NADOR).initProcess()
