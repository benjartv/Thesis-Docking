from Memetic import *
from molecule import *
from params import *
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
	__generations = 300
	__pocketSize = 5
	__treeNodes = 13
	__mutProbability = 0.2
	__isLocalSearch = False
	__typeCO = 0
	__typeLS = 0
	__distanceCriteria = None
	__nodeByTree = 3

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
						__nodeByTree)

	Memetic(parameters, ligand, NADOR).initProcess()
