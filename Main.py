from algorithm_py.Memetic import *
from objects_py.molecule import *
from utils_py.params import *
from utils_py.utils import *
import random
import sys
import getopt


__molecules = initConfig()
__ligand = __molecules[0]
__protein = __molecules[1]


if __name__ == "__main__":
	print "Preparing molecules..."
	originalLigand = Molecule(__ligand)

	print "Ligand: ", __ligand, " Complex: ", __protein
	print "Importing ligand...",
	moleculePath = setPathLig(__ligand)
	originalLigand.readPDBQT(moleculePath)
	originalLigand.calculateSegment()
	print "Complete."

	print "Importing protein...",
	impProtein(__protein)
	print "Complete."

	print "Preparing Ligand...",
	modligPath = setPathMLig(__ligand)
	modLigand = Molecule(__ligand)
	modLigand.readPDBQT(modligPath)
	modLigand.calculateSegment()
	print "Complete."

	print "Set vina configuration...",
	configParameters(__ligand)
	print "Complete."	

	spaceCenter = originalLigand.findCenter()
	__searchSpace = 2 #search space for translate center of ligand
	newSearchSpace = [random.uniform(-__searchSpace,__searchSpace)+spaceCenter[i] for i in range(3)]
	__generations = 10 #number of generations until the algorithm stop
						#1120 for 1.000.000 energy evaluation
	__pocketSize = 5 #size of the pocket of each agent
	__treeNodes = 13 #number of nodes of the hierarchical tree
	__mutProbability = 0.2 #probability of mutation
	__isLocalSearch = True 
	__typeCO = 0 #type of Crossover
	# 0: crossover Uniform
	# 1: crossover only center
	__typeLS = 1 #type of Local Search
	# 0: mutation Uniform
	# 1: mutation block
	# 2: mutation only Rotation (for crossover 1)
	__typeMut = 1 #type of Mutation (Memetic)
	# 0: mutation uniform
	# 1: mutation block
	__distanceCriteria = 0.5 #Acceptance criterion for each solution (gene)
	__nodeByTree = 3 #number of agent for tree-level
	__tempLS = 1000.0 #initial temperature for simulated annealing (LS)
	__minTemp = 1.0 #final temperature for simulated annealing (LS)
	__alphaTemp = 0.9 #alpha for simulated annealing (LS)
	__numberIteration = 1 #by Local Search loop
	__reset = -1 #number of generation between each reset (-1 for non reset)
	__typeReset = 0 #Type of reset
	# 0: generation reset
	# 1: molecule reset


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
						__numberIteration,
						__reset,
						__typeReset,
						__typeMut)
	print "Init memetic algorithm..."
	Memetic(parameters, modLigand, originalLigand).initProcess()
	print "Removing temporal data..."
	cleanTemp()
	print "Done."
