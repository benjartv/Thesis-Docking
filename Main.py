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
	
	__isKB = False
	__KBProb = 0.8 #Probability of using knowledge base

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
	if __isKB:
		print "Import Knowledge Base...",
		kbname = "kb-"+__ligand.lower()
		files = readKBfileName(kbname)
		anglesPath = setAnglesPath(__ligand)
		modLigand.importAngles(files, anglesPath)
		print "Complete."
	print "Set vina configuration...",
	configParameters(__ligand)
	print "Complete."	

	spaceCenter = originalLigand.findCenter()
	__searchSpace = 3 #search space for translate center of ligand
	newSearchSpace = [random.uniform(-__searchSpace,__searchSpace)+spaceCenter[i] for i in range(3)]
	__generations = 1220 #number of generations until the algorithm stop
						#1120 for 1.000.000 energy evaluation
	__pocketSize = 6 #size of the pocket of each agent
	__treeNodes = 13 #number of nodes of the hierarchical tree
	__mutProbability = 0.2 #probability of mutation
	__isLocalSearch = True

	__typeCO = 3 #type of Crossover
	# 0: crossover Uniform
	# 1: crossover Block
	# 2: crossover SPC
	# 3: crossover 50/50
	# 4: crossover only Center (ONLY for LS 2)

	__typeLS = 2 #type of Local Search
	# 0: mutation block
	# 1: mutation spacereduce by iteration
	# 2: mutation only Rotation (ONLY for CS 1)

	__typeMut = 1 #type of Mutation (Memetic)
	# 0: mutation uniform (not recomended)
	# 1: mutation block
	__distanceCriLVL = [1.5,1.0,0.5] #Aceptance criterion for each lvl of the tree
	__nodeByTree = 3 #number of agent for tree-level
	__tempLS = 1000.0 #initial temperature for simulated annealing (LS)
	__minTemp = 1.0 #final temperature for simulated annealing (LS)
	__alphaTemp = 0.9 #alpha for simulated annealing (LS)
	__numberIteration = 1 #by Local Search loop
	__reset = 100 #number of generation between each reset (-1 for non reset)
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
						__nodeByTree,
						__tempLS,
						__minTemp,
						__alphaTemp,
						__numberIteration,
						__reset,
						__typeReset,
						__typeMut,
						__isKB,
						__KBProb,
						__distanceCriLVL)
	print "Init memetic algorithm..."
	Memetic(parameters, modLigand, originalLigand).initProcess()
	print "Removing temporal data..."
	cleanTemp()
	print "Done."
