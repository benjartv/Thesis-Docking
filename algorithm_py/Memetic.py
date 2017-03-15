from objects_py.molecule import *
from objects_py.Gene import *
from objects_py.agent import *
from LocalSearch import *
import datetime
import copy
import os

class Memetic(object):
	def __init__(self, params, ligand, originalLig):
		self.__ligand = ligand
		self.__originalLigand = originalLig
		self.__searchSpace = params.searchSpace
		self.__centerSpace = params.centerSpace
		self.__generations = params.generations
		self.__pocketSize = params.pocketSize
		self.__treeNodes = params.treeNodes
		self.__mutProbability = params.mutProbability
		self.__isLocalSearch = params.isLocalSearch
		self.__LocalSearch = None
		self.__typeCO = params.typeCO
		self.__distanceCriteria = params.distanceCriteria
		self.__Time = datetime.datetime.now()
		self.__nodeByTree = params.nodeByTree
		self.__rootNode = None
		self.__fatherNode = [None]*self.__nodeByTree
		self.__leafNode = [None]*self.__nodeByTree*len(self.__fatherNode)
		self.__logPop = ""
		self.__logData = ""
		self.__numberScoring = 0
		self.__dirResult = ""
		self.__temporalDir = "temp/"
		self.__bestScore = 9999.9
		self.__bestGeneration = 0
		self.__isReset = False
		self.__rCount = 0
		self.__reset = params.reset
		if self.__reset != -1:
			self.__isReset = True

		if self.__isLocalSearch:
			self.__LocalSearch = LocalSearch(params.tempLS,
											params.minTemp,
											params.alphaTemp,
											self.__ligand,
											self.__searchSpace,
											self.__centerSpace,
											params.typeLS,
											params.numberIteration)
	
	def initProcess(self):
		startTime = datetime.datetime.now()
		time = convertTime(startTime)
		self.__dirResult = "results/Result_"+self.__ligand.recordName+"_"+time
		os.makedirs(self.__dirResult)
		self.initTree()
		self.initPopulation()
		for i in range(self.__generations):
			self.generation()
			bestNode = self.__rootNode.getBest().score
			self.updateBest(bestNode, i)
			if self.__isReset:
				self.resetCount(i)
			print "Generation: ", i, " Best: ", bestNode
			self.addlogPopulation(i)
			print "Root: ", self.__rootNode.getPocketScore()
		bestCell = self.__rootNode.getBest()
		print "Best ligand score: ", bestCell.score
		auxLigand = self.generateLigand(bestCell)
		auxRMSD = getRMSD(auxLigand, self.__originalLigand)
		print "RMSD: ", auxRMSD
		for i in range(self.__pocketSize):
			if self.__rootNode.pocket[i] != None:
				self.generateFinalBest(self.__rootNode.pocket[i], self.__dirResult+"/"+"best-lig-"+str(i)+".pdbqt")
		stopTime = datetime.datetime.now()
		self.__Time = stopTime - startTime
		print "Time: ", self.__Time
		if self.__isLocalSearch:
			self.__numberScoring += self.__LocalSearch.getNumberEvaluation()
		print "Number of Energy Evaluation: ", self.__numberScoring

		######### Add log info #########

		self.__logData += "Best ligand score: "+str(self.__bestScore)+"\n"
		self.__logData += "Reached on generation: "+str(self.__bestGeneration)+"\n"
		self.__logData += "RMSD: "+ str(auxRMSD)+"\n"
		self.__logData += "Time: "+ str(self.__Time)+"\n"
		self.__logData += "Number of Energy Evaluation: "+ str(self.__numberScoring)+"\n"

		print "Process complete"
		self.writeLog()

	def initTree(self):
		self.__rootNode = agent(self.__pocketSize,
							self.__distanceCriteria, 
							self.__ligand,
							self.__centerSpace)
		for i in range(self.__nodeByTree):
			self.__fatherNode[i] = agent(self.__pocketSize,
							self.__distanceCriteria, 
							self.__ligand,
							self.__centerSpace)
		for i in range((self.__nodeByTree*len(self.__fatherNode))):
			self.__leafNode[i] = agent(self.__pocketSize,
							self.__distanceCriteria, 
							self.__ligand,
							self.__centerSpace)

	def initPopulation(self):
		gene = Gene()
		gene.randomCell(len(self.__ligand.branchSegment), self.__searchSpace)
		gene = self.calculates(gene)
		self.__rootNode.addToPocket(copy.deepcopy(gene))
		for n in range(len(self.__fatherNode)):
			gene = Gene()
			gene.randomCell(len(self.__ligand.branchSegment), self.__searchSpace)
			gene = self.calculates(gene)
			self.__fatherNode[n].addToPocket(copy.deepcopy(gene))
		for n in range(len(self.__leafNode)):
			gene = Gene()
			gene.randomCell(len(self.__ligand.branchSegment), self.__searchSpace)
			gene = self.calculates(gene)
			self.__leafNode[n].addToPocket(copy.deepcopy(gene))
		self.initLog()

	def updateBest(self, score, generation):
		if score < self.__bestScore:
			self.__bestScore = score
			self.__bestGeneration = generation

	def calculates(self, cell):
		if self.__isLocalSearch:
			newcell = self.__LocalSearch.initLocalSearch(cell)
		else:
			newcell = copy.deepcopy(cell)
		auxLigand = copy.deepcopy(self.__ligand)
		auxLigand.translateToPoint([self.__centerSpace[0]+newcell.x, 
									self.__centerSpace[1]+newcell.y, 
									self.__centerSpace[2]+newcell.z])
		sphVect = spherePoint(1, newcell.sph_theta, newcell.sph_phi)
		auxLigand.rotateByVector(sphVect, newcell.theta)
		for i in range(len(auxLigand.branch)):
			auxLigand.rotateAtomsBranch(i, newcell.rotateBonds[i])
		auxLigand.writePDBQT(self.__temporalDir+"ligand.pdbqt")
		newcell.score = calculateFreeEnergy()
		self.__numberScoring += 1
		return copy.deepcopy(newcell)

	def generateLigand(self, cell):
		auxLigand = copy.deepcopy(self.__ligand)
		auxLigand.translateToPoint([self.__centerSpace[0]+cell.x, 
									self.__centerSpace[1]+cell.y, 
									self.__centerSpace[2]+cell.z])
		sphVect = spherePoint(1, cell.sph_theta, cell.sph_phi)
		auxLigand.rotateByVector(sphVect, cell.theta)
		for i in range(len(auxLigand.branch)):
			auxLigand.rotateAtomsBranch(i, cell.rotateBonds[i])
		return copy.deepcopy(auxLigand)

	def generateFinalBest(self, cell, name="best-ligand.pdbqt"):
		auxLigand = copy.deepcopy(self.__ligand)
		auxLigand.translateToPoint([self.__centerSpace[0]+cell.x, 
									self.__centerSpace[1]+cell.y, 
									self.__centerSpace[2]+cell.z])
		sphVect = spherePoint(1, cell.sph_theta, cell.sph_phi)
		auxLigand.rotateByVector(sphVect, cell.theta)
		for i in range(len(auxLigand.branch)):
			auxLigand.rotateAtomsBranch(i, cell.rotateBonds[i])
		auxLigand.writePDBQT(name)

	def generation(self):
		for i in range(len(self.__fatherNode)):
			for j in range(self.__nodeByTree):
				pop1 = self.__fatherNode[i].getRandom()
				pop2 = self.__leafNode[self.__nodeByTree*i+j].getRandom()
				npop = self.crossoverUniform(pop1, pop2)
				npop = self.mutation(npop)
				npop = self.calculates(npop)
				self.__leafNode[self.__nodeByTree*i+j].addToPocket(npop)
			pop1 = self.__rootNode.getRandom()
			pop2 = self.__fatherNode[i].getRandom()
			npop = self.crossoverUniform(pop1, pop2)
			npop = self.mutation(npop)
			npop = self.calculates(npop)
			self.__fatherNode[i].addToPocket(npop)
		self.updateTree()

	def updateTree(self):
		for i in range(len(self.__fatherNode)):
			for j in range(self.__nodeByTree):
				cell = self.__leafNode[self.__nodeByTree*i+j].getBest()
				#print "F: ",i,", L: ",j,", score: ",cell.score
				self.__fatherNode[i].addToPocket(cell)
			cell = self.__fatherNode[i].getBest()
			#print "F: ",i,", to Root, score: ",cell.score
			#print "P-F ",i,":",self.__fatherNode[i].getPocketScore()
			self.__rootNode.addToPocket(cell)

	def crossoverUniform(self, selectedCell1, selectedCell2):
		newCell = Gene()
		if random.randint(0,1) == 1:
			newCell.x = selectedCell1.x
		else:
			newCell.x = selectedCell2.x
		if random.randint(0,1) == 1:
			newCell.y = selectedCell1.y
		else:
			newCell.y = selectedCell2.y
		if random.randint(0,1) == 1:
			newCell.z = selectedCell1.z
		else:
			newCell.z = selectedCell2.z
		if random.randint(0,1) == 1:
			newCell.sph_theta = selectedCell1.sph_theta
		else:
			newCell.sph_theta = selectedCell2.sph_theta
		if random.randint(0,1) == 1:
			newCell.sph_phi = selectedCell1.sph_phi
		else:
			newCell.sph_phi = selectedCell2.sph_phi
		if random.randint(0,1) == 1:
			newCell.theta = selectedCell1.theta
		else:
			newCell.theta = selectedCell2.theta
		for i in range(len(selectedCell1.rotateBonds)):
			if random.randint(0,1)==1:
				newCell.rotateBonds.append(selectedCell1.rotateBonds[i])
			else:
				newCell.rotateBonds.append(selectedCell2.rotateBonds[i])
		return newCell

	def crossoverBlock(self, selectedCell1, selectedCell2):
		newCell = Gene()
		if random.randint(0,1) == 1:
			newCell.x = selectedCell1.x
			newCell.y = selectedCell1.y
			newCell.z = selectedCell1.z
		else:
			newCell.x = selectedCell2.x
			newCell.y = selectedCell2.y
			newCell.z = selectedCell2.z
		if random.randint(0,1) == 1:
			newCell.sph_theta = selectedCell1.sph_theta
			newCell.sph_phi = selectedCell1.sph_phi
			newCell.theta = selectedCell1.theta
		else:
			newCell.sph_theta = selectedCell2.sph_theta
			newCell.sph_phi = selectedCell2.sph_phi
			newCell.theta = selectedCell2.theta
		if random.randint(0,1) == 1:
			newcell.rotateBonds = selectedCell1.rotateBonds[:]
		else:
			newCell.rotateBonds = selectedCell2.rotateBonds[:]
		return newCell

	def crossover50(self, selectedCell1, selectedCell2):
		newCell = Gene()
		#switch center of ligand
		centrand = random.randint(0,5)
		if centrand == 0:
			newCell.x = selectedCell1.x
			newCell.y = selectedCell2.y
			newCell.z = selectedCell2.z
		elif centrand == 1:
			newCell.x = selectedCell1.x
			newCell.y = selectedCell1.y
			newCell.z = selectedCell2.z
		elif centrand == 2:
			newCell.x = selectedCell1.x
			newCell.y = selectedCell1.y
			newCell.z = selectedCell1.z
		elif centrand == 3:
			newCell.x = selectedCell2.x
			newCell.y = selectedCell2.y
			newCell.z = selectedCell2.z
		elif centrand == 4:
			newCell.x = selectedCell2.x
			newCell.y = selectedCell1.y
			newCell.z = selectedCell1.z
		else:
			newCell.x = selectedCell2.x
			newCell.y = selectedCell2.y
			newCell.z = selectedCell1.z
		#switch rotation of ligand
		rotrand = random.randint(0,5)
		if rotrand == 0:
			newCell.sph_theta = selectedCell1.sph_theta
			newCell.sph_phi = selectedCell2.sph_phi
			newCell.theta = selectedCell2.theta
		elif rotrand == 1:
			newCell.sph_theta = selectedCell1.sph_theta
			newCell.sph_phi = selectedCell1.sph_phi
			newCell.theta = selectedCell2.theta
		elif rotrand == 2:
			newCell.sph_theta = selectedCell1.sph_theta
			newCell.sph_phi = selectedCell1.sph_phi
			newCell.theta = selectedCell1.theta
		elif rotrand == 3:
			newCell.sph_theta = selectedCell2.sph_theta
			newCell.sph_phi = selectedCell2.sph_phi
			newCell.theta = selectedCell2.theta
		elif rotrand == 4:
			newCell.sph_theta = selectedCell2.sph_theta
			newCell.sph_phi = selectedCell1.sph_phi
			newCell.theta = selectedCell1.theta
		else:
			newCell.sph_theta = selectedCell2.sph_theta
			newCell.sph_phi = selectedCell2.sph_phi
			newCell.theta = selectedCell1.theta
		#switch rotational bonds
		bondrand = random.randint(0, (len(selectedCell1.rotateBonds)-1))
		if random.randint(0,1) == 1:
			newCell.rotateBonds = selectedCell1.rotateBonds[:bondrand]+selectedCell2.rotateBonds[bondrand:]
		else:
			newCell.rotateBonds = selectedCell2.rotateBonds[:bondrand]+selectedCell1.rotateBonds[bondrand:]
		return newCell

	def crossoverSPC(self, selectedCell1, selectedCell2):
		pop1 = []
		pop2 = []

		pop1.append(selectedCell1.x)
		pop1.append(selectedCell1.y)
		pop1.append(selectedCell1.z)
		pop1.append(selectedCell1.sph_theta)
		pop1.append(selectedCell1.sph_phi)
		pop1.append(selectedCell1.theta)

		for angle in selectedCell1.rotateBonds:
			pop1.append(angle)

		pop2.append(selectedCell2.x)
		pop2.append(selectedCell2.y)
		pop2.append(selectedCell2.z)
		pop2.append(selectedCell2.sph_theta)
		pop2.append(selectedCell2.sph_phi)
		pop2.append(selectedCell2.theta)

		for angle in selectedCell2.rotateBonds:
			pop2.append(angle)

		genSize = len(pop1)
		cutPoint = random.randint(0, genSize-1)
		length = random.randint(0, genSize-1)
		if length == 0:
			return selectedCell1
		elif length > (genSize - cutPoint):
			newPop = pop1[cutPoint:]
			turnTop = length - (genSize-cutPoint)
			for i in range(0, turnTop):
				newPop.insert(i, pop1[i])
			for i in range(turnTop, cutPoint):
				newPop.insert(i, pop2[i])
		else:
			newPop = pop1[cutPoint:cutPoint+length]
			for i in range(0, cutPoint):
				newPop.insert(i, pop2[i])
			for i in range(cutPoint+length, genSize):
				newPop.insert(i, pop2[i])

		newCell = Gene()
		newCell.x = newPop[0]
		newCell.y = newPop[1]
		newCell.z = newPop[2]
		newCell.sph_theta = newPop[3]
		newCell.sph_phi = newPop[4]
		newCell.theta = newPop[5]
		newCell.rotateBonds = newPop[6:]
		return newCell

	def initLog(self):
		node = 1
		self.__logData += "Molecule Ligand: "+self.__ligand.recordName+"\n"
		self.__logData += "Pocket: "+str(self.__pocketSize)+"\n"
		self.__logData += "Generations: "+str(self.__generations)+"\n"
		self.__logData += "Rotate bonds: "+str(self.__ligand.branch)+"\n"
		self.__logPop += "\n"
		#self.__logPop += "Best Score: "+str(self.__rootNode.getBest().score)+"\n"
		self.__logPop += "\n"
		self.__logPop += "Init Population\n"+"*"*74+"\n"
		self.__logPop += "Node 0 (Root)"+"\n"
		for i in range(self.__pocketSize):
			if self.__rootNode.pocket[i] == None:
				self.__logPop += str(i+1)+"-  "+str(self.__rootNode.pocket[i])+"\n"
			else:
				self.__logPop += str(i+1)+"-  "+str([self.__rootNode.pocket[i].x,
												self.__rootNode.pocket[i].y,
												self.__rootNode.pocket[i].z,
												self.__rootNode.pocket[i].sph_theta,
												self.__rootNode.pocket[i].sph_phi,
												self.__rootNode.pocket[i].theta,
												self.__rootNode.pocket[i].rotateBonds])
				self.__logPop += "score: "+str(self.__rootNode.pocket[i].score)+"\n"
		self.__logPop += "\n"
		for j in range(len(self.__fatherNode)):
			self.__logPop += "Node "+str(node)+"\n"
			for i in range(self.__pocketSize):
				if self.__fatherNode[j].pocket[i] == None:
					self.__logPop += str(i+1)+"-  "+str(self.__fatherNode[j].pocket[i])+"\n"
				else:	
					self.__logPop += str(i+1)+"-  "+str([self.__fatherNode[j].pocket[i].x,
													self.__fatherNode[j].pocket[i].y,
													self.__fatherNode[j].pocket[i].z,
													self.__fatherNode[j].pocket[i].sph_theta,
													self.__fatherNode[j].pocket[i].sph_phi,
													self.__fatherNode[j].pocket[i].theta,
													self.__fatherNode[j].pocket[i].rotateBonds])
					self.__logPop += "score: "+str(self.__fatherNode[j].pocket[i].score)+"\n"
			node += 1
			self.__logPop += "\n"
		for j in range(len(self.__leafNode)):
			self.__logPop += "Node "+str(node)+"\n"
			for i in range(self.__pocketSize):
				if self.__leafNode[j].pocket[i] == None:
					self.__logPop += str(i+1)+"-  "+str(self.__leafNode[j].pocket[i])+"\n"
				else:
					self.__logPop += str(i+1)+"-  "+str([self.__leafNode[j].pocket[i].x,
													self.__leafNode[j].pocket[i].y,
													self.__leafNode[j].pocket[i].z,
													self.__leafNode[j].pocket[i].sph_theta,
													self.__leafNode[j].pocket[i].sph_phi,
													self.__leafNode[j].pocket[i].theta,
													self.__leafNode[j].pocket[i].rotateBonds])
					self.__logPop += "score: "+str(self.__leafNode[j].pocket[i].score)+"\n"
			node += 1
			self.__logPop += "\n"

	def addlogPopulation(self, genera):
		node = 1
		self.__logPop += "Generation "+str(genera)+"\n"+"*"*74+"\n"
		self.__logPop += "Node 0 (Root)"+"\n"
		for i in range(self.__pocketSize):
			if self.__rootNode.pocket[i] == None:
				self.__logPop += str(i+1)+"-  "+str(self.__rootNode.pocket[i])+"\n"
			else:
				self.__logPop += str(i+1)+"-  "+str([self.__rootNode.pocket[i].x,
												self.__rootNode.pocket[i].y,
												self.__rootNode.pocket[i].z,
												self.__rootNode.pocket[i].sph_theta,
												self.__rootNode.pocket[i].sph_phi,
												self.__rootNode.pocket[i].theta,
												self.__rootNode.pocket[i].rotateBonds])
				self.__logPop += "score: "+str(self.__rootNode.pocket[i].score)+"\n"
		self.__logPop += "\n"
		for j in range(len(self.__fatherNode)):
			self.__logPop += "Node "+str(node)+"\n"
			for i in range(self.__pocketSize):
				if self.__fatherNode[j].pocket[i] == None:
					self.__logPop += str(i+1)+"-  "+str(self.__fatherNode[j].pocket[i])+"\n"
				else:	
					self.__logPop += str(i+1)+"-  "+str([self.__fatherNode[j].pocket[i].x,
													self.__fatherNode[j].pocket[i].y,
													self.__fatherNode[j].pocket[i].z,
													self.__fatherNode[j].pocket[i].sph_theta,
													self.__fatherNode[j].pocket[i].sph_phi,
													self.__fatherNode[j].pocket[i].theta,
													self.__fatherNode[j].pocket[i].rotateBonds])
					self.__logPop += "score: "+str(self.__fatherNode[j].pocket[i].score)+"\n"
			node += 1
			self.__logPop += "\n"
		for j in range(len(self.__leafNode)):
			self.__logPop += "Node "+str(node)+"\n"
			for i in range(self.__pocketSize):
				if self.__leafNode[j].pocket[i] == None:
					self.__logPop += str(i+1)+"-  "+str(self.__leafNode[j].pocket[i])+"\n"
				else:
					self.__logPop += str(i+1)+"-  "+str([self.__leafNode[j].pocket[i].x,
													self.__leafNode[j].pocket[i].y,
													self.__leafNode[j].pocket[i].z,
													self.__leafNode[j].pocket[i].sph_theta,
													self.__leafNode[j].pocket[i].sph_phi,
													self.__leafNode[j].pocket[i].theta,
													self.__leafNode[j].pocket[i].rotateBonds])
					self.__logPop += "score: "+str(self.__leafNode[j].pocket[i].score)+"\n"
			node += 1
			self.__logPop += "\n"

	def writeLog(self):
		file = open(self.__dirResult+"/iteration.log", "w")
		logWrite = self.__logData+self.__logPop
		file.write(logWrite)
		file.close()

	def resetCount(self, gen):
		if gen > self.__bestGeneration:
			self.__rCount += 1
			if self.__rCount == self.__reset:
				print "reset population..."
				#reset population
				self.resetPopulation()
				self.__rCount = 0
		else:
			self.__rCount == 0

	def resetPopulation(self):
		for node in self.__leafNode:
			bestCell = node.getBest()
			node.resetAgent()
			node.addToPocket(bestCell)
		for node in self.__fatherNode:
			bestCell = node.getBest()
			node.resetAgent()
			node.addToPocket(bestCell)
		bestCell = self.__rootNode.getBest()
		self.__rootNode.resetAgent()
		self.__rootNode.addToPocket(bestCell)
		#gene = Gene()
		#gene.randomCell(len(self.__ligand.branchSegment), self.__searchSpace)
		#gene = self.calculates(gene)
		#self.__rootNode.addToPocket(copy.deepcopy(gene))
		for n in range(len(self.__fatherNode)):
			gene = Gene()
			gene.randomCell(len(self.__ligand.branchSegment), self.__searchSpace)
			gene = self.calculates(gene)
			self.__fatherNode[n].addToPocket(copy.deepcopy(gene))
		for n in range(len(self.__leafNode)):
			gene = Gene()
			gene.randomCell(len(self.__ligand.branchSegment), self.__searchSpace)
			gene = self.calculates(gene)
			self.__leafNode[n].addToPocket(copy.deepcopy(gene))

	def mutation(self, cell):
		if random.randint(0,1) <= self.__mutProbability:
			select = random.randint(1, 7)
			if select == 1:
				cell.x = random.uniform(-self.__searchSpace, self.__searchSpace)
			elif select == 2:
				cell.y = random.uniform(-self.__searchSpace, self.__searchSpace)
			elif select == 3:
				cell.z = random.uniform(-self.__searchSpace, self.__searchSpace)
			elif select == 4:
				cell.sph_theta = random.uniform(0,2)*math.pi
			elif select == 5:
				cell.sph_phi = random.uniform(0,1)*math.pi
			elif select == 6:
				cell.theta = random.uniform(0,2)*math.pi
			elif select == 7:
				pos = random.randint(0, len(self.__ligand.branch)-1)
				cell.rotateBonds[pos] = random.uniform(0,2)*math.pi
		return cell



























