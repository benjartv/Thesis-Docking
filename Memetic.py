from molecule import *
from Gene import *
from agent import *
import datetime
import copy

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
		self.__typeLS = params.typeLS
		self.__typeCO = params.typeCO
		self.__distanceCriteria = params.distanceCriteria
		self.__Time = datetime.datetime.now()
		self.__nodeByTree = params.nodeByTree
		self.__rootNode = None
		self.__fatherNode = [None]*self.__nodeByTree
		self.__leafNode = [None]*self.__nodeByTree*len(self.__fatherNode)
		self.__logPop = ""
		self.__numberScoring = 0
	
	def initProcess(self):
		startTime = datetime.datetime.now()
		self.initTree()
		self.initPopulation()
		for i in range(self.__generations):
			print "Generation: ", i, " Best: ", self.__rootNode.getBest().score
			self.generation()
			self.addlogPopulation(i)
			print "Root: ", self.__rootNode.getPocketScore()
		bestCell = self.__rootNode.getBest()
		print "Best ligand score: ", bestCell.score
		auxLigand = self.generateLigand(bestCell)
		print "RMSD: ", getRMSD(auxLigand, self.__originalLigand)
		'''
		self.generateFinalBest(bestCell)
		for i in range(self.__pocketSize):
			self.generateFinalBest(self.__rootNode.pocket[i], "best-lig-"+str(i)+".pdbqt")
		'''
		stopTime = datetime.datetime.now()
		self.__Time = stopTime - startTime
		print "Time: ", self.__Time
		print "Number of Energy Evaluation: ", self.__numberScoring
		self.writeLog()

	def initTree(self):
		self.__rootNode = agent(self.__pocketSize)
		for i in range(self.__nodeByTree):
			self.__fatherNode[i] = agent(self.__pocketSize)
		for i in range((self.__nodeByTree*len(self.__fatherNode))):
			self.__leafNode[i] = agent(self.__pocketSize)

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

	def calculates(self, cell):
		auxLigand = copy.deepcopy(self.__ligand)
		if self.__isLocalSearch:
			pass
		else:
			auxLigand.translateToPoint([self.__centerSpace[0]+cell.x, 
										self.__centerSpace[1]+cell.y, 
										self.__centerSpace[2]+cell.z])
			sphVect = spherePoint(1, cell.sph_theta, cell.sph_phi)
			auxLigand.rotateByVector(sphVect, cell.theta)
			for i in range(len(auxLigand.branch)):
				auxLigand.rotatateBranch(i, cell.rotateBranch[i])
			auxLigand.writePDBQT("ligand.pdbqt")
			cell.score = calculateFreeEnergy()
			self.__numberScoring += 1
		return copy.deepcopy(cell)

	def generateLigand(self, cell):
		auxLigand = copy.deepcopy(self.__ligand)
		if self.__isLocalSearch:
			pass
		else:
			auxLigand.translateToPoint([self.__centerSpace[0]+cell.x, 
										self.__centerSpace[1]+cell.y, 
										self.__centerSpace[2]+cell.z])
			sphVect = spherePoint(1, cell.sph_theta, cell.sph_phi)
			auxLigand.rotateByVector(sphVect, cell.theta)
			for i in range(len(auxLigand.branch)):
				auxLigand.rotatateBranch(i, cell.rotateBranch[i])
		return copy.deepcopy(auxLigand)

	def generateFinalBest(self, cell, name="best-ligand.pdbqt"):
		auxLigand = copy.deepcopy(self.__ligand)
		auxLigand.translateToPoint([self.__centerSpace[0]+cell.x, 
									self.__centerSpace[1]+cell.y, 
									self.__centerSpace[2]+cell.z])
		sphVect = spherePoint(1, cell.sph_theta, cell.sph_phi)
		auxLigand.rotateByVector(sphVect, cell.theta)
		for i in range(len(auxLigand.branch)):
			auxLigand.rotatateBranch(i, cell.rotateBranch[i])
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
		for i in range(len(selectedCell1.rotateBranch)):
			if random.randint(0,1)==1:
				newCell.rotateBranch.append(selectedCell1.rotateBranch[i])
			else:
				newCell.rotateBranch.append(selectedCell2.rotateBranch[i])
		return newCell

	def initLog(self):
		node = 1
		self.__logPop += "Molecule Ligand: "+self.__ligand.recordName+"\n"
		self.__logPop += "Pocket: "+str(self.__pocketSize)+"\n"
		self.__logPop += "Generations: "+str(self.__generations)+"\n"
		self.__logPop += "Rotate bonds: "+str(self.__ligand.branch)+"\n"
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
												self.__rootNode.pocket[i].rotateBranch])
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
													self.__fatherNode[j].pocket[i].rotateBranch])
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
													self.__leafNode[j].pocket[i].rotateBranch])
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
												self.__rootNode.pocket[i].rotateBranch])
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
													self.__fatherNode[j].pocket[i].rotateBranch])
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
													self.__leafNode[j].pocket[i].rotateBranch])
					self.__logPop += "score: "+str(self.__leafNode[j].pocket[i].score)+"\n"
			node += 1
			self.__logPop += "\n"

	def writeLog(self):
		file = open("iteration.log", "w")
		file.write(self.__logPop)
		file.close()

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
				cell.rotateBranch[pos] = random.uniform(0,2)*math.pi
		return cell



























