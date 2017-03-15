from molecule import *
import random
import copy
from utils_py.utils import *

class agent(object):
	def __init__(self, pocketSize, distanceCriteria, ligand, newCenter):
		self.id = 0
		self.amountSolution = 0
		self.maxPocket = pocketSize
		self.current = None
		self.pocket = [None]*self.maxPocket
		self.distanceCriteria = distanceCriteria
		self.ligand = ligand
		self.centerSpace = newCenter

	def resetAgent(self):
		self.amountSolution = 0
		self.current = None
		self.pocket = [None]*self.maxPocket
		
	def getRandom(self):
		realPocket = []
		for p in self.pocket:
			if p != None:
				realPocket.append(p)
		select = random.randint(0,len(realPocket)-1)
		#self.current = realPocket[select]
		return copy.deepcopy(realPocket[select])

	def addToCurrent(self, cell):
		self.current = copy.deepcopy(cell)

	def generateStructure(self, cell):
		auxLigand = copy.deepcopy(self.ligand)
		auxLigand.translateToPoint([self.centerSpace[0]+cell.x, 
									self.centerSpace[1]+cell.y, 
									self.centerSpace[2]+cell.z])
		sphVect = spherePoint(1, cell.sph_theta, cell.sph_phi)
		auxLigand.rotateByVector(sphVect, cell.theta)
		for i in range(len(auxLigand.branch)):
			auxLigand.rotateAtomsBranch(i, cell.rotateBonds[i])
		return copy.deepcopy(auxLigand)

	def calculateDistance(self, cell1, cell2):
		ligand1 = self.generateStructure(cell1)
		ligand2 = self.generateStructure(cell2)
		distance = getRMSD(ligand1, ligand2)
		return distance

	def addToPocket(self, cell):
		minDist = 9999.9
		idDist = 0
		distCount = 0
		if self.amountSolution < self.maxPocket:
			auxPocket = []
			auxScore = []
			noneId = 0
			for p in range(self.maxPocket):
				if self.pocket[p] != None:
					auxPocket.append(self.pocket[p])
					auxScore.append(self.pocket[p].score)
					new_distance = self.calculateDistance(self.pocket[p], cell)
					if new_distance < minDist:
						minDist = new_distance
						idDist = p
					if new_distance >= self.distanceCriteria:
						distCount += 1
				else:
					noneId = p
			auxPocket.sort(key=lambda x: x.score, reverse=False)

			if (cell.score not in auxScore):
				if distCount == len(auxPocket):
					self.pocket[noneId] = copy.deepcopy(cell)
					self.amountSolution += 1
				elif cell.score < auxPocket[0].score:
					self.pocket[idDist] = copy.deepcopy(cell)
		else:
			self.pocket.sort(key=lambda x: x.score, reverse=False)
			for p in range(self.maxPocket):
				new_distance = self.calculateDistance(self.pocket[p], cell)
				if new_distance < minDist:
					minDist = new_distance
					idDist = p
				if new_distance >= self.distanceCriteria:
					distCount += 1
			#replace min dist
			if cell.score < self.pocket[0].score:
				self.pocket[idDist] = copy.deepcopy(cell)
			elif (cell.score >= self.pocket[0].score) and (distCount == self.maxPocket) and (cell.score < self.pocket[-1].score):
				self.pocket[-1] = copy.deepcopy(cell)

	def getPocketScore(self):
		newPocket = []
		for p in self.pocket:
			if p != None:
				newPocket.append(p.score)
		newPocket.sort()
		return newPocket[:]

	def getBest(self):
		if self.amountSolution < self.maxPocket:
			auxPocket = []
			for p in self.pocket:
				if p != None:
					auxPocket.append(p)
			auxPocket.sort(key=lambda x: x.score, reverse=False)
			return copy.deepcopy(auxPocket[0])
		else:
			self.pocket.sort(key=lambda x: x.score, reverse=False)
			return copy.deepcopy(self.pocket[0])








