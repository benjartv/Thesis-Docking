import copy
import math
from objects_py.molecule import *
from utils_py.utils import *

class LocalSearch(object):
	"""docstring for localsearch"""
	def __init__(self, temp, tempMin, tempAlpha, ligand, searchSpace, centerSpace, typeLS, numIteration):
		self.__temp = temp
		self.__tempMin = tempMin
		self.__tempAlpha = tempAlpha
		self.__ligand = ligand
		self.__searchSpace = searchSpace
		self.__centerSpace = centerSpace
		self.__typeLS = typeLS
		self.__numIteration = numIteration
		self.__numberScoring = 0
		self.__temporalDir = "temp/"

	def changeLigand(self, ligand):
		self.__ligand = ligand
		#self.__centerSpace = center

	def initLocalSearch(self, cell):
		#print "Init Local Search..."
		T = self.__temp
		oldCell = copy.deepcopy(cell)
		while T > self.__tempMin:
			j = 1
			while j <= self.__numIteration:
				newCell = self.getNeighbor(oldCell)
				criteria = self.accptance(oldCell, newCell, T)
				if criteria > random.random():
					oldCell = copy.deepcopy(newCell)
				j += 1
			T *= self.__tempAlpha
		return copy.deepcopy(oldCell)

	def getNeighbor(self, cell):
		newCell = copy.deepcopy(cell)
		if self.__typeLS == 0:
			newCell = self.mutation(newCell)
		elif self.__typeLS == 1:
			newCell = self.mutationBlock(newCell)
		elif self.__typeLS == 2:
			newCell = self.mutationRot(newCell)
		newCell = self.calculates(newCell)
		return copy.deepcopy(newCell)

	def accptance(self, oldCell, newCell, temp):
		oldScore = oldCell.score
		newScore = newCell.score
		if newScore < oldScore:
			return 1
		else:
			delta = abs(oldScore - newScore)
			k = (1 / (1 + delta / temp))
			eq = math.exp(-delta/(temp*k))
			return eq

	def calculates(self, cell):
		auxLigand = copy.deepcopy(self.__ligand)
		auxLigand.translateToPoint([self.__centerSpace[0]+cell.x, 
									self.__centerSpace[1]+cell.y, 
									self.__centerSpace[2]+cell.z])
		sphVect = spherePoint(1, cell.sph_theta, cell.sph_phi)
		auxLigand.rotateByVector(sphVect, cell.theta)
		auxLigand.writePDBQT(self.__temporalDir+"ligand.pdbqt")
		cell.score = calculateFreeEnergy()
		self.__numberScoring += 1
		return copy.deepcopy(cell)

	def mutation(self, cell):
		select = random.randint(1, 6)
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
		return copy.deepcopy(cell)

	def mutationBlock(self, cell):
		select = random.randint(0,1)
		if select == 0:
			sel2 = random.randint(0,2)
			if sel2 == 0:
				cell.x = random.uniform(-self.__searchSpace, self.__searchSpace)
			elif sel2 == 1:
				cell.y = random.uniform(-self.__searchSpace, self.__searchSpace)
			elif sel2 == 2:
				cell.z = random.uniform(-self.__searchSpace, self.__searchSpace)
		elif select == 1:
			sel2 = random.randint(0,2)
			if sel2 == 0:
				cell.sph_theta = random.uniform(0,2)*math.pi
			elif sel2 == 1:
				cell.sph_phi = random.uniform(0,1)*math.pi
			elif sel2 == 2:
				cell.theta = random.uniform(0,2)*math.pi
		return cell

	def mutationRot(self, cell):
		select = random.randint(1,3)
		if select == 1:
			cell.sph_theta = random.uniform(0,2)*math.pi
		elif select == 2:
			cell.sph_phi = random.uniform(0,1)*math.pi
		elif select == 3:
			cell.theta = random.uniform(0,2)*math.pi
		return copy.deepcopy(cell)

	def getNumberEvaluation(self):
		return self.__numberScoring

		