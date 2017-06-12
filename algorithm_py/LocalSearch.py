import copy
import math
from objects_py.molecule import *
from utils_py.utils import *

class LocalSearch(object):
	"""docstring for localsearch"""
	def __init__(self, temp, tempMin, tempAlpha, ligand, searchSpace, centerSpace, typeLS, numIteration, iskb, kbprob):
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
		self.__isKb = iskb
		self.__kbProb = kbprob
		self.__tempIterative = 1.0
		self.__alphaIt = None

		self.setAlphaIt()

	def changeLigand(self, ligand):
		self.__ligand = ligand
		#self.__centerSpace = center

	def initLocalSearch(self, cell):
		#print "Init Local Search..."
		self.__tempIterative = 1.0
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
			self.__tempIterative *= self.__alphaIt
		return copy.deepcopy(oldCell)

	def getNeighbor(self, cell):
		newCell = copy.deepcopy(cell)
		if self.__typeLS == 0:
			newCell = self.mutationBlock(newCell)
		elif self.__typeLS == 1:
			alpha = self.__tempIterative
			newCell = self.mutationReduce(newCell, alpha)
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
		if self.__isKb:
			for i in range(len(auxLigand.branch)):
				torAngle = auxLigand.rotateBranchKB(i, cell.rotateBonds[i])
				auxLigand.rotateAtomsBranch(i, torAngle)
		else:
			for i in range(len(auxLigand.branch)):
				auxLigand.rotateAtomsBranch(i, cell.rotateBonds[i])
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
			#cell.sph_phi = random.uniform(0,1)*math.pi
			cell.sph_phi = math.acos(2*random.uniform(0,1)-1)
		elif select == 6:
			cell.theta = random.uniform(0,2)*math.pi
		elif select == 7:
			pos = random.randint(0, len(self.__ligand.branch)-1)
			if self.__isKb:
				if random.uniform(0,1) <= self.__kbProb:
					ang = random.choice(self.__ligand.anglesArray[pos])
					cell.rotateBonds[pos] = np.radians(random.uniform(ang-1,ang+1))
				else:
					cell.rotateBonds[pos] = random.uniform(0,2)*math.pi
			else:
				cell.rotateBonds[pos] = random.uniform(0,2)*math.pi
		return copy.deepcopy(cell)

	def mutationBlock(self, cell):
		select = random.randint(0,2)
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
				#cell.sph_phi = random.uniform(0,1)*math.pi
				cell.sph_phi = math.acos(2*random.uniform(0,1)-1)
			elif sel2 == 2:
				cell.theta = random.uniform(0,2)*math.pi
		elif select == 2:
			pos = random.randint(0, len(self.__ligand.branch)-1)
			if self.__isKb:
				if random.uniform(0,1) <= self.__kbProb:
					ang = random.choice(self.__ligand.anglesArray[pos])
					cell.rotateBonds[pos] = np.radians(random.uniform(ang-1,ang+1))
				else:
					cell.rotateBonds[pos] = random.uniform(0,2)*math.pi
			else:
				cell.rotateBonds[pos] = random.uniform(0,2)*math.pi
		return cell


	def mutationReduce(self, cell, alpha):
		select = random.randint(0,2)
		if select == 0:
			sel2 = random.randint(0,2)
			if sel2 == 0:
				distright = self.__searchSpace - cell.x
				distleft = -self.__searchSpace + cell.x
				cell.x = random.uniform(distleft*alpha - cell.x, cell.x + distright*alpha)
			elif sel2 == 1:
				distright = self.__searchSpace - cell.y
				distleft = -self.__searchSpace + cell.y
				cell.y = random.uniform(distleft*alpha - cell.y, cell.y + distright*alpha)
			elif sel2 == 2:
				distright = self.__searchSpace - cell.z
				distleft = -self.__searchSpace + cell.z
				cell.z = random.uniform(distleft*alpha - cell.z, cell.z + distright*alpha)
		elif select == 1:
			sel2 = random.randint(0,2)
			if sel2 == 0:
				cell.sph_theta = random.uniform(0,2)*math.pi
			elif sel2 == 1:
				#cell.sph_phi = random.uniform(0,1)*math.pi
				cell.sph_phi = math.acos(2*random.uniform(0,1)-1)
			elif sel2 == 2:
				theta = cell.theta / math.pi
				if theta > 1:
					distleft = theta - 1
					distright = 1 + theta
				else:
					distleft = -1 + theta
					distright = 1 + theta
				angleFinal =  random.uniform(alpha*distleft,alpha*distright)
				if angleFinal > 2.0:
					angleFinal = 2.0 - angleFinal
				elif angleFinal < 0.0:
					angleFinal = 2.0 + angleFinal
				cell.theta = angleFinal * math.pi
		elif select == 2:
			pos = random.randint(0, len(self.__ligand.branch)-1)
			if self.__isKb:
				if random.uniform(0,1) <= self.__kbProb:
					ang = random.choice(self.__ligand.anglesArray[pos])
					cell.rotateBonds[pos] = np.radians(random.uniform(ang-1,ang+1))
				else:
					cell.rotateBonds[pos] = random.uniform(0,2)*math.pi
			else:
				angle = cell.rotateBonds[pos] / math.pi
				if angle > 1:
					distleft = angle - 1
					distright = 1 + angle
				else:
					distleft = -1 + angle
					distright = 1 + angle
				angleFinal =  random.uniform(alpha*distleft,alpha*distright)
				if angleFinal > 2.0:
					angleFinal = 2.0 - angleFinal
				elif angleFinal < 0.0:
					angleFinal = 2.0 + angleFinal
				cell.rotateBonds[pos] = angleFinal * math.pi
		return cell



	def mutationRot(self, cell):
		select = random.randint(1,4)
		if select == 1:
			cell.sph_theta = random.uniform(0,2)*math.pi
		elif select == 2:
			#cell.sph_phi = random.uniform(0,1)*math.pi
			cell.sph_phi = math.acos(2*random.uniform(0,1)-1)
		elif select == 3:
			cell.theta = random.uniform(0,2)*math.pi
		elif select == 4:
			pos = random.randint(0, len(self.__ligand.branch)-1)
			if self.__isKb:
				if random.uniform(0,1) <= self.__kbProb:
					ang = random.choice(self.__ligand.anglesArray[pos])
					cell.rotateBonds[pos] = np.radians(random.uniform(ang-1,ang+1))
				else:
					cell.rotateBonds[pos] = random.uniform(0,2)*math.pi
			else:
				cell.rotateBonds[pos] = random.uniform(0,2)*math.pi
		return copy.deepcopy(cell)

	def getNumberEvaluation(self):
		return self.__numberScoring

	def setAlphaIt(self):
		alpha = self.__tempAlpha
		InitTemp = self.__temp
		MinTemp = self.__tempMin
		count = 0
		while InitTemp > MinTemp:
			InitTemp *= alpha
			count += 1
		self.__alphaIt = 1 - (0.5 / count)














		