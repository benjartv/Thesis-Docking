import random
import copy


class agent(object):
	def __init__(self, pocketSize):
		self.id = 0
		self.amountSolution = 0
		self.maxPocket = pocketSize
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

	def addToPocket(self, cell):
		if self.amountSolution < self.maxPocket:
			for p in range(self.maxPocket):
				if self.pocket[p] == None:
					if cell.score not in self.getPocketScore():
						self.pocket[p] = copy.deepcopy(cell)
						self.amountSolution += 1
						break
		else:
			self.pocket.sort(key=lambda x: x.score, reverse=False)
			for i in range(self.maxPocket):
				if cell.score == self.pocket[i].score:
					break
				elif cell.score < self.pocket[i].score:
					self.pocket[-1] = copy.deepcopy(cell)
					#self.pocket.sort(key=lambda x: x.score, reverse=False)
					break

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








