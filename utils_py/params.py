class params(object):
	"""docstring for params"""
	def __init__(self, searchSpace,
					centerSpace,
					generations,
					pocketSize,
					treeNodes,
					mutProbability,
					isLocalSearch,
					typeLS,
					typeCO,
					distanceCriteria,
					nodeByTree,
					tempLS,
					minTemp,
					alphaTemp,
					numberIteration):
		self.searchSpace = searchSpace
		self.centerSpace = centerSpace
		self.generations = generations
		self.pocketSize = pocketSize
		self.treeNodes = treeNodes
		self.mutProbability = mutProbability
		self.isLocalSearch = isLocalSearch
		self.typeLS = typeLS
		self.typeCO = typeCO
		self.distanceCriteria = distanceCriteria
		self.nodeByTree = nodeByTree
		self.tempLS = tempLS
		self.minTemp = minTemp
		self.alphaTemp = alphaTemp
		self.numberIteration = numberIteration
		