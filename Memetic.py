from molecule import *


class Memetic(object):
	def __init__(self, params, ligand):
		self.__searchSpace = params.searchSpace
		self.__generations = params.generationNumber
        self.__pocketSize = params.pocketSize
        self.__treeNodes = params.treeNodes
        self.__mutProbability = params.mutProbability
        self.__typeLS = params.typeLS
        self.__typeCO = params.typeCO
        self.__distanceCriteria = params.distanceCriteria

