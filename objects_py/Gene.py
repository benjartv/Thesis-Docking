import random
import math

# --------------------------- Gene Object ---------------------------
class Gene:

    def __init__(self):
        self.id = 0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.sph_theta = 0.0
        self.sph_phi = 0.0
        self.theta = 0.0
        self.rotateBranch = []
        self.score = 0.0

    def randomCell(self, rotateBond, searchSpace):
    	self.x = random.uniform(-searchSpace, searchSpace)
    	self.y = random.uniform(-searchSpace, searchSpace)
    	self.z = random.uniform(-searchSpace, searchSpace)
    	self.sph_theta = random.uniform(0,2)*math.pi
    	self.sph_phi = random.uniform(0,1)*math.pi
    	self.theta = random.uniform(0,2)*math.pi
    	for i in range(rotateBond):
    		angle = random.uniform(0,2)*math.pi
    		self.rotateBranch.append(angle)
