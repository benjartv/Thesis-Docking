import random
import math
import numpy as np

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
        self.rotateBonds = []
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
    		self.rotateBonds.append(angle)

    def randomCellKB(self, rotateBond, searchSpace, kbase, prob):
        self.x = random.uniform(-searchSpace, searchSpace)
        self.y = random.uniform(-searchSpace, searchSpace)
        self.z = random.uniform(-searchSpace, searchSpace)
        self.sph_theta = random.uniform(0,2)*math.pi
        self.sph_phi = random.uniform(0,1)*math.pi
        self.theta = random.uniform(0,2)*math.pi
        for i in range(rotateBond):
            if random.uniform(0,1) <= prob:
                ang = random.choice(random.choice(kbase[i]))
                angle = np.radians(random.uniform(ang-1,ang+1))
            else:
                angle = random.uniform(0,2)*math.pi
            self.rotateBonds.append(angle)

