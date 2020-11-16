import random
import math
import standardize
import numpy as np
import CaLingam
import CCM

TEST_POINTS = 3
L_VALUES = 33
L_COUNT = 3
EXCLUDE_LINEAR_NEIGHBORS = False

class IC(CaLingam.IC):
	def __init__(self, dataFileName, limit=None, prune=True):
		print('CaCCM')
		super().__init__(dataFileName, limit=limit, prune=prune)
		return
		
	def direction2(self, s1, s2):
		c = CCM.ccm(s1,s2)
		return c.direction()
		

		