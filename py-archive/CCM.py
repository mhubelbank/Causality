import random
import math
import standardize
import numpy as np

TEST_POINTS = 3
L_VALUES = 33
L_COUNT = 3
EXCLUDE_LINEAR_NEIGHBORS = False

class ccm:
	def __init__(self, s1, s2, tau=2, dim=4):
		self.s1 = standardize.standardize(s1)
		self.s2 = standardize.standardize(s2)
		self.tau = tau
		self.dim = dim
		self.pointsToTest = []
		self.sLen = len(s1) - (self.dim-1)*self.tau
		self.pointsToTest = []
		self._generateShadowManifolds()
		return
	
	def _generateShadowManifolds(self):
		self.sm1 = self._getShadowManifold(self.s1)
		self.sm2 = self._getShadowManifold(self.s2)
		return
	
	def _getShadowManifold(self, series):
		SS = []
		tau = self.tau
		dim = self.dim
		for i in range(dim):
			if i == dim-1:
				ss = series[i*tau:]
			else:
				ss = series[i*tau:-(dim-1-i)*tau]
			SS.append(ss)
		SM = []
		for t in range(self.sLen):
			point = []
			for d in range(dim):
				#print('t = ', t, ', d = ', d, ', SS = ', SS)
				point.append(SS[d][t])
			point = tuple(point)
			SM.append(point)
		return SM

	def direction(self):
		forwardSeries = []
		reverseSeries = []
		neighborCount = self.dim + 1
		lValues = range(50, self.sLen, int((self.sLen-50)/L_VALUES))
		lValues = list(lValues[:L_COUNT]) + list(lValues[-L_COUNT:])
		for i in range(len(lValues)):
			L = lValues[i]
			# Check for convergence in forward and reverse directions
			distanceF = 0.0
			distanceR = 0.0
			pointsToTest = range(0, L, max([int(L/TEST_POINTS/3), 1]))
			testedPoints = 0
			for t in pointsToTest:
				neighbors = self.getClosestNeighbors(self.sm2[:L], t, neighborCount)
				neighbors2 = self.getClosestNeighbors(self.sm1[:L], t, neighborCount)
				if len(neighbors) != neighborCount or len(neighbors2) != neighborCount:
					continue
				# Forward
				dist = self.distanceToNeighbors(self.sm1[:L], t, neighbors)
				distanceF += dist
				# Reverse
				dist = self.distanceToNeighbors(self.sm2[:L], t, neighbors2)
				distanceR += dist
				testedPoints += 1
				if testedPoints == TEST_POINTS:
					break
			if testedPoints < TEST_POINTS:
				#print('skipping')
				continue
			forwardSeries.append(distanceF/len(pointsToTest))
			reverseSeries.append(distanceR/len(pointsToTest))
		fStren = self.convergenceStrength(forwardSeries)
		rStren = self.convergenceStrength(reverseSeries)
		print('fStren, rStren = ', fStren, rStren)
		if fStren > rStren:
			return 1
		return -1
		

	def getClosestNeighbors(self, manifold, t, count):
		tLoc = manifold[t]
		distances = []
		distSum = 0.0
		for t1 in range(len(manifold)):
			if t1 == t:
				continue
			loc = manifold[t1]
			dist = self.calcOneDist(tLoc, loc)
			distances.append((dist, t1))
		distances.sort()
		neighbors = []
		if EXCLUDE_LINEAR_NEIGHBORS:
			radius = 1
			for i in range(int(len(distances)/5)):
				d = distances[i]
				dist, t1 = d
				tDist = math.fabs(t-t1)
				if tDist > radius:
					neighbors.append(d)
					if len(neighbors) == count:
						break
				radius = tDist + 1
			if len(neighbors) < count:
				pass
				#print('***** Error.  Series too linear', len(neighbors))
		else:
			neighbors = distances[:count]
		return distances[:count]
		
	def distanceToNeighbors(self, manifold, t, neighbors):
		totalDist = 0.0 # Sum of weighted distances
		# Calculate Weights:
		weights = self.calcWeights(neighbors)
		#weights = [1] * len(neighbors)
		tPos = manifold[t]
		for i in range(len(neighbors)):
			nt = neighbors[i][1]
			weight = weights[i]
			pos = manifold[nt]
			dist = self.calcOneDist(tPos, pos)
			#print('dist = ', dist)
			weightedDist = dist * weight
			totalDist += weightedDist
		return totalDist
	
	def calcWeights(self, distances):
		rawWeights = []
		totalWeight = 0.0
		n1Dist = distances[0][0]
		for item in distances:
			iDist = item[0]
			rawWeight = math.e**(-iDist/n1Dist)
			totalWeight += rawWeight
			rawWeights.append(rawWeight)
		weights = []
		for rw in rawWeights:
			weight = rw / totalWeight
			weights.append(weight)
		return weights
		
	def calcOneDist(self, pos1, pos2):
		cumSquares = 0.0
		for i in range(self.dim):
			sq = (pos1[i] - pos2[i])**2
			cumSquares += sq
		dist = cumSquares**.5
		return dist
		
	def doesSeriesConverge(self, s):
		print('series = ', s)
		converges = False
		s = np.array(s)
		s1 = s[:int(len(s)*.6)]
		s2 = s[int(len(s)*.4):]
		print('s1.std, s2.std = ', s1.std(), s2.std(), s2.std()/s1.std())
		if s2.std() < .1 * s1.std():
			converges = True
		return converges

	def convergenceStrength(self, s):
		#print('s = ', s)
		s1 = np.array(s[:L_COUNT])
		s2 = np.array(s[-L_COUNT:])
		#print('s1, s2 = ', s1, s2)
		distance = s2.mean()
		convergence = s1.mean() / distance
		convergence = s1.std() / s2.std()
		if convergence > 1.2:
			pass
			print ('series are causally related -- convergence = ', convergence, ', dist = ', distance)
		else:
			distance = float('Inf')
			print ('series are independent -- convergence = ', convergence, ', dist = ', distance)
		#print('strength = ', strength, s1.std()/s1.mean(), s2.std()/s2.mean(), s3.std()/s3.mean(), 'distance = ', s3.mean())
		return 1/distance
			
		
		