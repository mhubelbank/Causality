# Implements the independence test of Garcia and Gonzales-Lopez (2009)
# Test is based on longest increasing sequence and longest decreasing sequence
import numpy as np
import math
from scipy.stats.distributions import norm
import random
from debug import *

maxTries = 20
maxData = 1000000


import calibration

cal = calibration.Calibration('IndLn')
DEPENDENCE_THRESHOLD = float(cal.get('DEPENDENCE_THRESHOLD', '.5'))

def isIndependent(X, Y):
	result = False
	score = scoreDependence(X, Y)
	#print('score = ', score, len(X))
	if score < DEPENDENCE_THRESHOLD:
		result = True
	return result

def scoreDependence(X, Y):
	XY = np.array([X[:maxData],Y[:maxData]]).T
	# Sort the array by its X value and create an index Xi
	Xi = XY[:, 0].argsort()
	# Rearrange the sample pairs by Xi
	XY = XY[Xi]
	# Get a new index Xi
	Xi = XY[:, 0].argsort()
	# Replace X with Xi
	XY[:,0] = Xi
	# Get index for Y = Yi
	Yi = XY[:, 1].argsort()
	# Replace Y with Yi
	XY[:,1] = Yi
	# Find the longest increasing subsegment
	#  and the longest decreasing subsegment
	lis = 0
	lds = 0
	sequence = Yi

	# Longest I / D sequence
	lastyI = -1
	lastyD = 10**100
	for i in range(0, len(sequence)):
		y = sequence[i]
		if y > lastyI:
			lis += 1
			lastyI = y
		if y <	lastyD:
			lds += 1
			lastyD = y
	
	lcis = 0
	lcds = 0
	tempLcis = 0
	tempLcds = 0
	lastI = -1
	lastD = 10**100
	# Longest contiguous sequence
	for i in range(len(sequence)):
		y = sequence[i]
		if y >= lastI:
			tempLcis += 1
		else:
			#if tempLcis > lcis:
			lcis += tempLcis
			tempLcis = 0
		lastI = y
		if y <= lastD:
			tempLcds += 1
		else:
			#if tempLcds > lcds:
			lcds += tempLcds
			tempLcds = 0
		lastD = y	
	#if tempLcis > lcis:
	lcis += tempLcis
	#if tempLcds > lcds:
	lcds += tempLcds
	n = len(X)
	#basis = max([lis,lds])
	basis = max([lis, lds]) + math.fabs(lcis - lcds)
	#print('lis, lds, lcis, lcds = ', lis, lds, lcis, lcds, basis)
	Dprint('IndLn: basis = ', basis)
	#score = (basis - 2*n**.5) / n**(1/6.0)
	score = basis / (2*n)
	#print('IndLn: rawScore = ', score)
	#score += 10
	#score /= 8.0
	#score = math.tanh(score * 20)
	#print('IndLn: score = ', score)
	return score
	
def scoreDependence_orig(X, Y):
	XY = np.array([X[:maxData],Y[:maxData]]).T
	#print('org XY = ', XY)
	# Sort the array by its X value and create an index Xi
	Xi = XY[:, 0].argsort()
	#print('org Xi = ', Xi)
	# Rearrange the sample pairs by Xi
	#XY = np.take(XY, Xi, 0)
	XY = XY[Xi]
	#print('XY1 = ', XY)
	# Get a new index Xi
	Xi = XY[:, 0].argsort()
	#print('Xi1 = ', Xi)
	# Replace X with Xi
	XY[:,0] = Xi
	# Get index for Y = Yi
	Yi = XY[:, 1].argsort()
	#print('Yi = ', Yi)
	# Replace Y with Yi
	XY[:,1] = Yi
	# Find the longest increasing subsegment
	#  and the longest decreasing subsegment
	lis = 0
	lds = 0
	sequence = Yi
	# for j in range(maxTries):
		# lasty = -1
		# templis = 0
		# templds = 0
		# for i in range(j, len(sequence)):
			# if lis > templis + len(sequence) - i:
				# break
			# if lasty == len(sequence) - 1:
				# break
			# y = sequence[i]
			# if y > lasty:
				# templis += 1
				# lasty = y
		# lasty = 10**10
		# for i in range(j, len(sequence)):
			# if lds > templds + len(sequence) - i:
				# break
			# if lasty == len(sequence) -1:
				# break
			# y = sequence[i]
			# if y < lasty:
				# templds += 1
				# lasty = y
		# if templis > lis:
			# lis = templis
		# if templds > lds:
			# lds = templds
		# if (lis > len(sequence) - j) and (lds > len(sequence) - j):
			# break

	lastyI = -1
	lastyD = 10**100
	for i in range(0, len(sequence)):
		y = sequence[i]
		if y >= lastyI:
			lis += 1
			lastyI = y
		if y <=		lastyD:
			lds += 1
			lastyD = y
	
	n = len(X)
	#basis = max([lis,lds])
	basis = lis + lds
	Dprint('IndLn: basis = ', basis)
	score = (basis - 2*n**.5) / n**(1/6.0)
	Dprint('IndLn: rawScore = ', score)
	score += 18
	score /= 8.0
	score = math.tanh(score)
	Dprint('IndLn: score = ', score)
	return score
	
def pi(XY,X):
	return XY[X,1]
	
X = [5, 4,2,9]
Y = [.1, .2, .9, .3]
#print (relIndependence(X,Y))
#print ('***Entropy (standard normal) = ', norm.entropy())
	