# Implements the independence test using Mutual Information (MI)

import numpy as np
import math
from scipy.stats.distributions import norm
import random
from debug import *
import entropy_estimators as ee

maxTries = 20
maxData = 1000

import calibration

cal = calibration.Calibration('IndMI')
DEPENDENCE_THRESHOLD = float(cal.get('DEPENDENCE_THRESHOLD', '.2'))

def isIndependent(X, Y):
	result = False
	score = scoreDependence(X, Y)
	#print('score = ', score)
	if score < DEPENDENCE_THRESHOLD:
		result = True
	return result
	
def scoreDependence(X, Y):
	dep = ee.mi(ee.vectorize(X), ee.vectorize(Y))
	return dep