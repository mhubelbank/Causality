import numpy as np
import math
import calibration

cal = calibration.Calibration('IndGauss')
DEPENDENCE_THRESHOLD = float(cal.get('DEPENDENCE_THRESHOLD', '.5'))

def isIndependent(X, Y):
	result = False # no ind
	# Test for conditional ind
	result = scoreDependence(X, Y) < DEPENDENCE_THRESHOLD
	return result

	
def scoreDependence(X, Y):
	corrCoefs = np.corrcoef([X,Y])
	cc = corrCoefs[1,0]
	n = len(X)
	corr = math.fabs(cc * ((n-2) / (1-cc**2))**.5)
	#print ('corr = ', corr, cc)
	return corr




	
			


