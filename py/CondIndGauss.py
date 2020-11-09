import numpy as np
import math
import calibration

cal = calibration.Calibration('CondIndGauss')
DEPENDENCE_THRESHOLD = float(cal.get('DEPENDENCE_THRESHOLD', '.5'))

def isIndependent(X, Y, Z):
	result = False # no cond ind
	# Test for conditional ind
	result = scoreDependence(X, Y, Z) < DEPENDENCE_THRESHOLD
	#print ('CondIndGauss: isIndependent ',Z, '_||_', Y, '|', Z, '=', result)
	return result

	
def scoreDependence(X, Y, Z):
	# Test for conditional ind
	corrData = [X, Y, Z]
	corrCoefs = np.corrcoef(corrData)
	Pxy = corrCoefs[0,1]
	Pxz = corrCoefs[0,2]
	Pyz = corrCoefs[1,2]
	Pxy_z = (Pxy - Pxz * Pyz) / ((1-Pxz**2)**.5 * (1-Pyz**2)**.5)
	fz = .5 * math.log((1+Pxy_z) / (1-Pxy_z)) # Fischer's Z(Pxy_z)
	dep = math.fabs(fz)
	#print('dep = ', dep)
	return dep



	
			


