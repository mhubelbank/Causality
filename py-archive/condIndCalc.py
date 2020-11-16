import numpy as np
import math

inFn = '..\seriesTest4.csv'


labelTokens = []
labelDict = {}
corrCoefs = None
indThreshold2 = .05 # Threshold against corr-coef**2
pearlCIThreshold2 = 0.01
discreteCIThreshold = 0.02
method = 'pearl'
#method = 'discretize'

def prepare():
	global labelTokens, corrCoefs
	f = open(inFn, 'r')
	lines = []
	while 1:
		l = f.readline()
		if len(l) <= 0:
			break
		lines.append(l)
	labelLine = lines[0]
	labelTokens = labelLine[:-1].split(',')
	for token in labelTokens:
		labelDict[token] = []
	for line in lines[1:]:
		dataTokens = line[:-1].split(',')
		for indx in range(len(dataTokens)):
			label = labelTokens[indx]
			val = float(dataTokens[indx])
			series = labelDict[label]
			series.append(val)
	#print('d = ', labelDict)
	corrCoefs = buildCorrCoefs()
	f.close()

def findEdges():
	edges = []
	varCount = len(labelTokens)
	for i in range(varCount):
		for j in range(i, varCount):
			if j == i:
				continue
			print ('Testing ', labelTokens[i], labelTokens[j])
			if isInd(labelTokens[i], labelTokens[j]): # Unconditionally independent.  Can't be adjacent.
				print ('Unconditionally Independent')
				continue
			found = False
			for k in range(varCount):
				if k == i or k == j:
					continue
				if isCI(labelTokens[i], labelTokens[j], labelTokens[k]):
					found = True
					break
			if not found: # Found no mediating node, must be adjacent
				edges.append(labelTokens[i] + '--' + labelTokens[j])
	print ('edges = ', edges)
	return edges

def isInd(x,y): # Returns true if X _||_ Y
	cc = getCorrCoef(x,y)
	cc2 = cc * cc
	if cc2 < indThreshold2:
		#print ('Independent')
		return True
	#print ('Not Independent')
	return False

def isCI(x, y, z):  # Returns true if X _||_ Y | Z (X is conditionally independent of Y given Z)
	# Pearl method
	if method == 'pearl':
		ci = pearlMethod(x,y,z)
	elif method == 'discretize':
		ci = discretizeMethod(x,y,z)
	#ci = pearlMethod(x,y,z)
	return ci

def getCorrCoef(varName1, varName2):
	indx1 = -1
	indx2 = -1
	for i in range(len(labelTokens)):
		if labelTokens[i] == varName1:
			indx1 = i
		if labelTokens[i] == varName2:
			indx2 = i
		if indx1 >= 0 and indx2 >= 0:
			break
	corrCoef = corrCoefs[indx1, indx2]
	return corrCoef
	
def pearlMethod(x,y,z):
	result = False # no cond ind
	# Test for conditional ind
	Pxy = getCorrCoef(x,y)
	Pxz = getCorrCoef(x,z)
	Pyz = getCorrCoef(y,z)
	Pxy_z = (Pxy - Pxz * Pyz) / ((1-Pxz**2)**.5 * (1-Pyz**2)**.5)

	fz = .5 * math.log((1+Pxy_z) / (1-Pxy_z)) # Fischer's Z(Pxy_z)
	print ('Fischers Z = ', fz, ' Pxy, Pxz, Pyz = ',Pxy, Pxz, Pyz)
	result = math.fabs(fz) < pearlCIThreshold2
	print ('Pearl Method: ',x, '_||_', y, '|', z, '=', result)
	return result

def discretizeMethod(x,y,z):
	bins = []
	xSamples = labelDict[x]
	ySamples = labelDict[y]
	zSamples = labelDict[z]
	#binCount = int(len(zSamples) / 20)
	binCount = int(len(zSamples)**.5)
	#print ('binCount = ', binCount)
	# Find the min and max for Z
	zMin = 1000000000
	zMax = -1000000000
	for sample in zSamples:
		if sample < zMin:
			zMin = sample
		if sample > zMax:
			zMax = sample
	zRange = zMax - zMin
	binSize = zRange / binCount
	for i in range(binCount):
		bin = [[],[]]
		bins.append(bin)
	for j in range(len(zSamples)):
		zSample = zSamples[j]
		binNum = int((zSample - zMin) / binSize)
		#print ('binNum = ', binNum)
		if binNum >= binCount:
			binNum = binCount-1
		bin = bins[binNum]
		bin[0].append(xSamples[j])
		bin[1].append(ySamples[j])
	corrsum = 0
	countedSamples = 0
	for i in range(binCount):
		bin = bins[i]
		if len(bin[0]) < 20:
			#print ('skipping bin')
			continue
		da = np.array(bin)
		dcorr = np.corrcoef(da)
		corrsum += dcorr[0,1] * len(bin[0])
		countedSamples += len(bin[0])
		#print ('bin[' + str(i)		+ ']: ', dcorr[0,1])
	corravg = corrsum / countedSamples
	print('corravg = ', str(corravg), str(corrsum), str(countedSamples))
	fz = .5 * math.log((1+corravg) / (1-corravg)) # Fischer's Z(Pxy_z)
	result = math.fabs(fz) < discreteCIThreshold
	print ('Fischers Z = ', fz)
	print ('Discretize Method: ',x, '_||_', y, '|', z, '=', result)
	#print ('bins = ', bins)
	return result
	
		
		
def buildCorrCoefs():
	columns = []
	arrayData = []

	for i in range(len(labelTokens)):
		label = labelTokens[i]
		series = labelDict[label]
		arrayData.append(series)
		
	a = np.array(arrayData)

	print ('a=', a.shape)
	corr = np.corrcoef(a)
	print ('corr = ', corr.shape)
	return corr
				


	
			
prepare()
findEdges()
# isCI('A','B','C')
# isCI('A','D','C')
# isCI('B','D','C')

