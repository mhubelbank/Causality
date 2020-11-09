from debug import *
import numpy as np
import getData
import math
import FlexDAG
import entropy_estimators as ee
import CaIC
import IndHSIC
#import CaPC
import statsmodels.nonparametric.kernel_regression as kr
from scipy.stats.distributions import norm
from scipy.linalg import lstsq

# Pruning Tuning Parameters:
# Eliminate any parents with relative total effect less than PRUNING_THRESHOLD
PRUNING_THRESHOLD = .03
# Void parent groups that include any parent with less than PRUNING_THRESHOLD2 coefficient
PRUNING_THRESHOLD2 = .001
# Penalize larger number of parents by multiplying the score by PARENT_COUNT_FACTOR**N;  where N is the number of parents
PARENT_COUNT_FACTOR = 1.001
# Maximum parent count to consider.  This primarily affects performance with large numbers of variables
MAX_PARENTAGE = 5

# Direction Test Method:  ('PL'--Pairwise Lingam, 'MI' -- Mutual Information, 'HSIC' -- Hilbert Space Independence Test)
DIR_TEST = 'PL'

class IC(CaIC.IC):
	# See above (DIR_TEST) for valid values for method
	# Set adjustData True to perform adjustments to the data series to make it stationary or to linearize relationships.
	def __init__(self, dataFileName, limit=None, prune=True, method='PL', adjustData=True):
		global DIR_TEST
		DIR_TEST = method
		self.entropies = {}
		super().__init__(dataFileName, limit=limit, prune=prune, adjustData=adjustData)
		self.doPrune = prune
		self.entropies = {}
		#self.computeEntropies()
		return
	
	def computeEntropies(self):
		eList = []
		eList2 = []
		for var in self.data.getSeriesNames():
			d = self.data.getSeries(var)
			#print('len = ', len(d))
			#print('p1')
			p = self.mi_prepare(d)
			#print('p2')
			entropy = ee.entropy(p, base=math.e)
			self.entropies[var] = entropy
			dentropy = self.dentropy(var)
			#print("Entropy, dentropy of ", var, "=", entropy, dentropy)
			eList.append((dentropy, var))
			eList2.append((entropy, var))
		eList.sort()
		#print('eList = ', eList)
		eList2.sort()
		#print('eList2 = ', eList2)
		return
		
	def getEntropy(self, varName):
		return self.entropies[varName]
		
	def getDAG(self):
		dag = FlexDAG.FlexDAG()
		# Order the vars in causal order
		oVars = self.getVarOrder()
		#print('Causal Order is: ', oVars)
		# Build the DAG based on the order with all possible edges
		for v in oVars:
			dag.addNode(v)
		for i in range(len(oVars)):
			v1 = oVars[i]
			for j in range(i+1, len(oVars)):
				v2 = oVars[j]
				#print(v, '<==>', v2, ' = ', self.scoreDependence(self.data.getSeries(v), self.data.getSeries(v2)))
				#if not self.isIndependent(v1, v2):
				dag.addEdge(v1, v2, 'D')
		Dprint('Raw Dag = ', str(dag))
		if self.doPrune:
			# Allow pruning to be turned off to increase perf when only order testing
			dag = self.prune(dag)
			Dprint('Pruned Dag = ', str(dag))
		return dag
	
	def prune_quick(self, dag):
		for node in dag.getVarOrder():
			#print('\nnode=', node)
			nodeD = np.array(self.getStandardSeries(node))
			nodeStd = nodeD.std()
			parents = dag.getParents(node)
			while len(parents) > 0:
				relStrengths = []
				coefs = []
				parentDs = []
				# Calculate coefs and relative strengths of affect for each parent
				for i in range(len(parents)):
					parentD = np.array(self.getStandardSeries(parents[i]))
					results = np.polyfit(parentD, nodeD.T, 1, full=True)
					fit, residuals = results[:2]
					coefs.append(fit)
					# Calculate relative strength of the effect for this parent
					x1coef = fit[0]
					absCoef = math.fabs(x1coef)
					nodeStd = nodeD.std()
					parentStd = parentD.std()
					relStrength = parentStd * absCoef / nodeStd
					relStrengths.append(relStrength)
					parentDs.append(parentD)
				# Sort parents in order of relative magnitude of effect
				zipped = list(zip(relStrengths, parentDs, parents, coefs))
				zipped.sort()
				zipped.reverse()
				strongestParentData = zipped[0]
				relStrength, parentD, parent, coef = strongestParentData
				x1coef, x0coef = coef
				dep = self.scoreDependence(parentD, nodeD)
				Dprint('Prune: relStrength (', parent, ' --> ', node, ') = ', relStrength, ', coefs = ', coef, ', dep = ', dep)
				if relStrength > PRUNING_THRESHOLD:
					# This parent is significant
					# Regress this parent's data from node data and then repeat the loop to find the next parent
					ycDSlope = parentD * x1coef
					ycD = ycDSlope + x0coef
					resD = nodeD - ycD
					nodeD = resD
					parents.remove(parent)
					Dprint('parent alllowed')
				else:
					# This parent is not significant.  This means that none that follow will be.  We've got
					# all significant parents now in outParents.  End the loop by settings parents to empty.
					# Remove edges from this parent and any subsequent parents to this node
					for parent in parents:
						dag.removeEdges(parent, node)
					parents = []
		return dag

	def prune(self, dag):
		deletes = []
		lowestError = 10**20
		lowestCombo = None
		# Remove independent edges
		for node in dag.getVarOrder():
			parentD = []
			origParents = dag.getParents(node)
			parents = origParents
			nodeD = self.getStandardSeries(node)
			for parent in parents:
				parentD.append(self.getStandardSeries(parent))
			if len(parents) > 0:
				coefs = []
				rawError = 0
				tempNodeD = np.array(nodeD)
				# Calculate coefs for each parent
				for i in range(len(parents)):
					results = np.polyfit(np.array(parentD[i]), tempNodeD, 1, full=True)
					fit, residuals = results[:2]
					coefs.append(fit)
					#print('Original coef for ', parents[i], ' --> ', node, ' = ', fit[0])
					#print('Original resid = ', residuals[0])
				# Sort parents in order of magnitude of coefs
				relStrengths = []
				outParents = []
				for i in range(len(parents)):
					coef = coefs[i][0]
					absCoef = math.fabs(coef)
					# Calculate relative strength of the affect
					nodeStd = tempNodeD.std()
					parentStd = np.array(self.getStandardSeries(parents[i])).std()
					relStrength = parentStd * absCoef / nodeStd
					#print('relStrength (', parents[i], ' --> ', node, ') = ', relStrength)
					if relStrength >= PRUNING_THRESHOLD:
					#if not self.isIndependent(node, parents[i]) and relStrength >= PRUNING_THRESHOLD:
						outParents.append(parents[i])
						relStrengths.append(relStrength)
						#Dprint('coef = ', absCoef)
					else:
						Dprint('CaLingam.prune: Deleting edge(1)', parents[i], ' --> ', node, '.  Coef = ', coef)
						deletes.append((parents[i], node))
				zipped = list(zip(relStrengths, outParents))
				zipped.sort()
				zipped.reverse()
				parents = [parent for (coef, parent) in zipped]
				#print('Significant Parents = ', parents)
				parentIter = self.expandParentList(parents)
				#print('parentIter = ', parentIter)
				lowestError = 10**20
				lowestCombo = parents
				for parentCombo in parentIter:
					#print('parentCombo = ', parentCombo)
					tempParents = list(parentCombo)
					tempParentD = []
					for tempParent in tempParents:
						tempParentD.append(self.getStandardSeries(tempParent))
					coefs = []
					rawError = 0
					tempNodeD = np.array(nodeD)
					for i in range(len(tempParents)):
						tempParentArray = np.array(tempParentD[i])
						results = np.polyfit(tempParentArray, tempNodeD.T, 1, full=True)
						fit, residuals = results[:2]
						#print('fit, residuals for ', tempParents[i], ' = ', fit, residuals)
						x1coef = fit[0]
						x0coef = fit[1]
						# If any parent has an X coef that is less than the threshold, then bail out.  It means that one of the parents is redundant.
						if x1coef < PRUNING_THRESHOLD2:
							# Don't consider this group
							rawError = 10**20 # Infinity
							break
						tempNodeD = tempNodeD - (tempParentArray * x1coef)
						coefs.append(fit)
						rawError = residuals[0] # Resids after regressing all parents.  The residuals are inherently cumulative, so we only care about the last one.
						#print('coef for ', tempParents[i], ' --> ', node, ' = ', fit[0])
						#print('resid = ', residuals[0])
					error = self.scaleError(rawError, len(tempParents), self.data.sampleCount)
					#error = rawError
					
					Dprint('Evaluating Edges = ', tempParents, '-->', node, '. Error = ', error)

					#print('error = ', error)
					#print('coefs = ', coefs)
					if error < lowestError:
						lowestError = error
						lowestCombo = tempParents
						Dprint('Best combo so far')
						
					# if (error2 - error)/error < .01:
						# for parent in parentCombo:
							# deletes.append((node, parent))
							# print('**** Print Deleting Edge = ', parent, '-->', node)
						# parents = tempParents
						# parentD = tempParentD
						# coefs = coefs2
						# error = error2
						# parentIter = self.expandParentList(parents)
						# print('parentIter = ', parentIter)
			#print('Best fit is: ', lowestCombo)
			# Prune out the parents not in this combo
			for parent in parents:
				if parent not in lowestCombo:
					deletes.append((parent, node))
					Dprint('CaLingam.prune: Deleting Edge(2) = ', parent, '-->', node)
			# if len(parents) == 1:
				# pNode = parents[0]
				# print('Evaluating Edge = ', pNode, '-->', node, coefs[0])

				# if math.fabs(coefs[0]) <= .2:
					#Special case for the last parent
					# deletes.append((node, pNode))
					# print('**** Print Deleting Edge = ', pNode, '-->', node)

					
		for delete in deletes:
			node, child = delete
			dag.removeEdges(node, child)
		return dag		
		
	def parameterizeLinearModel(self, dag):
		# First set node values
		for node in dag.getNodes():
			nodeD = np.array(self.data.getSeries(node))
			mean = nodeD.mean()
			#print('mean = ', mean)
			std = nodeD.std(ddof=1)
			dag.setNodeValue(node, (mean,std))
			# Set the weights for each parent edge
			parentEdges = dag.getParentEdges(node)
			for edge in parentEdges:
				n1, n2, eT, eW = edge
				parent = n1
				parentD = np.array(self.data.getSeries(parent))
				results = np.polyfit(parentD, nodeD, 1, full=True)
				fit = results[0]
				x1coef, x0coef = fit
				dag.setEdgeWeight(parent, node, fit)
				nodeD = nodeD - (parentD * x1coef + x0coef)
			dag.nodeD[node] = nodeD
		model = dag
		return model

				
			
	def getStandardSeries(self, seriesName):
		#return self.standardize(self.data.getSeries(seriesName))
		return self.getDataSeries(seriesName)
		
	def scaleError(self, raw, parentCount, dataCount):
		error = raw * PARENT_COUNT_FACTOR**(parentCount-1)
		#error = raw / (parentCount * dataCount) / (1+(parentCount * PARENT_COUNT_FACTOR))
		return error
		
	def expandParentList(self, parents):
		import itertools
		out = []
		for i in range(1,min([len(parents),MAX_PARENTAGE])):
			combinations = itertools.combinations(parents, i)
			out.extend(combinations)
		out.append(tuple(parents))
		return out
		
	# Return vars in causal order
	def getVarOrder(self):
		oVars = []
		R={}
		D = {}
		vNames = self.data.getSeriesNames()
		for vName in vNames:
			#D[vName] = self.standardize(self.data.getSeries(vName))
			D[vName] = self.data.getSeries(vName)
		while len(oVars) < len(vNames) - 1:
			dirM = {} # Direction matrix
			R = {}
			miCum = {}
			mostIndependent = None
			miMin = 10.0**100
			for j in range(len(vNames)):
				yD = []
				yNames = []
				xName = vNames[j]
				if xName in oVars:
					continue
				xD = D[xName]
				for vName in vNames:
					if vName == xName:
						continue
					if vName in oVars:
						continue
					d = D[vName]
					yNames.append(vName)
					yD.append(d)
				R[xName] = {}
				if xName not in miCum:
					miCum[xName] = 0
				for i in range(0, len(yNames)):
					yName = yNames[i]
					#if yName not in R:
					#	R[yName] = {}
					ayD = yD[i]
					axD = xD
					#depXY = self.scoreDependence(axD, ayD)
					#Dprint('depXY = ', depXY)
					dirKey = xName + '|' + yName
					if dirKey not in dirM:
						direction = self.direction2(axD, ayD)
						dirM[dirKey] = direction
					else:
						direction = dirM[dirKey]
					
					Dprint('xName, yName, direction = ', xName, yName, direction)
					
					if direction > 0:
						# First direction is correct
						# Penalize the more dependent direction (Y)
						loserName = yName
						winnerName = xName
						loserScore = direction
						winnerScore = 0
					elif direction < 0 :
						loserName = xName
						winnerName = yName
						loserScore = -direction
						winnerScore = 0
					else:
						loserName = None # If equal (unlikely), don't penalize anyone, or we could create some bias
						winnerName = None
						print('Warning: Direction score is equal.  This is an unusual situation', winnerScore, loserScore)
					Dprint('Predecessor = ', winnerName)
					if loserName is not None:
						if loserName not in miCum:
							miCum[loserName]  = 0
						cum = miCum[loserName]
						cum += loserScore
						miCum[loserName] = cum
						if winnerName not in miCum:
							miCum[winnerName] = 0
						
						cum = miCum[winnerName]
						cum += winnerScore
						miCum[winnerName] = cum
					else:
						5/0
					R[xName][yName] = ayD # residV1
					if yName not in R:
						R[yName] = {}
					R[yName][xName] = axD # residV2
					Dprint('R = ', R.keys())

			lowestCorr = 10**20
			mostIndependent = None
			#print ('miCum = ', miCum)
			zipped = list(zip(list(miCum.values()), list(miCum.keys())))
			zipped.sort()
			scores = [key for (value, key) in zipped]
			Dprint('score order = ', scores)
			keys = list(miCum.keys())
			keys.sort()
			for k in keys:
				cum = miCum[k]
				Dprint('cum(', k, ')=', cum)
			for k in miCum.keys():
				cum = miCum[k]
				if cum < lowestCorr:
					lowestCorr = cum
					mostIndependent = k
			Dprint('mostIndependent2 = ', mostIndependent, lowestCorr)
			Dprint()
			oVars.append(mostIndependent)
			# Regress all variables on the mostIndependent
			varReg = R[mostIndependent]
			D = varReg 
			#Dprint('D = ', keys)
			
		for v in D.keys():
			if v != mostIndependent:
				oVars.append(v)
		Dprint('oVars = ', oVars)
		return oVars

	def getVarOrder_2_9(self):
#	def getVarOrder(self):
		oVars = []
		R={}
		D = {}
		vNames = self.data.getSeriesNames()
		for vName in vNames:
			#D[vName] = self.standardize(self.data.getSeries(vName))
			D[vName] = self.data.getSeries(vName)
		while len(oVars) < len(vNames) - 1:
			R = {}
			miCum = {}
			mostIndependent = None
			miMin = 10.0**100
			for j in range(len(vNames)):
				yD = []
				yNames = []
				xName = vNames[j]
				if xName in oVars:
					continue
				xD = D[xName]
				for vName in vNames:
					if vName == xName:
						continue
					if vName in oVars:
						continue
					d = D[vName]
					yNames.append(vName)
					yD.append(d)
				R[xName] = {}
				if xName not in miCum:
					miCum[xName] = 0
				for i in range(0, len(yNames)):
					yName = yNames[i]
					#if yName not in R:
					#	R[yName] = {}
					ayD = np.array(yD[i])
					axD = np.array(xD)
					depXY = self.scoreDependence(axD, ayD)
					Dprint('depXY = ', depXY)
					
					#coefs = lstsq(np.array([axD]).T, ayD)[0]
					coefs = np.polyfit(axD, ayD.T,1)
					x1coef = coefs[0]
					x0coef = coefs[1]
					Dprint('xname, yname, x1coef = ', xName, yName, x1coef)
					ycD = axD * x1coef + x0coef
					residV1 = ayD - ycD
					dep1 = self.scoreDependence(axD, residV1)
					
					#coefs2 = lstsq(np.array([ayD]).T, axD)[0]
					coefs2 = np.polyfit(ayD, axD.T,1)
					x1coef2 = coefs2[0]
					x0coef2 = coefs2[1]
					Dprint('xname2, yname2, x1coef2 = ', yName, xName, x1coef2)
					xcD2 = ayD * x1coef2 + x0coef
					residV2 = axD - xcD2
					dep2 = self.scoreDependence(ayD, residV2)
					Dprint(xName, '-->', yName, ': dep1 = ', dep1)
					Dprint(yName, '-->', xName, ': dep2 = ', dep2)
					direction = self.direction(xName, yName)
					Dprint('direction = ', direction)
					if dep1 < dep2:
						# First direction is correct
						# Penalize the more dependent direction (Y)
						loserName = yName
						winnerName = xName
						loserScore = math.fabs((dep1-dep2)/(dep1 + dep2)/2 ) + direction
						winnerScore = dep1
					elif dep2 < dep1:
						loserName = xName
						winnerName = yName
						loserScore = math.fabs((dep1-dep2)/(dep1 + dep2)/2) - direction
						winnerScore = dep2
					else:
						loserName = None # If equal (unlikely), don't penalize anyone, or we could create some bias
						winnerName = None
						5/0
					Dprint('Predecessor = ', winnerName, 'dep1, dep2 = ', dep1, dep2)
					if loserName is not None:
						if loserName not in miCum:
							miCum[loserName]  = 0
						cum = miCum[loserName]
						if loserScore > cum:
							cum = loserScore
							miCum[loserName] = cum
						if winnerName not in miCum:
							miCum[winnerName] = 0
						#cum = miCum[winnerName]
						#if winnerScore > cum:
						#	cum = winnerScore
						#	miCum[winnerName] = cum
					else:
						5/0
					R[xName][yName] = ayD # residV1
					if yName not in R:
						R[yName] = {}
					R[yName][xName] = axD # residV2
					Dprint('R = ', R.keys())

			lowestCorr = 10**20
			mostIndependent = None
			#print ('miCum = ', miCum)
			zipped = list(zip(list(miCum.values()), list(miCum.keys())))
			zipped.sort()
			scores = [key for (value, key) in zipped]
			Dprint('score order = ', scores)
			keys = list(miCum.keys())
			keys.sort()
			for k in keys:
				cum = miCum[k]
				Dprint('cum(', k, ')=', cum)
			for k in miCum.keys():
				cum = miCum[k]
				if cum < lowestCorr:
					lowestCorr = cum
					mostIndependent = k
			Dprint('mostIndependent2 = ', mostIndependent, lowestCorr)
			Dprint()
			oVars.append(mostIndependent)
			# Regress all variables on the mostIndependent
			varReg = R[mostIndependent]
			D = varReg 
			Dprint('D = ', D.keys())
		for v in D.keys():
			if v != mostIndependent:
				oVars.append(v)
		Dprint('oVars = ', oVars)
		return oVars

	# Return vars in causal order
	def getVarOrder_2_5(self):
		oVars = []
		R={}
		D = {}
		vNames = self.data.getSeriesNames()
		for vName in vNames:
			#D[vName] = self.standardize(self.data.getSeries(vName))
			D[vName] = self.data.getSeries(vName)
		while len(oVars) < len(vNames) - 1:
			R = {}
			miCum = {}
			mostIndependent = None
			miMin = 10.0**100
			for j in range(len(vNames)):
				yD = []
				yNames = []
				xName = vNames[j]
				if xName in oVars:
					continue
				xD = D[xName]
				for vName in vNames:
					if vName == xName:
						continue
					if vName in oVars:
						continue
					d = D[vName]
					yNames.append(vName)
					yD.append(d)
				R[xName] = {}
				miCum[xName] = 0
				for i in range(len(yNames)):
					yName = yNames[i]
					#if yName not in R:
					#	R[yName] = {}
					ayD = np.array(yD[i])
					axD = np.array(xD)
					depXY = self.scoreDependence(axD, ayD)
					#coefs = lstsq(np.array([axD]).T, ayD)[0]
					coefs = np.polyfit(axD, ayD.T,1)
					x1coef = coefs[0]
					x0coef = coefs[1]
					Dprint('xname, yname, x1coef = ', xName, yName, x1coef, x0coef)
					Dprint('depXY = ', depXY)
					ycD = axD * x1coef + x0coef
					residV = ayD - ycD

					dep = self.scoreDependence(axD, residV)
					cum = miCum[xName]
					if dep > cum:
						cum = dep
						miCum[xName] = cum
					#R[xName][yName] = ayD
					R[xName][yName] = ayD
					Dprint(xName, '-->', yName, ': dep = ', dep)
					#print()
			lowestCorr = 10**20
			mostIndependent = None
			#Dprint ('miCum = ', miCum)
			keys = list(miCum.keys())
			keys.sort()
			for k in keys:
				cum = miCum[k]
				Dprint('cum(', k, ')=', cum)
			for k in miCum.keys():
				cum = miCum[k]
				if cum < lowestCorr:
					lowestCorr = cum
					mostIndependent = k
			Dprint('mostIndependent2 = ', mostIndependent, lowestCorr)
			Dprint()
			oVars.append(mostIndependent)
			# Regress all variables on the mostIndependent
			varReg = R[mostIndependent]
			D = varReg 
			#print('D = ', len(D), D.keys())
		for v in D.keys():
			if v != mostIndependent:
				oVars.append(v)
		Dprint('oVars = ', oVars)
		return oVars

		
	def isIndependent(self, v1, v2):
		d1 = self.getDataSeries(v1)
		d2 = self.getDataSeries(v2)
		result = self.isIndependent2(d1, d2)		
		Dprint('isIndependent ', v1, '<==>', v2, ' = ', result)
		return result
		
	def isIndependent2(self, d1, d2):
		import IndGauss
		result = IndGauss.isIndependent(d1,d2)
		return result
	
	def scoreDependencexxx(self, d1, d2):
		#d = norm.cdf(self.direction2(d1, d2) * 1000.0, scale=1)
		
		#d = self.direction2(d1, d2)
		#print ('d =', d)
		d1s = self.standardize(d1)
		d2s = self.standardize(d2)
		
		#score = self.scoreDependence4(d1, d2)
		score = self.scoreDependence1(d1s, d2s)
		
		#score += 10
		#if score < 0:
		#	score = 0
		#score = score -min([-d, 0])**2
		#print('master score = ', score)
		return score
		
	def scoreDependence3(self, d1, d2):
		#self.genCorrMatrix()
		d1 = self.standardize(d1)
		d2 = self.standardize(d2)
		arrayData = [d1, d2]
		a = np.array(arrayData)
		corrMatrix = np.corrcoef(a)
		#print ('corrM = ', corrMatrix)
		cc = corrMatrix[1,0]
		n = len(d1)

		corr = math.fabs(cc * ((n-2) / (1-cc**2))**.5)
		return corr	
	
	def getVarOrderOrg(self):
		oVars = []
		R={}
		D = {}
		vNames = self.data.getSeriesNames()
		for vName in vNames:
			#D[vName] = self.standardize(self.data.getSeries(vName))
			D[vName] = self.data.getSeries(vName)
		while len(oVars) < len(vNames) - 1:
			R = {}
			miCum = {}
			mostIndependent = None
			miMin = 10.0**20
			for j in range(len(vNames)):
				yD = []
				yNames = []
				xName = vNames[j]
				if xName in oVars:
					continue
				xD = D[xName]
				for vName in vNames:
					if vName == xName:
						continue
					if vName in oVars:
						continue
					d = D[vName]
					yNames.append(vName)
					yD.append(d)
				#coefs, resids = np.polyfit(xD, np.array(yD).T, 1, full=True)[:2]
				#print('coefs = ', coefs)
				#xcoefs = coefs[0]
				#kcoefs = coefs[1]
				R[xName] = {}
				miCum[xName] = 0
				for i in range(len(yNames)):
					yName = yNames[i]
					ayD = np.array(yD[i])
					axD = np.array(xD)
					#coefs, resid = np.polyfit(yD[i], xD, 1, full=True)[:2]
					xcoef = lstsq(np.array([yD[i]]).T, xD)[0][0]
					#txcoef = coefs[0]
					#print('coefs = ', coefs)
					#xcoef =	xcoefs[i]
					#kcoef = kcoefs[i]
					#resid = resids[i]
					#print('xcoef = ', xcoef,1/xcoef )
					
					# if not yName in R:
						# R[yName] = {}
					# if not yName in miCum:
						# miCum[yName] = 0
					ycD = axD * xcoef
					#print('ycD = ', len(ycD))
					residV = ayD - ycD
					#temp
					# testxcoef = 10**20
					# best = 10**20
					# bestR = None
					# count = 0
					# while math.fabs(testxcoef) > .01 and count < 5:
						# coefs2, resid2 = np.polyfit(residV, xD, 1, full=True)[:2]
						# testxcoef = coefs2[0]
						# axD2 = np.array(xD)
						# ayD2 = np.array(yD)
						# ycD2 = axD2 * testxcoef
						# r = ayD2 - ycD2
						# print ('testxcoef = ', testxcoef)
						# count += 1
						# if math.fabs(testxcoef) < best:
							# residV = r
					
					#print ('residV = ', residV.sum(), residV.mean(), residV.std())
					#print ('xName,yName, xcoef= ', xName, yName,  xcoef)
					#print('xD, residV, yD[i] = ', xD, residV, yD[i])
					corr = self.scoreDependence(xD, residV)
					cum = miCum[xName]
					if corr > cum:
						cum = corr
					miCum[xName] = cum
					R[xName][yName] = residV
					#print(xName, '-->', yName, ': corr = ', corr)
					#print()
			lowestCorr = 10**20
			mostIndependent = None
			#print ('miCum = ', miCum)
			keys = list(miCum.keys())
			keys.sort()
			for k in keys:
				cum = miCum[k]
				#print('cum(', k, ')=', cum)
			for k in miCum.keys():
				cum = miCum[k]
				if cum < lowestCorr:
					lowestCorr = cum
					mostIndependent = k
			#print('mostIndependent2 = ', mostIndependent, lowestCorr)
			#print()	
			oVars.append(mostIndependent)
			# Regress all variables on the mostIndependent
			varReg = R[mostIndependent]
			D = varReg 
			#print('D = ', len(D), D.keys())
		for v in D.keys():
			if v != mostIndependent:
				oVars.append(v)
		#print('oVars = ', oVars)
		return oVars

	def getVarOrder_org(self):
		oVars = self.data.getSeriesNames()
		scores = []
		for i in range(len(oVars)):
			score = 0
			var = oVars[i]
			d1 = self.data.getSeries(var)
			for j in range(len(oVars)):
				if i == j:
					continue
				d2 = self.data.getSeries(oVars[j])
				direction = self.direction(var, oVars[j])
				if direction < 0:
					score += 1
			scores.append((score, var))
				# nw = kr.KernelReg([d1], [d2], 'c', reg_type='ll', bw=[.002])
				# means, marginals = nw.fit()
				# r2 = nw.r_squared()
				# print('r2 = ', r2)
				# print ('means, marginals = ', means, marginals, var, oVars[j])
				# print ('mean = ', means.mean(), var, oVars[j])
				# break
			#break
		scores.sort()
		#print('scoreArray = ', scores)
		return oVars
	
	
	
	# Calculate causal direction of between two vars using PairwiseLiNGAM tanh method
	# Returns: direction := > 0 if Var1 causes Var2 or < 0 if Var2 causes Var1.  
	#                      Undefined if Var1 and Var2 are not causally related.
	
	def direction(self, var1, var2):
		s1 = self.getDataSeries(var1)
		s2 = self.getDataSeries(var2)
		result = self.direction2(s1, s2)
		#print ('R(', var1,'->', var2, ') =', result)
		return result
	
	def direction2(self, s1, s2):
		s1 = self.standardize(s1)
		s2 = self.standardize(s2)
		if DIR_TEST == 'PL':
			# Use pairwise direction determination
			dir = self.direction_PW(s1, s2)
		else:
			# Use the Darmois-Skitovich Theorem to determine pairwise direction based on a dependence test
			dir = self.direction_DS(s1, s2)
		return dir

	def direction_PW(self, s1, s2):
		# Pairwise Lingam Algorithm
		cum = 0
		for i in range(len(s1)):
			v1 = s1[i]
			v2 = s2[i]
			cumulant = v1*math.tanh(v2) - v2*math.tanh(v1)
			cum += cumulant
		avg = cum / float(len(s1))
		rho = self.getCorrCoef2(s1, s2)
		R = rho * avg
		Dprint('R = ', R)
		return R

	def direction_DS(self, s1, s2):
		# Darmois-Skitovich method for determining pairwise causal direction based on independence of residuals
		res1 = self.calcResidual(s1, s2)
		res2 = self.calcResidual(s2, s1)
		if DIR_TEST == 'MI':
			import IndMI
			indMod = IndMI
		elif DIR_TEST == 'HSIC':
			import IndHSIC
			indMod = IndHSIC
		dep1 = indMod.scoreDependence(s1, res1)
		dep2 = indMod.scoreDependence(s2, res2)
		# Direction is indicated by the order with the least dependence
		dir = dep2 - dep1
		return dir

	def calcResidual(self, s1, s2):
		coefs = np.polyfit(s1, s2.T, 1)
		x1coef, x0coef = coefs
		#Dprint('coefs = ', coefs)
		cs2 = s1 * x1coef + x0coef
		res = s2 - cs2
		return res

	def getCorrCoef2(self, s1, s2):
		cc = np.corrcoef([s1,s2])
		return cc[1,0]

	def mi(self, d1, d2):
		h1 = self.dentropy2(d1)
		h2 = self.dentropy2(d2)
		h3 = self.dentropy2(d1 + d2)
		result = h1 + h2 - h3
		#print('dmi = ', result)
		return result
		
	def dentropy(self, var):
		d = self.data.getSeries(var)
		return self.dentropy2(d)
		
	def dentropy2(self, d):
		d = self.standardize(d)
		t1 = (math.log(math.pi*2)+1) / 2
		t2a = 0
		t3a = 0
		for i in range(len(d)):
			t2a += math.log(math.cosh(d[i]))
			t3a += d[i]*math.e**((-d[i]**2)/2.0)
		t2 = 79.047*(t2a/float(len(d)) - .37457)**2
		t3 = 7.4129 * (t3a / float(len(d)))**2
		ent = t1 - t2 -t3
		#print ('t1, t2, t3 = ', t1, t2, t3, t3a)
		#print ('ent = ', ent)
		return ent
		
#inputFileName = '..\data\synthData1.csv'	
#cal = IC(inputFileName)
#print('dir(F,E) = ',cal.direction('F', 'E'))
#cal.plotData()
#print('dep = ', cal.scoreDependence(cal.data.getSeries('A'),cal.data.getSeries('E')))
#dag = cal.getDAG()
#print(dag)

# cal.direction('138.42.94.173|Fa3/1/0|bytesIn', '138.42.94.173|Gi0/0/0|bytesOut')

# cal.direction('A', 'E')
# cal.direction('E', 'F')
# cal.direction('G', 'F')
# cal.direction('A', 'I')
# cal.direction('E', 'D')
# cal.direction('F', 'E')
