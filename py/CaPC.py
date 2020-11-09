from debug import *
import numpy as np
import getData
import math
import FlexDAG
import entropy_estimators as ee
import CaIC
import IndGauss
import CondIndGauss



indThreshold = 2
pearlCIThreshold2 = .01
miIndThreshold = .02
miCondIndThreshold = .002



# indMethod and condIndMethod have the following values:
# 	- 'Gauss' -- Use Gaussian methods based on correlation matrix
#   - 'MI' -- Use estimated Mutual Information
class IC(CaIC.IC):
	# def __init__(self, dataFile, limit=None, prune=True, indMethod='MI', condIndMethod='MI'):
	def __init__(self, dataFile, limit=None, prune=True, indMethod='Gauss', condIndMethod='Gauss'):
		self.indMethod = indMethod
		self.condIndMethod = condIndMethod
		super().__init__(dataFile, limit=limit, prune=prune)
		self.varIndexes = {}
		self.genCorrMatrix()
		self.entropies = {}
		self.ciCache = {}
		self.indCache = {}
		if self.indMethod == 'MI' or self.condIndMethod == 'MI':
			self.computeEntropies()
		return
	
	def computeEntropies(self):
		for var in self.data.getSeriesNames():
			d = self.data.getSeries(var)
			#print('p1')
			p = self.mi_prepare(d)
			#print('p2')
			entropy = ee.entropy(p)
			self.entropies[var] = entropy
			Dprint("Entropy of ", var, "=", entropy)
		return
	
	def getEntropy(self, varName):
		return self.entropies[varName]
	
	def temp(self):
		d1 = self.data.getSeries('A')
		d2 = self.data.getSeries('B')
		p1 = self.mi_prepare(d1)
		p2 = self.mi_prepare(d2)
		mi_ab = ee.mi(p1, p2, k=3)
		mi_aa = ee.mi(p1, p1, k=3)
		Dprint('mi ab = ', mi_ab)
		Dprint('mi aa = ', mi_aa)
		return
				
	def getDAG(self):
		dag = FlexDAG.FlexDAG()
		vars = self.data.getSeriesNames()
		for var in vars:
			dag.addNode(var)
		for i in range(len(vars)):
			for j in range(i+1, len(vars)):
				var1 = vars[i]
				var2 = vars[j]
			
				if self.isIndependent(var1, var2):
					continue
				if self.isAdjacent(var1, var2):
					if self.indMethod != 'Gauss':
						# Pairwise direction test only works for non-Gaussian.
						if self.direction(var1, var2) < 0:
							dag.addEdge(var2, var1, 'D')
						else:
							dag.addEdge(var1, var2, 'D')
					else:
						# For gaussian, add undirected edges for adjacent variables.  Then we'll come back and 
						# orient the edges once all adjacencies are found
						dag.addEdge(var1, var2, 'U')
		Dprint('Dag before orientation = ', dag)
		if self.indMethod == 'Gauss':
			dag = self.orientEdges(dag)
		if self.doPrune:
			dag = self.prune(dag)
		Dprint('Dag = ', dag)
		return dag
		
	def isAdjacent(self, var1, var2):
		adj = True
		ind1 = self.data.getIndexForSeries(var1)
		ind2 = self.data.getIndexForSeries(var2)
		Dprint('ind1, ind2 = ', ind1, ind2)
		vars = self.data.getSeriesNames()
		for i in range(len(vars)):
			varName = vars[i]
			indZ = self.data.getIndexForSeries(varName)
			Dprint('i = ', i)
			if (indZ != ind1) and (indZ != ind2):
				if self.isCondIndependent(var1, var2, varName):
					adj = False
					break
		if adj:
			Dprint(var1, 'is adjacent to', var2)
		else:
			Dprint(var1, 'not adjacent to', var2)
		return adj
	
	def isIndependent(self, var1, var2): # Returns true if var1 _||_ var2
		# Use cache if there
		keyTerms = [var1, var2]
		keyTerms.sort()
		key = tuple(keyTerms)
		if key in self.indCache:
			return self.indCache[key]
		result = True
		if self.indMethod == 'Gauss':
			result = self.isIndependentGauss(var1, var2)
		elif self.indMethod == 'MI':
			# Mutual Information
			result = self.isIndependentMI(var1, var2)
		self.indCache[key] = result
		return result
		
	def isIndependentGauss(self, var1, var2):
		X = self.data.getSeries(var1)
		Y = self.data.getSeries(var2)
		result = IndGauss.isIndependent(X,Y)
		return result

	def isIndependentMI(self, var1, var2):
		#print ('var1 = ', var2)
		d1 = self.data.getSeries(var1)
		d2 = self.data.getSeries(var2)
		cum = 0
		count = 0
		for i in range(0, len(d1), 100):
			td1 = d1[i:i+100]
			td2 = d2[i:i+100]
			tmi = self.isIndependentMI2(td1, td2)
			Dprint('tmi = ', tmi)
			cum += tmi
			count += 1
		mi = cum / float(count)
		Dprint('mi = ', mi)
		Dprint('org_mi = ', ee.mi(self.mi_prepare(d1), self.mi_prepare(d2)))
		# Normalize the MI using:  corr = I(X;Y)/(H(X)*H(Y))**.5
		h1 = self.getEntropy(var1)
		h2 = self.getEntropy(var2)
		corr = mi / ((h1*h2)**.5)
		if corr < miIndThreshold:
			Dprint (var1, 'is independent of', var2, corr)
			return True
		Dprint (var1, 'not independent of', var2, corr)
		return False
		
	def isIndependentMI2(self, d1, d2):
		p1 = self.mi_prepare(d1)
		p2 = self.mi_prepare(d2)
		#print ('d1 = ', d2)
		mi = ee.mi(p1, p2)
		return mi

	def isIndependentMI_Org(self, var1, var2):
		#print ('var1 = ', var2)
		d1 = self.data.getSeries(var1)
		d2 = self.data.getSeries(var2)
		p1 = self.mi_prepare(d1)
		p2 = self.mi_prepare(d2)
		#print ('d1 = ', d2)
		mi = ee.mi(p1, p2)
		# Normalize the MI using:  corr = I(X;Y)/(H(X)*H(Y))**.5
		h1 = self.getEntropy(var1)
		h2 = self.getEntropy(var2)
		corr = mi / ((h1*h2)**.5)
		if corr < miIndThreshold:
			Dprint (var1, 'is independent of', var2, corr)
			return True
		print (var1, 'not independent of', var2, corr)
		return False
	
		
	# Returns True iff var1 _||_ var2 | var3
	def isCondIndependent(self, var1, var2, var3):
		# check cache first
		keyTerms = [var1, var2]
		keyTerms.sort()
		keyTerms.append(var3)
		key = tuple(keyTerms)
		if key in self.ciCache:
			#print('returning ', key, ' from ciCache')
			return self.ciCache[key]
		result = False
		if self.condIndMethod == 'Gauss':
			result = self.isCondIndGauss(var1, var2, var3)
		elif self.condIndMethod == 'MI':
			# Conditional Mutual Information
			result =  self.isCondIndMI(var1, var2, var3)
		#print('Encaching ', key, 'in ciCache')
		self.ciCache[key] = result
		return result
		
	def isCondIndMI(self, var1, var2, var3):
		d1 = self.data.getSeries(var1)
		d2 = self.data.getSeries(var2)
		d3 = self.data.getSeries(var3)
		# d1 = self.standardize(d1)
		# d2 = self.standardize(d2)
		# d3 = self.standardize(d3)
		p1 = self.mi_prepare(d1)
		p2 = self.mi_prepare(d2)
		p3 = self.mi_prepare(d3)
		cmi = ee.cmi(p1,p2,p3,k=3)
		e1 = self.getEntropy(var1)
		e2 = self.getEntropy(var2)
		e3 = self.getEntropy(var3)
		normCmi = cmi / (e1 * e2 * e3)**(1/3.0)
		result = normCmi < miCondIndThreshold
		Dprint ('CI: ',var1, '_||_', var2, '|', var3, '=', result, normCmi)
		return result
		

	def isCondIndGauss(self, x,y,z):
		X = self.data.getSeries(x)
		Y = self.data.getSeries(y)
		Z = self.data.getSeries(z)
		result = CondIndGauss.isIndependent(X,Y,Z)
		#result = False # no cond ind
		# Test for conditional ind
		#Pxy = self.getCorrCoef(x,y)
		#Pxz = self.getCorrCoef(x,z)
		#Pyz = self.getCorrCoef(y,z)
		#Pxy_z = (Pxy - Pxz * Pyz) / ((1-Pxz**2)**.5 * (1-Pyz**2)**.5)

		#fz = .5 * math.log((1+Pxy_z) / (1-Pxy_z)) # Fischer's Z(Pxy_z)
		#print ('Fischers Z = ', fz, ' Pxy, Pxz, Pyz = ',Pxy, Pxz, Pyz)
		#result = math.fabs(fz) < pearlCIThreshold2
		#result = fz < pearlCIThreshold2
		Dprint ('Gaussian Method: ',x, '_||_', y, '|', z, '=', result)
		return result
		
	# Calculate causal direction of between two vars using PairwiseLiNGAM tanh method
	# Returns: direction := > 0 if Var1 causes Var2 or < 0 if Var2 causes Var1.  
	#                      Undefined if Var1 and Var2 are not causally related.
	
	def direction(self, var1, var2):
		s1 = self.data.getSeries(var1)
		s2 = self.data.getSeries(var2)
		s1 = self.standardize(s1)
		s2 = self.standardize(s2)
		cum = 0
		cum1 = 0
		cum2 = 0
	
		for i in range(len(s1)):
			v1 = s1[i]
			v2 = s2[i]
			cumulant = v1*math.tanh(v2) - v2*math.tanh(v1)
			cumulant1 = (v1**3)*v2
			cumulant2 = v1*(v2**3)
			cum1 += cumulant1
			cum2 += cumulant2
			cum += cumulant
		avg1 = cum1 / float(len(s1))
		avg2 = cum2 / float(len(s1))
		avg = cum / float(len(s1))
		rho = self.getCorrCoef(var1, var2)
		#print ('rho = ', rho)
		R = rho * avg
		#R = avg1 - avg2
		Dprint ('R(',var1, '->', var2, ') =', R)
		return R

	def getCorrCoef2(self, series1, series2):
		arrayData = [series1, series2]
		a = np.array(arrayData)
		corr = np.corrcoef(a)
		rho = corr[1,0]
		return rho
		
	def prune(self, dag):
		return dag

	def orientEdges(self, dag):	
		Dprint('CaPC.orientEdges')
		# First find any potential v-structures, and orient them based on independence
		vertexNodes = []
		for node in dag.getNodes():
			Dprint('node = ', node)
			related = dag.getEdgesForNode(node, 'A')
			Dprint('Related = ', related)
			if len(related) > 1:
				for i in range(len(related)):
					n1, n2, eT, eW = related[i]
					if n1 == node:
						a = n2
					else:
						a = n1
					for j in range(i+1, len(related)):
						n1, n2, eT, eW = related[j]
						if n1 == node:
							b = n2
						else:
							b = n1
						Dprint('Testing cond ind for ', a, b, node)
						if not self.isCondIndependent(a,b, node):
							# It must be a v-structure a-> node <- b
							# Orient the arrows accordingly.  Note that addEdge removes previous edges between the nodes, and is directed
							# by default
							dag.addEdge(a, node)
							Dprint('Added Directed Edge: ', a, node)
							dag.addEdge(b, node)
							Dprint('Added Directed Edge: ', b, node)

							vertexNodes.append(node)
		# Now orient any undirected edges such that 1) No new v-structures are created and 2) no cycles are created
		nodes = dag.getNodes()
		for i in range(len(nodes)):
			addedVertexNodes = 0
			for node in nodes:
				related = dag.getEdgesForNode(node, 'A')
				for edge in related:
					n1, n2, eT, eW = edge
					if eT != 'U':
						continue
					if node == n1:
						other = n2
					else:
						other = n1
					if node in vertexNodes:
						# We are a vertex node.  This rel should go away from us.
						dag.addEdge(node, other)
						Dprint('Added Directed Edge: ', node, other)
						# Other is now a potential vertex node (it has at least on incoming arrow)
						vertexNodes.append(other)
						addedVertexNodes += 1
					elif other in vertexNodes:
						# Other is a vertex node, but we are not.  The rel should point toward us, and we now become a potential vertex node
						dag.addEdge(other, node)
						Dprint('Added Directed Edge: ', other, node)
						vertexNodes.append(node)
						addedVertexNodes += 1
			if addedVertexNodes == 0:
				# We didn't add anything this time around.  We're done.
				break				
		return dag
		
if __name__ == '__main__':
	inputFileName = '..\data\synthData4.csv'	
	cal = IC(inputFileName, indMethod='MI', condIndMethod='MI')
	#cal = IC(inputFileName, indMethod='Gauss', condIndMethod='Gauss')
	dag = cal.getDAG()
	print('DAG = ', dag)


# cal.direction('138.42.94.173|Fa3/1/0|bytesIn', '138.42.94.173|Gi0/0/0|bytesOut')

# cal.direction('A', 'B')
# cal.direction('A', 'C')
# cal.direction('B', 'C')
# cal.direction('C', 'D')
# cal.direction('D', 'C')
# cal.direction('A', 'D')
