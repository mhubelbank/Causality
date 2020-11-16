from scipy.linalg import lstsq
import getData
import IndHSIC
import IndLn
import IndMI
import numpy as np
import synthDataGen
from debug import *
import math
import time


class SemAnalyzer:
	def __init__(self, filePath, dataCounts=None):
		dot_pos = filePath.rfind('.')
		basePath = filePath[:dot_pos]
		ext = filePath[dot_pos+1:]
		self.genFile = basePath + '.py'
		self.testFile = basePath + '.csv'

		try:
			f = open(basePath + '.py', 'r')
			exec(f.read(), globals())
			self.validations = validation
		except:
			self.validations = []
		if dataCounts is None:
			dataCounts = [10000, 5000, 2000, 1000, 100]
		self.dataCounts = dataCounts
		return
		
	def analyzeLingam(self):
		for dataCount in self.dataCounts:
			stats = self.analyzeOneLingam(dataCount)
			#print('DataCount = ', dataCount)
			#print('   stats = ', stats)
		return
	
	def analyzeOneLingam(self, dataCount, gen=True, valData=None):
		if gen:
			synthDataGen.run(self.genFile, dataCount)
		dr = getData.DataReader(self.testFile, dataCount)
		vars = dr.getSeriesNames()
		if valData is None:
			valData = synthDataGen.getValidation()
		lowestSNR = 10**100
		cumSNR = 0
		pairCount = 0
		lowestPair = None
		stats = {}
		nodeDiffs = {}
		for tuple in valData:
			successor, predecessors = tuple
			if len(predecessors) == 0:
				continue
			for predecessor in predecessors:
				nodeDiffs[predecessor] = {}
				pD = np.array(dr.getSeries(predecessor))
				sD = np.array(dr.getSeries(successor))
				coefs = np.polyfit(pD, sD.T, 1)
				x1coef, x0coef = coefs
				#Dprint('coefs = ', coefs)
				cs2D = pD * x1coef + x0coef
				res = sD - cs2D
				# Note that signal and noise are reversed from what you might expect.
				# In this case, the variance of the error term is the signal we are looking for,
				# while the linear component is actually the noise
				noise = cs2D.var(ddof=1)
				signal = res.var(ddof=1)
				#print('pair = ', predecessor, '--->', successor)
				#print('noise = ', noise)
				#print('signal = ', signal)
				snr = signal / noise
				#print('snr = ', snr)
				nodeDiffs[predecessor][successor] = 10* math.log(1.0/snr, 10)
				if snr < lowestSNR:
					lowestSNR = snr
					lowestPair = (predecessor, successor)
				cumSNR += snr
				pairCount += 1
		stats['minSnr'] = 10*math.log(lowestSNR,10)
		avgSNR = cumSNR/float(pairCount)
		stats['avgSnr'] = 10*math.log(avgSNR, 10)
		stats['weakestPair'] = lowestPair
		difficulty = max([10 * math.log(1.0 / lowestSNR, 10),0.0])
		stats['difficulty'] = difficulty
		stats['weakestPairing'] = lowestPair
		stats['variableDifficulty'] = nodeDiffs
		stats['normDifficulty'] = 100.0 * difficulty / (dataCount**.5)		
		return stats
		
	def getDependencies(self):
		dependencies = []
		independencies = []
		if len(self.validations) > 1:
			exogenousVars = []
			for v in self.validations:
				successor, predecessors = v
				if len(predecessors) == 0:
					# successor is exogenous
					exogenousVars.append(successor)
				else:
					for parent in predecessors:
						dependencies.append((successor, parent))
			for i in range(len(exogenousVars)):
				for j in range(i+1, len(exogenousVars)):
					independencies.append((exogenousVars[i], exogenousVars[j]))
		return (dependencies, independencies)
	
	def getCondDependencies(self):
		dependencies = []
		independencies = []
		if len(self.validations) > 1:
			nodeParents = {}
			nodeChildren = {}
			for v in self.validations:
				successor, predecessors = v
				nodeParents[successor] = predecessors
				for parent in predecessors:
					if parent not in nodeChildren:
						nodeChildren[parent] = []
					nodeChildren[parent].append(successor)
			for node in nodeChildren.keys():
				children = nodeChildren[node]
				# Two types of Cond Independencies:
				#	Type1: A -> B -> C ==>  A _||_ C | B, B not _||_ C | A, A not _||_ B | C
				#	Type2: A <- B -> C ==> A _||_ C | B, A not _||_ B | C, B not _||_ C | A
				#	Type3: A -> B <- C ==>  A not _||_ C | B, A not _||_ B | C, B not _||_ C | A
				# Start Type1.  
				# All children of this node should be independent of this node's parents given this node.
				# Children of this node should not be independent of this node, given its parents.
				# Parents of this node should not be independent of this node, given its children
				parents = nodeParents[node]
				for i in range(len(children)):
					for j in range(len(parents)):
						independencies.append((children[i], parents[j], node))
						dependencies.append((children[i], node, parents[j]))
						dependencies.append((parents[j], node, children[i]))
				# Note that we should technically include every combination of node's ancestors and descendants, but this should do for now.
				# Now Type2.  
				# All children of this node should be independent of one another given this node.
				# Children should not be independent of this node, given one another
				for i in range(len(children)):
					for j in range(i+1, len(children)):
						independencies.append((children[i], children[j], node))
						dependencies.append((children[i], node, children[j]))
				# Note that we should technically include every combination of child's descendants, but this should do for now.
				# Now Type3
				# Parents should not be independent of one another given this node
				# This node should not be independent of parent given the other parents
				for i in range(len(parents)):
					for j in range(i+1, len(parents)):
						dependencies.append((parents[i], parents[j], node ))
						dependencies.append((parents[i], node, parents[j]))
		return (dependencies, independencies)
				
