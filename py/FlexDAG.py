import math
import numpy as np
from debug import *

EdgeTypes = ['D', 'U']
class FlexDAG():
	def __init__(self):
		self.nodes = []
		self.nodeValues = {}
		self.edges = []
		self.noises = {}
		self.nodeD = {}
		return
		
	def addNode(self, node):
		if not node in self.nodes:
			self.nodes.append(node)
		self.nodeValues[node] = (0,0)
		return
	
	def addEdge(self, node1, node2, edgeType = 'D', edgeWeight=(0,0)):
		self.removeEdges(node1, node2) # Remove any existing edges between the 2 nodes
		edge = (node1, node2, edgeType, edgeWeight)
		self.edges.append(edge)
		return
	
	# Removes any edges between the two nodes
	def removeEdges(self, node1, node2):
		deletes = []
		for i in range(len(self.edges)):
			n1, n2, eT, eW = self.edges[i]
			if ((n1 == node1) and (n2 == node2)) or ((n1 == node2) and (n2 == node1)):
				deletes.append(i)
		deletes.reverse()
		for i in deletes:
			self.edges.pop(i)
		return
	
	def getEdgesForNode(self, node, direction='A'):
		# direction := 'C' (childward), 'P', parentward, or 'A' all.  
		#  Note that A returns undirected edges as well as directed. C and P return only directed edges
		#print('direction = ', direction)
		#print('self.edges = ', self.edges)
		nodeEdges = []
		for i in range(len(self.edges)):
			n1, n2, eT, eW = self.edges[i]
			if (direction == 'A' and (n1==node or n2==node)) or \
				(eT == 'D' and ((direction == 'C' and n1 == node) or (direction == 'P' and n2 == node))):
				nodeEdges.append(self.edges[i])
		#print('Edges for node ', node, ' = ', nodeEdges)
		return nodeEdges
	

	def getEdges(self):
		return self.edges

	def getNodes(self):
		return self.nodes[:]
		
	def getExogenousNodes(self):
		xNodes = []
		for node in self.getNodes():
			parents = self.getParents(node)
			if len(parents) == 0:
				xNodes.append(node)
		return xNodes

	def _getEdge(self, node1, node2):
		outEdge = None
		for i in range(len(self.edges)):
			edge = self.edges[i]
			n1, n2, eT, eW = edge
			if n1 == node1 and n2 == node2:
				outEdge = i
				break
		return outEdge
		
	def setEdgeWeight(self, node1, node2, weight):
		edgeInd = self._getEdge(node1, node2)
		if edgeInd != None:
			n1, n2, eT, eW = self.edges[edgeInd]
			self.edges[edgeInd] = (n1, n2, eT, weight)
		return
	
	def getNodeValue(self, node):
		return self.nodeValues[node]
	
	def setNodeValue(self, node, value):
		self.nodeValues[node] = value
		return
	
	def calculate(self):
		Dprint('Calculate: Dag = ')
		Dprint(str(self))
		for i in range(len(self.nodes)):
			node = self.nodes[i]
			cumMean = 0
			cumVar = 0
			origMean, origStd = self.getNodeValue(node)
			parentEdges = self.getParentEdges(node)
			if len(parentEdges) > 0:
				for edge in parentEdges:
					n1, n2, eT, eW = edge
					parentVal = self.getNodeValue(n1)
					pMean, pStd = parentVal
					#print('parentVal = ', parentVal, eW)
					x1coef, x0coef = eW
					Dprint('   parentVal = ', parentVal)
					covTerm = 0
					for edge2 in parentEdges:
						nn1, nn2, neT, neW = edge2
						if edge2 == edge:
							continue
						nx1coef, nx0coef = neW
						
						tempcovTerm = self.getCov(n1, nn1) * 2 * x1coef * nx1coef
						covTerm += tempcovTerm
					cumMean += pMean * x1coef + x0coef
					cumVar += pStd**2 * x1coef**2 + covTerm
				newMean = cumMean
				Dprint( 'cumvar = ', cumVar)
				newVar = cumVar
				newStd = newVar**.5
				if node not in self.noises:
					noiseMean = origMean - newMean
					noiseVar = origStd**2 - newVar
					noiseStd = noiseVar**.5
					Dprint('node, noiseMean, noiseStd = ', node, noiseMean, noiseStd)
					self.noises[node] = (noiseMean, noiseStd)
				else:
					noiseMean, noiseStd = self.noises[node]
					noiseVar = noiseStd**2
				#print('prev, cum = ', prev, cum)
				newMean += noiseMean
				newVar += noiseVar
				newStd = newVar**.5
				
				if math.fabs(newMean - origMean) > 10**-10 or math.fabs(newStd - origStd) > 10**-10:
					Dprint('changing node value for ', node, ' from (', origMean,',',origStd,')' , ' to (', newMean, ',', newStd, ')')
					self.setNodeValue(node, (newMean, newStd ))
				else:
					Dprint('node value for ', node, ' is unchanged')
		return

	
	def getNodeNoises(self):
		# returns dict := {<nodeName>:<noiseLevel>}
		return self.noises.copy()
		
	def getChildEdges(self, node):
		edges = self.getEdgesForNode(node, direction='C')
		return edges
		
	def getParentEdges(self, node):
		edges = self.getEdgesForNode(node, direction='P')
		return edges
		
	def hasAncestor(self, node, ancestorNode):
		if ancestorNode in self.getAncestors(node):
			return True
		return False
	
	def getParents(self, node):
		parents = []
		pEdges = self.getParentEdges(node)
		for pEdge in pEdges:
			n1, n2, eT, eW = pEdge
			parents.append(n1)
		return parents
	
	def getChildren(self, node):
		children = []
		cEdges = self.getChildEdges(node)
		for cEdge in cEdges:
			n1, n2, eT, eW = cEdge
			children.append(n2)
		return children
		
	def getAncestors(self, node):
		#print('G = ', str(self))
		ancestors = []
		parents = self.getParents(node)
		#print('parents = ', parents)
		for parent in parents:
			ancestors.append(parent)
			pancestors = self.getAncestors(parent)
			for pa in pancestors:
				if pa not in ancestors:
					ancestors.append(pa)
		#print('Ancestors for: ', node, ' = ', ancestors)
		return ancestors
	
	def getVarOrder(self):
		return self.nodes[:]
		
	def __str__(self):
		outlines = []
		outlines.append('FlexDAG.FlexDag:')
		nodes = self.nodes[:]
		nodes.sort()
		for node in nodes:
			edgeCount = 0
			nodeV = self.getNodeValue(node)
			#print('node = ', node, ', nodeV = ', nodeV)
			nMean, nStd = nodeV
			outlines.append('  ' + node + ' (' + str(nMean) + ', ' + str(nStd) + '):')
			edges = self.getEdgesForNode(node, 'A')
			for edge in edges:
				n1, n2, eT, eW = edge
				if eT == 'U':
					if n1 == node:
						other = n2
					else:
						other = n1
					outlines.append('     <-->' + other)
					edgeCount += 1
				else:
					if n1 == node:
						outlines.append('     -->' + n2 + ' ([' + str(eW[0]) + ',' + str(eW[1]) + '])')
						edgeCount += 1
			if edgeCount == 0:
				outlines.append('     <Empty>')
		outstr = str.join('\n', outlines)
		return outstr
	
	def implications(self):
		outStrs = []
		for node in self.getNodes():
			impStr = self.implicationsForNode(node)
			outStrs.append(impStr)
		outStr = str.join('\n\n', outStrs)
		return outStr
		
	def implicationsForNode(self, node):
		outLines = []
		pEdges = self._sortEdges(self.getParentEdges(node))
		outLines.append('Observations for metric ' + node + ' are:')
		if len(pEdges) == 0:
			outLines.append('   This metric is not affected by any other metrics in the dataset')
		else:
			outLines.append('   This metric is affected by the following metrics:')
			for pEdge in pEdges:
				n1, n2, eT, eW = pEdge
				strength = self._getStrength(pEdge)
				outLines.append('      ' + n1 + ' -- ' + strength)
		cEdges = self._sortEdges(self.getChildEdges(node))
		if len(cEdges) == 0:
			outLines.append('   This metric does not affect any other metrics in the dataset')
		else:
			outLines.append('   This metric affects the following other metrics:')
			for cEdge in cEdges:
				n1, n2, eT, eW = cEdge
				strength = self._getStrength(cEdge)
				outLines.append('      ' + n2 + ' -- ' + strength)
		outStr = str.join('\n', outLines)
		return outStr
	
	def _getStrength(self, edge):
		n1, n2, eT, eW = edge
		pVal = self.getNodeValue(n1)
		pMean, pStd = pVal
		cVal = self.getNodeValue(n2)
		cMean, cStd = cVal
		influence = math.fabs(eW[0] * pStd**2/cStd**2)
		strength = 'moderately'
		x1coef, x0coef = eW
		if influence >= .75:
			strength = 'very strongly'
		elif influence >= .5:
			strength = 'strongly'
		elif influence <= .05:
			strength = 'faintly'
		elif influence <= .1:
			strength = 'very weakly'
		elif influence <= .2:
			strength = 'weakly'
		if x1coef < 0:
			strength += ' inverse'
			
		return strength
			
	def _sortEdges(self, edges):
		influences = []
		for edge in edges:
			n1, n2, eT, eW = edge
			parent = n1
			parentVal = self.getNodeValue(parent)
			pMean, pStd = parentVal
			childVal = self.getNodeValue(n2)
			cMean, cStd = childVal
			relInfluence = math.fabs(pStd * eW[0] / cStd)
			influences.append(relInfluence)
		zipped = list(zip(influences, edges))
		zipped.sort()
		zipped.reverse()
		edges = [edge for (influence, edge) in zipped]
		return edges

	def getCov(self, v1, v2):
		s1 = self.nodeD[v1]
		s2 = self.nodeD[v2]
		cc = np.cov([s1,s2], ddof=1)
		#print('cov = ', cc[1,0])
		return cc[1,0]
		
	def getEquivalencyKey(self):
		equivalent = True
		edges = self.getEdges()
		edges.sort()
		outEdges = []
		for edge in edges:
			n1, n2 = edge[:2]
			outEdges.append((n1,n2))
		equivKey = str(outEdges)
		#print('EquivKey = ', equivKey)
		return equivKey
		