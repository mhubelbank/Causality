from debug import *
import numpy as np
import getData
import math
import FlexDAG
import CaIC
import IndHSIC
import CaLingam
import CaPC


LINGAM_EFFECTIVENESS_THRESHOLD = .002



class IC(CaIC.IC):
	def __init__(self, dataFile, limit=None, prune=True):
		super().__init__(dataFile, limit=limit, prune=prune)
		self.myIC = None
		self.myIndTest = IndHSIC
		self.cal = CaLingam.IC(dataFile, limit=limit) # Temporary use for direction algo
		return
					
	def method(self):
		if self.myIC is None:
			self.chooseICMethod()
		return self.myIC.method()

	def getDAG(self):
		vars = self.data.getSeriesNames()
		if self.myIC is None:
			self.chooseICMethod()
		dag = self.myIC.getDAG()
		return dag

	def parameterizeLinearModel(self, dag):
		# Modifies the dag passed in and returns it as the model
		if self.myIC is None:
			self.chooseICMethod()
		dag = self.myIC.parameterizeLinearModel(self, dag)
		return dag
		
	def chooseICMethod(self):
		vars = self.data.getSeriesNames()
		score = 0.0
		tests = 0
		for i in range(len(vars)):
			for j in range(len(vars)):
				v1 = vars[i]
				v2 = vars[j]
				if self.isIndependent(v1, v2):
					continue					
				tempScore = self.scoreLingam(v1, v2)
				score += tempScore
				# if tempScore > score:
					# score = tempScore
				tests += 1
		score = score / tests
		Dprint('CaICAny.chooseICMethod: score, threshold = ', score, LINGAM_EFFECTIVENESS_THRESHOLD)
		if score > LINGAM_EFFECTIVENESS_THRESHOLD:
			icType = CaLingam
			Dprint('CaICAny.chooseICMethod: Using LiNGAM Method')
		else:
			Dprint('CaICAny.chooseICMethod: Using PC Method', score)
			icType = CaPC
			
		self.myIC = icType.IC(self.dataFile, limit=self.limit, prune=self.doPrune)
		return
		
	def scoreLingam(self, v1, v2):
		direction1 = self.cal.direction(v1, v2)
		direction2 = self.cal.direction(v2, v1)
		#score = math.fabs(direction1)
		score = math.fabs(direction2 - direction1)
		return score
		
	def scoreLingam_alt(self, v1, v2):
		s1D = np.array(self.data.getSeries(v1))
		s2D = np.array(self.data.getSeries(v2))
		coefs = np.polyfit(s1D, s2D.T, 1)
		x1coef, x0coef = coefs
		Dprint('coefs = ', coefs)
		cs2D = s1D * x1coef + x0coef
		res = s2D - cs2D
		dep = IndHSIC.scoreDependence(s1D, res)
		return dep

	def isIndependent(self, v1, v2):
		d1 = self.data.getSeries(v1)
		d2 = self.data.getSeries(v2)
		result = self.myIndTest.isIndependent(d1, d2)
		return result
		
		
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
