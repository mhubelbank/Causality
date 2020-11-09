import numpy as np
import getData
import math
import FlexDAG
import entropy_estimators as ee
import PlotN
import statsmodels.api as sm
from statsmodels.tsa.stattools import adfuller


indThreshold = 2
pearlCIThreshold2 = .01
miIndThreshold = .02
miCondIndThreshold = -.002

indMethod = 'Gauss'
condIndMethod = 'Gauss'
dataLimit = 100

class IC():

	def __init__(self, dataFile, limit=None, prune=True, adjustData=True):
		if limit is None:
			limit = dataLimit
		self.limit = limit
		self.doPrune = prune
		self.adjustData = adjustData
		self.dataFile = dataFile
		self.dataCache = {}
		self.data = getData.DataReader(dataFile, limit=limit)
		self.corrMatrix = None
		self.varIndexes = {}

		return
	
	def method(self):
		return self.__class__.__module__
		
	# Public function to generate parameterized model of causal network.  Returned as a FlexDAG
	def getModel(self, mtype='linear'):
		dag = self.getDAG()
		model = None
		if mtype=='linear':
			model = self.parameterizeLinearModel(dag)
		else:
			print('CaIC.getModel: Unsupported model type = ', mtype)

		return model	
	
	# Virtual functions
	def getDAG(self):
		dag = FlexDAG.FlexDAG()
		return dag
	
	def prune(self, dag):
		return dag
	
	def parameterizeLinearModel(self, dag):
		# Modifies the dag passed in and returns it as the model
		return dag
		
	# Helper Function
	
	def getDataSeries(self, sName):
		if sName in self.dataCache:
			X = self.dataCache[sName]
		else:
			X = self.data.getSeries(sName)
			if self.adjustData:
				X = self.makeStationary(X)
			self.dataCache[sName] = X
		return X
		
	def standardize(self, series):
		# Recenter at zero mean, and rescale to unit variance
		a = np.array(series)
		mu = a.mean()
		aCentered = a - mu
		sigma = aCentered.std()
		#print ('sigma= ', sigma)
		aScaled = aCentered / sigma
		#print('aScaled = ', aScaled.mean(), aScaled.std())
		#print('ac=', aScaled)
		return aScaled
		
	def getCorrCoef(self, var1, var2):
		self.genCorrMatrix()
		indx1 = -1
		indx2 = -1
		nodeNames = self.data.getSeriesNames()
		for i in range(len(nodeNames)):
			nodeName = nodeNames[i]
			if nodeName == var1:
				indx1 = i
			if nodeName == var2:
				indx2 = i
			if indx1 >= 0 and indx2 >= 0:
				break
		corrCoef = self.corrMatrix[indx1, indx2]
		return corrCoef
	
	def genCorrMatrix(self):
		if self.corrMatrix is None:
			arrayData = []
			for var in self.data.getSeriesNames():
				sData = self.standardize(self.getDataSeries(var))
				#sData = self.data.getSeries(var)
				arrayData.append(sData)
			a = np.array(arrayData)
			self.corrMatrix = np.corrcoef(a)
		return
		
	def mi_prepare(self, ds):
		out = []
		for dp in ds:
			out.append([dp])
		return out
		
	def plotData(self):
		slices = 1000
		x_grid = np.linspace(-2.0, 2.0, slices)
		D = []
		PDF = []
		vNames = self.data.getSeriesNames()
		
		for i in range(len(vNames)):
			D.append(self.standardize(self.data.getSeries(vNames[i])))
			k = sm.nonparametric.KDEMultivariate(data=D[i], var_type='c', bw='normal_reference')
			bw = k.bw[0]
			PDF.append(k.pdf(data_predict=x_grid))
		PlotN.plot(PDF, x_grid)
	
	def makeStationary(self, X):
		# X is a series.
		# Run Augmented Dickey-Fuller Test to detect non-stationarity.  Then difference the series until it is stationary.
		Xa = np.array(X)
		result = adfuller(Xa)
		pvalue = result[1]
		#print('result = ', result)
		while pvalue > .05:
			# Non stationary -- Difference the series and try again
			print('differencing')
			X2 = [0]
			X2.extend(X[:-1])
			X2a = np.array(X2)
			X3a = Xa - X2a
			result = adfuller(X3a)
			pvalue = result[1]
			Xa = X3a
			#print('result2 = ', result)
		return list(Xa)

# cal.direction('138.42.94.173|Fa3/1/0|bytesIn', '138.42.94.173|Gi0/0/0|bytesOut')

