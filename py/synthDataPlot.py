import synthDataGen
import getData
import numpy as np
import statsmodels.api as sm
from scipy.stats.distributions import norm
import matplotlib.pyplot as plt
import sys

PlotColors = ['green', 'red', 'blue', 'orange']

def run(dataDefFile, datapoints):
	print('Generating data using: ', dataDefFile)
	fileExt = dataDefFile.split('.')[-1]
	if fileExt == 'py':
		# Data Definition File
		outFile = synthDataGen.run(dataDefFile, datapoints)
	elif fileExt == 'csv':
		# Data file
		outFile = dataDefFile
	else:
		print('*** Invalid file type = ', fileExt)
		return
	d = getData.DataReader(outFile, limit=datapoints)
	print('Vars = ', d.getSeriesNames())
	plotData(d, datapoints)


def plotData(d, datapoints):
	if d.sampleCount < datapoints:
		datapoints = d.sampleCount
	slices = datapoints
	x_grid = np.linspace(0, 1, slices)
	vNames = d.getSeriesNames()
	D = []
	PDF = []
	vNames.sort()
	for i in range(min(len(vNames),100)):
		dS = standardize(d.getSeries(vNames[i]))
		D.append(dS)
		k = sm.nonparametric.KDEMultivariate(data=dS, var_type='c', bw='normal_reference')
		#bw = k.bw[0]
		PDF.append(scalePdf(k.pdf(data_predict=x_grid)))
	plot(vNames, D, PDF, x_grid)

def plot(vNames, D, PDF, x_grid=None):
	numplots = len(D)
	fillData  = None
	# Plot the  estimates
	rows = max([int(numplots / 5.0 + .99),2])
	cols = min([numplots, 5])
	#print ('rows, cols = ', rows, cols)
	fig, ax = plt.subplots(rows, cols, sharey=True,
						   figsize=(13, 13))
	#fig.subplots_adjust(wspace=0)
	
	row = 0
	col = 0
	if x_grid is None:
		x_grid = range(Dcount)

	for i in range(numplots):
		lineData = [D[i], PDF[i]]
		vName = vNames[i]
		for j in range(len(lineData)):
			if j == 0:
				#print('ax: ', ax, ', x_grid: ', x_grid, 'lineData[j]: ', lineData[j], ', PlotColors[j]', PlotColors[j])
				ax[row, col].plot(x_grid, lineData[j], color=PlotColors[j], alpha=.5, lw=1)
			else:
				#ax[row, col].plot(x_grid, lineData[j], color='black', alpha=.5, lw=(j+1)/2)
				lineData[j][0] = 0.0
				lineData[j][-1] = 0.0
				ax[row, col].fill(x_grid, lineData[j], ec='black', fc='gray', alpha=0.4)
		ax[row, col].set_title(vName)
		ax[row, col].set_ylim(0,1)
		col += 1
		if col >= 5:
			col = 0
			row += 1
		
		#ax[i].set_xlim(0, Dcount)
		
	
	from IPython.display import HTML
	HTML("<font color='#666666'>Gray = True underlying distribution</font><br>"
     "<font color='6666ff'>Blue = Estimated distribution (500 pts)</font>")
	plt.show()
	return
	
def standardize(series):
	# Recenter and scale to interval [0,1] 
	a = np.array(series)
	maxD = a.max()
	minD = a.min()
	rangeD = maxD - minD
	dZeroed = a - minD
	dScaled = dZeroed / rangeD
	#print('minD, maxD = ', minD, maxD)
	#print('dScaled = ', dScaled)
	return dScaled

def scalePdf(series):
	a = np.array(series)
	max = a.max()
	aScaled = a/max
	return aScaled
	
#print('args = ', sys.argv)	
dataGenFile = sys.argv[1]
#print('File = ', dataGenFile)
datacount = 1000
if len(sys.argv) > 2:
	datacount = int(sys.argv[2])
run(dataGenFile, datacount)
