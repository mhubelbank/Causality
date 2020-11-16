import synthDataGen
import getData
import numpy as np
import statsmodels.api as sm
from scipy.stats.distributions import norm
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d.axes3d import Axes3D
import standardize

TAU = 5

standard = True

PlotColors = ['green', 'red', 'blue', 'orange']

def run(dataDefFile, var1, var2, datapoints):
	tau = TAU
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
	fig = plt.figure()
	ax = fig.add_axes([.1,.1,.8,.8],projection='3d')
	
	# If var1 and var2 are not specified, then build a master manifold out of three vars at a time.  Otherwise, build shadow manifolds
	# from the two specified vars
	if var1 is None:
		vars = d.getSeriesNames()
		vars.sort()
		print('Vars = ', vars)
		colors = ['b', 'g', 'r','o'] 
		for i in range(4):
			if len(vars) < i*3+3:
				break
			X = d.getSeries(vars[i*3])
			Y = d.getSeries(vars[i*3+1])
			Z = d.getSeries(vars[i*3+2])
			if standard:
				 X = standardize.standardize(X)
				 Y = standardize.standardize(Y)
				 Z = standardize.standardize(Z)
			color = colors[i]
			ax.plot(X, Y, Z)
			
	else:
		var1D = d.getSeries(var1)
		if standard:
			var1D = standardize.standardize(var1D)
		X1 = var1D[:-2*tau]
		Y1 = var1D[tau:-tau]
		Z1 = var1D[2*tau:]

		var2D = d.getSeries(var2)
		if standard:
			var2D = standardize.standardize(var2D)
		X2 = var2D[:-2*tau]
		Y2 = var2D[tau:-tau]
		Z2 = var2D[2*tau:]
		
		#plotData(d, datapoints)
		ax.plot(X1,Y1,Z1)
		ax.plot(X2,Y2,Z2, 'r')
		#ax.plot(X1, Y1, Z2, 'g')
	
	# `ax` is a 3D-aware axis instance because of the projection='3d' keyword argument to add_subplot
	#ax = fig.add_subplot(1, 2, 1, projection='3d')
	plt.show()
	
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
	


def scalePdf(series):
	a = np.array(series)
	max = a.max()
	aScaled = a/max
	return aScaled
	
#print('args = ', sys.argv)	
dataGenFile = sys.argv[1]
datacount = 1000

if len(sys.argv) > 3:
	var1 = sys.argv[2]
	var2 = sys.argv[3]
	if len(sys.argv) > 4:
		datacount = int(sys.argv[4])
else:
	var1 = None
	var2 = None
	if len(sys.argv) > 2:
		datacount = int(sys.argv[2])
#print('File = ', dataGenFile)

run(dataGenFile, var1, var2, datacount)
