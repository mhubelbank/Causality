import numpy as np
from scipy.stats.distributions import norm
import matplotlib.pyplot as plt
PlotColors = ['blue', 'green', 'red', 'orange']

# data = [(est, true)]
def plot(lineData, x_grid=None):
	numplots = 2
	fillData  = None
	# Plot the  estimates
	fig, ax = plt.subplots(1, numplots, sharey=False,
						   figsize=(13, 3))
	fig.subplots_adjust(wspace=0)
	

	for i in range(numplots):
		if x_grid is None:
			x_grid = range(Dcount)
		for j in range(len(lineData)):
			ax[i].plot(x_grid, lineData[j], color=PlotColors[j], alpha=0.5, lw=3)
		
		#ax[i].fill(x_grid, trueD, ec='gray', fc='gray', alpha=0.4)
		ax[i].set_title('Test')
		#ax[i].set_xlim(0, Dcount)
		#ax[i].set_ylim(0,2)
		
	
	from IPython.display import HTML
	HTML("<font color='#666666'>Gray = True underlying distribution</font><br>"
     "<font color='6666ff'>Blue = Estimated distribution (500 pts)</font>")
	plt.show()
	return
	

