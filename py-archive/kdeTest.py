from statsmodels.nonparametric.kernel_density import KDEMultivariate
import numpy as np
from scipy.stats.distributions import norm
from scipy.stats.distributions import beta
import matplotlib.pyplot as plt

def kde_statsmodels_m(x, x_grid, bandwidth=0.5, **kwargs):
    """Multivariate Kernel Density Estimation with Statsmodels"""
    kde = KDEMultivariate(x, bw=bandwidth * np.ones_like(x),
                          var_type='c', **kwargs)
    return kde.pdf(x_grid)
	
def plot(dist_est, dist_true, x_grid=None):
	if x_grid is None:
		x_grid = np.linspace(-4.5, 3.5, 1000)
	
	pdf_true = dist_true.pdf(x_grid)
	pdf_est = dist_est.pdf(x_grid)
	print ('pdf_true = ', pdf_true)
	print ('pdf_est = ', pdf_est)
	#pdf_true = (0.8 * norm(-1, 1).pdf(x_grid) +
	#			0.2 * norm(1, 0.3).pdf(x_grid))

	# Plot the  estimates
	fig, ax = plt.subplots(1, 2, sharey=True,
						   figsize=(13, 3))
	fig.subplots_adjust(wspace=0)

	for i in range(2):
		
		#pdf = kde_statsmodels_m(x, x_grid, bandwidth=0.2)
		ax[i].plot(x_grid, pdf_est, color='blue', alpha=0.5, lw=3)
		ax[i].fill(x_grid, pdf_true, ec='gray', fc='gray', alpha=0.4)
		ax[i].set_title('Test')
		ax[i].set_xlim(-4.5, 3.5)
		ax[i].set_ylim(0,2)
	from IPython.display import HTML
	HTML("<font color='#666666'>Gray = True underlying distribution</font><br>"
     "<font color='6666ff'>Blue = KDE model distribution (500 pts)</font>")
	plt.show()
	return
	

import statsmodels.api as sm
nobs = 50
#np.random.seed(1234)  # Seed random generator
c1 = np.random.beta(.2,.3, size=(nobs,1))
c2 = np.random.normal(50, 25, size=(nobs,1))
#print('c1 = ', c1)

# Estimate a bivariate distribution and display the bandwidth found:

#dens_u = sm.nonparametric.KDEMultivariate(data=[c1], var_type='c', bw='normal_reference')
dens_u = sm.nonparametric.KDEMultivariate(data=[c1], var_type='c', bw=[.1])


cum = 0
max = -100
maxx = -100


# for i in range(len(pdf)):
	# p = pdf[i]
	# x = tests[i]
	# if p > max:
		# print ('i = ', i)
		# max = p
		# maxx = x
	# cum += p

#mean = cum / len(pdf)


plot(dens_u, beta(.2, .3))