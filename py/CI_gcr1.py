import math
import numpy as np
import statsmodels.api as sm
from scipy.stats import norm

Zero_threshold = .01
GAMMA_size = 10

def isCondIndependent(ser1, ser2, ser3):
	ser1 = standardize(ser1)
	ser2 = standardize(ser2)
	ser3 = standardize(ser3)
	h = .5
	result = False
	k = sm.nonparametric.KDEMultivariate(data=ser3, var_type='c', bw='normal_reference')
	GAMMA = chooseGAMMA(GAMMA_size)
	cum = 0.0
	for gamma in GAMMA:
		dnh = delta_n_h(ser1, ser2, ser3, gamma, k, h)
		print('dnh = ', dnh)
		cum += dnh
	mu = cum / GAMMA_size
	print ('mu = ', mu)
	if mu <= Zero_threshold:
		result = True
	return result

def chooseGAMMA(n):
	GAMMA = []
	for i in range(n):
		gamma = [np.random.rand(), np.random.rand(), np.random.rand(), np.random.rand()]
		GAMMA.append(gamma)
	return GAMMA
	
def delta_n_h(s1, s2, s3, gamma, k, h):
	n = len(s1)
	factor = 1.0 / (n * (n-1))
	sum = 0
	for i in range(n):
		for j in range(n):
			if i >= j:
				continue
			Wi = (s1[i], s2[i], s3[i])
			Wj = (s1[j], s2[j], s3[j])
			term =  Kh2(gamma, h, k, Wi, Wj)
			#print ('term = ', term, i, j)
			sum += term
	result = factor * sum
	return result
	
def isIndependent(ser1, ser2):

	return

	
def Kh2(gamma, h, k, Wi, Wj):
	Zi = Wi[2]
	Zj = Wj[2]
	diff1 = Zi - Zj
	diff2 = Zj - Zi
	factor1 = Kh(h, k, diff1) * .5
	factor2 = Kh(h, k, diff2) * .5
	term1 = (phi_i(Wi, gamma) - phi_i_j(Wi, Wj, gamma)) * factor1
	term2 = (phi_i(Wj, gamma) - phi_i_j(Wj, Wi, gamma)) * factor2
	result = term1 + term2
	return result

def Kh(h, k, u):
	est = k.pdf([u])
	#print ('est = ', est)
	return est
	
def phi_i(Wi, gamma):
	result = GCRF(gamma[0] + gamma[1] * Wi[0] + gamma[2] * Wi[1] + gamma[3] * Wi[2])
	return result
	
def phi_i_j(Wi, Wj, gamma):
	result = GCRF(gamma[0] + gamma[1] * Wi[0] + gamma[2] * Wj[1] + gamma[3] * Wi[2])
	return result
	

	
def normcdf(x):
	result = norm.cdf(x)
	return result
	
def sine(x):
	return math.sin(x)

def standardize(series):
	# Recenter at zero mean, and rescale to unit variance
	a = np.array(series)
	mu = a.mean()
	muArray = np.full_like(a, mu)
	aCentered = a - mu
	sigma2 = a.std()
	#print ('sigma2= ', sigma2)
	stdArray = np.full_like(a, sigma2)
	aScaled = aCentered / stdArray
	#print('ac=', aScaled)
	return aScaled

GCRF = normcdf

import PlotN
x = np.random.normal(1,1, 100)
y = np.random.normal(1,1, 100) + x * .5
#z = np.random.normal(1,1, 100) + y * 5
z = np.random.normal(1,1, 100)

print ('x _|_ z | y = ', isCondIndependent(x,z,y))
print ('x _|_ y | y = ', isCondIndependent(x,y,y))
print ('y _|_ z | x = ', isCondIndependent(y,z,x))
print ('z _|_ y | x = ', isCondIndependent(z,y,x))
x = standardize(x)
y = standardize(y)
z = standardize(z)
x_grid = np.linspace(-4.5, 3.5, 1000)
k = sm.nonparametric.KDEMultivariate(data=x, var_type='c', bw='normal_reference')
xpdf = k.pdf(data_predict=x_grid)
k = sm.nonparametric.KDEMultivariate(data=y, var_type='c', bw='normal_reference')
ypdf = k.pdf(data_predict=x_grid)
k = sm.nonparametric.KDEMultivariate(data=z, var_type='c', bw='normal_reference')
zpdf = k.pdf(data_predict=x_grid)
PlotN.plot([(xpdf,xpdf), (ypdf,ypdf), (zpdf,zpdf)])

