import math
import numpy as np
import statsmodels.api as sm
from scipy.stats import norm
import statsmodels.nonparametric.kernel_regression as kr
from scipy.stats.distributions import norm

Zero_threshold = .3
GAMMA_size = 3
GAMMAS = [
		[np.random.rand(), np.random.rand(),np.random.rand(),np.random.rand()],
		[np.random.rand(), np.random.rand(),np.random.rand(),np.random.rand()],
		[np.random.rand(), np.random.rand(),np.random.rand(),np.random.rand()],
		[np.random.rand(), np.random.rand(),np.random.rand(),np.random.rand()],
		[np.random.rand(), np.random.rand(),np.random.rand(),np.random.rand()]
		]
gIndex = 0

def isCondIndependent(ser1, ser2, ser3):
	ser1 = standardize(ser1)
	ser2 = standardize(ser2)
	ser3 = standardize(ser3)
	#print ('x= ', ser1)
	h = .5
	result = False
	k = sm.nonparametric.KDEMultivariate(data=ser3, var_type='c', bw='normal_reference') # cv_ls
	GAMMA = chooseGAMMA(GAMMA_size)
	GAMMA = GAMMAS[:GAMMA_size]
	cum = 0.0
	for gamma in GAMMA:
		dnh = delta_n_h(ser1, ser2, ser3, gamma, k, h)
		#print('dnh = ', dnh)
		cum += math.fabs(dnh)
	mu = cum / GAMMA_size
	print ('mu = ', mu)
	if mu <= Zero_threshold:
		result = True
	return result

def chooseGAMMA(n):
	GAMMA = []
	for i in range(n):
		#gamma = [np.random.rand(), np.random.rand(), np.random.rand(), np.random.rand()]
		gamma = [.76642, .98837, .013333, .54321]
		GAMMA.append(gamma)
	return GAMMA
	
def delta_n_h(s1, s2, s3, gamma, kz, h):

	#print ('gamma = ', gamma)
	t1 = term1(s1, s2, s3, gamma, kz, h)
	t2 = term2(s1, s2, s3, gamma, kz, h)
	result = t1 - t2
	#print ('delta_n_h = ', result)
	# n = len(s1)
	# factor = 1.0 / (n * (n-1))
	# sum = 0
	# for i in range(n):
		# for j in range(n):
			# if i >= j:
				# continue
			# Wi = (s1[i], s2[i], s3[i])
			# Wj = (s1[j], s2[j], s3[j])
			# term =  Kh2(gamma, h, k, Wi, Wj)
			# print ('term = ', term, i, j)
			# sum += term
	# result = factor * sum
	return result
	
def isIndependent(ser1, ser2):

	return

def term1(x, y, z, gamma, kz, h):
	n = len(x)
	cum = 0.0
	for i in range(n):
		Wi = (x[i], y[i], z[i])
		xp1 = phi_i(Wi, gamma)
		xp2 = Kh(h, kz, z[i])
		xp =  xp1 * xp2
		cum += xp
	result = cum / n
	print ('E(phi(x)*pdf(z)) = ', result)
	return result
	
def term2(x, y, z, gamma, kz, h):
	n = len(x)
	xps = []
	for i in range(n):
		Wi = (x[i], y[i], z[i])
		xp1 = phi_i(Wi, gamma)
		xps.append(xp1)
	nw = kr.KernelReg(xps, [z], 'c', reg_type='lc', bw='aic')
	means, marginals = nw.fit([z])
	for i in range(n):
		xp2 = Kh(h, kz, z[i])
		means[i] *= xp2
	mean = means[i].mean()
	result = mean
	print('E(phi(X|Z)*pdf(Z)) = ', result)
	return result
	
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
	

def tanh(x):
	return math.tanh(x)
	
def normcdf(x):
	result = norm.cdf(x)
	return result
	
def sine(x):
	return math.sin(x)

def cossin(x):
	return math.sin(x) + math.cos(x)
	
def standardize(series):
	# Map to [0,1]
	a = np.array(series)
	aMin = a.min()
	aMax = a.max()
	aRange = aMax - aMin
	aTrans = a - aMin
	aScaled = aTrans / aRange
	return aScaled

GCRF = cossin


import PlotN
numpoints = 1000
x = np.random.normal(.5,.5, numpoints)
#x = np.random.random(numpoints) * 47
y = np.random.normal(-10,.5, numpoints) + x * .5
z = np.random.normal(0,1, numpoints) + y * -5
#z = np.random.normal(0,1, numpoints)
print ('x _|_ y | z) = ', isCondIndependent(x,y,z))
# print ('x _|_ y | y) = ', isCondIndependent(x,y,y))
# print ('x _|_ z | y) = ', isCondIndependent(x,z,y), True)
# print ('x _|_ y | y) = ', isCondIndependent(x,y,y), True)
# print ('x _|_ z | z) = ', isCondIndependent(x,z,z), True)
# print ('y _|_ z | x) = ', isCondIndependent(y,z,x), False)
# print ('z _|_ y | x) = ', isCondIndependent(z,y,x), False)
# print ('y _|_ y | x) = ', isCondIndependent(y,y,x), False)
x = standardize(x)
y = standardize(y)
z = standardize(z)
slices = 1000

x_grid = np.linspace(0.0, 1.0, slices)
k = sm.nonparametric.KDEMultivariate(data=x, var_type='c', bw='normal_reference')
bw = k.bw[0]
xpdf = k.pdf(data_predict=x_grid)
k = sm.nonparametric.KDEMultivariate(data=y, var_type='c', bw='normal_reference')
ypdf = k.pdf(data_predict=x_grid)
k = sm.nonparametric.KDEMultivariate(data=z, var_type='c', bw='normal_reference')
zpdf = k.pdf(data_predict=x_grid)
PlotN.plot([xpdf, ypdf,zpdf], x_grid)

