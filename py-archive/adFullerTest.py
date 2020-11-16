from statsmodels.tsa.stattools import *
import getData
import numpy as np

def run():
	#d = getData.DataReader(input='..\data\experiment1.csv', limit=2000)
	# d = getData.DataReader(input='..\Tests\IS1c.csv', limit=2000)
	# X = d.getSeries('R11')
	d = getData.DataReader(input='..\data\experiment1.csv', limit=2000)
	X = d.getSeries('Io')
	print('acf = ', acf(X), ', pacf = ', pacf(X))
	print('arma = ', arma_order_select_ic(X))
	result = adfuller(X)
	pvalue = result[1]
	print('result = ', result)
	if pvalue > .05:
		# Non stationary -- Difference the series and try again
		Xa = np.array(X)
		X2 = [0]
		X2.extend(X[:-1])
		X2a = np.array(X2)
		X3a = Xa - X2a
		result = adfuller(X3a)
		print('result2 = ', result)
	
	
run()