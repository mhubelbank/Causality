from math import sin, cos, tan, pi, tanh, cosh, sinh, log
# ----------------------------------------------------------------------------
# SEM Definitions
#
# Simulated Time Series

varNames = ['A', 'B', 'C','D', 'E', 'F', 'G', 'H', 'I']

trend1 = 0
trend2 = 10
trend3 = 0

varEquations = ['A= logistic(10,1)', 
				'trend1 = trend1 + .0075',
				'B = logistic(20,1) + trend1', 
				'trend2 = trend2 - .02',
				'C = logistic(-100,1)', 
				'trend3 = trend3 + pi / 47',
				'D = logistic(11,1)', 
				'E = A - .2*B + .3 * C - .5 * D + logistic(-55,1)',
				'F = -1.1* E - logistic(100, 1)',
				'G = -.3*F + logistic(-1000,1)',
				'H = 1.2*F + logistic(0,1)',
				'I = -1.1* F + logistic(20,1)',
				]
				
validation =  [('A', []), ('B', []), ('C', []), ('D', []), ('E', ['A','B','C','D']),
				('F', ['E']), ('G', ['F']), ('H',['F']), ('I', ['F']), 
				]
# -----------------------------------------------------------------------------