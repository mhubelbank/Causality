# ----------------------------------------------------------------------------
# SEM Definitions
#

from math import sin, pi

testDescript = '3-Chain Test with exogenous sine wave input'

varNames = ['A', 'B', 'C']

t = 0

varEquations = [
				'A = 1 * (2+sin(t))',
				'B = A + logistic(0,1)',
				'C = B + logistic(0,1)',
				't = t + pi / 5000',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['B']), 
				]
# -----------------------------------------------------------------------------