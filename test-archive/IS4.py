# ----------------------------------------------------------------------------
# SEM Definitions
#
from math import sin, pi

testDescript = 'Sine Wave Dual Input -- V Structure'

varNames = ['A', 'B', 'C']

trendFactor = 10
t1 = 0
t2 = 0


varEquations = [
				'A = trendFactor * sin(t1)',
				'B = trendFactor * sin(t2)',
				'C = .5 * A + .5 * B + logistic(0,1)',
				't1 += pi / 500',
				't2 += pi / 290',
				]
				
validation =  [
				('A', []), ('B', []), ('C', ['A','B']), 
				]
# -----------------------------------------------------------------------------