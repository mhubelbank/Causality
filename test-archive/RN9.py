# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Non-linear Hourglass Pattern Test (3 level)'

from math import sin, cos, tanh

varNames = ['R11', 'R12',
			'R21',
			'R31', 'R32'
			]

t1 = 0.0
t2 = 0.0

varEquations = [
				'R11 = sin(t1) + logistic(0,.01)',
				'R12 = sin(t2) + logistic(0,.01)',
				'R21 = sin(.5 * R11 + .5 * R12) + logistic(0,.1)',
				'R31 = tanh(R21) + logistic(0,.1)',
				'R32 = cos(R21) + logistic(0,.1)',
				't1 = t1 + .1',
				't2 = t1 + .05',
				]
				
validation =  [
				('R11', []), ('R12', []), ('R21', ['R11','R12']), ('R31', ['R21']), ('R32', ['R21']),
				]
# -----------------------------------------------------------------------------