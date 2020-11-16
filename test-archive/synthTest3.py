from math import sin, cos, tan, pi, tanh, cosh, sinh, log
# ----------------------------------------------------------------------------
# SEM Definition
#
# Non-linear relationships


varNames = ['A', 'B', 'C']
t = 0

varEquations = ['A = cos(t) + normal(0,.1)',
				'B = A + cos(A) + normal(0,.1)',
				'C = pi**B + normal(0,.1)',
				]

validation =  [('A', []), ('B', ['A']), ('C', ['B']), 
				]
# -----------------------------------------------------------------------------