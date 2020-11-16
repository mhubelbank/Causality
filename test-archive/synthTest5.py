# ----------------------------------------------------------------------------
# SEM Definitions
#
# Experimental set
from math import pi, sin, cos, tan, log, tanh, fabs, e
varNames = [ 'A', 'B', 'C']

density = 1.0
t = 0
t2 = 0

varEquations = [
				'A = sin(t) + noise()/100',
				'B = A + sin(t-.5) + noise()/100',
				'C = sin(t+.2) + noise()/100',
				't = t + .1',
				't2 = t2 - (logistic(.22, .0000001)/density)',
				]
				
validation = []
# -----------------------------------------------------------------------------