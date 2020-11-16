# ----------------------------------------------------------------------------
# SEM Definitions
#
# Experimental set
from math import pi, sin, cos, tan, log, tanh, fabs, e, log
varNames = [ 'TestF', 'MainF', 'HighF', 'LowF']

density = 1.0
t = .01
t2 = 0

varEquations = [
				'TestF = logistic(100,10)',
				'MainF = logistic(200, 8)',
				'HighF = logistic(.6,.01) * TestF + logistic(.6, .01) * MainF',
				'LowF  = TestF + MainF - HighF',
				]
			
validation =  [('A', []), ('B', ['A']), ('C', ['B']), 
				]
# -----------------------------------------------------------------------------