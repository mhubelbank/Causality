# ----------------------------------------------------------------------------
# SEM Definitions
#
# Experimental set
from math import pi, sin, cos, tan, log, tanh, fabs, e, log
testDescript = 'attractor'
varNames = ['s1', 's2', 's3']

density = 1.0
incr = .2345
t = .0
t1 = t + incr
t2 = t - 2 * incr

varEquations = [ 
			's1 = sin(t) + noise()/200',
			's2 = sin(t/3) + noise()/200',
			#'s2 = sin(t/3) + noise()/100',
			's3 = s2 + s1',
			#'s3 = logistic(0,1)',
			't = t+ incr',
			]
			
validation =  [
				]
# -----------------------------------------------------------------------------