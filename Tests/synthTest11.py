# ----------------------------------------------------------------------------
# SEM Definitions
#
# Experimental set
from math import pi, sin, cos, tan, log, tanh, fabs, e, log
varNames = [ 'Linear', 'Erlang', 'NatLog', 'XSquared', 'XCubed', 'SqrRoot', 'NLogN']

density = 1.0
t = .01
t2 = 0

varEquations = [
				'Linear = t',
				'Erlang = 1 - e**-t',
				'NatLog = log(t,e)',
				'NLogN = t * log(t,e)',
				'XSquared = t**2',
				'XCubed =t**3',
				'SqrRoot = t**.5',
				't=t+.01',
				]
				
validation =  [('A', []), ('B', []), ('C', []), ('D', []), ('E', ['A','B','C','D']),
				('F', ['E']), ('G', ['F']), ('H',['F']), ('I', ['F']), 
				]
# -----------------------------------------------------------------------------