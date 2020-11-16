# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Gaussian Diamond Pattern Test'

varNames = ['A', 'B', 'C', 'D']

varEquations = [
				'A = normal(0,1)',
				'B = A + normal(0,1)',
				'C = A + normal(0,1)',
				'D = .5 * B + .5 * C + normal(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A']), ('D', ['B','C']), 
				]
# -----------------------------------------------------------------------------