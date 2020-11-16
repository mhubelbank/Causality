# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Gaussian V-Structure Test'

varNames = ['A', 'B', 'C']

varEquations = [
				'A= normal(0,1)',
				'C = normal(0,1)',
				'B = .5 * A + .5 * C + normal(0,1)',

				]
				
validation =  [
				('A', []), ('B', ['A','C']), ('C', []), 
				]
# -----------------------------------------------------------------------------