# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Gaussian Inverted V-Structure (IV) 1->2 Test'

varNames = ['A', 'B', 'C']

varEquations = [
				'A = normal(0,1)',
				'B = A + normal(0,1)',
				'C = A + normal(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A']), 
				]
# -----------------------------------------------------------------------------