# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Inverted V-Structure (IV) 1->2 Test'

varNames = ['A', 'B', 'C']

varEquations = [
				'A = logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = A + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A']), 
				]
# -----------------------------------------------------------------------------