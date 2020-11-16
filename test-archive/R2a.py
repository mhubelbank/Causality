# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Inverted V-Structure (IV) 1->3 Test'

varNames = ['A', 'B', 'C', 'D']

varEquations = [
				'A = logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = A + logistic(0,1)',
				'D = A + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A']), ('D', ['A']), 
				]
# -----------------------------------------------------------------------------