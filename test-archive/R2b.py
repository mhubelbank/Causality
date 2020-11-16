# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Inverted V-Structure (IV) 1->5 Test'

varNames = ['A', 'B', 'C', 'D', 'E', 'F']

varEquations = [
				'A = logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = A + logistic(0,1)',
				'D = A + logistic(0,1)',
				'E = A + logistic(0,1)',
				'F = A + logistic(0,1)',				
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A']), ('D', ['A']), ('E', ['A']), ('F', ['A']), 
				]
# -----------------------------------------------------------------------------