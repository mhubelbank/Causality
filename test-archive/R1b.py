# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'V-Structure 5->1'

varNames = ['A', 'B', 'C', 'D', 'E', 'F']

varEquations = [
				'A = logistic(0,1)',
				'B = logistic(0,1)',
				'C = logistic(0,1)',
				'D = logistic(0,1)',
				'E = logistic(0,1)',
				'F = .2 * A + .2 * B + .2 * C + .2 * D + .2 * E + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', []), ('C', []), ('D', []), ('E', []), ('F', ['A', 'B','C','D','E']), 
				]
# -----------------------------------------------------------------------------