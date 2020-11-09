# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Diamond Pattern Test'

varNames = ['A', 'B', 'C', 'D']

varEquations = [
				'A = logistic(0,1)',
				'B = logistic(0,1)',
				'C = logistic(0,1)',
				'D = .33 * A + .33 * B + .33 * C + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', []), ('C', []), ('D', ['A', 'B','C']), 
				]
# -----------------------------------------------------------------------------