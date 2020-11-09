# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'W - Pattern Test'

varNames = ['A', 'B', 'C', 'D', 'E']

varEquations = [
				'A = logistic(0,1)',
				'B = logistic(0,1)',
				'C = logistic(0,1)',
				'D = .5 * A + .5 * B + logistic(0,1)',
				'E = .5 * B + .5 * C + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', []), ('C', []), ('D', ['A','B']), ('E', ['B', 'C'])
				]
# -----------------------------------------------------------------------------