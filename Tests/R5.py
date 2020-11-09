# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Inverted W - Pattern Test'

varNames = ['A', 'B', 'C', 'D', 'E']

varEquations = [
				'A = logistic(0,1)',
				'B = logistic(0,1)',
				'C = A + logistic(0,1)',
				'D = .5 * A + .5 * B + logistic(0,1)',
				'E = B + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', []), ('C', ['A']), ('D', ['A','B']), ('E', ['B'])
				]
# -----------------------------------------------------------------------------