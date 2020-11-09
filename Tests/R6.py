# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Triangle Pattern Test (3 level)'

varNames = ['A', 
			'B', 'C', 
			'D', 'E', 'F']

varEquations = [
				'A = logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = A + logistic(0,1)',
				'D = B + logistic(0,1)',
				'E = .5 * B + .5 * C + logistic(0,1)',
				'F = C + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A']), ('D', ['B']), ('E', ['B','C']), ('F', ['C']),
				]
# -----------------------------------------------------------------------------