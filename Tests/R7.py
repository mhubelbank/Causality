# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Inverted Triangle Pattern Test (3 level)'

varNames = ['A', 'B', 'C',
			'D', 'E', 
			'F',
			]

varEquations = [
				'A = logistic(0,1)',
				'B = logistic(0,1)',
				'C = logistic(0,1)',
				'D = .5 * A + .5 * B + logistic(0,1)',
				'E = .5 * B + .5 * C + logistic(0,1)',
				'F = .5 * D + .5 * E + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', []), ('C', []), ('D', ['A','B']), ('E', ['B','C']), ('F', ['D', 'E']),
				]
# -----------------------------------------------------------------------------