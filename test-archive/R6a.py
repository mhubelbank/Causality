# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Triangle Pattern Test (4 level)'

varNames = ['A', 
			'B', 'C', 
			'D', 'E', 'F', 
			'G', 'H', 'I', 'J',
			]

varEquations = [
				'A = logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = A + logistic(0,1)',
				'D = B + logistic(0,1)',
				'E = .5 * B + .5 * C + logistic(0,1)',
				'F = C + logistic(0,1)',
				'G = D + logistic(0,1)',
				'H = .5 * D + .5 * E + logistic(0,1)',
				'I = .5 * E + .5 * F + logistic(0,1)',
				'J = F + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A']), ('D', ['B']), ('E', ['B','C']), ('F', ['C']), ('G',['D']), ('H',['D','E']), ('I',['E','F']), ('J', ['F']),
				]
# -----------------------------------------------------------------------------