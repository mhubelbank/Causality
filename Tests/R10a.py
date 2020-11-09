# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Cascade Pattern Test (4 Node)'

varNames = ['A', 'B',  'C', 'D',
			]

varEquations = [
				'A = logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = .5 * A + .5 * B + logistic(0,1)',
				'D = .33 * A + .33 * B + .33 * C + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A','B']), ('D', ['A', 'B', 'C']),
				]
# -----------------------------------------------------------------------------