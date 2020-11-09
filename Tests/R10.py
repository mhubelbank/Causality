# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Cascade Pattern Test (3 Node)'

varNames = ['A', 'B',  'C',
			]

varEquations = [
				'A = logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = .5 * A + .5 * B + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A','B']), 
				]
# -----------------------------------------------------------------------------