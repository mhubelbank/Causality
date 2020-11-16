# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = '7-Chain Test'

varNames = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

varEquations = [
				'A = logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = B + logistic(0,1)',
				'D = C + logistic(0,1)',
				'E = D + logistic(0,1)',
				'F = E + logistic(0,1)',
				'G = F + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['B']), ('D', ['C']), ('E', ['D']), ('F', ['E']), ('G', ['F']), 
				]
# -----------------------------------------------------------------------------