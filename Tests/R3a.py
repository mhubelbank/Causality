# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = '4-Chain Test'

varNames = ['A', 'B', 'C', 'D']

varEquations = [
				'A = logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = B + logistic(0,1)',
				'D = C + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['B']), ('D', ['C']), 
				]
# -----------------------------------------------------------------------------