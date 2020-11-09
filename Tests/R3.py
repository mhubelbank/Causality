# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = '3-Chain Test'

varNames = ['A', 'B', 'C']

varEquations = [
				'A = logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = B + logistic(0,1)',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['B']), 
				]
# -----------------------------------------------------------------------------