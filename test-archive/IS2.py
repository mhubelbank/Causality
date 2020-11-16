# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Random Walk Single Input -- IV Structure'

varNames = ['A', 'B', 'C']

last = 0

varEquations = [
				'A = last + logistic(0,1)',
				'B = A + logistic(0,1)',
				'C = A + logistic(0,1)',
				'last = A',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A']), 
				]
# -----------------------------------------------------------------------------