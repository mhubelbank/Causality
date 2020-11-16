# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Random Walk Single Input -- IV Structure (Gaussian)'

varNames = ['A', 'B', 'C']

last = 0

varEquations = [
				'A = last + normal(0,1)',
				'B = A + normal(0,1)',
				'C = A + normal(0,1)',
				'last = A',
				]
				
validation =  [
				('A', []), ('B', ['A']), ('C', ['A']), 
				]
# -----------------------------------------------------------------------------