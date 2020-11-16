# ----------------------------------------------------------------------------
# SEM Definitions
#

testDescript = 'Basic V-Structure Test'

varNames = ['A', 'B']

varEquations = [
				'A= logistic(0,1)',
				#'C = logistic(0,1)',
				'B = .5 * A + .5 + logistic(0,1)',

				]
				
validation =  [
				('A', []), ('B', ['A','C']), ('C', []), 
				]
# -----------------------------------------------------------------------------